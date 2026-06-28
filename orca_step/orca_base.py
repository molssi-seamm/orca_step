# -*- coding: utf-8 -*-

"""Shared base class for the ORCA step and its sub-steps.

Holds the command-line/seamm.ini parser, the geometry block, and the routine
that writes the ORCA input, runs ORCA through the flowchart executor, and parses
the output. The main ``ORCA`` node and the ``Energy`` sub-step both inherit from
this (``Optimization`` inherits from ``Energy``), mirroring ``MOPACBase`` in the
MOPAC step.
"""

import configparser
import importlib
import logging
import os
from pathlib import Path
import shutil

import seamm
import seamm_exec
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("ORCA")

# Hartree -> the energy unit ORCA reports (E_h); kept explicit for clarity.


class ORCABase(seamm.Node):
    """Common functionality for ORCA nodes."""

    def create_parser(self):
        """Set up the command-line / seamm.ini parser for the ORCA step.

        All ORCA nodes share the ``[orca-step]`` section (step_type), so the
        options are registered once here.
        """
        parser_name = self.step_type  # 'orca-step'
        parser = self.flowchart.parser

        parser_exists = parser.exists(parser_name)

        result = super().create_parser(name=parser_name)

        if parser_exists:
            return result

        parser.add_argument(
            parser_name,
            "--orca-path",
            default="",
            help="the path to the ORCA executable, if not on the PATH",
        )
        parser.add_argument(
            parser_name,
            "--ncores",
            default="1",
            help=(
                "How many cores/processes ORCA may use (via %%pal). Default is "
                "serial (1); parallel ORCA needs a matching OpenMPI runtime, set "
                "via --library-path."
            ),
        )
        parser.add_argument(
            parser_name,
            "--library-path",
            default="",
            help=(
                "Directory containing ORCA's required OpenMPI libraries "
                "(e.g. libmpi.40.dylib). Prepended to the (DY)LD_LIBRARY_PATH "
                "when running ORCA in parallel."
            ),
        )
        parser.add_argument(
            parser_name,
            "--max-atoms-to-print",
            default=25,
            help="Maximum number of atoms for which to print detailed results.",
        )

        return result

    @property
    def is_runable(self):
        """Whether this node actually runs (vs. only contributing input)."""
        return True

    @property
    def version(self):
        """The semantic version of this module."""
        import orca_step

        return orca_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        import orca_step

        return orca_step.__git_revision__

    # ------------------------------------------------------------------
    # Input generation
    # ------------------------------------------------------------------
    def geometry_block(self, configuration, charge, multiplicity):
        """Return the ORCA coordinate block for the current configuration."""
        lines = [f"* xyz {charge} {multiplicity}"]
        symbols = configuration.atoms.symbols
        xyzs = configuration.atoms.get_coordinates(fractionals=False, in_cell=True)
        for symbol, (x, y, z) in zip(symbols, xyzs):
            lines.append(f"{symbol:2s} {x:15.8f} {y:15.8f} {z:15.8f}")
        lines.append("*")
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------
    def run_orca(self, keyword_line, extra_blocks="", extra_files=None, make_wfx=False):
        """Write the ORCA input, run ORCA, and return the parsed results.

        Parameters
        ----------
        keyword_line : str
            The contents of the ORCA "!" simple-input line (without the "!").
        extra_blocks : str
            Any additional ``%`` blocks to place before the geometry.
        extra_files : dict | None
            Extra input files to write into the run directory (e.g. an external
            basis file referenced by ``%basis GTOName ... end``).
        make_wfx : bool
            After ORCA finishes, run ``orca_2aim`` (shipped alongside ORCA) to
            convert the density (retained by the ``keepdensity`` keyword) into an
            AIMPAC ``orca.wfx`` wavefunction file for a following Atomic Charges
            step. Requires the ORCA input to include ``keepdensity``.

        Returns
        -------
        dict
            Parsed results, at least ``{"energy": <E_h>, "success": bool}``.
        """
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        _, configuration = self.get_system_configuration(None)
        charge = configuration.charge
        multiplicity = configuration.spin_multiplicity

        # Resources. run_orca is called by a sub-step, whose options live on its
        # parent (the main ORCA node), mirroring the MOPAC step.
        ce = seamm_exec.computational_environment()
        options = self.parent.options
        n_cores = ce["NTASKS"]
        if options.get("ncores", "available") not in ("available", "default"):
            try:
                n_cores = min(n_cores, int(options["ncores"]))
            except ValueError:
                pass
        n_cores = max(1, n_cores)
        memory_mb = 2000  # %maxcore is per-process MB; conservative default

        lines = [f"! {keyword_line.strip()}"]
        if n_cores > 1:
            lines.append(f"%pal nprocs {n_cores} end")
        lines.append(f"%maxcore {memory_mb}")
        if extra_blocks.strip() != "":
            lines.append(extra_blocks.rstrip())
        lines.append(self.geometry_block(configuration, charge, multiplicity))
        lines.append("")
        input_text = "\n".join(lines)

        files = {"orca.inp": input_text}
        if extra_files:
            files.update(extra_files)
        logger.debug("orca.inp:\n" + input_text)

        config = self._orca_config()

        # For parallel runs ORCA's MPI launcher needs the matching OpenMPI
        # libraries on the dynamic-library path.
        env = {}
        library_path = options.get("library-path", "") or ""
        if n_cores > 1 and library_path != "":
            for var in ("DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
                existing = os.environ.get(var, "")
                env[var] = library_path + (os.pathsep + existing if existing else "")

        # ORCA must be invoked by its full path so it can find its sub-programs.
        # When a wavefunction file is wanted, chain orca_2aim (which lives next
        # to the orca binary) in the same shell so it runs in this directory
        # right after ORCA, reading the just-written orca.gbw/orca.densities.
        cmd = ["{code}", "orca.inp", ">", "orca.out", "2>", "orca.err"]
        return_files = [
            "orca.out",
            "orca.err",
            "orca.gbw",
            "*.bibtex",
            "*.txt",
            "*.engrad",
        ]
        if make_wfx:
            orca_2aim = str(Path(config["code"]).parent / "orca_2aim")
            cmd += ["&&", orca_2aim, "orca", ">", "orca_2aim.out", "2>&1"]
            return_files += ["orca.wfx", "orca_2aim.out"]

        result = self.flowchart.executor.run(
            cmd=cmd,
            config=config,
            directory=self.directory,
            files=files,
            # Wildcards: ORCA writes orca.bibtex (suggested citations) and
            # orca.property.txt (and other *.txt depending on options); the
            # executor discards files not requested. orca.engrad appears only
            # when gradients are requested (! EnGrad).
            return_files=return_files,
            in_situ=True,
            shell=True,
            env=env,
        )
        if not result:
            raise RuntimeError("There was an error running ORCA.")

        return self._parse_output(directory / "orca.out")

    def _orca_config(self):
        """Resolve the ORCA executable configuration.

        Reads ``<root>/orca.ini`` if present (like mopac.ini); otherwise locates
        ORCA on the PATH (or via --orca-path) and uses its full path, which ORCA
        requires to launch its sub-programs.
        """
        executor_type = self.flowchart.executor.name
        seamm_options = self.parent.global_options
        ini_dir = Path(seamm_options["root"]).expanduser()
        ini_path = ini_dir / "orca.ini"

        full_config = configparser.ConfigParser()
        if ini_path.exists():
            full_config.read(ini_path)

        if (
            executor_type in full_config
            and full_config[executor_type].get("code", "") != ""
        ):
            return dict(full_config.items(executor_type))

        # Fall back to locating ORCA ourselves.
        options = self.parent.options
        code = options.get("orca-path", "") or ""
        if code != "":
            code = str(Path(code).expanduser() / "orca")
        else:
            code = shutil.which("orca")
        if code is None:
            raise RuntimeError(
                "Could not find the 'orca' executable. Put it on your PATH, set "
                "--orca-path, or add it to orca.ini."
            )
        return {"code": code, "installation": "local"}

    # ------------------------------------------------------------------
    # Output parsing
    # ------------------------------------------------------------------
    def _parse_output(self, path):
        """Parse ORCA output: final single-point energy and success."""
        data = {"success": False, "energy": None}
        if not path.exists():
            raise RuntimeError(f"ORCA produced no output ({path}).")
        text = path.read_text()
        for line in text.splitlines():
            if "FINAL SINGLE POINT ENERGY" in line:
                try:
                    data["energy"] = float(line.split()[-1])
                except (ValueError, IndexError):
                    pass
        data["success"] = "ORCA TERMINATED NORMALLY" in text
        if not data["success"]:
            raise RuntimeError(
                f"ORCA did not terminate normally; see {path} and orca.err."
            )
        return data

    @staticmethod
    def _add_properties():
        """Register this package's properties with molsystem, if present."""
        import molsystem

        path = importlib.resources.files("orca_step") / "data"
        csv_file = path / "properties.csv"
        if path.exists() and csv_file.exists():
            molsystem.add_properties_from_file(csv_file)


# Run once at import to register any custom properties.
try:
    ORCABase._add_properties()
except Exception:  # pragma: no cover - non-fatal
    pass

# Re-export for convenience
__all__ = ["ORCABase", "printer", "job", "os", "__"]
