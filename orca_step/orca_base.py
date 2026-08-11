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
import re
import shlex
import shutil

import seamm
import seamm_exec
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("ORCA")

# Byte multipliers for parsing a memory string (SI 'GB' and binary 'GiB').
_MEMORY_UNITS = {
    "": 1,
    "b": 1,
    "k": 1000,
    "kb": 1000,
    "ki": 1024,
    "kib": 1024,
    "m": 1000**2,
    "mb": 1000**2,
    "mi": 1024**2,
    "mib": 1024**2,
    "g": 1000**3,
    "gb": 1000**3,
    "gi": 1024**3,
    "gib": 1024**3,
    "t": 1000**4,
    "tb": 1000**4,
    "ti": 1024**4,
    "tib": 1024**4,
}


def _dehumanize_bytes(text):
    """Parse a memory string ('3 GB', '512MB', '2Gi', or a bare number of bytes)
    into an integer number of bytes. Raises ValueError on anything unparseable."""
    m = re.fullmatch(r"\s*([0-9]*\.?[0-9]+)\s*([a-zA-Z]*)\s*", str(text))
    if not m:
        raise ValueError(f"cannot parse memory '{text}'")
    unit = m.group(2).lower()
    if unit not in _MEMORY_UNITS:
        raise ValueError(f"unknown memory unit '{m.group(2)}'")
    return int(float(m.group(1)) * _MEMORY_UNITS[unit])


def _library_path_vars(n_cores, library_path, environ=None):
    """The dynamic-library search variables to set for a parallel ORCA run.

    Returns a list of ``(name, value)`` pairs prepending ``library_path`` to
    ``DYLD_LIBRARY_PATH`` and ``LD_LIBRARY_PATH`` (preserving any existing
    value). Empty for a serial run or when no path is given. ORCA's OpenMPI
    launcher needs these to find e.g. ``libmpi.40.dylib``.
    """
    if n_cores <= 1 or not library_path:
        return []
    environ = os.environ if environ is None else environ
    pairs = []
    for var in ("DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
        existing = environ.get(var, "")
        pairs.append((var, library_path + (os.pathsep + existing if existing else "")))
    return pairs


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

        # User-configurable run options (the [orca-step] section of the main
        # seamm.ini). How to find/launch ORCA -- its executable path and the
        # OpenMPI library directory -- lives in ~/SEAMM/orca.ini instead (see
        # data/orca.ini and _orca_config).
        parser.add_argument(
            parser_name,
            "--ncores",
            default="available",
            help=(
                "How many cores/processes ORCA may use (via %%pal). 'available' "
                "(the default) uses all cores the job/machine provides; give an "
                "integer to cap it, or '1' to force serial. Parallel ORCA needs a "
                "matching OpenMPI runtime, whose location is set with "
                "'library-path' in orca.ini."
            ),
        )
        parser.add_argument(
            parser_name,
            "--memory",
            default="available",
            help=(
                "Memory ORCA may use per process (its %%maxcore). 'available' (the "
                "default) scales to the memory-per-core of the machine; 'all' uses "
                "the whole node divided among the processes; or give an explicit "
                "amount such as '3 GB' (per process)."
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
    def geometry_block(
        self, configuration, charge, multiplicity, atom_indices=None, ghost_atoms=None
    ):
        """Return the ORCA coordinate block for `configuration`.

        Parameters
        ----------
        atom_indices : Sequence[int] | None
            0-based atom indices to include, in this order. ``None`` (the
            default) is every atom, in the configuration's own order -- a BSSE
            sub-job passes a fragment's (sub-)list instead.
        ghost_atoms : Container[int] | None
            The subset of `atom_indices` to write as ORCA ghost centres (the
            element symbol with a trailing ``:``, e.g. ``O:`` -- basis
            functions only, no nucleus/electrons). ``None``/empty writes every
            atom as real.
        """
        lines = [f"* xyz {charge} {multiplicity}"]
        symbols = configuration.atoms.symbols
        xyzs = configuration.atoms.get_coordinates(fractionals=False, in_cell=True)
        indices = range(len(symbols)) if atom_indices is None else atom_indices
        ghosts = set(ghost_atoms) if ghost_atoms else set()
        for i in indices:
            symbol = symbols[i]
            x, y, z = xyzs[i]
            label = f"{symbol}:" if i in ghosts else symbol
            lines.append(f"{label:4s} {x:15.8f} {y:15.8f} {z:15.8f}")
        lines.append("*")
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------
    def run_orca(self, keyword_line, extra_blocks="", extra_files=None, make_wfx=False):
        """Write the ORCA input for the node's own configuration, run ORCA,
        and return the parsed results.

        A thin wrapper around :meth:`run_orca_job` for the common case: the
        whole of the node's own system, at its own charge/multiplicity, in
        the node's own directory. See :meth:`run_orca_job` for a sub-job
        against an explicit geometry/charge/multiplicity/directory (e.g. one
        of a BSSE correction's counterpoise sub-jobs).

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
        _, configuration = self.get_system_configuration(None)
        return self.run_orca_job(
            keyword_line,
            configuration,
            configuration.charge,
            configuration.spin_multiplicity,
            extra_blocks=extra_blocks,
            extra_files=extra_files,
            make_wfx=make_wfx,
        )

    def run_orca_job(
        self,
        keyword_line,
        configuration,
        charge,
        multiplicity,
        atom_indices=None,
        ghost_atoms=None,
        directory=None,
        extra_blocks="",
        extra_files=None,
        make_wfx=False,
    ):
        """Write an ORCA input for an explicit geometry/charge/multiplicity,
        run ORCA, and return the parsed results.

        The primitive ``run_orca`` and a BSSE counterpoise correction's
        per-fragment sub-jobs both drive: geometry/charge/multiplicity come
        from arguments instead of ``self.get_system_configuration()``/
        ``self.directory``, since a BSSE sub-job is a different atom subset
        (and, for a ghost-augmented fragment, a different charge) than the
        node's own system, run in its own sub-directory so the 2N + 1 jobs of
        one BSSE correction do not collide.

        Parameters
        ----------
        keyword_line : str
            The contents of the ORCA "!" simple-input line (without the "!").
        configuration : molsystem.Configuration
            The system to take atoms/coordinates from.
        charge, multiplicity : int
            The charge/multiplicity for *this job* -- not necessarily
            ``configuration.charge``/``configuration.spin_multiplicity`` (a
            BSSE fragment sub-job has its own).
        atom_indices : Sequence[int] | None
            0-based atom indices to include, in this order. ``None`` (the
            default) is every atom, in the configuration's own order.
        ghost_atoms : Container[int] | None
            The subset of `atom_indices` to write as ORCA ghost centres. See
            :meth:`geometry_block`.
        directory : str | Path | None
            Where to run this job. ``None`` (the default) is ``self.directory``
            -- the node's own job directory, as ``run_orca`` uses.
        extra_blocks : str
            Any additional ``%`` blocks to place before the geometry.
        extra_files : dict | None
            Extra input files to write into the run directory.
        make_wfx : bool
            After ORCA finishes, run ``orca_2aim`` to convert the retained
            density into an AIMPAC ``orca.wfx`` in *this job's* directory.

        Returns
        -------
        dict
            Parsed results, at least ``{"energy": <E_h>, "success": bool}``.
        """
        run_directory = (
            Path(directory) if directory is not None else Path(self.directory)
        )
        run_directory.mkdir(parents=True, exist_ok=True)

        # Resources (cores + per-process memory), shared with run_orca_compound.
        n_cores, memory_mb = self._resources()

        lines = [f"! {keyword_line.strip()}"]
        if n_cores > 1:
            lines.append(f"%pal nprocs {n_cores} end")
        lines.append(f"%maxcore {memory_mb}")
        if extra_blocks.strip() != "":
            lines.append(extra_blocks.rstrip())
        lines.append(
            self.geometry_block(
                configuration,
                charge,
                multiplicity,
                atom_indices=atom_indices,
                ghost_atoms=ghost_atoms,
            )
        )
        lines.append("")
        input_text = "\n".join(lines)

        files = {"orca.inp": input_text}
        if extra_files:
            files.update(extra_files)
        logger.debug(f"orca.inp ({run_directory}):\n" + input_text)

        config = self._orca_config()
        env, lib_prefix = self._mpi_env(n_cores, config)

        # ORCA must be invoked by its full path so it can find its sub-programs.
        # When a wavefunction file is wanted, chain orca_2aim (which lives next
        # to the orca binary) in the same shell so it runs in this directory
        # right after ORCA, reading the just-written orca.gbw/orca.densities.
        cmd = lib_prefix + ["{code}", "orca.inp", ">", "orca.out", "2>", "orca.err"]
        return_files = [
            "orca.out",
            "orca.err",
            "orca.gbw",
            "orca.xyz",
            "orca.hess",
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
            directory=str(run_directory),
            files=files,
            # Wildcards: ORCA writes orca.bibtex (suggested citations) and
            # orca.property.txt (and other *.txt depending on options); the
            # executor discards files not requested. orca.engrad appears only
            # when gradients are requested (! EnGrad).
            return_files=return_files,
            # None: run in node-local scratch under a scheduler (SLURM's job
            # dir is NFS, unsafe for ORCA's MPI scratch I/O), in place
            # otherwise. See molssi-seamm/orca_step#20.
            in_situ=None,
            shell=True,
            env=env,
        )
        if not result:
            raise RuntimeError(f"There was an error running ORCA in {run_directory}.")
        self._report_run_location(result, run_directory)

        return self._parse_output(run_directory / "orca.out")

    def _report_run_location(self, result, directory):
        """Note in step.out where ORCA actually ran -- the job directory, or
        (under a scheduler) node-local scratch, with only results copied back.
        See molssi-seamm/orca_step#20: this distinction matters for diagnosing
        NFS-related failures and for knowing where any leftover scratch
        (e.g. after a crash) actually landed.
        """
        if result.get("in_situ", True):
            printer.normal(
                __(
                    f"Ran ORCA directly in the job directory, {directory}.",
                    indent=self.indent + 4 * " ",
                )
            )
        else:
            printer.normal(
                __(
                    "Ran ORCA in node-local scratch "
                    f"({result.get('directory')}), not the job directory, "
                    "because this job is running under a batch scheduler and "
                    "ORCA's MPI scratch I/O is not safe on NFS; only the "
                    "requested result files were copied back here.",
                    indent=self.indent + 4 * " ",
                )
            )

    def _resources(self):
        """Resolve (n_cores, memory_mb) for an ORCA run from the executor's
        computational environment and the [orca-step] user options.

        run_orca (and run_orca_compound) are called by a sub-step, whose options
        live on its parent (the main ORCA node), mirroring the MOPAC step.
        """
        ce = seamm_exec.computational_environment()
        options = self.parent.options
        seamm_options = self.parent.global_options

        # --- Cores/processes (ORCA %pal) ---
        available_cores = max(1, int(ce.get("NTASKS", 1) or 1))
        ncores_opt = str(options.get("ncores", "available")).strip().lower()
        if ncores_opt in ("available", "default", "all", ""):
            n_cores = available_cores
        else:
            try:
                n_cores = min(available_cores, int(ncores_opt))
            except ValueError:
                n_cores = available_cores
        # Respect a global cap from the SEAMM options, if one is set.
        global_ncores = str(seamm_options.get("ncores", "available")).strip().lower()
        if global_ncores not in ("available", "default", "all", ""):
            try:
                n_cores = min(n_cores, int(global_ncores))
            except ValueError:
                pass
        n_cores = max(1, n_cores)

        # --- Memory per process (ORCA %maxcore, in MB) ---
        mem_per_node = int(ce.get("MEM_PER_NODE", 0) or 0)  # bytes
        mem_per_cpu = int(ce.get("MEM_PER_CPU", 0) or 0)  # bytes
        memory_opt = str(options.get("memory", "available")).strip().lower()
        if memory_opt in ("available", "default", ""):
            # ~80% of the per-core memory, leaving headroom for the OS/driver.
            per_proc_bytes = int(0.8 * mem_per_cpu) if mem_per_cpu else 0
        elif memory_opt == "all":
            per_proc_bytes = int(mem_per_node / n_cores) if mem_per_node else 0
        else:
            try:
                per_proc_bytes = _dehumanize_bytes(options["memory"])
            except ValueError:
                per_proc_bytes = 0
        # ORCA's %maxcore is per-process MB; fall back to a conservative default
        # if the machine's memory could not be determined.
        memory_mb = int(per_proc_bytes / 1_000_000) if per_proc_bytes else 2000
        memory_mb = max(256, memory_mb)
        return n_cores, memory_mb

    def _mpi_env(self, n_cores, config):
        """Return ``(env, lib_prefix)`` for launching a (possibly parallel) ORCA.

        Parallel ORCA needs a matching OpenMPI runtime. Two things must line up
        and they are handled differently:

          1. mpirun (PATH): ORCA launches its workers with whatever ``mpirun`` it
             finds on PATH; it MUST be the same OpenMPI as the libraries the
             workers link, or the ranks corrupt each other's data (ORCA aborts
             with a "BLAS-ERROR: incompatible matrices"). We prepend the OpenMPI
             bin (the sibling of the lib dir given by library-path) so a
             different mpirun on PATH (e.g. a newer Homebrew OpenMPI) is not used.
             PATH is an ordinary variable, so it reaches ORCA's children.

          2. libmpi (dynamic-loader path): on Linux, LD_LIBRARY_PATH is honored
             and inherited, so exporting it is enough. On macOS, ORCA does NOT
             pass DYLD_* to the MPI sub-processes it spawns (and SIP strips DYLD_*
             through /bin/sh anyway), so the OpenMPI libraries must instead be on
             dyld's default search path -- e.g. symlink libmpi.*.dylib into
             /usr/local/lib. That is a one-time machine setup, documented in the
             User Guide; nothing here can substitute for it. We still export the
             loader variables for Linux.

          3. Core binding: ORCA calls its own ``mpirun`` internally once it
             parses ``%pal``, so we never see that command line -- but its
             environment is ours to set. On an unmanaged, interactive host
             (no SLURM etc. coordinating which cores belong to which job),
             every independent ``mpirun`` invocation applies OpenMPI's
             default binding policy (bind each rank to a core, starting from
             a low core number) with no knowledge of other concurrently-
             running ORCA jobs, so several jobs launched at once all pile
             onto the same one or two cores instead of spreading across the
             machine. Disabling binding (``OMPI_MCA_hwloc_base_binding_policy
             =none``, the environment-variable spelling of ``mpirun
             --bind-to none``, same on Linux and macOS -- no ``DYLD_*``/SIP
             wrinkle here) lets the OS scheduler load-balance ranks from
             concurrent jobs across all cores instead. Under a real
             scheduler this is unnecessary and can be counter-productive: a
             SLURM allocation already restricts the job to specific cores
             (cgroups), so we leave OpenMPI's binding on in that case, the
             same ``SLURM_JOB_ID`` check seamm_exec's ``in_situ`` auto-
             detection and ``computational_environment()`` use.

        The OpenMPI library directory is part of *how to run ORCA*, so it comes
        from the executor config (~/SEAMM/orca.ini), not the user-facing
        [orca-step] options. configparser lower-cases keys.
        """
        env = {}
        lib_prefix = []
        library_path = config.get("library-path", "") or ""
        quoted_library_path = shlex.quote(library_path)
        for var, value in _library_path_vars(n_cores, library_path):
            env[var] = value
            # Prepend at shell *runtime*, not with this Python-computed
            # snapshot of os.environ: an installation=modules code
            # (seamm_exec.Local.exec()) runs `module load ...` earlier in
            # this same generated script, which Python's os.environ can't
            # see yet at this point -- a static value here would silently
            # discard whatever that module load just added (confirmed for
            # real: ORCA's own liborca_tools_*.so disappeared from a
            # parallel run's LD_LIBRARY_PATH this way).
            #
            # Deliberately brace-free (an if/then/else, not
            # ``${VAR:+...}`` parameter expansion): this text is later run
            # through seamm_exec.local.Local.exec()'s own
            # ``command.format(**config, **ce)`` (for ITS "{code}"-style
            # placeholders), called in a loop *until the output stops
            # changing* -- a single-braced ``${VAR:+...}`` here survives
            # exactly one pass (mistaken for one of its own fields on the
            # next), and even doubling the braces only buys one extra
            # pass before the same KeyError (confirmed for real, both
            # ways). A construct with no braces at all is immune
            # regardless of how many passes run.
            lib_prefix.append(
                f'if [ -n "${var}" ]; then '
                f"export {var}={quoted_library_path}:${var}; "
                f"else export {var}={quoted_library_path}; fi;"
            )
        if n_cores > 1 and library_path:
            bindir = Path(library_path).expanduser().parent / "bin"
            if bindir.is_dir():
                lib_prefix.insert(0, f"export PATH={shlex.quote(str(bindir))}:$PATH;")
        if n_cores > 1 and "SLURM_JOB_ID" not in os.environ:
            env["OMPI_MCA_hwloc_base_binding_policy"] = "none"
        return env, lib_prefix

    def run_orca_compound(
        self,
        compound_block,
        extra_files=None,
        engrad="result.engrad",
        make_wfx=False,
        wfx_step="orca_Compound_5",
    ):
        """Write and run an ORCA *Compound* job, returning ``(energy, gradient)``.

        Unlike ``run_orca`` (a single ``! keywords`` + inline geometry job), a
        Compound job runs several calculations from a script and writes its own
        result. The ``%pal``/``%maxcore`` preamble is emitted here (applying to
        every sub-calculation); the caller supplies the ``%Compound ... end``
        block and any files it references (the ``.cmp`` script, the geometry).

        Parameters
        ----------
        compound_block : str
            The ``%Compound "..." ... end`` block for ``orca.inp``.
        extra_files : dict | None
            Files the Compound job reads (e.g. ``{"bssegradient.cmp": ...,
            "bsse.xyz": ...}``), written into the run directory.
        engrad : str | None
            The EnGrad file the Compound writes with the final energy and
            gradient. Pass ``None`` for an energy-only job that writes no EnGrad
            (e.g. a method with no analytic gradient); then ``(None, None)`` is
            returned and the caller gets the energy elsewhere.

        Returns
        -------
        tuple(float | None, list | None)
            The energy (E_h) and gradient (``[n_atoms][3]``, E_h/bohr) from the
            Compound's EnGrad file, or ``(None, None)`` when ``engrad`` is None.
        """
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        n_cores, memory_mb = self._resources()
        config = self._orca_config()
        env, lib_prefix = self._mpi_env(n_cores, config)

        lines = []
        if n_cores > 1:
            lines.append(f"%pal nprocs {n_cores} end")
        lines.append(f"%maxcore {memory_mb}")
        lines.append(compound_block.rstrip())
        lines.append("")
        input_text = "\n".join(lines)

        files = {"orca.inp": input_text}
        if extra_files:
            files.update(extra_files)
        logger.debug("orca.inp (Compound):\n" + input_text)

        return_files = ["orca.out", "orca.err", "*.bibtex", "*.txt", "*.engrad"]
        if engrad:
            return_files.append(engrad)
        cmd = lib_prefix + ["{code}", "orca.inp", ">", "orca.out", "2>", "orca.err"]
        # Convert the dimer step's retained density into a .wfx (as run_orca does
        # for the Energy substep), for a following Atomic Charges (DDEC6) step.
        # wfx_step is the Compound sub-job that kept its density (the dimer).
        if make_wfx:
            orca_2aim = str(Path(config["code"]).parent / "orca_2aim")
            cmd += [
                "&&",
                orca_2aim,
                wfx_step,
                ">",
                "orca_2aim.out",
                "2>&1",
                "&&",
                "cp",
                f"{wfx_step}.wfx",
                "orca.wfx",
            ]
            return_files += ["orca.wfx", "orca_2aim.out"]
        result = self.flowchart.executor.run(
            cmd=cmd,
            config=config,
            directory=self.directory,
            files=files,
            return_files=return_files,
            # See the comment in run_orca -- molssi-seamm/orca_step#20.
            in_situ=None,
            shell=True,
            env=env,
        )
        if not result:
            raise RuntimeError("There was an error running the ORCA Compound job.")
        self._report_run_location(result, directory)

        # Read the corrected gradient from the EnGrad file the Compound wrote,
        # when one was requested (energy-only jobs pass engrad=None).
        if engrad:
            return self._read_engrad(directory / engrad)
        return None, None

    @staticmethod
    def _read_engrad(path):
        """Read ``(energy, gradient)`` from an ORCA EnGrad file.

        Format (comment lines start with ``#``; blank lines are ignored): the
        atom count, then the energy (E_h), then ``3*n`` gradient components
        (E_h/bohr). Returns ``(None, None)`` if it cannot be parsed.
        """
        path = Path(path)
        if not path.exists():
            return None, None
        lines = [
            ln.strip()
            for ln in path.read_text().splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        try:
            n = int(lines[0])
            energy = float(lines[1])
            vals = [float(x) for x in lines[2 : 2 + 3 * n]]  # noqa: E203
        except (ValueError, IndexError):
            return None, None
        if len(vals) != 3 * n:
            return energy, None
        gradient = [vals[3 * i : 3 * i + 3] for i in range(n)]  # noqa: E203
        return energy, gradient

    def _orca_config(self):
        """Resolve the ORCA executable configuration.

        Reads the executor's section of ``<root>/orca.ini`` (the ``code`` = full
        path to ORCA, and the optional ``library-path`` for parallel runs);
        falls back to locating ``orca`` on the PATH. ORCA must be invoked by its
        full path so it can find its sub-programs.
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

        # Fall back to finding ORCA on the PATH.
        code = shutil.which("orca")
        if code is None:
            raise RuntimeError(
                "Could not find the 'orca' executable. Put it on your PATH, or "
                "set 'code' (the full path) in the relevant section of "
                f"{ini_path}."
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
