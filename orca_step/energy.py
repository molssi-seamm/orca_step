# -*- coding: utf-8 -*-

"""An ORCA single-point energy sub-step."""

import csv
import json
import logging
from pathlib import Path
import pprint  # noqa: F401
import re
import textwrap

import bibtexparser
from bibtexparser.bwriter import BibTexWriter
from tabulate import tabulate

import orca_step
import seamm
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("ORCA")


class Energy(orca_step.ORCABase):
    """A single-point energy with ORCA.

    See Also
    --------
    TkEnergy, EnergyParameters, Optimization
    """

    def __init__(self, flowchart=None, title="Energy", extension=None, logger=logger):
        """Initialize the Energy sub-step."""
        logger.debug(f"Creating ORCA Energy {self}")

        super().__init__(
            flowchart=flowchart,
            title=title,
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._calculation = "energy"
        self._metadata = orca_step.metadata
        self.parameters = orca_step.EnergyParameters()

    def description_text(self, P=None):
        """Describe what this sub-step will do."""
        if not P:
            P = self.parameters.values_to_dict()

        use_mc = P["use model chemistry"]
        if not isinstance(use_mc, bool):
            use_mc = use_mc == "yes"

        if use_mc:
            text = (
                "Single-point energy with ORCA using the model chemistry from a "
                "preceding Model Chemistry step."
            )
        else:
            method, basis = self._resolve_method_basis(P)
            text = f"Single-point energy with ORCA at {method}/{basis}."
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    @staticmethod
    def _basis_name(value):
        """The basis name from the 'basis' parameter, which is a
        ``{"name", "elements"}`` dict (from the shared BasisSetField) or, for
        older flowcharts, a plain string."""
        if isinstance(value, dict):
            return value.get("name", "") or ""
        return value or ""

    def _resolve_method_basis(self, P):
        """Resolve (method, basis name) from the model chemistry (if used) or the
        explicit method/basis parameters."""
        use_mc = P["use model chemistry"]
        if not isinstance(use_mc, bool):
            use_mc = use_mc == "yes"
        if use_mc:
            return self._method_basis_from_model_chemistry(P)
        method = P["method"]
        # For DFT the '!' keyword is the chosen functional, not the word "DFT".
        if method == "DFT":
            method = P["functional"]
        return method, self._basis_name(P["basis"])

    def _using_bse(self, P):
        """Whether the orbital basis comes from the Basis Set Exchange -- either
        the 'basis source' toggle is set to it, or the basis name uses the
        'bse:NAME' shorthand (which forces the BSE regardless of the toggle)."""
        _, basis = self._resolve_method_basis(P)
        if P.get("basis source") == "Basis Set Exchange":
            return True
        return basis.strip().lower().startswith("bse:")

    @staticmethod
    def _strip_bse(basis):
        """The basis name without any 'bse:' shorthand prefix."""
        if basis.strip().lower().startswith("bse:"):
            return basis.split(":", 1)[1].strip()
        return basis

    def keyword_line(self, P):
        """Build the ORCA '!' keyword line.

        Method and basis come from the model chemistry (if used) or the explicit
        parameters. When the basis source is the Basis Set Exchange the basis is
        embedded as a file instead, so its name is omitted from the '!' line. The
        auxiliary basis and extra keywords always come from this node.
        """
        method, basis = self._resolve_method_basis(P)

        keywords = [method]
        # When the basis comes from the Basis Set Exchange (either via the
        # 'basis source' toggle or a 'bse:NAME' basis), it is embedded as a file
        # by extra_input, so its name is left off the '!' line.
        if not self._using_bse(P):
            keywords.append(basis)
        aux = P["auxiliary basis"]
        if aux and aux.lower() != "none":
            keywords.append(aux)
        # Compute the Cartesian gradient when the gradients result is requested
        # (e.g. by a driver step such as Reaction Path or Thermochemistry). Use
        # the analytic gradient (EnGrad) when ORCA has one for this method, else
        # fall back to the numerical gradient (NumGrad).
        if self._wants_gradients(P):
            keywords.append(self._gradient_keyword(P))
        # Retain the density so orca_2aim can write a .wfx for a following
        # Atomic Charges step.
        if P.get("save wavefunction", "no") == "yes":
            keywords.append("keepdensity")
        extra = P["extra keywords"].strip()
        if extra:
            keywords.append(extra)
        return " ".join(k for k in keywords if k)

    @staticmethod
    def _wants_gradients(P):
        """Whether the energy gradient was requested in the results."""
        results = P.get("results") or {}
        return "gradients" in results

    def _gradient_availability(self, P):
        """Whether ORCA has an ANALYTIC nuclear gradient for the resolved
        method/functional ('analytic') or only a numerical one ('numeric').

        For DFT this is a property of the functional (e.g. wB97M(2)/wB97X-2 are
        non-self-consistent -> 'numeric'); for the other methods it comes from
        metadata['methods'] (the (T) in (DLPNO-)CCSD(T) has no analytic
        gradient). Unknown methods default to 'analytic'.
        """
        md = orca_step.metadata
        if self._is_dft(P):
            functional, _ = self._resolve_method_basis(P)
            rec = md["functionals"].get(functional)
            if rec is None:
                rec = next(
                    (
                        r
                        for k, r in md["functionals"].items()
                        if k.upper() == functional.upper()
                    ),
                    None,
                )
            return (rec or {}).get("gradients", "analytic")
        method, _ = self._resolve_method_basis(P)
        return md["methods"].get(method, {}).get("gradients", "analytic")

    def _gradient_keyword(self, P):
        """The ORCA keyword requesting the gradient: 'EnGrad' when an analytic
        gradient exists, otherwise 'NumGrad' (numerical, much more expensive)."""
        if self._gradient_availability(P) == "analytic":
            return "EnGrad"
        return "NumGrad"

    def extra_input(self, P):
        """Return ``(extra_blocks, extra_files)`` for the ORCA input: the BSE
        basis block (if used), and directives for the optional properties
        (Hirshfeld charges, polarizability)."""
        blocks = []
        files = {}

        if self._using_bse(P):
            _, basis = self._resolve_method_basis(P)
            basis = self._strip_bse(basis)
            _, configuration = self.get_system_configuration(None)
            znums = sorted(set(configuration.atoms.atomic_numbers))
            files["basis.bas"] = self._bse_basis_file(basis, znums)
            blocks.append('%basis GTOName "basis.bas" end')

        if P["Hirshfeld charges"] != "no":
            blocks.append("%output Print[ P_Hirshfeld ] 1 end")
        if P["polarizability"] == "yes":
            blocks.append("%elprop Polar 1 end")

        return "\n".join(blocks), files

    @staticmethod
    def _bse_basis_file(basis_name, atomic_numbers):
        """ORCA-format basis-file contents for `basis_name` and `atomic_numbers`,
        fetched from the Basis Set Exchange (header omitted so ORCA does not read
        the comment lines as keywords)."""
        import basis_set_exchange as bse

        try:
            return bse.get_basis(
                basis_name,
                elements=list(atomic_numbers),
                fmt="orca",
                header=False,
            )
        except Exception as e:
            raise RuntimeError(
                f"Could not get basis set '{basis_name}' from the Basis Set "
                f"Exchange for elements {list(atomic_numbers)}: {e}. Use a basis "
                "available in the BSE, or set the basis source to 'ORCA internal'."
            )

    def _method_basis_from_model_chemistry(self, P):
        """Resolve (method, basis) from the global model chemistry, erroring
        clearly if ORCA cannot provide it.
        """
        if not self.variable_exists("_model_chemistry"):
            raise RuntimeError(
                "'Use the global model chemistry' is selected, but no Model "
                "Chemistry step has defined one (the '_model_chemistry' variable "
                "is not set). Add a Model Chemistry step before this ORCA step, "
                "or turn off 'use model chemistry' and set the method and basis "
                "explicitly."
            )
        mc = self.get_variable("_model_chemistry")
        level = mc.get("level", "") or mc.get("method", "")
        owner = (mc.get("owner") or "").strip()
        mtype = (mc.get("type") or "").strip()
        method = (mc.get("method") or "").strip()
        basis = mc.get("basis") or ""

        # ORCA can evaluate an ORCA-owned or program-agnostic level of theory,
        # but not another program's (e.g. MOPAC/xTB/LAMMPS) or a non-QC type.
        if owner not in ("", "ORCA"):
            raise RuntimeError(
                f"The model chemistry '{level}' is owned by '{owner}', which "
                "ORCA cannot provide. Choose an ORCA (or program-agnostic) model "
                "chemistry, or set the method explicitly in the ORCA step."
            )
        if mtype.upper() in ("SQM", "FF", "MLFF"):
            raise RuntimeError(
                f"The model chemistry '{level}' is of type '{mtype}', which ORCA "
                "does not provide. ORCA handles HF, DFT, MP2, and coupled-cluster "
                "(QC) methods."
            )
        if method == "":
            raise RuntimeError(
                f"The model chemistry '{level}' does not name a method ORCA can " "use."
            )
        # If the model chemistry omits a basis, fall back to this node's basis.
        if basis == "":
            basis = self._basis_name(P["basis"])
        return method, basis

    def run(self, keywords=None):
        """Run the single-point energy.

        Note: this is a sub-step, driven by the main ORCA node, which has already
        set up printing and cited the plug-in. So (like the Gaussian step) it does
        NOT call ``super().run()`` -- doing so would cite the plug-in a second
        time. The working directory is created by ``run_orca``.
        """
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        printer.important(__(self.description_text(P), indent=self.indent))

        keyword_line = self.keyword_line(P)
        if keywords:
            keyword_line += " " + " ".join(keywords)

        # Warn when we had to fall back to a numerical gradient -- it is much
        # more costly (a displaced single point per degree of freedom).
        if self._wants_gradients(P) and self._gradient_availability(P) == "numeric":
            method, _ = self._resolve_method_basis(P)
            printer.important(
                __(
                    f"Note: ORCA has no analytic gradient for {method}, so the "
                    "gradient is computed numerically (NumGrad). This is much "
                    "more expensive than an analytic gradient.",
                    indent=self.indent + 4 * " ",
                )
            )

        extra_blocks, extra_files = self.extra_input(P)
        data = self.run_orca(
            keyword_line,
            extra_blocks=extra_blocks,
            extra_files=extra_files,
            make_wfx=P.get("save wavefunction", "no") == "yes",
        )

        self._data = data
        self._cite_references(P)
        self.analyze(P=P, data=data)

        return self.next()

    # ------------------------------------------------------------------
    # Citations
    # ------------------------------------------------------------------
    def _cite_references(self, P):
        """Cite ORCA's suggested papers (from orca.bibtex), the supporting
        libraries (libint, libXC), and the basis set (via the Basis Set
        Exchange). Citation problems are logged, never fatal.
        """
        for what, fn in (
            ("ORCA suggested references", self._cite_orca_bibtex),
            ("ORCA support libraries", lambda: self._cite_libraries(P)),
            ("the DFT functional", lambda: self._cite_functional(P)),
            ("the basis set", lambda: self._cite_basis(P)),
        ):
            try:
                fn()
            except Exception as e:
                logger.warning(f"Could not cite {what}: {e}")

    def _essential_dois(self):
        """DOIs ORCA lists under 'List of essential papers' in orca.out -- the
        minimum-necessary (main) ORCA citation(s)."""
        out = Path(self.directory) / "orca.out"
        if not out.exists():
            return set()
        text = out.read_text()
        m = re.search(
            r"List of essential papers.*?(?=List of papers to cite|\Z)",
            text,
            re.DOTALL,
        )
        section = m.group(0) if m else ""
        return {
            d.strip().lower().rstrip(".")
            for d in re.findall(r"doi\.org/(\S+)", section)
        }

    def _cite_orca_bibtex(self):
        """Ingest ORCA's orca.bibtex. The essential (main) paper is cited at
        level 1; the rest at level 2, to keep the primary list short."""
        path = Path(self.directory) / "orca.bibtex"
        if not path.exists():
            return
        essential = self._essential_dois()
        writer = BibTexWriter()
        for key, entry in bibtexparser.loads(path.read_text()).entries_dict.items():
            doi = (entry.get("doi", "") or "").strip().lower().rstrip(".")
            level = 1 if (doi and doi in essential) else 2
            self.references.cite(
                raw=writer._entry_to_bibtex(entry),
                alias=f"orca-{key}",
                module="orca_step",
                level=level,
                note="Suggested by ORCA for this calculation.",
            )

    def _cite_libraries(self, P):
        """Cite the supporting libraries. libint is used for the 2-electron
        integrals in every run; libXC only for DFT."""
        aliases = ["libint2"]
        if self._is_dft(P):
            aliases.append("libxc")
        for alias in aliases:
            if alias in self._bibliography:
                self.references.cite(
                    raw=self._bibliography[alias],
                    alias=alias,
                    module="orca_step",
                    level=2,
                    note=f"ORCA uses the {alias} library.",
                )

    def _cite_functional(self, P):
        """For a DFT run, cite the functional's references (mined from the ORCA
        manual; keys in metadata['dft functionals'], bibtex in references.bib)."""
        if not self._is_dft(P):
            return
        method, _ = self._resolve_method_basis(P)
        funcs = orca_step.metadata.get("dft functionals", {})
        # Case-insensitive match on the ORCA functional keyword.
        keys = funcs.get(method)
        if keys is None:
            keys = next(
                (v for k, v in funcs.items() if k.upper() == method.upper()), None
            )
        for key in keys or []:
            if key in self._bibliography:
                self.references.cite(
                    raw=self._bibliography[key],
                    alias=key,
                    module="orca_step",
                    level=2,
                    note=f"Reference for the {method} density functional.",
                )

    def _is_dft(self, P):
        """Whether the requested method is a DFT functional (so libXC is used)."""
        use_mc = P["use model chemistry"]
        if not isinstance(use_mc, bool):
            use_mc = use_mc == "yes"
        if use_mc:
            if self.variable_exists("_model_chemistry"):
                mc = self.get_variable("_model_chemistry")
                return (mc.get("type", "") or "").upper() == "DFT"
            return False
        method = P["method"]
        if method == "DFT":
            return True
        # Legacy flowcharts (pre-catalog) stored the functional keyword directly
        # as the method; recognize those too.
        if method in orca_step.metadata["functionals"]:
            return True
        return orca_step.metadata["methods"].get(method, {}).get("type") == "DFT"

    def _cite_basis(self, P):
        """Cite the orbital basis set's own references (always -- the basis is
        used regardless of source). The Basis Set Exchange *tool* citation is
        added only when the basis was actually fetched from the BSE."""
        _, basis = self._resolve_method_basis(P)
        using_bse = self._using_bse(P)
        basis = self._strip_bse(basis)
        if not basis:
            return
        _, configuration = self.get_system_configuration(None)
        znums = sorted(set(configuration.atoms.atomic_numbers))

        import basis_set_exchange as bse

        # fmt="json" returns only the basis's own references -- NOT the BSE tool
        # citation (which fmt="bib" would prepend).
        data = json.loads(bse.get_references(basis, elements=znums, fmt="json"))
        seen = set()
        for block in data:
            for info in block.get("reference_info", []):
                for pair in info.get("reference_data", []):
                    key, entry = pair[0], pair[1]
                    if key in seen:
                        continue
                    seen.add(key)
                    self.references.cite(
                        raw=self._json_to_bibtex(key, entry),
                        alias=f"basis-{key}",
                        module="orca_step",
                        level=2,
                        note=f"Reference for the {basis} basis set.",
                    )

        # The BSE tool itself -- only when we actually used the BSE.
        if using_bse and "bse" in self._bibliography:
            self.references.cite(
                raw=self._bibliography["bse"],
                alias="bse",
                module="orca_step",
                level=2,
                note="Basis set obtained from the Basis Set Exchange.",
            )

    @staticmethod
    def _json_to_bibtex(key, e):
        """Build a bibtex entry from a BSE json reference dict."""
        lines = [f"@{e.get('_entry_type', 'article')}{{{key},"]
        if e.get("authors"):
            lines.append(f"    author = {{{' and '.join(e['authors'])}}},")
        for field in ("title", "journal", "volume", "pages", "year", "doi"):
            if e.get(field):
                lines.append(f"    {field} = {{{e[field]}}},")
        return "\n".join(lines).rstrip(",") + "\n}"

    def analyze(self, indent="", P=None, data=None, **kwargs):
        """Parse the properties, store them, print a summary, write CSV files,
        and optionally apply bond orders / Hirshfeld charges to the structure."""
        if P is None:
            P = self.parameters.current_values_to_dict(
                context=seamm.flowchart_variables._data
            )
        directory = Path(self.directory)
        props = self._parse_properties(directory)
        # Keep the run_orca energy if parsing missed it.
        if "energy" not in props and data and data.get("energy") is not None:
            props["energy"] = data["energy"]

        _, configuration = self.get_system_configuration(None)
        n_atoms = configuration.n_atoms
        try:
            max_print = int(self.parent.options.get("max_atoms_to_print", 25))
        except (TypeError, ValueError):
            max_print = 25

        # Store the scalar/array results (writes the storable properties to the
        # configuration and any results the user asked to save).
        try:
            self.store_results(configuration=configuration, data=props)
        except Exception as e:  # pragma: no cover
            logger.warning(f"Could not store results: {e}")

        self._print_scalar_summary(props)
        self._report_gradients(props)
        self._report_charges(props, configuration, directory, n_atoms, max_print)
        self._report_mayer(P, props, configuration, directory, n_atoms, max_print)
        self._report_hirshfeld(P, props, configuration, directory, n_atoms, max_print)

    # ------------------------------------------------------------------
    # Output parsing
    # ------------------------------------------------------------------
    def _parse_properties(self, directory):
        """Parse the available properties from orca.out."""
        path = Path(directory) / "orca.out"
        text = path.read_text() if path.exists() else ""
        p = {}

        def last_float(pattern):
            m = re.findall(pattern, text)
            return float(m[-1]) if m else None

        for key, pat in (
            ("energy", r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"),
            ("ccsd energy", r"E\(CCSD\)\s*\.*\s*(-?\d+\.\d+)"),
            ("ccsd(t) energy", r"E\(CCSD\(T\)\)\s*\.*\s*(-?\d+\.\d+)"),
            ("mp2 energy", r"MP2 TOTAL ENERGY:?\s*(-?\d+\.\d+)"),
            ("dipole magnitude", r"Magnitude \(Debye\)\s*:\s*(-?\d+\.\d+)"),
            (
                "isotropic polarizability",
                r"Isotropic polarizability\s*:?\s*(-?\d+\.\d+)",
            ),
        ):
            v = last_float(pat)
            if v is not None:
                p[key] = v

        m = re.findall(
            r"TOTAL SCF ENERGY.*?Total Energy\s*:\s*(-?\d+\.\d+)\s*Eh", text, re.DOTALL
        )
        if m:
            p["scf energy"] = float(m[-1])

        m = self._last(r"<S\*\*2>\s*[:=]\s*(-?\d+\.\d+)", text)
        if m:
            p["S^2"] = float(m.group(1))

        m = self._last(
            r"Total Dipole Moment\s*:\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",
            text,
        )
        if m:  # a.u. -> Debye
            au2d = 2.5417464519
            p["dipole moment"] = [round(float(m.group(i)) * au2d, 6) for i in (1, 2, 3)]

        m = self._last(
            r"Rotational constants in MHz\s*:\s*"
            r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",
            text,
        )
        if m:  # MHz -> GHz
            p["rotational constants"] = [
                round(float(m.group(i)) / 1000, 6) for i in (1, 2, 3)
            ]

        self._parse_orbitals(text, p)
        p["mulliken charges"] = self._parse_atom_charges(
            text, "MULLIKEN ATOMIC CHARGES"
        )
        p["löwdin charges"] = self._parse_atom_charges(text, "LÖWDIN ATOMIC CHARGES")
        p["hirshfeld charges"] = self._parse_hirshfeld(text)
        p["mayer valences"], p["mayer bonds"] = self._parse_mayer(text)
        p["gradients"] = self._parse_gradients(directory, text)
        return {k: v for k, v in p.items() if v not in (None, [], {})}

    def _parse_gradients(self, directory, text=""):
        """The Cartesian energy gradient (E_h/bohr) as ``[n_atoms][3]``.

        Prefer the structured ``orca.engrad`` file (written when ``! EnGrad`` is
        requested); fall back to the ``CARTESIAN GRADIENT`` block in orca.out.
        """
        path = Path(directory) / "orca.engrad"
        if path.exists():
            lines = [
                ln.strip()
                for ln in path.read_text().splitlines()
                if ln.strip() and not ln.strip().startswith("#")
            ]
            try:
                n = int(lines[0])  # number of atoms
                vals = [float(x) for x in lines[2 : 2 + 3 * n]]  # noqa: E203
                if len(vals) == 3 * n:
                    return [vals[3 * i : 3 * i + 3] for i in range(n)]  # noqa: E203
            except (ValueError, IndexError):
                pass

        # Fall back to the orca.out CARTESIAN GRADIENT block.
        m = self._last(r"CARTESIAN GRADIENT\s*\n-+\n(.*?)\n\s*\n", text, re.DOTALL)
        if not m:
            return None
        grad = []
        for row in m.group(1).splitlines():
            nums = re.findall(r"-?\d+\.\d+", row)
            if len(nums) >= 3:
                grad.append([float(x) for x in nums[-3:]])
        return grad or None

    @staticmethod
    def _last(pattern, text, flags=0):
        """The last match of `pattern` in `text` (for optimizations, the property
        blocks are printed per geometry -- we want the final one)."""
        matches = list(re.finditer(pattern, text, flags))
        return matches[-1] if matches else None

    def _parse_orbitals(self, text, p):
        """HOMO/LUMO (and nHOMO/nLUMO) energies + the gap from the orbital list."""
        m = self._last(r"ORBITAL ENERGIES\s*\n-+\n(.*?)\n\s*\n", text, re.DOTALL)
        if not m:
            return
        rows = re.findall(
            r"^\s*\d+\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", m.group(1), re.M
        )
        occ = [i for i, r in enumerate(rows) if float(r[0]) > 0.5]
        if not occ or occ[-1] + 1 >= len(rows):
            return
        h, lu = occ[-1], occ[-1] + 1
        p["HOMO energy"] = float(rows[h][1])
        p["LUMO energy"] = float(rows[lu][1])
        p["HOMO-LUMO gap"] = round(float(rows[lu][2]) - float(rows[h][2]), 4)  # eV
        p["HOMO index"], p["LUMO index"] = h, lu
        if h - 1 >= 0:
            p["nHOMO energy"] = float(rows[h - 1][1])
        if lu + 1 < len(rows):
            p["nLUMO energy"] = float(rows[lu + 1][1])

    def _parse_atom_charges(self, text, header):
        # Block runs from the header to the next blank line. (Mulliken ends with
        # a "Sum of atomic charges" line, Löwdin does not -- so terminate on the
        # blank line and pull out only the per-atom "<i> <El> : <q>" lines.) Use
        # the last block, so an optimization reports the final geometry.
        m = self._last(re.escape(header) + r"\s*\n-+\n(.*?)\n\s*\n", text, re.DOTALL)
        if not m:
            return []
        return [
            float(x)
            for x in re.findall(r"^\s*\d+\s+\w+\s*:\s*(-?\d+\.\d+)", m.group(1), re.M)
        ]

    def _parse_hirshfeld(self, text):
        m = self._last(
            r"HIRSHFELD ANALYSIS.*?CHARGE\s+SPIN\s*\n(.*?)\n\s*\n", text, re.DOTALL
        )
        if not m:
            return []
        return [
            float(r)
            for r in re.findall(r"^\s*\d+\s+\w+\s+(-?\d+\.\d+)", m.group(1), re.M)
        ]

    def _parse_mayer(self, text):
        # Use the LAST occurrence of each block (final geometry of an optimization).
        m = self._last(
            r"ATOM\s+NA\s+ZA\s+QA\s+VA\s+BVA\s+FA\s*\n(.*?)\n\s*\n", text, re.DOTALL
        )
        valences = []
        if m:
            for line in m.group(1).splitlines():
                f = line.split()
                if len(f) >= 7:
                    valences.append(float(f[5]))  # VA column
        bonds = []
        m = self._last(
            r"Mayer bond orders larger than.*?\n(.*?)\n\s*\n", text, re.DOTALL
        )
        if m:
            bonds = [
                (int(i), int(j), float(o))
                for i, j, o in re.findall(
                    r"B\(\s*(\d+)-\w+\s*,\s*(\d+)-\w+\s*\)\s*:\s*(-?\d+\.\d+)",
                    m.group(1),
                )
            ]
        return valences, bonds

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------
    def _print_scalar_summary(self, p):
        rows = []

        def add(prop, value, units="", fmt="{:.6f}"):
            if value is not None:
                rows.append([prop, fmt.format(value), units])

        add("Total energy", p.get("energy"), "E_h", "{:.8f}")
        add("SCF energy", p.get("scf energy"), "E_h", "{:.8f}")
        add("MP2 energy", p.get("mp2 energy"), "E_h", "{:.8f}")
        add("CCSD energy", p.get("ccsd energy"), "E_h", "{:.8f}")
        add("CCSD(T) energy", p.get("ccsd(t) energy"), "E_h", "{:.8f}")
        if "HOMO energy" in p:
            add(
                f"HOMO energy (MO {p['HOMO index']})", p["HOMO energy"], "E_h", "{:.5f}"
            )
            add(
                f"LUMO energy (MO {p['LUMO index']})", p["LUMO energy"], "E_h", "{:.5f}"
            )
        add("HOMO-1 energy", p.get("nHOMO energy"), "E_h", "{:.5f}")
        add("LUMO+1 energy", p.get("nLUMO energy"), "E_h", "{:.5f}")
        add("HOMO-LUMO gap", p.get("HOMO-LUMO gap"), "eV", "{:.3f}")
        if "dipole moment" in p:  # 3-vector -> one row per component
            for comp, v in zip("XYZ", p["dipole moment"]):
                add(f"Dipole moment {comp}", v, "debye", "{:.4f}")
        add("Dipole magnitude", p.get("dipole magnitude"), "debye", "{:.4f}")
        if "rotational constants" in p:
            for comp, v in zip("ABC", p["rotational constants"]):
                add(f"Rotational constant {comp}", v, "GHz", "{:.4f}")
        add(
            "Isotropic polarizability",
            p.get("isotropic polarizability"),
            "a.u.",
            "{:.4f}",
        )
        add("<S^2>", p.get("S^2"), "", "{:.4f}")

        if rows:
            tmp = tabulate(
                rows,
                headers=["Property", "Value", "Units"],
                tablefmt="rounded_outline",
                colalign=("left", "right", "left"),
                disable_numparse=True,
            )
            printer.normal("")
            printer.normal(textwrap.indent(tmp, self.indent + 7 * " "))

    def _report_gradients(self, p):
        """Note that the gradient was computed and the largest component, which
        is the useful at-a-glance number; the full array goes to Results.json."""
        grad = p.get("gradients")
        if not grad:
            return
        gmax = max(abs(g) for row in grad for g in row)
        printer.normal(
            __(
                f"Computed the energy gradient on {len(grad)} atoms "
                f"(largest component {gmax:.6f} E_h/bohr).",
                indent=self.indent + 4 * " ",
            )
        )

    def _report_charges(self, p, configuration, directory, n_atoms, max_print):
        symbols = list(configuration.atoms.symbols)
        for name, key in (
            ("Mulliken", "mulliken charges"),
            ("Löwdin", "löwdin charges"),
        ):
            charges = p.get(key)
            if not charges or len(charges) != n_atoms:
                continue
            self._write_charge_csv(directory, key.replace(" ", "_"), symbols, charges)
        if n_atoms <= max_print and p.get("mulliken charges"):
            tbl = {"Atom": list(range(1, n_atoms + 1)), "Element": symbols}
            if len(p.get("mulliken charges", [])) == n_atoms:
                tbl["Mulliken"] = [f"{q:.4f}" for q in p["mulliken charges"]]
            if len(p.get("löwdin charges", [])) == n_atoms:
                tbl["Löwdin"] = [f"{q:.4f}" for q in p["löwdin charges"]]
            self._print_table("Atomic charges (e)", tbl)

    @staticmethod
    def _write_charge_csv(directory, name, symbols, charges):
        with open(Path(directory) / f"{name}.csv", "w", newline="") as fd:
            w = csv.writer(fd)
            w.writerow(["Atom", "Element", "Charge"])
            for i, (s, q) in enumerate(zip(symbols, charges), start=1):
                w.writerow([i, s, f"{q:.6f}"])

    def _report_mayer(self, P, p, configuration, directory, n_atoms, max_print):
        control = P["bond orders"]
        if control == "no":
            return
        bonds = p.get("mayer bonds", [])
        # Always write the bond orders (above threshold) to a CSV.
        with open(Path(directory) / "mayer_bond_orders.csv", "w", newline="") as fd:
            w = csv.writer(fd)
            w.writerow(["i", "j", "bond order", "multiplicity"])
            for i, j, order in bonds:
                if order > 0.5:
                    w.writerow(
                        [i + 1, j + 1, f"{order:.4f}", self._multiplicity(order)[1]]
                    )

        names = configuration.atoms.names
        signif = [(i, j, o) for i, j, o in bonds if o > 0.5]
        if n_atoms <= max_print and signif:
            tbl = {
                "i": [names[i] for i, _, _ in signif],
                "j": [names[j] for _, j, _ in signif],
                "bond order": [f"{o:.3f}" for _, _, o in signif],
                "multiplicity": [self._multiplicity(o)[1] for _, _, o in signif],
            }
            self._print_table("Mayer bond orders", tbl)

        if control == "yes, and apply to structure" and signif:
            ids = configuration.atoms.ids
            iatoms = [ids[i] for i, _, _ in signif]
            jatoms = [ids[j] for _, j, _ in signif]
            orders = [self._multiplicity(o)[0] for _, _, o in signif]
            configuration.new_bondset()
            configuration.bonds.append(i=iatoms, j=jatoms, bondorder=orders)
            printer.normal(
                __(
                    "Replaced the bonds in the structure with those from the Mayer "
                    "bond orders.",
                    indent=self.indent + 4 * " ",
                )
            )

    @staticmethod
    def _multiplicity(order):
        """(bondorder code, label) from a bond order, MOPAC-style thresholds."""
        if 1.3 < order < 1.7:
            return 5, "aromatic"
        n = round(order)
        return n, {1: "single", 2: "double", 3: "triple"}.get(n, str(n))

    def _report_hirshfeld(self, P, p, configuration, directory, n_atoms, max_print):
        control = P["Hirshfeld charges"]
        if control == "no":
            return
        charges = p.get("hirshfeld charges", [])
        symbols = list(configuration.atoms.symbols)
        if not charges:
            printer.normal(
                __(
                    "Hirshfeld charges were requested but not found in the output.",
                    indent=self.indent + 4 * " ",
                )
            )
            return
        self._write_charge_csv(directory, "hirshfeld_charges", symbols, charges)
        if n_atoms <= max_print and len(charges) == n_atoms:
            self._print_table(
                "Hirshfeld charges (e)",
                {
                    "Atom": list(range(1, n_atoms + 1)),
                    "Element": symbols,
                    "Charge": [f"{q:.4f}" for q in charges],
                },
            )
        if control == "yes, and apply to structure" and len(charges) == n_atoms:
            atoms = configuration.atoms
            if "charge" not in atoms:
                atoms.add_attribute(
                    "charge", coltype="float", configuration_dependent=True
                )
            atoms["charge"][0:] = list(charges)
            printer.normal(
                __(
                    "Stored the Hirshfeld charges as the atomic charges on the "
                    "structure.",
                    indent=self.indent + 4 * " ",
                )
            )

    def _print_table(self, title, table):
        tmp = tabulate(
            table, headers="keys", tablefmt="rounded_outline", disable_numparse=True
        )
        printer.normal("")
        printer.normal(textwrap.indent(title + "\n" + tmp, self.indent + 7 * " "))

    def cleanup(self):
        """Nothing to clean up for a single point."""
        pass
