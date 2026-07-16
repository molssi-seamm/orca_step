# -*- coding: utf-8 -*-

"""An ORCA counterpoise (BSSE) sub-step.

Computes the counterpoise-corrected (Boys--Bernardi) energy and gradient of a
two-fragment complex in a single ORCA run, by driving the ORCA *Compound* script
``bssegradient.cmp`` (D. G. Liakos & F. Neese). SEAMM prepares the ghost-flagged
geometry, injects the level of theory and options, runs ORCA once, and reads the
corrected energy and gradient from the EnGrad file the script writes.

See the campaign design note
(``docs/developer_guide/campaigns/2026-07-09/bsse_scope.rst``) for the physics
and the Phase-1 scope boundaries.
"""

import importlib.resources
import logging
from pathlib import Path
import re
import textwrap

from tabulate import tabulate

import orca_step
from .energy import Energy
import seamm
from seamm_util.printing import FormattedText as __
import seamm_util.printing as printing

#: Hartree -> kcal/mol, for reporting the (small) correction in familiar units.
_HARTREE_TO_KCAL = 627.509474

logger = logging.getLogger(__name__)
printer = printing.getPrinter("ORCA")

#: Boys & Bernardi counterpoise method -- the correction this step applies.
_BSSE_CITATION = """\
@article{Boys1970,
    author = {Boys, S. F. and Bernardi, F.},
    title = {The calculation of small molecular interactions by the differences
             of separate total energies. Some procedures with reduced errors},
    journal = {Molecular Physics},
    volume = {19},
    number = {4},
    pages = {553--566},
    year = {1970},
    doi = {10.1080/00268977000101561}
}
"""


class BSSE(Energy):
    """A counterpoise (BSSE) correction with ORCA.

    See Also
    --------
    TkBSSE, BSSEParameters, Energy
    """

    def __init__(self, flowchart=None, title="BSSE", extension=None, logger=logger):
        logger.debug(f"Creating ORCA BSSE {self}")
        super().__init__(
            flowchart=flowchart, title=title, extension=extension, logger=logger
        )
        self._calculation = "bsse"
        self.parameters = orca_step.BSSEParameters()

    def _extrapolating(self, P):
        """Never extrapolate: a CBS-extrapolated energy has no gradient, so it
        cannot drive the counterpoise gradient. Hidden in the GUI; also ignored
        here for a hand-edited flowchart."""
        return False

    def description_text(self, P=None):
        if not P:
            P = self.parameters.values_to_dict()
        what = (
            "energy and gradient"
            if P.get("compute gradient", "yes") == "yes"
            else "energy"
        )
        text = (
            f"Counterpoise (BSSE) corrected {what} with ORCA at "
            f"{self._level_of_theory_text(P)}."
        )
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    # ------------------------------------------------------------------
    # Fragments
    # ------------------------------------------------------------------
    def _fragment_atoms(self, P, configuration):
        """Return ``(fragmentA, fragmentB)`` as 0-based atom-index lists.

        Fragment B is written as ghost atoms in the input; the correction is
        symmetric, so which is A vs B is immaterial.
        """
        n_atoms = configuration.n_atoms
        mode = P["fragments"]
        if mode == "specified":
            fragA = self._parse_indices(P["fragment A atoms"], n_atoms)
            if not fragA:
                raise RuntimeError(
                    "BSSE: 'Fragment A atoms' is empty; list the atoms of "
                    "fragment A (1-based), e.g. '1-3, 5'."
                )
            fragB = [i for i in range(n_atoms) if i not in set(fragA)]
            if not fragB:
                raise RuntimeError(
                    "BSSE: fragment A is the whole system; nothing is left for "
                    "fragment B."
                )
            return sorted(fragA), fragB

        # "auto (2 molecules)"
        molecules = configuration.find_molecules(as_indices=True)
        if len(molecules) != 2:
            raise RuntimeError(
                f"BSSE 'auto' fragments require exactly two separate molecules, "
                f"but the structure has {len(molecules)}. Use 'specified' to "
                "define the fragments by atom, or check the bonding."
            )
        return sorted(molecules[0]), sorted(molecules[1])

    @staticmethod
    def _parse_indices(text, n_atoms):
        """Parse a 1-based index/range list (``'1-3, 5 7'``) to 0-based indices."""
        indices = []
        for token in str(text).replace(",", " ").split():
            if "-" in token[1:]:  # a range like 1-3 (not a leading minus)
                lo, hi = token.split("-", 1)
                indices.extend(range(int(lo), int(hi) + 1))
            elif token:
                indices.append(int(token))
        out = []
        for i in indices:
            j = i - 1
            if j < 0 or j >= n_atoms:
                raise RuntimeError(
                    f"BSSE: atom index {i} is out of range (1..{n_atoms})."
                )
            if j not in out:
                out.append(j)
        return out

    # ------------------------------------------------------------------
    # Input generation
    # ------------------------------------------------------------------
    def _ghost_xyz(self, configuration, ghost_atoms):
        """An .xyz file of the whole complex with `ghost_atoms` (0-based) flagged
        as ORCA ghost centres (element symbol + ':')."""
        ghosts = set(ghost_atoms)
        symbols = configuration.atoms.symbols
        xyzs = configuration.atoms.get_coordinates(fractionals=False, in_cell=True)
        lines = [str(len(symbols)), ""]
        for i, (symbol, (x, y, z)) in enumerate(zip(symbols, xyzs)):
            label = f"{symbol}:" if i in ghosts else symbol
            lines.append(f"{label:4s} {x:15.8f} {y:15.8f} {z:15.8f}")
        return "\n".join(lines) + "\n"

    def _compound_input(self, P, xyz_filename, script_name):
        """Return ``(compound_block, method, basis)`` -- the ``%Compound`` block
        that runs ``script_name`` with the resolved level of theory injected via
        a ``with`` clause."""
        method, basis = self._resolve_method_basis(P)

        # "rest of input": the ! keywords the Compound sub-steps also use. SThresh
        # (a %scf block) and the property blocks are NOT plumbed through Compound
        # in Phase 1.
        rest = []
        aux = P["auxiliary basis"]
        if aux and aux.lower() != "none":
            rest.append(aux)
        grid = P.get("grid", "default")
        if grid == "default":
            grid = self._auto_grid(P)
        if grid and grid != "default":
            rest.append(grid)
        scf = P.get("scf convergence", "default")
        if scf and scf != "default":
            rest.append(scf)
        extra = P["extra keywords"].strip()
        if extra:
            rest.append(extra)
        # Explicitly-correlated (F12) methods need the matching CABS basis.
        cabs = self._cabs_keyword(P)
        if cabs:
            rest.append(cabs)
        rest_of_input = " ".join(rest)

        do_opt = "true" if P.get("optimize monomers", "no") == "yes" else "false"
        make_wfx = "true" if P.get("save wavefunction", "no") == "yes" else "false"

        block = "\n".join(
            [
                f'%Compound "{script_name}"',
                "  with",
                f'    molecule       = "{xyz_filename}";',
                f'    method         = "{method}";',
                f'    basis          = "{basis}";',
                f'    restOfInput    = "{rest_of_input}";',
                # Neutral closed-shell only in Phase 1 (guarded in run()); the
                # same charge/mult is applied to the monomer sub-calculations.
                "    charge         = 0;",
                "    mult           = 1;",
                f"    DoOptimization = {do_opt};",
                f"    ProduceWavefunction = {make_wfx};",
                "  end",
            ]
        )
        return block, method, basis

    def _compound_script(self, script_name):
        """The shipped Compound script `script_name` (as text)."""
        path = importlib.resources.files("orca_step") / "data" / script_name
        return path.read_text()

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------
    def run(self, keywords=None):
        """Run the counterpoise correction as an ORCA Compound job.

        Like the other ORCA sub-steps, this is driven by the main ORCA node
        (which set up printing and cited the plug-in), so it does not call
        ``super().run()``.
        """
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        printer.important(__(self.description_text(P), indent=self.indent))

        _, configuration = self.get_system_configuration(None)
        self._check_supported(P, configuration)

        fragmentA, fragmentB = self._fragment_atoms(P, configuration)
        printer.important(
            __(
                f"Fragment A has {len(fragmentA)} atoms and fragment B has "
                f"{len(fragmentB)} (written as ghost atoms).",
                indent=self.indent + 4 * " ",
            )
        )

        # Energy-only mode uses a gradient-free Compound script, so a method with
        # no analytic gradient (e.g. CCSD(T)) still runs; otherwise the gradient
        # script is used and the corrected gradient comes from its EnGrad file.
        want_gradient = P.get("compute gradient", "yes") == "yes"
        script_name = "bssegradient.cmp" if want_gradient else "bssenergy.cmp"

        xyz_filename = "bsse.xyz"
        xyz_text = self._ghost_xyz(configuration, ghost_atoms=fragmentB)
        compound_block, method, basis = self._compound_input(
            P, xyz_filename, script_name
        )
        extra_files = {
            script_name: self._compound_script(script_name),
            xyz_filename: xyz_text,
        }

        # The counterpoise-corrected GRADIENT comes from the script's EnGrad file
        # (its ghost-atom bookkeeping is done on the full-method Nuclear_Gradient,
        # so it is correct for every method). Its ENERGY is not used -- see below.
        # Optionally write a .wfx from the dimer (the last COMPOUND JOB, step 5)
        # for a following Atomic Charges (DDEC6) step, mirroring the Energy step.
        make_wfx = P.get("save wavefunction", "no") == "yes"
        gradient = None
        if want_gradient:
            _, gradient = self.run_orca_compound(
                compound_block, extra_files=extra_files, make_wfx=make_wfx
            )
            if gradient is None:
                raise RuntimeError(
                    "The ORCA BSSE Compound job produced no result.engrad "
                    "gradient; see orca.out and orca.err."
                )
        else:
            self.run_orca_compound(
                compound_block, extra_files=extra_files, engrad=None, make_wfx=make_wfx
            )

        # Compute the corrected ENERGY here from the five sub-calculations' total
        # energies (each step's FINAL SINGLE POINT ENERGY), NOT from the script's
        # result.engrad. The script reads ORCA's SCF_Energy, which for a double
        # hybrid omits the MP2 correlation (and would be wrong by that amount);
        # FINAL SINGLE POINT ENERGY is the full-method total that matches the
        # gradient. Dispersion has no BSSE (ghosts have no nuclei) and cancels in
        # the correction terms, so this is consistent for -D methods too.
        energies = self._parse_compound_energies(self.directory)
        if energies is None:
            raise RuntimeError(
                "Could not read the five BSSE sub-calculation energies from "
                "orca.out; see orca.out and orca.err."
            )
        e_fragA, e_monA, e_fragB, e_monB, e_total = energies
        corrected = e_total - (e_fragA - e_monA) - (e_fragB - e_monB)

        # Tag stored properties with the level of theory so BSSE-corrected data
        # is distinguishable in the database.
        self.model = f"{method}/{basis}"

        data = {
            "success": True,
            "energy": corrected,
            "uncorrected energy": e_total,
            "bsse correction": corrected - e_total,
        }
        if gradient is not None:
            data["gradients"] = gradient
        self._data = data
        self._cite_references(P)
        self._cite_bsse()
        self.analyze(P=P, data=data)

        return self.next()

    def _check_supported(self, P, configuration):
        """Refuse the cases the Phase-1 Compound path cannot do correctly, with a
        clear message, rather than returning wrong numbers."""
        if self._using_bse(P):
            raise RuntimeError(
                "The ORCA BSSE sub-step does not yet support Basis Set Exchange "
                "bases; choose an ORCA-internal basis set."
            )
        # F12 methods need an F12 basis (matching CABS) -- fail early if not.
        self._check_f12(P)
        # The energy is taken from each step's FINAL SINGLE POINT ENERGY (the
        # full-method total), so double hybrids and MP2 are fine. When the
        # gradient is also requested, the method must have an analytic gradient
        # (the gradient script requests EnGrad); energy-only lifts that.
        method, _ = self._resolve_method_basis(P)
        if (
            P.get("compute gradient", "yes") == "yes"
            and self._gradient_availability(P) != "analytic"
        ):
            raise RuntimeError(
                f"The ORCA BSSE sub-step needs an analytic gradient, but {method} "
                "has only a numerical one (e.g. (DLPNO-)CCSD(T)). Set 'Compute "
                "the gradient' to 'no' for an energy-only correction, or choose a "
                "method/functional with an analytic gradient."
            )
        # The Compound script applies the complex's charge/multiplicity to the
        # monomer sub-calculations too, so it is only valid when each neutral
        # closed-shell fragment shares them -- i.e. a neutral singlet complex.
        if configuration.charge != 0 or configuration.spin_multiplicity != 1:
            raise RuntimeError(
                "The ORCA BSSE sub-step (Phase 1) supports only a neutral, "
                "closed-shell (charge 0, multiplicity 1) complex, because the "
                "same charge/multiplicity is applied to each monomer. Charged or "
                "open-shell fragments need per-fragment charge/multiplicity "
                "(the general BSSE step)."
            )

    def _cite_bsse(self):
        """Cite the counterpoise method (Boys & Bernardi). Best-effort."""
        try:
            self.references.cite(
                raw=_BSSE_CITATION,
                alias="boys-bernardi-1970",
                module="orca_step",
                level=1,
                note="The counterpoise correction for BSSE.",
            )
        except Exception as e:  # pragma: no cover
            logger.warning(f"Could not cite the counterpoise method: {e}")

    def _parse_compound_energies(self, directory):
        """The five sub-calculation total energies, in the script's order:
        ``[E_fragA(AB), E_monA(A), E_fragB(AB), E_monB(B), E_total(AB)]`` (E_h).

        Each value is the *last* ``FINAL SINGLE POINT ENERGY`` inside its
        ``COMPOUND JOB N`` block -- the true full-method total (SCF + any MP2
        correlation + dispersion), and the last one so it is the converged value
        when the free monomers are optimized. Returns ``None`` if the five blocks
        or their energies cannot be found.
        """
        path = Path(directory) / "orca.out"
        if not path.exists():
            return None
        text = path.read_text()
        jobs = list(re.finditer(r"COMPOUND JOB\s+(\d+)", text))
        if len(jobs) < 5:
            return None
        energies = []
        for i in range(5):
            start = jobs[i].end()
            end = jobs[i + 1].start() if i + 1 < len(jobs) else len(text)
            fspe = re.findall(
                r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text[start:end]
            )
            if not fspe:
                return None
            energies.append(float(fspe[-1]))
        return energies

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------
    def analyze(self, indent="", P=None, data=None, **kwargs):
        """Store and report the corrected energy, the uncorrected energy, the
        BSSE correction, and the corrected gradient.

        Unlike Energy, this does NOT parse orca.out for the usual properties:
        with a Compound job orca.out holds the sub-calculations, whose
        per-fragment properties would be misleading. Only the corrected
        energy/gradient (from the EnGrad file), the raw complex energy, and the
        correction are meaningful.
        """
        if data is None:
            data = self._data or {}
        props = {
            key: data[key]
            for key in ("energy", "gradients", "uncorrected energy", "bsse correction")
            if data.get(key) is not None
        }

        _, configuration = self.get_system_configuration(None)
        try:
            self.store_results(configuration=configuration, data=props)
        except Exception as e:  # pragma: no cover
            logger.warning(f"Could not store results: {e}")

        # Energy breakdown: uncorrected -> correction -> corrected. The
        # correction is small, so also show it in kcal/mol.
        rows = []
        if "uncorrected energy" in props:
            rows.append(
                ["Uncorrected energy", f"{props['uncorrected energy']:.8f}", "E_h"]
            )
        if "energy" in props:
            rows.append(["BSSE-corrected energy", f"{props['energy']:.8f}", "E_h"])
        if "bsse correction" in props:
            corr = props["bsse correction"]
            rows.append(["BSSE correction", f"{corr:.8f}", "E_h"])
            rows.append(
                ["BSSE correction", f"{corr * _HARTREE_TO_KCAL:.4f}", "kcal/mol"]
            )
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

        self._report_gradients(props)
