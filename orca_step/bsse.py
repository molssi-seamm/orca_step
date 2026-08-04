# -*- coding: utf-8 -*-

"""An ORCA counterpoise (BSSE) sub-step.

Computes the counterpoise-corrected (Boys--Bernardi) energy and gradient of an
N-fragment complex, with independent per-fragment charge. SEAMM generates the
2N + 1 job specs (via ``seamm_bsse``), runs each as an ordinary ORCA job (real
atoms + ghost-flagged atoms of the other fragments, or ghost-free for a
fragment alone in its own basis), and combines the results into the corrected
total energy and gradient. See the campaign design note
(``docs/developer_guide/campaigns/2026-08-03/bsse_architecture.rst``) for the
physics, and ``seamm_bsse`` itself for the job-spec/combine algebra.

The original ORCA *Compound* script path (``bssegradient.cmp``/
``bssenergy.cmp``, D. G. Liakos & F. Neese -- one ORCA process running all 5
sub-calculations of a 2-fragment, neutral-singlet complex internally) is kept
in this module and in ``data/`` as the N = 2 regression oracle for the new
path (see the architecture doc's milestone M2), not as a second production
code path. ``_compound_input``/``run_orca_compound``/
``_parse_compound_energies`` below are unused by :meth:`run` now.
"""

import importlib.resources
import logging
from pathlib import Path
import re
import shutil
import textwrap

from tabulate import tabulate

import orca_step
from .energy import Energy
import seamm
import seamm_bsse
from seamm_util import Q_
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
    def _fragments(self, P, configuration):
        """Return the N ``seamm_bsse.Fragment`` for this complex, per the
        'fragments' parameter, with per-fragment charge threaded in from
        'fragment charges'. Validates against the complex's own
        charge/multiplicity (``seamm_bsse.validate_fragments``) before
        returning -- catches a fragment-charge typo, or an open-shell/
        multiplicity>1 fragment (not yet supported), with a clear message.
        """
        n_atoms = configuration.n_atoms
        mode = P["fragments"]
        if mode == "specified":
            groups = self._parse_fragment_groups(P["fragment atoms"], n_atoms)
        else:
            # "auto (molecules)"
            molecules = configuration.find_molecules(as_indices=True)
            if len(molecules) < 2:
                raise RuntimeError(
                    f"BSSE 'auto' fragments require at least two separate "
                    f"molecules, but the structure has {len(molecules)}. Use "
                    "'specified' to define the fragments by atom, or check "
                    "the bonding."
                )
            groups = [sorted(molecule) for molecule in molecules]

        charges = self._parse_charges(P["fragment charges"])
        if not charges:
            charges = [0] * len(groups)
        elif len(charges) != len(groups):
            raise RuntimeError(
                f"BSSE: {len(charges)} 'fragment charges' given but "
                f"{len(groups)} fragments found/specified; give one charge "
                "per fragment (in the same order), or leave 'fragment "
                "charges' empty for all-neutral."
            )

        fragments = [
            seamm_bsse.Fragment(label=str(i + 1), atom_indices=group, charge=charge)
            for i, (group, charge) in enumerate(zip(groups, charges))
        ]
        try:
            seamm_bsse.validate_fragments(
                fragments,
                cluster_charge=configuration.charge,
                cluster_multiplicity=configuration.spin_multiplicity,
            )
        except ValueError as e:
            raise RuntimeError(f"BSSE: {e}") from e
        return fragments

    @staticmethod
    def _parse_fragment_groups(text, n_atoms):
        """Parse 'specified'-mode fragment atoms: semicolon-separated
        1-based index/range groups, one per fragment, e.g. ``'1-3; 4-6; 7'``
        for three fragments -- to 0-based index lists."""
        groups = [group.strip() for group in str(text).split(";") if group.strip()]
        if len(groups) < 2:
            raise RuntimeError(
                "BSSE: 'Fragment atoms' (specified mode) needs at least two "
                "semicolon-separated groups of atoms, e.g. '1-3; 4-6' for two "
                "fragments."
            )
        return [BSSE._parse_indices(group, n_atoms) for group in groups]

    @staticmethod
    def _parse_charges(text):
        """Parse 'fragment charges': a comma/space separated list of integers,
        in fragment order. '' (the default, all-neutral) -> []."""
        text = str(text).strip()
        if not text:
            return []
        return [int(token) for token in text.replace(",", " ").split()]

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
    # Compound-script input generation -- N = 2 regression oracle only.
    # Not called by run() (see the module docstring); kept for the M2
    # validation gate and any future re-check against it.
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
    def _wants_gradients(self, P):
        """Whether this run wants the gradient. Overrides Energy's version
        (which checks the driver-supplied 'results' entry) with BSSE's own
        'compute gradient' toggle, so the inherited ``keyword_line`` adds
        EnGrad/NumGrad correctly for every counterpoise sub-job."""
        return P.get("compute gradient", "yes") == "yes"

    def run(self, keywords=None):
        """Run the N-fragment counterpoise correction: generate the 2N + 1
        job specs (``seamm_bsse``), run each as an ordinary ORCA job in its
        own sub-directory (real atoms + the other fragments' atoms ghosted,
        or a fragment alone in its own basis), and combine the results into
        the corrected total energy and gradient.

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
        fragments = self._fragments(P, configuration)
        for fragment in fragments:
            printer.important(
                __(
                    f"Fragment {fragment.label} has {len(fragment.atom_indices)} "
                    f"atom(s), charge {fragment.charge:+d}.",
                    indent=self.indent + 4 * " ",
                )
            )

        specs = seamm_bsse.generate_job_specs(fragments)
        printer.important(
            __(
                f"Running {len(specs)} ORCA job(s) for the "
                f"{len(fragments)}-fragment counterpoise correction.",
                indent=self.indent + 4 * " ",
            )
        )

        want_gradient = P.get("compute gradient", "yes") == "yes"
        optimize_monomers = P.get("optimize monomers", "no") == "yes"
        make_wfx = P.get("save wavefunction", "no") == "yes"
        # 'keepdensity'/wavefunction conversion is only meaningful for the
        # full cluster (what a following Atomic Charges step wants), not
        # every fragment sub-job -- build the shared keyword line without it.
        base_keyword_line = self.keyword_line({**P, "save wavefunction": "no"})

        results = {}
        for spec in specs:
            job_directory = Path(self.directory) / spec.label
            keyword_line = base_keyword_line
            # 'optimize monomers' (the old Compound script's DoOptimization)
            # relaxes only the free (fragment-alone) sub-jobs.
            if optimize_monomers and spec.kind == seamm_bsse.FRAGMENT_ALONE:
                keyword_line = f"Opt {keyword_line}"
            job_make_wfx = make_wfx and spec.kind == seamm_bsse.CLUSTER
            if job_make_wfx:
                keyword_line = f"{keyword_line} keepdensity"

            outcome = self.run_orca_job(
                keyword_line,
                configuration,
                spec.charge,
                spec.multiplicity,
                atom_indices=spec.atom_indices,
                ghost_atoms=spec.ghost_indices,
                directory=job_directory,
                make_wfx=job_make_wfx,
            )
            gradient = self._parse_gradients(job_directory) if want_gradient else None
            if want_gradient and gradient is None:
                raise RuntimeError(
                    f"The ORCA BSSE job '{spec.label}' produced no gradient; "
                    f"see {job_directory}/orca.out and orca.err."
                )
            results[spec.label] = seamm_bsse.JobResult(
                energy=outcome["energy"], gradient=gradient
            )

        cp = seamm_bsse.combine(specs, results, n_atoms=configuration.n_atoms)

        # For a following Atomic Charges (DDEC6) step: the full cluster's
        # wavefunction, copied up to this node's own directory (where a
        # downstream step expects it), mirroring the Energy step.
        if make_wfx:
            cluster_wfx = Path(self.directory) / seamm_bsse.CLUSTER / "orca.wfx"
            if cluster_wfx.exists():
                shutil.copy(cluster_wfx, Path(self.directory) / "orca.wfx")

        # Tag stored properties with the level of theory so BSSE-corrected data
        # is distinguishable in the database.
        self.model = self._model_string(P)

        data = {
            "success": True,
            "energy": cp.energy,
            "uncorrected energy": results[seamm_bsse.CLUSTER].energy,
            "bsse correction": cp.bsse_correction,
            "interaction energy": Q_(cp.interaction_energy, "E_h").m_as("kJ/mol"),
            "uncorrected interaction energy": Q_(
                cp.uncorrected_interaction_energy, "E_h"
            ).m_as("kJ/mol"),
        }
        if cp.gradient is not None:
            data["gradients"] = cp.gradient
        self._data = data
        self._cite_references(P)
        self._cite_bsse()
        self.analyze(P=P, data=data)

        return self.next()

    def _check_supported(self, P, configuration):
        """Refuse the cases this sub-step cannot do correctly, with a clear
        message, rather than returning wrong numbers. Per-fragment
        charge/multiplicity consistency (vs. the complex's own) is checked
        separately, in ``_fragments`` -- it needs the fragment definition
        first."""
        if self._using_bse(P):
            raise RuntimeError(
                "The ORCA BSSE sub-step does not yet support Basis Set Exchange "
                "bases; choose an ORCA-internal basis set."
            )
        # F12 methods need an F12 basis (matching CABS) -- fail early if not.
        self._check_f12(P)
        # When the gradient is also requested, the method must have an
        # analytic gradient (every job requests EnGrad); energy-only lifts
        # that (e.g. for gold-standard (DLPNO-)CCSD(T) interaction energies).
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
        BSSE correction, the interaction (binding) energy relative to the
        separated monomers (corrected and uncorrected, in kJ/mol), the
        corrected gradient, and the formation-referenced DfE0 computed from
        the corrected energy.

        Unlike Energy, this does NOT parse orca.out for the usual properties:
        with a Compound job orca.out holds the sub-calculations, whose
        per-fragment properties would be misleading. Only the corrected
        energy/gradient (from the EnGrad file), the raw complex energy, and the
        correction are meaningful.
        """
        if data is None:
            data = self._data or {}
        if P is None:
            P = self.parameters.current_values_to_dict(
                context=seamm.flowchart_variables._data
            )
        props = {
            key: data[key]
            for key in (
                "energy",
                "gradients",
                "uncorrected energy",
                "bsse correction",
                "interaction energy",
                "uncorrected interaction energy",
            )
            if data.get(key) is not None
        }

        _, configuration = self.get_system_configuration(None)

        # DfE0 (and the atomization energy it needs), referenced to the
        # BSSE-corrected complex energy in `props["energy"]` -- NOT the
        # uncorrected energy, and not any of the five intermediate
        # sub-calculation energies. Mutates `props` in place, so DfE0/E
        # atomization flow into store_results/the printed summary below
        # exactly as the plain Energy sub-step's does.
        if "energy" in props:
            formation_text = self.calculate_energy_of_formation(P, props, configuration)
            if formation_text:
                (Path(self.directory) / "Thermochemistry.txt").write_text(
                    formation_text
                )

        try:
            self.store_results(configuration=configuration, data=props)
        except Exception as e:  # pragma: no cover
            logger.warning(f"Could not store results: {e}")

        # Energy breakdown: formation/atomization energies first, then the
        # interaction (binding) energy relative to the separated monomers
        # (the headline number for a dimer/cluster), then the absolute
        # uncorrected -> correction -> corrected complex energies. The BSSE
        # correction is small, so also show it in kcal/mol.
        rows = []
        if "DfE0" in props:
            rows.append(["Energy of formation (0 K)", f"{props['DfE0']:.2f}", "kJ/mol"])
        if "E atomization" in props:
            rows.append(
                ["Atomization energy", f"{props['E atomization']:.2f}", "kJ/mol"]
            )
        if "interaction energy" in props:
            rows.append(
                [
                    "Interaction energy (CP-corrected)",
                    f"{props['interaction energy']:.4f}",
                    "kJ/mol",
                ]
            )
        if "uncorrected interaction energy" in props:
            rows.append(
                [
                    "Interaction energy (uncorrected)",
                    f"{props['uncorrected interaction energy']:.4f}",
                    "kJ/mol",
                ]
            )
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
