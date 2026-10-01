# -*- coding: utf-8 -*-
"""Control parameters for an ORCA single-point energy."""

import logging

import orca_step
import seamm

logger = logging.getLogger(__name__)


class EnergyParameters(seamm.Parameters):
    """The control parameters for an ORCA energy calculation.

    The default path takes the method from a preceding Model Chemistry step; the
    user can turn that off and choose the method and basis explicitly (similar to
    the Gaussian step). The explicit choices are driven by the step's metadata.
    """

    parameters = {
        "use model chemistry": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Use the global model chemistry:",
            "help_text": (
                "Use the model chemistry defined by a preceding Model Chemistry "
                "step (the '_model_chemistry' variable). If ORCA cannot provide "
                "it, an error is raised. Turn this off to set the method and basis "
                "explicitly below."
            ),
        },
        "method": {
            "applies_when": {"use model chemistry": "no"},
            "default": "DLPNO-CCSD(T)",
            "kind": "enum",
            "default_units": "",
            "enumeration": tuple(orca_step.metadata["methods"].keys()),
            "format_string": "",
            "description": "Method:",
            "help_text": (
                "The ORCA method. For 'DFT' the exchange-correlation functional "
                "is chosen with the two controls below; every other choice is an "
                "ORCA '!' keyword on its own."
            ),
        },
        "functional type": {
            "applies_when": {"method": "DFT"},
            "default": "global hybrid",
            "kind": "enum",
            "default_units": "",
            "enumeration": tuple(orca_step.metadata["functional categories"]),
            "format_string": "",
            "description": "Functional type:",
            "help_text": (
                "ORCA's classification of the density functional (local, GGA, "
                "meta-GGA, (range-separated) hybrid, (range-separated) "
                "double-hybrid). Picking a type filters the functional list. Only "
                "used when the method is 'DFT'."
            ),
        },
        "functional": {
            "applies_when": {"method": "DFT"},
            "default": "B3LYP",
            "kind": "enum",
            "default_units": "",
            "enumeration": tuple(orca_step.metadata["functionals"].keys()),
            "format_string": "",
            "description": "Functional:",
            "help_text": (
                "The exchange-correlation functional (an ORCA '!' keyword), "
                "restricted to the chosen functional type. Only used when the "
                "method is 'DFT'. Double-hybrid functionals need an auxiliary "
                "'/C' basis, which 'AutoAux' provides."
            ),
        },
        "basis": {
            "applies_when": {"use model chemistry": "no"},
            "default": "def2-TZVP",
            "kind": "special",
            "widget": "seamm_widgets.BasisSetField",
            "default_units": "",
            "enumeration": tuple(orca_step.metadata["basis sets"]),
            "format_string": "",
            "description": "Basis set:",
            "help_text": (
                "The orbital basis set. Type a name, pick a common one from the "
                "list, or press '...' to choose any basis from the Basis Set "
                "Exchange (filtered to the elements you need); a Basis Set "
                "Exchange choice is stored as 'bse:NAME'. How a typed name is "
                "resolved is set by 'Basis set source' below."
            ),
        },
        "basis source": {
            "applies_when": {"use model chemistry": "no"},
            "default": "ORCA internal",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("ORCA internal", "Basis Set Exchange"),
            "format_string": "",
            "description": "Basis set source:",
            "help_text": (
                "Where the orbital basis comes from. 'ORCA internal' uses ORCA's "
                "built-in definition (the basis name goes on the '!' line). "
                "'Basis Set Exchange' fetches the named basis from the Basis Set "
                "Exchange and embeds it, for cross-code-identical definitions or "
                "bases ORCA does not ship. The auxiliary basis (AutoAux) is "
                "unaffected."
            ),
        },
        "basis set extrapolation": {
            "applies_when": {"use model chemistry": "no"},
            "default": "none",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("none", "2/3", "3/4", "4/5"),
            "format_string": "",
            "description": "Basis-set extrapolation (CBS):",
            "help_text": (
                "Extrapolate the energy to the complete-basis-set limit from two "
                "successive cardinal numbers (ORCA's Extrapolate keyword): '2/3' "
                "uses double- and triple-zeta, '3/4' triple/quadruple, '4/5' "
                "quadruple/quintuple. This is a single ORCA job that runs both "
                "basis sets and extrapolates the SCF and correlation parts "
                "separately. When on, the fixed basis set above is ignored, and "
                "gradients are not available (ORCA has no gradient for an "
                "extrapolated energy)."
            ),
        },
        "extrapolation family": {
            "applies_when": {"basis set extrapolation": {"not": "none"}},
            "default": "cc",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("cc", "aug-cc", "def2", "ANO"),
            "format_string": "",
            "description": "Extrapolation family:",
            "help_text": (
                "The basis-set family for the extrapolation: 'cc' (cc-pVnZ), "
                "'aug-cc' (aug-cc-pVnZ, with diffuse functions), 'def2' (the "
                "Karlsruhe def2 sets), or 'ANO'. Only used when basis-set "
                "extrapolation is on."
            ),
        },
        "auxiliary basis": {
            "default": "AutoAux",
            "kind": "enum",
            "default_units": "",
            "enumeration": tuple(orca_step.metadata["auxiliary basis sets"]),
            "format_string": "",
            "description": "Auxiliary (fitting) basis:",
            "help_text": (
                "The auxiliary/fitting basis. 'AutoAux' generates a fitting basis "
                "automatically and is the robust choice for correlated methods "
                "(DLPNO, MP2). 'none' omits it."
            ),
        },
        "grid": {
            "default": "default",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("default", "DEFGRID1", "DEFGRID2", "DEFGRID3"),
            "format_string": "",
            "description": "Integration grid:",
            "help_text": (
                "ORCA's numerical integration grid preset (for the DFT "
                "exchange-correlation and RIJCOSX/COSX grids). 'default' leaves "
                "ORCA's own default (DEFGRID2). DEFGRID1 is coarser/faster, "
                "DEFGRID3 is finer/more accurate. Only affects methods that use a "
                "grid (DFT, RIJCOSX); ignored otherwise."
            ),
        },
        "scf convergence": {
            "default": "TIGHTSCF",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "default",
                "SLOPPYSCF",
                "LOOSESCF",
                "NORMALSCF",
                "STRONGSCF",
                "TIGHTSCF",
                "VERYTIGHTSCF",
                "EXTREMESCF",
            ),
            "format_string": "",
            "description": "SCF convergence:",
            "help_text": (
                "ORCA's SCF convergence-tolerance preset. 'default' leaves ORCA's "
                "own default (NORMALSCF for a single point; ORCA tightens it to "
                "TIGHTSCF for optimizations). The presets run SLOPPYSCF (loosest) "
                "-> LOOSESCF -> NORMALSCF -> STRONGSCF -> TIGHTSCF -> VERYTIGHTSCF "
                "-> EXTREMESCF (tightest). TIGHTSCF (the default here) is a good "
                "choice for smooth energies/forces."
            ),
        },
        "sthresh": {
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": ("default", "1.0e-07"),
            "format_string": "",
            "description": "SCF SThresh:",
            "help_text": (
                "ORCA's SCF convergence threshold ('SThresh' in the '%scf' "
                "block), in E_h. 'default' emits nothing, so the "
                "SCF-convergence preset above (or ORCA's own default) governs "
                "SThresh. Any explicit value -- e.g. ORCA's nominal 1.0e-07 --"
                " is written to '%scf SThresh', overriding whatever the preset "
                "would otherwise set. Lower it for a tighter SCF (smoother "
                "energies/forces), raise it to converge more loosely."
            ),
        },
        "initial guess": {
            "default": "default",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "default",
                "Hueckel",
                "HCore",
                "PAtom",
                "PModel",
                "SAD",
                "SADNO",
                "Previous wavefunction",
                "Specified orbitals",
            ),
            "format_string": "",
            "description": "Initial guess:",
            "help_text": (
                "ORCA's SCF starting guess ('Guess' in the '%scf' block). "
                "'default' leaves ORCA's own default (SAD for most systems). "
                "'PModel' is often much more stable than SAD for a single atom "
                "or ion, where a superposition of atomic densities has little "
                "meaning. 'Previous wavefunction' seeds from the nearest "
                "earlier ORCA step in this flowchart (its 'orca.gbw'); "
                "'Specified orbitals' seeds from the file named by "
                "'Specified orbitals' below. Either way ORCA projects the "
                "orbitals onto this job's basis when it differs (the same "
                "trick its own CBS extrapolation uses internally), so this "
                "also works across a basis-set escalation. 'If wavefunction "
                "not found' below controls what happens when there is "
                "nothing to read."
            ),
        },
        "if wavefunction not found": {
            "applies_when": {
                "initial guess": ["Previous wavefunction", "Specified orbitals"]
            },
            "default": "Throw an error",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "Throw an error",
                "Use default guess",
                "Use Hueckel guess",
                "Use HCore guess",
                "Use PAtom guess",
                "Use PModel guess",
                "Use SAD guess",
                "Use SADNO guess",
            ),
            "format_string": "",
            "description": "If wavefunction not found:",
            "help_text": (
                "What to do when 'Initial guess' above is 'Previous "
                "wavefunction' or 'Specified orbitals' but nothing is "
                "there to read (e.g. no earlier ORCA step yet, or the first "
                "basis set of a loop where 'Specified orbitals' has not been "
                "written yet). 'Throw an error' (the default) fails loudly, "
                "matching what an explicit wavefunction request implies. "
                "Choose one of the 'Use ... guess' options instead to make "
                "this safe to leave on for every iteration of a loop, e.g. "
                "the first pass through a basis-set escalation."
            ),
        },
        "save orbital checkpoint": {
            "default": "no",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("no", "yes"),
            "format_string": "",
            "description": "Save orbital checkpoint:",
            "help_text": (
                "After a successful run, copy the converged orbitals "
                "('orca.gbw') to the file named by 'Checkpoint name' below, "
                "for a later step to read back -- e.g. with 'Initial guess' "
                "= 'Specified orbitals' there, naming the same "
                "checkpoint via 'Specified orbitals'. This reaches across a "
                "later iteration of an enclosing loop (e.g. over basis sets "
                "for the same atom), which 'Previous wavefunction' (a "
                "graph walk) cannot: each loop iteration gets its own "
                "directory, so the previous iteration is not a preceding "
                "node in the flowchart graph."
            ),
        },
        "checkpoint name": {
            "applies_when": {"save orbital checkpoint": "yes"},
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Checkpoint name:",
            "help_text": (
                "Where 'Save orbital checkpoint' writes this run's orbitals. "
                "'default' (the default) derives a label automatically from "
                "the system's composition, charge, and multiplicity (e.g. "
                "'Co_q0_m4') -- the usual choice, since it naturally gives "
                "one checkpoint per atom/electronic state, shared across an "
                "inner loop (e.g. over basis sets) for that atom, and reset "
                "automatically when an outer loop moves to a new atom. Any "
                "other bare name is stored in a 'checkpoints' folder inside "
                "this job; an absolute path is used as-is, e.g. to keep the "
                "checkpoint outside this job and reuse it across separate "
                "flowchart runs."
            ),
        },
        "specified orbitals": {
            "applies_when": {"initial guess": "Specified orbitals"},
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Specified orbitals:",
            "help_text": (
                "Which file 'Initial guess' = 'Specified orbitals' reads "
                "from. Same rules as 'Checkpoint name': 'default' derives "
                "the label automatically from this system's composition, "
                "charge, and multiplicity, matching what a step upstream "
                "saved with 'Save orbital checkpoint' left on 'default'; a "
                "bare name looks in this job's 'checkpoints' folder; an "
                "absolute path is used as-is. This is the one field that "
                "may also point at ANOTHER job (reading only -- a job "
                "cannot write into another job): 'job://<job number>/<name>' "
                "-- or 'job://<job number>/default' to pick up that other "
                "job's auto-derived name for this same system, when you do "
                "not know what it resolved to."
            ),
        },
        "extra keywords": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Extra keywords:",
            "help_text": (
                "Any additional ORCA '!' keywords to append, e.g. 'RIJCOSX', "
                "'NoFrozenCore', 'SlowConv'. The SCF tolerance and integration "
                "grid have their own controls above."
            ),
        },
        "extra blocks": {
            "default": "",
            "kind": "special",
            "widget": "seamm_widgets.LabeledText",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Extra ORCA blocks:",
            "help_text": (
                "Any additional literal ORCA input, inserted verbatim right "
                "before the geometry -- after the blocks generated by the "
                "controls above. Use this for one-off SCF-stabilization "
                "tricks with no dedicated control, e.g. for a difficult atom:"
                "\n%scf\n  MaxIter 400\n"
                "  Shift Shift 0.3 ErrStart 0.05 end\n  DIISBfac 1.1\nend"
            ),
        },
        "bond orders": {
            "default": "yes",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("no", "yes", "yes, and apply to structure"),
            "format_string": "",
            "description": "Mayer bond orders:",
            "help_text": (
                "Analyze the Mayer bond orders (always written to a CSV file; "
                "printed for small systems). 'apply to structure' also replaces "
                "the bonds in the structure with single/aromatic/double/triple "
                "bonds based on the bond orders."
            ),
        },
        "Hirshfeld charges": {
            "default": "no",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("no", "yes", "yes, and apply to structure"),
            "format_string": "",
            "description": "Hirshfeld charges:",
            "help_text": (
                "Compute Hirshfeld atomic charges (written to a CSV file; printed "
                "for small systems). 'apply to structure' also stores them as the "
                "atomic charges on the structure."
            ),
        },
        "polarizability": {
            "default": "no",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("no", "yes"),
            "format_string": "",
            "description": "Polarizability:",
            "help_text": (
                "Compute the dipole polarizability (analytic for HF and DFT). "
                "This adds to the cost of the calculation."
            ),
        },
        "save wavefunction": {
            "default": "no",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("no", "yes"),
            "format_string": "",
            "description": "Write the wavefunction (wfx) file:",
            "help_text": (
                "Retain the electron density ('keepdensity') and convert it to an "
                "AIMPAC wavefunction (.wfx) file with orca_2aim. This analytic "
                "wavefunction is read by a following Atomic Charges step "
                "(DDEC6 via Chargemol), mirroring the Gaussian wfx path."
            ),
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    # Rules shared by the dialog and the flowchart builder (see seamm.Parameters).
    # The simple conditions are the "applies_when" entries above.

    extrapolation = True
    """Whether basis-set extrapolation can be used. False for steps that need a
    gradient or Hessian, which an extrapolated energy does not have."""

    unused = ()
    """Parameters that this kind of step never uses."""

    @staticmethod
    def _is_f12(values):
        return "F12" in str(values.get("method", "")).upper()

    def _extrapolating(self, values):
        return self.applies("basis set extrapolation", values) and values.get(
            "basis set extrapolation"
        ) not in (None, "none")

    def applies(self, key, values=None, _seen=None):
        """As seamm.Parameters.applies, plus: no extrapolation for F12 methods (they
        need their own F12 basis) or for steps that need a gradient; the basis and
        its source are replaced by the extrapolation when it is used; and the basis
        source is always ORCA's own for F12 (the Basis Set Exchange has no CABS)."""
        if values is None:
            values = self.current_values()
        if key in self.unused:
            return False
        if not super().applies(key, values, _seen):
            return False
        if key == "basis set extrapolation":
            return self.extrapolation and not self._is_f12(values)
        if key == "basis":
            return not self._extrapolating(values)
        if key == "basis source":
            return not self._extrapolating(values) and not self._is_f12(values)
        return True

    def not_applicable_reason(self, key, values=None):
        """Why a parameter does not apply, for the builder's messages."""
        if values is None:
            values = self.current_values()
        if key in self.unused:
            return "this kind of step does not use it"
        reason = super().not_applicable_reason(key, values)
        if reason:
            return reason
        if key == "basis set extrapolation":
            if not self.extrapolation:
                return (
                    "this step needs a gradient (or Hessian), which an extrapolated "
                    "energy does not have"
                )
            if self._is_f12(values):
                return "F12 methods use their own F12 basis, not an extrapolation"
        if key in ("basis", "basis source") and self._extrapolating(values):
            return "the basis-set extrapolation replaces it"
        if key == "basis source" and self._is_f12(values):
            return "F12 methods use ORCA's own basis sets (they need a CABS)"
        return ""

    def choices(self, key, values=None):
        """F12 methods take only the F12 bases; the functionals are those of the
        functional type."""
        if values is None:
            values = self.current_values()
        if key == "basis" and self._is_f12(values):
            return tuple(
                b
                for b in orca_step.metadata["basis sets"]
                if b.upper().endswith("-F12")
            )
        if key == "functional":
            ftype = values.get("functional type")
            if not self._is_expr(ftype):
                return tuple(
                    name
                    for name, record in orca_step.metadata["functionals"].items()
                    if record["category"] == ftype
                )
        return super().choices(key, values)

    def implied(self, values=None):
        """An F12 method needs an F12 basis from ORCA itself; a functional implies
        its functional type."""
        if values is None:
            values = self.current_values()
        result = {}
        if self.applies("method", values) and self._is_f12(values):
            f12 = self.choices("basis", values)
            basis = values.get("basis")
            name = basis.get("name") if isinstance(basis, dict) else basis
            if name not in f12 and not self._is_expr(name):
                result["basis"] = "cc-pVTZ-F12" if "cc-pVTZ-F12" in f12 else f12[0]
            result["basis source"] = "ORCA internal"
        if self.applies("functional", values):
            functional = values.get("functional")
            record = orca_step.metadata["functionals"].get(functional)
            if record is not None:
                result["functional type"] = record["category"]
        return result

    # Methods renamed since flowcharts were saved: old name -> current name
    renamed_methods = {
        # 2026-07-13: ORCA rejects the bare keyword; the canonical one has /RI
        "CCSD(T)-F12D": "CCSD(T)-F12D/RI",
    }

    def __init__(self, defaults={}, data=None):
        """Initialize with the parameters above plus any overrides."""
        logger.debug("EnergyParameters.__init__")
        if data is not None:
            entry = data.get("method")
            if isinstance(entry, dict) and entry.get("value") in self.renamed_methods:
                entry["value"] = self.renamed_methods[entry["value"]]
        super().__init__(
            defaults={**EnergyParameters.parameters, **defaults}, data=data
        )
