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
        "extra keywords": {
            "default": "TightSCF",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Extra keywords:",
            "help_text": (
                "Any additional ORCA '!' keywords to append, e.g. 'TightSCF', "
                "'RIJCOSX', 'Grid5'."
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

    def __init__(self, defaults={}, data=None):
        """Initialize with the parameters above plus any overrides."""
        logger.debug("EnergyParameters.__init__")
        super().__init__(
            defaults={**EnergyParameters.parameters, **defaults}, data=data
        )
