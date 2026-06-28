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
            "help_text": "The ORCA method (an ORCA '!' keyword).",
        },
        "basis": {
            "default": "def2-TZVP",
            "kind": "enum",
            "default_units": "",
            "enumeration": tuple(orca_step.metadata["basis sets"]),
            "format_string": "",
            "description": "Basis set:",
            "help_text": (
                "The orbital basis set. ORCA's built-in name is used; a valid "
                "name may also be typed in. (Basis Set Exchange is a planned "
                "opt-in alternative.)"
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
