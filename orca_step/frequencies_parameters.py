# -*- coding: utf-8 -*-

"""Control parameters for the ORCA Frequencies (Hessian) sub-step.

Extends the Energy parameters (level of theory, basis, ...) with the choice of
analytic vs numerical second derivatives and the thermochemistry temperature.
"""

import logging

from .energy_parameters import EnergyParameters

logger = logging.getLogger(__name__)


class FrequenciesParameters(EnergyParameters):
    """The Frequencies parameters: the Energy parameters plus the second-
    derivative method and the thermochemistry temperature."""

    parameters = {
        "second derivatives": {
            "default": "analytic",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("analytic", "numerical"),
            "format_string": "",
            "description": "Second derivatives:",
            "help_text": (
                "How to compute the Hessian. 'analytic' (ORCA's AnFreq) is much "
                "faster but requires an analytic second derivative for the method "
                "(available for HF, most DFT functionals, and MP2). 'numerical' "
                "(NumFreq) finite-differences the gradient -- it works for any "
                "method that has a gradient, but is considerably more expensive."
            ),
        },
        "temperature": {
            "default": "298.15",
            "kind": "float",
            "default_units": "K",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Temperature:",
            "help_text": (
                "The temperature for the thermochemistry (zero-point energy, "
                "thermal corrections, entropy, and Gibbs free energy). The "
                "pressure is ORCA's default of 1 atm."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        logger.debug("FrequenciesParameters.__init__")
        super().__init__(
            defaults={**FrequenciesParameters.parameters, **defaults}, data=data
        )
