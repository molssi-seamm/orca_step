# -*- coding: utf-8 -*-

"""Control parameters for the ORCA Frequencies (Hessian) sub-step.

Extends the Energy parameters (level of theory, basis, ...) with the choice of
analytic vs numerical second derivatives and the thermochemistry temperature.
"""

import logging

import seamm

from .energy_parameters import EnergyParameters

logger = logging.getLogger(__name__)


class FrequenciesParameters(EnergyParameters):
    """The Frequencies parameters: the Energy parameters plus the second-
    derivative method and the thermochemistry temperature."""

    parameters = {
        "second derivatives": {
            "default": "default",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("default", "analytic", "numerical"),
            "format_string": "",
            "description": "Second derivatives:",
            "help_text": (
                "How to compute the Hessian. 'default' uses ORCA's analytic "
                "second derivative (AnFreq) when one is available for the "
                "method (HF, most DFT functionals, MP2), and falls back to "
                "the numerical one (NumFreq) otherwise -- the right choice "
                "for almost every case. 'analytic' and 'numerical' force one "
                "or the other regardless of whether the method has an "
                "analytic Hessian: 'analytic' fails for a method without one "
                "(e.g. a double hybrid or (DLPNO-)CCSD(T)); 'numerical' "
                "(finite-differencing the gradient) works for any method "
                "that has a gradient, but is considerably more expensive."
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
            defaults={
                **FrequenciesParameters.parameters,
                **seamm.standard_parameters.structure_handling_parameters,
                **defaults,
            },
            data=data,
        )

        # A frequency calculation does not change the geometry, so default to
        # overwriting the current configuration (storing the results/properties
        # there) and keeping its name.
        self["structure handling"].description = "Structure handling:"
        self["configuration name"].default = "keep current name"
