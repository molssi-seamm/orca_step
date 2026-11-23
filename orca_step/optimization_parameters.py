# -*- coding: utf-8 -*-
"""Control parameters for an ORCA geometry optimization."""

import logging

import seamm

from .energy_parameters import EnergyParameters

logger = logging.getLogger(__name__)


class OptimizationParameters(EnergyParameters):
    """Optimization parameters: the energy parameters plus optimization controls
    and the standard structure-handling options (where to put the optimized
    geometry)."""

    parameters = {
        "optimization convergence": {
            "default": "NormalOpt",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("LooseOpt", "NormalOpt", "TightOpt", "VeryTightOpt"),
            "format_string": "",
            "description": "Convergence:",
            "help_text": "The ORCA geometry-optimization convergence preset.",
        },
        "if not converged": {
            # What happens after the calculation, not what it computes: left out
            # of the fingerprint a resumed job must match.
            "resource": True,
            "default": "warn and continue",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("warn and continue", "fail the step"),
            "format_string": "",
            "description": "If not converged:",
            "help_text": (
                "What to do when ORCA reports that the optimization did not "
                "converge. 'warn and continue' prints a warning and keeps the last "
                "geometry, as before. 'fail the step' stops it with an error before "
                "anything is stored, so a Loop counts the iteration as failed and "
                "the steps after it never use an unconverged structure. Either way "
                "the properties 'optimization converged' and 'optimization cycles' "
                "can be stored with the results."
            ),
        },
    }

    # An optimization needs a gradient, which an extrapolated energy lacks.
    extrapolation = False

    def __init__(self, defaults={}, data=None):
        logger.debug("OptimizationParameters.__init__")
        super().__init__(
            defaults={
                **OptimizationParameters.parameters,
                **seamm.standard_parameters.structure_handling_parameters,
                **defaults,
            },
            data=data,
        )

        # Default to overwriting the current configuration so the optimized
        # geometry flows into any following sub-step (e.g. a frequency
        # calculation) unless the user asks otherwise.
        self["structure handling"].description = "Structure handling:"

        tmp = self["system name"]
        tmp._data["enumeration"] = (*tmp.enumeration, "optimized")
        tmp.default = "keep current name"

        tmp = self["configuration name"]
        tmp._data["enumeration"] = (
            "{current name} optimized with {model}",
            "optimized",
            *tmp.enumeration,
        )
        tmp.default = "{current name} optimized with {model}"
