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
    }

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
        tmp._data["enumeration"] = ("optimized", *tmp.enumeration)
        tmp.default = "optimized"
