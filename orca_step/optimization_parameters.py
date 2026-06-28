# -*- coding: utf-8 -*-
"""Control parameters for an ORCA geometry optimization."""

import logging

from .energy_parameters import EnergyParameters

logger = logging.getLogger(__name__)


class OptimizationParameters(EnergyParameters):
    """Optimization parameters: the energy parameters plus optimization controls."""

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
            defaults={**OptimizationParameters.parameters, **defaults}, data=data
        )
