# -*- coding: utf-8 -*-

"""Control parameters for the ORCA BSSE (counterpoise) sub-step.

Extends the Energy parameters (level of theory, basis, grid, SCF controls) with
the fragment definition and the monomer-relaxation switch that the counterpoise
correction needs.
"""

import logging

from .energy_parameters import EnergyParameters

logger = logging.getLogger(__name__)


class BSSEParameters(EnergyParameters):
    """The counterpoise (BSSE) parameters: the Energy parameters plus the
    fragment definition and the monomer-optimization switch."""

    parameters = {
        "fragments": {
            "default": "auto (2 molecules)",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("auto (2 molecules)", "specified"),
            "format_string": "",
            "description": "Fragments:",
            "help_text": (
                "How to split the complex into the two fragments for the "
                "counterpoise correction. 'auto (2 molecules)' uses the two "
                "separate molecules found in the structure (an error if there "
                "are not exactly two). 'specified' takes the atoms of fragment A "
                "from the field below; the rest of the atoms are fragment B."
            ),
        },
        "fragment A atoms": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Fragment A atoms:",
            "help_text": (
                "The atoms making up fragment A when 'Fragments' is 'specified', "
                "given as atom numbers (1-based, as shown in the structure) -- a "
                "comma/space separated list and/or ranges, e.g. '1-3, 5 7'. The "
                "remaining atoms form fragment B (written as ghosts). Ignored "
                "when the fragments are found automatically."
            ),
        },
        "optimize monomers": {
            "default": "no",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("no", "yes"),
            "format_string": "",
            "description": "Optimize free monomers:",
            "help_text": (
                "Whether to relax each free monomer before taking the "
                "correction (the ORCA script's 'DoOptimization'). 'no' (the "
                "default) is correct for a fixed-geometry potential-energy "
                "surface or MLFF training target; 'yes' gives the counterpoise "
                "correction relative to the relaxed monomers."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        logger.debug("BSSEParameters.__init__")
        super().__init__(defaults={**BSSEParameters.parameters, **defaults}, data=data)
