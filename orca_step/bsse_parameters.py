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
            "default": "auto (molecules)",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("auto (molecules)", "specified"),
            "format_string": "",
            "description": "Fragments:",
            "help_text": (
                "How to split the complex into fragments for the counterpoise "
                "correction. 'auto (molecules)' uses every separate molecule "
                "found in the structure (an error if there are fewer than "
                "two -- e.g. a Na+/Cl- pair is two molecules). 'specified' "
                "takes the fragments from 'Fragment atoms' below."
            ),
        },
        "fragment atoms": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Fragment atoms:",
            "help_text": (
                "The atoms making up each fragment when 'Fragments' is "
                "'specified': one semicolon-separated group per fragment, each "
                "a comma/space list and/or ranges of 1-based atom numbers, "
                "e.g. '1-3; 4-6' for two fragments or '1-3; 4-6; 7' for three. "
                "Ignored when the fragments are found automatically."
            ),
        },
        "fragment charges": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Fragment charges:",
            "help_text": (
                "The formal charge of each fragment, in the same order as the "
                "fragments (the order 'auto' finds molecules, or the "
                "semicolon-separated groups in 'Fragment atoms') -- a "
                "comma/space separated list of integers, e.g. '1, -1' for a "
                "Na+/Cl- pair. Leave empty for all-neutral fragments (the usual "
                "case for a neutral H-bonded complex). Every fragment is "
                "closed-shell (multiplicity 1); the complex's own charge and "
                "multiplicity are checked against these for consistency."
            ),
        },
        "compute gradient": {
            "default": "yes",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Compute the gradient:",
            "help_text": (
                "Whether to compute the counterpoise-corrected gradient (forces) "
                "as well as the energy. 'yes' (the default) is what MLFF training "
                "needs. 'no' (energy only) is cheaper and, importantly, allows "
                "methods that have no analytic gradient in ORCA -- notably "
                "CCSD(T) / DLPNO-CCSD(T) -- for gold-standard counterpoise "
                "interaction energies."
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
