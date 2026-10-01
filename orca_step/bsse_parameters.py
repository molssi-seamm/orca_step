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
            "applies_when": {"fragments": "specified"},
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
                "The last group may be 'rest' for the atoms in no other group, "
                "e.g. '1-3; rest'. Ignored when the fragments are found "
                "automatically."
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
                "Na+/Cl- pair. Leave empty to use each fragment's net formal "
                "charge from the input structure, if the structure format "
                "carries one (e.g. an ion marked with an SDF/MOL 'M  CHG' "
                "record) -- otherwise all-neutral (the usual case for a "
                "neutral H-bonded complex). Every fragment is closed-shell "
                "(multiplicity 1); the complex's own charge and multiplicity "
                "are checked against these for consistency."
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

    # The counterpoise gradient needs a real gradient, which an extrapolated energy
    # lacks; and these Energy settings are not used by the BSSE step.
    extrapolation = False
    unused = (
        "sthresh",
        "initial guess",
        "save orbital checkpoint",
        "checkpoint name",
        "extra blocks",
        "bond orders",
        "Hirshfeld charges",
        "polarizability",
    )

    def __init__(self, defaults={}, data=None):
        logger.debug("BSSEParameters.__init__")
        if data is not None:
            self._translate_old(data)
        super().__init__(defaults={**BSSEParameters.parameters, **defaults}, data=data)

    @staticmethod
    def _translate_old(data):
        """Flowcharts saved before BSSE took N fragments (2026-08-04) had two:
        'auto (2 molecules)', or 'fragment A atoms' with the rest as fragment B."""
        entry = data.get("fragments")
        if isinstance(entry, dict) and entry.get("value") == "auto (2 molecules)":
            entry["value"] = "auto (molecules)"
        if "fragment A atoms" in data:
            old = data.pop("fragment A atoms")
            atoms = old.get("value") if isinstance(old, dict) else old
            if atoms and "fragment atoms" not in data:
                data["fragment atoms"] = {"value": f"{atoms}; rest", "units": None}
