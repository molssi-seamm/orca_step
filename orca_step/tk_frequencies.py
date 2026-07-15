# -*- coding: utf-8 -*-

"""The graphical part of an ORCA Frequencies sub-step."""

import logging

import orca_step  # noqa: F401
from .tk_energy import TkEnergy
import seamm_widgets as sw  # noqa: F401

logger = logging.getLogger(__name__)


class TkFrequencies(TkEnergy):
    """Graphical ORCA Frequencies sub-step: the energy dialog's level-of-theory
    controls plus the second-derivative choice and the thermochemistry
    temperature. CBS extrapolation is hidden (an extrapolated energy has no
    Hessian); the Energy property toggles are not shown.
    """

    def _show_cbs(self):
        """Hide CBS extrapolation: a frequency calculation needs a Hessian,
        which an extrapolated energy does not have."""
        return False

    def _run_detail_keys(self):
        return (
            "auxiliary basis",
            "grid",
            "scf convergence",
            "extra keywords",
            "second derivatives",
            "temperature",
        )
