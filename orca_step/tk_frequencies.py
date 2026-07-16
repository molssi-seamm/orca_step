# -*- coding: utf-8 -*-

"""The graphical part of an ORCA Frequencies sub-step."""

import logging
import tkinter as tk

import orca_step  # noqa: F401
from .tk_energy import TkEnergy
import seamm_widgets as sw

logger = logging.getLogger(__name__)


class TkFrequencies(TkEnergy):
    """Graphical ORCA Frequencies sub-step: the energy dialog's level-of-theory
    controls plus the second-derivative choice, the thermochemistry temperature,
    and the standard structure-handling options (where to store the structure and
    its properties). CBS extrapolation is hidden (an extrapolated energy has no
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

    def reset_dialog(self, widget=None):
        """Lay out the energy controls, then the structure-handling options."""
        row = super().reset_dialog(widget)
        widgets = []
        for key in ("structure handling", "system name", "configuration name"):
            if key in self:
                self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                widgets.append(self[key])
                row += 1
        sw.align_labels(widgets, sticky=tk.E)
        return row
