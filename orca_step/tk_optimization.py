# -*- coding: utf-8 -*-

"""The graphical part of an ORCA Optimization sub-step."""

import logging
import tkinter as tk

import orca_step  # noqa: F401
from .tk_energy import TkEnergy
import seamm_widgets as sw  # noqa: F401

logger = logging.getLogger(__name__)


class TkOptimization(TkEnergy):
    """Graphical ORCA Optimization sub-step: the energy dialog plus the
    optimization-convergence control.
    """

    def _show_cbs(self):
        """Hide the CBS basis-set-extrapolation controls: a geometry
        optimization needs a gradient, which an extrapolated energy does not
        have."""
        return False

    def create_dialog(self, title="ORCA Optimization"):
        frame = super().create_dialog(title=title)
        P = self.node.parameters
        if "optimization convergence" not in self:
            self["optimization convergence"] = P["optimization convergence"].widget(
                frame
            )
        for key in ("structure handling", "system name", "configuration name"):
            if key not in self:
                self[key] = P[key].widget(frame)
        self.reset_dialog()
        return frame

    def reset_dialog(self, widget=None):
        row = super().reset_dialog(widget)
        if "optimization convergence" in self:
            self["optimization convergence"].grid(row=row, column=0, sticky=tk.EW)
            row += 1
        for key in ("structure handling", "system name", "configuration name"):
            if key in self:
                self[key].grid(row=row, column=0, sticky=tk.EW)
                row += 1
        return row
