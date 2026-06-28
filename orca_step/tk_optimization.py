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

    def create_dialog(self, title="ORCA Optimization"):
        frame = super().create_dialog(title=title)
        P = self.node.parameters
        if "optimization convergence" not in self:
            self["optimization convergence"] = P["optimization convergence"].widget(
                frame
            )
        self.reset_dialog()
        return frame

    def reset_dialog(self, widget=None):
        row = super().reset_dialog(widget)
        if "optimization convergence" in self:
            self["optimization convergence"].grid(row=row, column=0, sticky=tk.EW)
            row += 1
        return row
