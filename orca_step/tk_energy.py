# -*- coding: utf-8 -*-

"""The graphical part of an ORCA Energy sub-step."""

import logging
import tkinter as tk

import orca_step  # noqa: F401
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw  # noqa: F401

logger = logging.getLogger(__name__)


class TkEnergy(seamm.TkNode):
    """The graphical part of an ORCA Energy sub-step.

    See Also
    --------
    Energy, EnergyParameters
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        self.dialog = None
        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )

    def create_dialog(self, title="ORCA Energy"):
        """Create the dialog and its widgets."""
        frame = super().create_dialog(title=title)
        P = self.node.parameters

        for key in P:
            if key != "results":
                self[key] = P[key].widget(frame)

        # React to the model-chemistry toggle.
        for item in ("use model chemistry",):
            w = self[item]
            w.combobox.bind("<<ComboboxSelected>>", self.reset_dialog)
            w.combobox.bind("<Return>", self.reset_dialog)
            w.combobox.bind("<FocusOut>", self.reset_dialog)

        self.reset_dialog()
        return frame

    def reset_dialog(self, widget=None):
        """Lay out the widgets, hiding the explicit method/basis controls when
        the model chemistry is used.
        """
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        use_mc = self["use model chemistry"].get() == "yes"

        row = 0
        self["use model chemistry"].grid(row=row, column=0, sticky=tk.EW)
        row += 1

        widgets = [self["use model chemistry"]]
        # Method/basis come from the model chemistry when it is used; the
        # auxiliary basis and extra keywords are ORCA run details either way.
        keys = ["auxiliary basis", "extra keywords"]
        if not use_mc:
            keys = ["method", "basis"] + keys
        for key in keys:
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        sw.align_labels(widgets, sticky=tk.E)
        return row
