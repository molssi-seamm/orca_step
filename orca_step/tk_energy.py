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
        frame = super().create_dialog(title=title, widget="notebook", results_tab=True)
        P = self.node.parameters

        for key in P:
            if key != "results":
                self[key] = P[key].widget(frame)

        # The basis field is a seamm_widgets.BasisSetField (entry/combobox + a
        # '...' button to the Basis Set Exchange). Give it the curated quick-pick
        # names and a way to preselect the current system's elements.
        self["basis"].config(values=list(orca_step.metadata["basis sets"]))
        self["basis"].elements_callback = self._current_elements

        # React to the model-chemistry toggle.
        for item in ("use model chemistry",):
            w = self[item]
            w.combobox.bind("<<ComboboxSelected>>", self.reset_dialog)
            w.combobox.bind("<Return>", self.reset_dialog)
            w.combobox.bind("<FocusOut>", self.reset_dialog)

        self.reset_dialog()
        return frame

    def _current_elements(self):
        """Element symbols in the current configuration, to preselect in the
        Basis Set Exchange dialog. Best-effort: empty if there is none yet."""
        try:
            _, configuration = self.node.get_system_configuration(None)
            return sorted(set(configuration.atoms.symbols))
        except Exception:
            return []

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
        # Method, basis, and the basis source come from the model chemistry when
        # it is used (a 'bse:' basis there forces the Basis Set Exchange), so they
        # are hidden in that mode. The auxiliary basis and the rest are ORCA run
        # details that apply either way.
        keys = [
            "auxiliary basis",
            "extra keywords",
            "bond orders",
            "Hirshfeld charges",
            "polarizability",
            "save wavefunction",
        ]
        if not use_mc:
            keys = ["method", "basis", "basis source"] + keys
        for key in keys:
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        sw.align_labels(widgets, sticky=tk.E)

        # Lay out the Results tab from the metadata (energy, gradients, charges,
        # ...). Without this the Results tab is created but stays empty.
        self.setup_results()
        return row
