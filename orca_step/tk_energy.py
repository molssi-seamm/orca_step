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

        # React to the model-chemistry toggle and the method choice (the
        # method controls whether the DFT functional pulldowns are shown),
        # and to the wavefunction-guess controls that reveal/hide their own
        # indented sub-controls (see reset_dialog).
        for item in (
            "use model chemistry",
            "method",
            "initial guess",
            "save orbital checkpoint",
        ):
            w = self[item]
            w.combobox.bind("<<ComboboxSelected>>", self.reset_dialog)
            w.combobox.bind("<Return>", self.reset_dialog)
            w.combobox.bind("<FocusOut>", self.reset_dialog)

        # Changing the functional type re-filters the functional list in place.
        w = self["functional type"]
        w.combobox.bind("<<ComboboxSelected>>", self.reset_functionals)
        w.combobox.bind("<Return>", self.reset_functionals)
        w.combobox.bind("<FocusOut>", self.reset_functionals)

        # Choosing the Basis Set Exchange as the source opens the periodic-table
        # picker so the user can specify which BSE basis (also on the '...' btn).
        self["basis source"].combobox.bind(
            "<<ComboboxSelected>>", self._on_basis_source
        )

        # Turning CBS extrapolation on/off swaps the basis controls in and out.
        w = self["basis set extrapolation"]
        w.combobox.bind("<<ComboboxSelected>>", self.reset_dialog)
        w.combobox.bind("<Return>", self.reset_dialog)
        w.combobox.bind("<FocusOut>", self.reset_dialog)

        self.reset_dialog()
        return frame

    def right_click(self, event):
        """Post the node's popup menu. The base class builds the menu (with the
        Delete command) but only shows it for the bare TkNode, so each sub-step
        must add its own items and pop it up -- without this, right-clicking a
        sub-step in the ORCA sub-flowchart shows no menu (so it can't be
        deleted). Inherited by the Optimization and BSSE sub-steps."""
        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)
        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def _current_elements(self):
        """Element symbols in the current configuration, to preselect in the
        Basis Set Exchange dialog. Best-effort: empty if there is none yet."""
        try:
            _, configuration = self.node.get_system_configuration(None)
            return sorted(set(configuration.atoms.symbols))
        except Exception:
            return []

    def reset_dialog(self, widget=None):
        """Lay out the widgets. Hide the explicit method/basis controls when the
        model chemistry is used; for DFT, show the functional-type and functional
        pulldowns indented one and two levels under Method, mirroring Gaussian.
        """
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()
        # Clear any indentation left from a previous (DFT) layout.
        frame.columnconfigure(0, minsize=0)
        frame.columnconfigure(1, minsize=0)

        # Which controls to show, and what they offer, come from the parameters'
        # rules (orca_step.EnergyParameters), which the flowchart builder uses too.
        P = self.node.parameters
        values = self._widget_values()

        def applies(key):
            return P.applies(key, values)

        row = 0
        widgets = []  # full-width (column 0) controls
        type_widgets = []  # indented one level (column 1): functional type
        func_widgets = []  # indented two levels (column 2): functional

        def add_full(key):
            nonlocal row
            self[key].grid(row=row, column=0, columnspan=3, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        def add_indented(key):
            nonlocal row
            self[key].grid(row=row, column=1, columnspan=2, sticky=tk.EW)
            type_widgets.append(self[key])
            row += 1

        add_full("use model chemistry")

        # Method, basis, and basis source come from the model chemistry when it
        # is used, so they do not apply then. The auxiliary basis and the rest
        # are ORCA run details that apply either way.
        if applies("method"):
            add_full("method")
            # Narrow the basis list to what the method allows (F12 -> F12 bases).
            self._filter_basis_sets(values)
            values = self._widget_values()
            if applies("functional type"):
                # Functional type indented one level, functional two levels.
                self._filter_functionals(values)
                add_indented("functional type")
                self["functional"].grid(row=row, column=2, columnspan=1, sticky=tk.EW)
                func_widgets.append(self["functional"])
                row += 1
            # CBS extrapolation replaces the fixed basis. It does not apply to
            # sub-steps that need a gradient (e.g. Optimization) or to F12 methods.
            if applies("basis set extrapolation"):
                add_full("basis set extrapolation")
            if applies("extrapolation family"):
                add_full("extrapolation family")
            if applies("basis"):
                add_full("basis")
            if applies("basis source"):
                add_full("basis source")

        for key in self._run_detail_keys():
            if not applies(key):
                continue
            # 'checkpoint name' is shown indented under 'save orbital checkpoint'.
            if key == "checkpoint name":
                add_indented(key)
                continue

            add_full(key)

            # 'initial guess' choosing a wavefunction reveals its own indented
            # sub-controls right below it (see energy.extra_input).
            if key == "initial guess":
                for sub in ("specified orbitals", "if wavefunction not found"):
                    if applies(sub):
                        add_indented(sub)

        # Align the full-width labels; indent the nested widgets by the leftover
        # label width plus a fixed gap, so each nested combobox sits ~30 px to the
        # right of its parent's (the Gaussian idiom).
        width0 = sw.align_labels(widgets, sticky=tk.E)
        if type_widgets:
            width1 = sw.align_labels(type_widgets, sticky=tk.E)
            width2 = sw.align_labels(func_widgets, sticky=tk.E)
            frame.columnconfigure(0, minsize=max(0, width0 - width1) + 30)
            frame.columnconfigure(1, minsize=max(0, width1 - width2) + 30)

        # Lay out the Results tab from the metadata (energy, gradients, charges,
        # ...). Without this the Results tab is created but stays empty.
        self.setup_results()
        return row

    def _widget_values(self):
        """The dialog's current values, {name: value}, for the parameters' rules."""
        values = {}
        for key in self.node.parameters:
            if key == "results" or key not in self:
                continue
            try:
                value = self[key].get()
            except Exception:
                continue
            values[key] = value[0] if isinstance(value, tuple) else value
        return values

    def _run_detail_keys(self):
        """The full-width 'run detail' controls laid out below the level of
        theory, in order. Sub-steps override this to add or drop controls."""
        return (
            "auxiliary basis",
            "grid",
            "scf convergence",
            "sthresh",
            "initial guess",
            "save orbital checkpoint",
            "checkpoint name",
            "extra keywords",
            "extra blocks",
            "bond orders",
            "Hirshfeld charges",
            "polarizability",
            "save wavefunction",
        )

    def _on_basis_source(self, widget=None):
        """When the user selects the Basis Set Exchange as the source, open the
        picker so they can specify which basis (writes back as 'bse:NAME')."""
        if self["basis source"].get() == "Basis Set Exchange":
            browse = getattr(self["basis"], "_browse", None)
            if callable(browse):
                browse()

    def _filter_basis_sets(self, values=None):
        """Offer only the bases the method allows, and set what it requires.

        From the parameters' rules: explicitly-correlated F12 methods take only the
        F12-optimized bases (they need a matching CABS), from ORCA itself (the Basis
        Set Exchange has no CABS). Otherwise the full curated list is offered.
        """
        P = self.node.parameters
        values = self._widget_values() if values is None else values
        allowed = P.choices("basis", values)
        if allowed is None:
            allowed = orca_step.metadata["basis sets"]
        self["basis"].config(values=list(allowed))
        for key, value in P.implied(values).items():
            if key in ("basis", "basis source"):
                self[key].set(value)

    def _filter_functionals(self, values=None):
        """Offer only the functionals of the chosen functional type, keeping the
        selection valid."""
        P = self.node.parameters
        values = self._widget_values() if values is None else values
        funcs = list(P.choices("functional", values) or ())
        self["functional"].combobox.configure(values=funcs)
        if funcs and self["functional"].get() not in funcs:
            self["functional"].set(funcs[0])

    def reset_functionals(self, widget=None):
        """Re-filter the functional list when the functional type changes."""
        self._filter_functionals()
