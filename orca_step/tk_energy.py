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

        # React to the model-chemistry toggle and the method choice (the method
        # controls whether the DFT functional pulldowns are shown).
        for item in ("use model chemistry", "method"):
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

        use_mc = self["use model chemistry"].get() == "yes"
        is_dft = (not use_mc) and self["method"].get() == "DFT"
        is_f12 = (not use_mc) and "F12" in self["method"].get().upper()

        row = 0
        widgets = []  # full-width (column 0) controls
        type_widgets = []  # indented one level (column 1): functional type
        func_widgets = []  # indented two levels (column 2): functional

        def add_full(key):
            nonlocal row
            self[key].grid(row=row, column=0, columnspan=3, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        add_full("use model chemistry")

        # Method, basis, and basis source come from the model chemistry when it
        # is used, so they are hidden in that mode. The auxiliary basis and the
        # rest are ORCA run details that apply either way.
        if not use_mc:
            add_full("method")
            # Narrow the basis list to what the method allows (F12 -> F12 bases).
            self._filter_basis_sets()
            if is_dft:
                # Functional type indented one level, functional two levels.
                self._filter_functionals()
                self["functional type"].grid(
                    row=row, column=1, columnspan=2, sticky=tk.EW
                )
                type_widgets.append(self["functional type"])
                row += 1
                self["functional"].grid(row=row, column=2, columnspan=1, sticky=tk.EW)
                func_widgets.append(self["functional"])
                row += 1
            # CBS extrapolation replaces the fixed basis: show the family instead
            # of the basis/basis-source controls when it is on. It is hidden for
            # sub-steps that need a gradient (e.g. Optimization), where an
            # extrapolated energy is unusable, and for F12 methods (which fix a
            # specific F12 basis).
            show_cbs = self._show_cbs() and not is_f12
            if show_cbs:
                add_full("basis set extrapolation")
            if show_cbs and self["basis set extrapolation"].get() != "none":
                add_full("extrapolation family")
            else:
                add_full("basis")
                add_full("basis source")

        for key in self._run_detail_keys():
            add_full(key)

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

    def _show_cbs(self):
        """Whether to show the CBS basis-set-extrapolation controls. Off for
        sub-steps that need a gradient (Optimization overrides this): an
        extrapolated energy has no gradient, so it cannot drive them."""
        return True

    def _run_detail_keys(self):
        """The full-width 'run detail' controls laid out below the level of
        theory, in order. Sub-steps override this to add or drop controls."""
        return (
            "auxiliary basis",
            "grid",
            "scf convergence",
            "sthresh",
            "extra keywords",
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

    def _filter_basis_sets(self):
        """Restrict the basis dropdown to the choices valid for the current
        method, so the user can only pick a usable basis.

        Explicitly-correlated F12 methods work only with the F12-optimized
        orbital bases (they need a matching CABS), so when one is chosen the list
        is narrowed to those and the source is forced to ORCA-internal (the Basis
        Set Exchange has no CABS). Otherwise the full curated list is offered.
        """
        all_bases = list(orca_step.metadata["basis sets"])
        if "F12" in self["method"].get().upper():
            f12 = [b for b in all_bases if b.upper().endswith("-F12")]
            self["basis"].config(values=f12)
            if self["basis"].get() not in f12:
                self["basis"].set("cc-pVTZ-F12" if "cc-pVTZ-F12" in f12 else f12[0])
            self["basis source"].set("ORCA internal")
        else:
            self["basis"].config(values=all_bases)

    def _filter_functionals(self):
        """Restrict the functional pulldown to the functionals of the currently
        selected functional type, keeping the selection valid."""
        ftype = self["functional type"].get()
        funcs = [
            name
            for name, rec in orca_step.metadata["functionals"].items()
            if rec["category"] == ftype
        ]
        self["functional"].combobox.configure(values=funcs)
        if funcs and self["functional"].get() not in funcs:
            self["functional"].set(funcs[0])

    def reset_functionals(self, widget=None):
        """Re-filter the functional list when the functional type changes."""
        self._filter_functionals()
