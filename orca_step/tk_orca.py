# -*- coding: utf-8 -*-

"""The graphical part of the main ORCA node: a sub-flowchart editor."""

import pprint  # noqa: F401

import orca_step  # noqa: F401
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401


class TkORCA(seamm.TkNode):
    """Graphical ORCA node. Opens a sub-flowchart canvas where ORCA capabilities
    (Energy, Optimization, ...) are added, mirroring the Gaussian/MOPAC steps.
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        namespace="org.molssi.seamm.orca.tk",
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        self.namespace = namespace
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
        self.create_dialog()

    def create_dialog(self):
        """Create the dialog holding the sub-flowchart editor."""
        frame = super().create_dialog(title="ORCA")

        screen_w = self.dialog.winfo_screenwidth()
        screen_h = self.dialog.winfo_screenheight()
        w = int(0.9 * screen_w)
        h = int(0.8 * screen_h)
        x = int(0.05 * screen_w / 2)
        y = int(0.1 * screen_h / 2)
        self.dialog.geometry(f"{w}x{h}+{x}+{y}")

        self.tk_subflowchart = seamm.TkFlowchart(
            master=frame, flowchart=self.node.subflowchart, namespace=self.namespace
        )
        self.tk_subflowchart.draw()

    def right_click(self, event):
        """Add an 'Edit...' entry to the right-click menu."""
        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)
        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
