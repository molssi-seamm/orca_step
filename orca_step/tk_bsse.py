# -*- coding: utf-8 -*-

"""The graphical part of an ORCA BSSE (counterpoise) sub-step."""

import logging

import orca_step  # noqa: F401
from .tk_energy import TkEnergy
import seamm_widgets as sw  # noqa: F401

logger = logging.getLogger(__name__)


class TkBSSE(TkEnergy):
    """Graphical ORCA BSSE sub-step: the energy dialog's level-of-theory
    controls plus the fragment definition and the monomer-optimization switch.
    The Energy property toggles (bond orders, Hirshfeld, polarizability, saved
    wavefunction) and SThresh are not shown -- they are not plumbed through the
    ORCA Compound path in Phase 1.
    """

    def _show_cbs(self):
        """Hide CBS extrapolation: the counterpoise gradient needs a real
        gradient, which an extrapolated energy does not have."""
        return False

    def create_dialog(self, title="ORCA BSSE"):
        """Build the dialog and make the fragment controls reactive: the
        'Fragment A atoms' field is shown only when the fragments are
        'specified'."""
        frame = super().create_dialog(title=title)
        w = self["fragments"]
        for sequence in ("<<ComboboxSelected>>", "<Return>", "<FocusOut>"):
            w.combobox.bind(sequence, self.reset_dialog)
        self.reset_dialog()
        return frame

    def _run_detail_keys(self):
        keys = [
            "auxiliary basis",
            "grid",
            "scf convergence",
            "extra keywords",
            "fragments",
        ]
        # 'Fragment A atoms' only applies when defining the fragments by hand.
        if self["fragments"].get() == "specified":
            keys.append("fragment A atoms")
        keys.append("optimize monomers")
        return tuple(keys)
