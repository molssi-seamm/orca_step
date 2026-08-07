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
    The Energy property toggles (bond orders, Hirshfeld, polarizability) and
    SThresh are not shown -- they are not plumbed through the per-fragment
    BSSE job set. 'Write the wavefunction (wfx) file' (for a following Atomic
    Charges step) IS shown -- see ``_run_detail_keys`` below.
    """

    def _show_cbs(self):
        """Hide CBS extrapolation: the counterpoise gradient needs a real
        gradient, which an extrapolated energy does not have."""
        return False

    def create_dialog(self, title="ORCA BSSE"):
        """Build the dialog and make the fragment controls reactive: the
        'Fragment atoms' field is shown only when the fragments are
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
        # 'Fragment atoms' only applies when defining the fragments by hand;
        # 'Fragment charges' applies either way (auto-detected fragments
        # still need their charges, e.g. a Na+/Cl- pair).
        if self["fragments"].get() == "specified":
            keys.append("fragment atoms")
        keys.append("fragment charges")
        keys.append("compute gradient")
        keys.append("optimize monomers")
        # For a following Atomic Charges (DDEC6) step: write the cluster's .wfx.
        keys.append("save wavefunction")
        return tuple(keys)
