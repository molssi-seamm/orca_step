# -*- coding: utf-8 -*-

"""Stevedore helper for the main ORCA node."""

import orca_step

# Basis sets advertised with each method to the Model Chemistry step. ORCA can
# use any basis it knows; this is a curated, bounded set so the model-chemistry
# list does not explode. A user wanting another basis can set it explicitly in
# the ORCA dialog. (Broadening what the Model Chemistry step offers is a
# model_chemistry_step concern, not ORCA's.)
_ADVERTISED_BASES = ("def2-SVP", "def2-TZVP", "def2-QZVP", "cc-pVTZ")


def mc_method_alias(functional):
    """A model-chemistry-safe spelling of a DFT functional keyword.

    The model-chemistry grammar reserves ``/`` (it separates method from basis),
    so functionals whose ORCA keyword contains ``/`` -- e.g.
    ``REVDSD-PBEP86-D4/2021`` -- cannot appear literally in a model-chemistry
    string. They are advertised with ``/`` replaced by ``_`` (no ORCA functional
    keyword contains ``_``, so this round-trips) and translated back to the real
    keyword when the ORCA step consumes the model chemistry (see
    ``Energy._method_basis_from_model_chemistry``).
    """
    return functional.replace("/", "_")


def mc_method_unalias(method):
    """Inverse of :func:`mc_method_alias`: the real ORCA functional keyword for a
    (possibly aliased) model-chemistry method, or the method unchanged if it is
    not an aliased functional."""
    if method in orca_step.metadata["functionals"]:
        return method
    candidate = method.replace("_", "/")
    if candidate in orca_step.metadata["functionals"]:
        return candidate
    return method


class ORCAStep(object):
    """Helper class for the stevedore integration of the ORCA step."""

    my_description = {
        "description": "An interface for ORCA",
        "group": "Simulations",
        "name": "ORCA",
    }

    def __init__(self, flowchart=None, gui=None):
        pass

    @classmethod
    def get_model_chemistry_options(cls, periodic_only=False, mdi_only=False):
        """Return the model chemistries ORCA can provide.

        ORCA is a molecular quantum-chemistry program: not periodic and not
        (currently) an MDI engine, so it returns nothing when either filter is
        set. Otherwise it advertises ``ORCA:<type>@<method>/<basis>`` for each
        method in the metadata paired with a curated set of basis sets.
        """
        if periodic_only or mdi_only:
            return {}

        options = {}
        for method, info in orca_step.metadata["methods"].items():
            mtype = info.get("type", "QC")
            # For DFT the specific functional is the "method"; advertise each
            # one (aliasing any '/' in the keyword). Other methods are literal.
            if method == "DFT":
                entries = [
                    (mc_method_alias(name), name, rec.get("note") or "DFT functional")
                    for name, rec in orca_step.metadata["functionals"].items()
                ]
            else:
                entries = [(method, method, info.get("description", method))]
            for adv_method, real, desc in entries:
                for basis in _ADVERTISED_BASES:
                    key = f"ORCA:{mtype}@{adv_method}/{basis}"
                    options[key] = {
                        "model_chemistry": key,
                        "type": mtype,
                        "method": adv_method,
                        "basis": basis,
                        "description": f"{real} / {basis}",
                        "periodic": False,
                        "mdi_capable": False,
                    }
        return options

    def description(self):
        """Return a description of what this step does."""
        return ORCAStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        """Return a new ORCA node."""
        return orca_step.ORCA(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        """Return a new graphical ORCA node."""
        return orca_step.TkORCA(canvas=canvas, **kwargs)
