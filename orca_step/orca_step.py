# -*- coding: utf-8 -*-

"""Stevedore helper for the main ORCA node."""

import configparser
import importlib.resources
from pathlib import Path
import shutil
import sys

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

        Advertises ``ORCA:<type>@<method>/<basis>`` for each method in the
        metadata paired with a curated set of basis sets. ORCA is molecular, so
        ``periodic_only`` returns nothing. ``mdi_only`` keeps only the methods
        drivable through the ``orca_mdi.py`` engine -- those with an analytic
        gradient (the engine always requests EnGrad); each such option carries
        the real ORCA keyword and basis in ``mdi_method_arg``/``mdi_basis_arg``.
        """
        # ORCA is molecular here, so nothing for a periodic-only request.
        if periodic_only:
            return {}

        options = {}
        for method, info in orca_step.metadata["methods"].items():
            mtype = info.get("type", "QC")
            # For DFT the specific functional is the "method"; advertise each one
            # (aliasing any '/' in the keyword) with its gradient availability.
            # Other methods are literal.
            if method == "DFT":
                entries = [
                    (
                        mc_method_alias(name),
                        name,
                        rec.get("gradients", "analytic"),
                    )
                    for name, rec in orca_step.metadata["functionals"].items()
                ]
            else:
                entries = [(method, method, info.get("gradients", "analytic"))]
            for adv_method, real, gradients in entries:
                # The MDI engine (orca_mdi.py) always requests EnGrad, so a
                # method is MDI-capable only if ORCA has an analytic gradient for
                # it -- excludes DLPNO-CCSD(T)/CCSD(T) and the non-self-consistent
                # wB97M(2)/wB97X-2. (Those remain available through the ordinary
                # file-based ORCA step, just not as a persistent MDI engine.)
                mdi_capable = gradients == "analytic"
                if mdi_only and not mdi_capable:
                    continue
                for basis in _ADVERTISED_BASES:
                    key = f"ORCA:{mtype}@{adv_method}/{basis}"
                    options[key] = {
                        "model_chemistry": key,
                        "type": mtype,
                        "method": adv_method,
                        "basis": basis,
                        "description": f"{real} / {basis}",
                        "periodic": False,
                        "mdi_capable": mdi_capable,
                        # The real (un-aliased) ORCA keyword + basis the engine
                        # needs, or None when not MDI-capable.
                        "mdi_method_arg": real if mdi_capable else None,
                        "mdi_basis_arg": basis if mdi_capable else None,
                    }
        return options

    @classmethod
    def get_executor_config(cls, executor, seamm_options):
        """How to launch ORCA (and its MDI engine) on this machine.

        Reads ``<root>/orca.ini`` for the current executor to find the ``orca``
        binary (falling back to the PATH), and adds ``mdi_script`` -- the
        absolute path to the bundled ``data/orca_mdi.py`` engine. That script
        imports only packages present in the SEAMM environment (pymdi, numpy,
        seamm_util) and runs the ``orca`` binary itself, so no separate conda
        environment is needed.
        """
        executor_type = executor.name
        ini_path = Path(seamm_options["root"]).expanduser() / "orca.ini"
        resources = importlib.resources.files("orca_step") / "data"

        full_config = configparser.ConfigParser()
        if ini_path.exists():
            full_config.read(ini_path)

        config = (
            dict(full_config.items(executor_type))
            if executor_type in full_config
            else {}
        )
        code = config.get("code", "") or ""
        if code == "":
            code = shutil.which("orca") or ""
        if code == "":
            raise RuntimeError(
                "Could not find the 'orca' executable for the MDI engine. Set "
                "'code' in the relevant section of orca.ini, or put orca on the "
                "PATH."
            )
        config["code"] = code
        config["version"] = orca_step.__version__
        config["mdi_script"] = str(resources / "orca_mdi.py")
        return config

    @classmethod
    def get_mdi_engine_command(
        cls,
        executor,
        seamm_options,
        *,
        method,
        basis="def2-SVP",
        port,
        hostname="localhost",
        charge=0,
        multiplicity=1,
        n_atoms=None,
        ncores=1,
        engine_name="ORCA",
        extra_args=None,
    ):
        """Build the argv that launches the ORCA MDI *engine* over TCP.

        The transport (TCP, ``port``, ``hostname``) is decided by the driver and
        passed in; everything ORCA-specific -- the bundled ``orca_mdi.py``, the
        orca binary, and the method/basis/charge/multiplicity flags -- is
        supplied here so the driver hardwires no ORCA knowledge. ``method``
        should be the real ORCA keyword (the ``mdi_method_arg`` from
        :meth:`get_model_chemistry_options`, not an aliased functional name).
        """
        config = cls.get_executor_config(executor, seamm_options)
        mdi_init = (
            f"-role ENGINE -name {engine_name} -method TCP "
            f"-port {port} -hostname {hostname}"
        )
        argv = [
            sys.executable,
            config["mdi_script"],
            "-mdi",
            mdi_init,
            "--orca",
            config["code"],
            "--method",
            method,
            "--basis",
            basis,
            "--charge",
            str(charge),
            "--multiplicity",
            str(multiplicity),
            "--ncores",
            str(ncores),
        ]
        if extra_args:
            argv.extend(extra_args)
        return argv

    def description(self):
        """Return a description of what this step does."""
        return ORCAStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        """Return a new ORCA node."""
        return orca_step.ORCA(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        """Return a new graphical ORCA node."""
        return orca_step.TkORCA(canvas=canvas, **kwargs)
