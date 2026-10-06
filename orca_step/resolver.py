# -*- coding: utf-8 -*-

"""Resolve ORCA on the machine that runs it (seamm_exec's resolver hook).

A task names the program ``orca`` and a bare command (``{code} orca.inp > ...``);
where it runs, ``seamm_exec`` reads the ``[local]`` section of that machine's
``<root>/orca.ini`` and calls :func:`resolve` with the task's share of the
machine. This does there what the ORCA step used to do in the evaluator: ORCA
must be invoked by its full path, and a parallel run needs the OpenMPI it was
built against on the paths (see :func:`orca_step.orca_base.mpi_env`).

Registered as the entry point ``orca`` in ``org.molssi.seamm.exec.resolvers``.
"""

import shutil

from .orca_base import _orca_2aim, check_mpirun, mpi_env
from .orca_step import full_orca_path


def resolve(config, cmd, env, ce, root):
    """``(config, cmd, env)`` for running ORCA here.

    Parameters
    ----------
    config : dict
        This machine's ``orca.ini`` section (empty if there is none).
    cmd : [str]
        The task's command template.
    env : dict
        The task's extra environment.
    ce : dict
        The task's computational environment (``NTASKS``, ...).
    root : str or Path
        The SEAMM root holding ``orca.ini``.
    """
    code = (config.get("code") or "").strip()
    if code == "":
        # No orca.ini, or no code in it: ORCA on the PATH, as the ORCA step does.
        found = full_orca_path("")
        if not found:
            raise RuntimeError(
                "Could not find the 'orca' executable. Put it on your PATH, or set "
                f"'code' (the full path) in the [local] section of {root}/orca.ini."
            )
        config = {"code": found, "installation": "local"}
    else:
        config = dict(config)
        # A bare name such as 'orca' must become the full path.
        config["code"] = full_orca_path(code)

    n_cores = max(1, int(ce.get("NTASKS", 1) or 1))
    check_mpirun(n_cores, config)
    extra_env, lib_prefix = mpi_env(n_cores, config)
    env = dict(env)
    for name, value in extra_env.items():
        env.setdefault(name, value)
    # For commands that run ORCA's companion orca_2aim
    config["orca_2aim"] = _orca_2aim(config)
    return config, lib_prefix + list(cmd), env


def _available(root):
    return shutil.which("orca") is not None


resolve.available = _available
