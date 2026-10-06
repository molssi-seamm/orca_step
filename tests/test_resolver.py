# -*- coding: utf-8 -*-

"""The 'orca' resolver reproduces what the ORCA step did in the evaluator.

Before phase 3 the ORCA step resolved its configuration itself (``_orca_config``,
``_mpi_env``) and passed it on the task; now the task names only the program and
the local pool resolves it with ``orca_step.resolver``. The command, the
environment and the configuration must be exactly the same.
"""

import os
import types

import pytest

import seamm_exec
from seamm_exec.local_pool import LocalPool
from seamm_exec.resolve import register

from orca_step import orca_base
from orca_step.orca_base import _orca_2aim, mpi_env
from orca_step.orca_step import full_orca_path
from orca_step.resolver import resolve


def _old_way(root, n_cores, make_wfx):
    """What run_orca_job built before phase 3."""
    me = types.SimpleNamespace(
        flowchart=types.SimpleNamespace(executor=seamm_exec.Local()),
        parent=types.SimpleNamespace(global_options={"root": str(root)}),
    )
    config = orca_base.ORCABase._orca_config(me)
    env, lib_prefix = mpi_env(n_cores, config)
    cmd = lib_prefix + ["{code}", "orca.inp", ">", "orca.out", "2>", "orca.err"]
    if make_wfx:
        cmd += ["&&", _orca_2aim(config), "orca", ">", "orca_2aim.out", "2>&1"]
    return config, cmd, env


def _new_way(root, n_cores, make_wfx):
    """What the local pool now runs for the bare task."""
    register("orca", resolve)
    cmd = ["{code}", "orca.inp", ">", "orca.out", "2>", "orca.err"]
    if make_wfx:
        cmd += ["&&", "{orca_2aim}", "orca", ">", "orca_2aim.out", "2>&1"]
    task = seamm_exec.Task(key="orca", program="orca", cmd=cmd, shell=True)
    pool = LocalPool(seamm_exec.Local(), root=root, ce={"NTASKS": 8})
    config, cmd, env = pool._configure(task, {"NTASKS": n_cores})
    return config, cmd, env


def _format(config, cmd):
    """The command as the executor formats it."""
    text = " ".join(cmd)
    while True:
        new = text.format(**config)
        if new == text:
            return text
        text = new


@pytest.mark.parametrize("n_cores", [1, 4])
@pytest.mark.parametrize("make_wfx", [False, True])
def test_resolver_gives_the_old_command(tmp_path, n_cores, make_wfx, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    orca = tmp_path / "orca_6" / "orca"
    orca.parent.mkdir()
    orca.write_text("#!/bin/sh\n")
    orca.chmod(0o755)
    lib = tmp_path / "openmpi" / "lib"
    (tmp_path / "openmpi" / "bin").mkdir(parents=True)
    lib.mkdir()
    # The OpenMPI 4.1 beside library-path, as a real installation has
    mpirun = tmp_path / "openmpi" / "bin" / "mpirun"
    mpirun.write_text("#!/bin/sh\necho 'mpirun (Open MPI) 4.1.6'\n")
    mpirun.chmod(0o755)
    (tmp_path / "orca.ini").write_text(
        f"[local]\ninstallation = local\ncode = {orca}\nlibrary-path = {lib}\n"
    )
    old_config, old_cmd, old_env = _old_way(tmp_path, n_cores, make_wfx)
    new_config, new_cmd, new_env = _new_way(tmp_path, n_cores, make_wfx)
    assert _format(new_config, new_cmd) == _format(
        {**old_config, "code_dir": str(orca.parent)}, old_cmd
    )
    assert new_env == old_env
    for key, value in old_config.items():
        assert new_config[key] == value


def test_resolver_falls_back_to_the_path(tmp_path, monkeypatch):
    """No orca.ini: ORCA from the PATH, as _orca_config did."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    orca = bindir / "orca"
    orca.write_text("#!/bin/sh\n")
    orca.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    config, cmd, env = resolve({}, ["{code}", "orca.inp"], {}, {"NTASKS": 1}, tmp_path)
    assert config["code"] == str(orca) == full_orca_path("")
    assert config["installation"] == "local"
    assert resolve.available(tmp_path)


def test_resolver_without_orca_says_so(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(RuntimeError, match="Could not find the 'orca' executable"):
        resolve({}, ["{code}"], {}, {"NTASKS": 1}, tmp_path)
