# -*- coding: utf-8 -*-

"""ORCA's batch contract: the batch path gives the MDI path's numbers.

Runs real ORCA (skipped where ORCA or pymdi is missing) through
``seamm_exec.Evaluator`` both ways for two water geometries.

Stated tolerances: the first geometry, where both paths start the SCF from
scratch, agrees to 1e-6 kJ/mol and 1e-6 kJ/mol/Å (the same input, the same
numbers). Later geometries differ by the SCF convergence: the MDI engine starts
each SCF from the previous geometry's orbitals, a batch task from scratch, and
ORCA's default SCF converges the energy to 1e-6 Eh -- so 3e-3 kJ/mol for
energies and 1e-2 kJ/mol/Å for gradients. (Observed: 2.5e-4 kJ/mol.)
"""

import logging
import shutil
import types

import numpy as np
import pytest

import seamm_exec
from seamm_exec import Evaluator, Geometry

from orca_step import ORCAStep

pytest.importorskip("mdi")
pytest.importorskip("seamm_mdi")

EXACT = 1e-6  # kJ/mol, kJ/mol/Å: the first geometry
E_TOL = 3e-3  # kJ/mol, ORCA's default SCF energy convergence (1e-6 Eh)
G_TOL = 1e-2  # kJ/mol/Å

WATERS = [
    [[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
    [[0.0, 0.0, 0.120], [0.0, 0.770, -0.470], [0.0, -0.750, -0.460]],
]


def _root():
    import pathlib

    root = pathlib.Path("~/SEAMM").expanduser()
    try:
        config = ORCAStep.get_executor_config(seamm_exec.Local(), {"root": str(root)})
    except Exception:
        return None
    code = config.get("code")
    return root if code and shutil.which(code) else None


def _node(directory, root):
    plugin_manager = types.SimpleNamespace(get=lambda name: ORCAStep)
    return types.SimpleNamespace(
        directory=str(directory),
        global_options={"root": str(root)},
        logger=logging.getLogger("test"),
        variable_exists=lambda name: False,
        flowchart=types.SimpleNamespace(
            executor=seamm_exec.Local(),
            plugin_manager=plugin_manager,
            root_directory=str(directory),
        ),
    )


def _mc(basis="def2-SVP"):
    return {
        "level": f"ORCA:DFT@B3LYP/{basis}",
        "owner": "ORCA",
        "type": "DFT",
        "method": "B3LYP",
        "basis": basis,
        "cutoff": None,
        "step": "ORCA",
        "options": {
            "mdi_capable": True,
            "mdi_method_arg": "B3LYP",
            "mdi_basis_arg": basis,
            "prefers_batch": True,
        },
    }


def _evaluate(tmp_path, root, mc, path):
    node = _node(tmp_path / path, root)
    with Evaluator(node, mc, path=path) as evaluator:
        for i, xyz in enumerate(WATERS):
            evaluator.submit(Geometry([8, 1, 1], xyz), key=f"w{i}")
        return {r.key: r for r in evaluator.results()}


@pytest.mark.parametrize("basis", ["def2-SVP", "bse:def2-SVP"])
def test_batch_equals_mdi(tmp_path, basis):
    root = _root()
    if root is None:
        pytest.skip("ORCA is not configured")
    mc = _mc(basis)
    mdi = _evaluate(tmp_path, root, mc, "mdi")
    batch = _evaluate(tmp_path, root, mc, "batch")
    for key, e_tol, g_tol in (("w0", EXACT, EXACT), ("w1", E_TOL, G_TOL)):
        assert mdi[key].ok and batch[key].ok, batch[key].reason
        assert batch[key].path == "batch" and mdi[key].path == "mdi"
        assert batch[key].energy == pytest.approx(mdi[key].energy, abs=e_tol)
        assert np.allclose(batch[key].gradients, mdi[key].gradients, atol=g_tol)
    # The batch tasks restart: a second run restores them
    again = _evaluate(tmp_path, root, mc, "batch")
    assert all(r.restored for r in again.values())
    assert again["w0"].energy == batch["w0"].energy


def test_ghost_atoms_and_fragment_charge(tmp_path):
    """A counterpoise sub-job: a fragment in the full basis (ghosts), with the
    fragment's own charge and multiplicity."""
    from orca_step.batch import get_task

    task = get_task(
        Geometry(
            [8, 1, 1, 8, 1, 1], WATERS[0] + [[2.9, 0, 0], [3.5, 0.8, 0], [3.5, -0.8, 0]]
        ),
        _mc(),
        key="cp",
        options={"atom_indices": [0, 1, 2, 3, 4, 5], "ghost_atoms": [3, 4, 5]},
    )
    text = task.files["orca.inp"]
    assert "O: " in text and text.count("H: ") == 2
    assert "* xyz 0 1" in text
    task = get_task(
        Geometry([3], [[0, 0, 0]], charge=1),
        _mc(),
        key="li",
        options={"charge": 1, "multiplicity": 1},
    )
    assert "* xyz 1 1" in task.files["orca.inp"]
    assert task.program == "orca" and task.config is None
    assert task.cmd[0] == "{code}"


def test_analyze_task_refuses_partial_results(tmp_path):
    from orca_step.batch import analyze_task
    from seamm_exec import AnalysisError, TaskResult

    result = TaskResult(key="x", state="finished", files={"orca.out": "no energy"})
    with pytest.raises(AnalysisError, match="has no energy, gradients"):
        analyze_task(result, _mc(), Geometry([1], [[0, 0, 0]]))


def test_periodic_open_shell_dlpno_and_guess_are_not_tasks():
    from orca_step.batch import can_run_task, get_task

    box = Geometry([8, 1, 1], WATERS[0], cell=np.eye(3) * 10)
    assert not can_run_task(box, _mc())
    with pytest.raises(ValueError, match="periodic"):
        get_task(box, _mc(), key="box")
    water = Geometry([8, 1, 1], WATERS[0])
    assert can_run_task(water, _mc())
    dlpno = dict(_mc(), method="DLPNO-REVDSD-PBEP86-D4/2021")
    dlpno["options"] = dict(
        dlpno["options"], mdi_method_arg="DLPNO-REVDSD-PBEP86-D4/2021"
    )
    assert not can_run_task(Geometry([8], [[0, 0, 0]], multiplicity=3), dlpno)
    with pytest.raises(ValueError, match="initial guess"):
        get_task(water, _mc(), key="g", options={"guess": "orca.gbw"})


def test_single_center_gets_nocosx_on_both_paths():
    from orca_step.batch import engine_helpers, get_task

    task = get_task(Geometry([11], [[0, 0, 0]], charge=1), _mc(), key="na")
    assert task.files["orca.inp"].splitlines()[0] == (
        "! B3LYP AutoAux NoCOSX def2-SVP EnGrad"
    )
    # A ghost atom is a center: a fragment in the cluster basis keeps RIJCOSX
    task = get_task(
        Geometry([11, 17], [[0, 0, 0], [2.4, 0, 0]]),
        _mc(),
        key="na-in-cluster",
        options={"ghost_atoms": [1], "charge": 1},
    )
    assert "NoCOSX" not in task.files["orca.inp"]
    helpers = engine_helpers()
    assert helpers.single_center_method("B3LYP AutoAux", 1) == "B3LYP AutoAux NoCOSX"
    assert helpers.single_center_method("B3LYP RIJK", 1) == "B3LYP RIJK"
    assert helpers.single_center_method("B3LYP", 2) == "B3LYP"


def test_lone_ion_batch_equals_mdi(tmp_path):
    root = _root()
    if root is None:
        pytest.skip("ORCA is not configured")
    results = {}
    for path in ("mdi", "batch"):
        node = _node(tmp_path / path, root)
        with Evaluator(node, _mc(), path=path) as evaluator:
            evaluator.submit(Geometry([11], [[0, 0, 0]], charge=1), key="na")
            results[path] = list(evaluator.results())[0]
    assert results["batch"].ok and results["mdi"].ok
    assert results["batch"].energy == pytest.approx(results["mdi"].energy, abs=EXACT)
