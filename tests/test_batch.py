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


def test_exact_exchange_on_both_paths():
    from orca_step.batch import engine_helpers, get_task

    task = get_task(Geometry([11], [[0, 0, 0]], charge=1), _mc(), key="na")
    assert task.files["orca.inp"].splitlines()[0] == (
        "! B3LYP AutoAux NoCOSX def2-SVP TIGHTSCF DEFGRID3 EnGrad"
    )
    # Na in a molecule, or as a ghost, also gets exact exchange (orca_step#44)
    nacl = Geometry([11, 17], [[0, 0, 0], [2.4, 0, 0]])
    task = get_task(nacl, _mc(), key="nacl")
    assert "NoCOSX" in task.files["orca.inp"].splitlines()[0]
    task = get_task(
        nacl, _mc(), key="cl-in-cluster", options={"ghost_atoms": [0], "charge": -1}
    )
    assert "NoCOSX" in task.files["orca.inp"].splitlines()[0]
    # ... unless the element set says otherwise, from the task or the model
    # chemistry
    task = get_task(nacl, _mc(), key="n", options={"exact_exchange_elements": "none"})
    assert "NoCOSX" not in task.files["orca.inp"]
    mc = _mc()
    mc["options"] = {**(mc.get("options") or {}), "exact_exchange_elements": "none"}
    assert "NoCOSX" not in get_task(nacl, mc, key="m").files["orca.inp"]
    # A ghost atom of another element is a center: Cl- in the cluster basis of
    # water keeps RIJCOSX
    task = get_task(
        Geometry(
            [17, 8, 1, 1], [[0, 0, 0], [3.1, 0, 0], [3.7, 0.76, 0], [3.7, -0.76, 0]]
        ),
        _mc(),
        key="cl-in-cluster",
        options={"ghost_atoms": [1, 2, 3], "charge": -1},
    )
    assert "NoCOSX" not in task.files["orca.inp"]
    helpers = engine_helpers()
    method = helpers.exact_exchange_method
    assert method("B3LYP AutoAux", ["Na"]) == "B3LYP AutoAux NoCOSX"
    assert method("B3LYP RIJK", ["Na"]) == "B3LYP RIJK"
    assert method("B3LYP", ["O", "H"]) == "B3LYP"
    assert method("B3LYP", ["Na", "O"]) == "B3LYP NoCOSX"
    assert method("B3LYP", ["Na", "O"], "none") == "B3LYP"
    assert method("B3LYP", ["Li", "O"]) == "B3LYP"
    assert method("B3LYP", ["Li", "O"], "Li") == "B3LYP NoCOSX"
    assert helpers.exact_exchange_elements(None) == helpers.EXACT_EXCHANGE_ELEMENTS
    assert helpers.exact_exchange_elements("li, NA") == {"Li", "Na"}
    assert helpers.exact_exchange_elements("") == frozenset()


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


def test_thermochemistry_reference_uses_the_dlpno_parent():
    """The MBE step's DfE0 offsets look up the atomic references through this
    hook: a DLPNO double hybrid uses its canonical parent's atoms, as the
    Energy sub-step's DfE0 does; other methods are unchanged."""
    from orca_step import ORCAStep

    def mc(method, basis="def2-QZVPPD"):
        return {**_mc(basis), "method": method}

    assert ORCAStep.thermochemistry_reference(mc("DLPNO-REVDSD-PBEP86-D4_2021")) == (
        "orca",
        "REVDSD-PBEP86-D4_2021",
        "def2-QZVPPD",
    )
    assert ORCAStep.thermochemistry_reference(mc("REVDSD-PBEP86-D4_2021")) == (
        "orca",
        "REVDSD-PBEP86-D4_2021",
        "def2-QZVPPD",
    )
    assert ORCAStep.thermochemistry_reference(mc("B3LYP", "bse:def2-TZVP")) == (
        "orca",
        "B3LYP",
        "def2-TZVP",
    )


def _engrad(gradient, energy=-76.3):
    """An orca.engrad file for a gradient given in E_h/bohr."""
    lines = ["#", "# Number of atoms", "#", f" {len(gradient)}", "#"]
    lines += ["# The current total energy in Eh", "#", f" {energy}", "#"]
    lines += ["# The current gradient in Eh/bohr", "#"]
    lines += [f" {g:.9f}" for row in gradient for g in row]
    return "\n".join(lines) + "\n"


# A sound water gradient (E_h/bohr); _shifted adds dx E_h/bohr (51.4 meV/Å per
# 1e-3) to one atom's x.
_SOUND = [[0.0, 0.0, 0.02], [0.0, 0.011, -0.01], [0.0, -0.011, -0.01]]


def _shifted(dx):
    return [[_SOUND[0][0] + dx, *_SOUND[0][1:]], _SOUND[1], _SOUND[2]]


def test_net_force_limits_depend_on_cosx():
    """Without RIJCOSX sound double-hybrid gradients stay within 2.6 meV/Å, so
    the limit is 5; RIJCOSX's grid leaves up to 18 meV/Å, so there it is 50
    (FEC pilot and NaCl cells, 2026-10-09)."""
    from orca_step.orca_base import engrad_gradient, net_force, net_force_problem

    cosx = "...\nCOSX GRID GENERATION\n..."
    assert engrad_gradient(_engrad(_SOUND)) == _SOUND
    assert net_force(_SOUND) < 1e-6
    assert net_force_problem(_SOUND, "") is None
    assert abs(net_force(_shifted(2e-4)) - 10.28) < 0.01
    # 2.6 meV/Å: sound either way
    assert net_force_problem(_shifted(5e-5), "") is None
    # 10.3 meV/Å: broken without RIJCOSX, grid noise with it
    assert "10.3 meV/Å" in net_force_problem(_shifted(2e-4), "")
    assert "without RIJCOSX" in net_force_problem(_shifted(2e-4), "")
    assert net_force_problem(_shifted(2e-4), cosx) is None
    assert "with RIJCOSX" in net_force_problem(_shifted(1.2e-3), cosx)
    assert net_force_problem(None, "") is None
    assert engrad_gradient("garbage") is None


def test_analyze_task_refuses_a_gradient_that_is_not_translation_invariant():
    """The batch path fails such a task rather than return its forces; the
    check sums every centre in orca.engrad, ghosts included."""
    from orca_step.batch import analyze_task
    from seamm_exec import AnalysisError, TaskResult

    out = "FINAL SINGLE POINT ENERGY       -76.300000000000\n"
    geometry = Geometry(
        [8, 1, 1], [[0, 0, 0.117], [0, 0.757, -0.467], [0, -0.757, -0.467]]
    )

    def result(gradient):
        return TaskResult(
            key="w",
            state="finished",
            files={"orca.out": out, "orca.engrad": _engrad(gradient)},
        )

    data = analyze_task(result(_SOUND), _mc(), geometry)
    assert data["gradients"].shape == (3, 3)
    assert analyze_task(result(_shifted(5e-5)), _mc(), geometry)  # 2.6 meV/Å
    with pytest.raises(AnalysisError, match="forces do not sum to zero"):
        analyze_task(result(_shifted(2e-4)), _mc(), geometry)
    # A ghost centre beyond the atoms counts towards the sum
    ghosted = _SOUND + [[0.0, 0.0, 1e-3]]
    with pytest.raises(AnalysisError, match="net force is 51.4"):
        analyze_task(result(ghosted), _mc(), geometry)


def test_energy_warns_about_a_gradient_that_is_not_translation_invariant(tmp_path):
    """The Energy sub-step's run path keeps the forces but warns loudly."""
    from types import SimpleNamespace

    from orca_step.orca_base import ORCABase

    node = SimpleNamespace(indent="")
    (tmp_path / "orca.out").write_text("FINAL SINGLE POINT ENERGY -76.3\n")
    (tmp_path / "orca.engrad").write_text(_engrad(_SOUND))
    assert ORCABase._check_net_force(node, tmp_path) is None
    (tmp_path / "orca.engrad").write_text(_engrad(_shifted(2e-4)))
    assert "10.3 meV/Å" in ORCABase._check_net_force(node, tmp_path)
    (tmp_path / "orca.out").write_text("COSX GRID GENERATION\n")
    assert ORCABase._check_net_force(node, tmp_path) is None
    (tmp_path / "orca.engrad").unlink()
    assert ORCABase._check_net_force(node, tmp_path) is None


def test_double_hybrid_gradients_run_with_verytightscf():
    """With TightSCF ORCA 6.1.1 gets some double-hybrid gradients wrong, so every
    double-hybrid gradient runs with VeryTightSCF -- on the batch path, the MDI
    engine and the steps' own runs -- and other methods keep TightSCF."""
    from orca_step.batch import engine_helpers, get_task
    from orca_step.orca_base import double_hybrid_scf
    from orca_step.orca_step import _engine_method_words

    water = Geometry(
        [8, 1, 1], [[0, 0, 0.117], [0, 0.757, -0.467], [0, -0.757, -0.467]]
    )
    mc = {
        **_mc(),
        "method": "REVDSD-PBEP86-D4_2021",
        "options": {
            **_mc()["options"],
            "mdi_method_arg": "REVDSD-PBEP86-D4/2021",
        },
    }
    line = get_task(water, mc, key="w").files["orca.inp"].splitlines()[0]
    assert "VERYTIGHTSCF" in line.split() and "TIGHTSCF" not in line.split()
    # An energy alone keeps TightSCF; so does a hybrid's gradient
    line = get_task(water, mc, key="e", properties=("energy",)).files["orca.inp"]
    assert "TIGHTSCF" in line.splitlines()[0].split()
    line = get_task(water, _mc(), key="b").files["orca.inp"].splitlines()[0]
    assert "TIGHTSCF" in line.split() and "VERYTIGHTSCF" not in line.split()

    # The MDI engine always computes the gradient
    assert _engine_method_words("REVDSD-PBEP86-D4/2021") == (
        "REVDSD-PBEP86-D4/2021 VERYTIGHTSCF"
    )
    assert _engine_method_words("DLPNO-REVDSD-PBEP86-D4/2021") == (
        "REVDSD-PBEP86-D4/2021 VERYTIGHTSCF"
    )
    assert _engine_method_words("B3LYP") == "B3LYP"
    text = engine_helpers().orca_input(
        "REVDSD-PBEP86-D4/2021 VERYTIGHTSCF AutoAux",
        "def2-SVP",
        0,
        1,
        ["O", "H", "H"],
        [[0, 0, 0.117], [0, 0.757, -0.467], [0, -0.757, -0.467]],
    )
    words = text.splitlines()[0].split()
    assert "VERYTIGHTSCF" in words and "TIGHTSCF" not in words

    # The steps' own runs: gradients, optimizations and numerical frequencies
    for line in (
        "REVDSD-PBEP86-D4/2021 def2-SVP TIGHTSCF EnGrad",
        "REVDSD-PBEP86-D4/2021 def2-SVP Opt",
        "B2PLYP def2-SVP NumFreq",
    ):
        new, note = double_hybrid_scf(line)
        assert "VERYTIGHTSCF" in new.split() and "TIGHTSCF" not in new.split()
        assert "double-hybrid gradient runs with VERYTIGHTSCF" in note
    for line in (
        "REVDSD-PBEP86-D4/2021 def2-SVP TIGHTSCF",  # an energy
        "B3LYP def2-SVP TIGHTSCF EnGrad",  # not a double hybrid
        "REVDSD-PBEP86-D4/2021 def2-SVP EXTREMESCF EnGrad",  # already tighter
    ):
        assert double_hybrid_scf(line) == (line, None)


def test_timing_flags_are_the_options_of_the_line():
    """The cost model learns a factor per option of the '!' line; the task,
    basis sets and method are not options (they are in the model already)."""
    timing_model = pytest.importorskip("seamm_exec.timing_model")
    if not hasattr(timing_model, "flag_tokens"):
        pytest.skip("seamm-exec without flags")
    from orca_step.orca_base import TIMING_SPEC

    flags = TIMING_SPEC["flags"]
    record = {
        "keywords": "REVDSD-PBEP86-D4/2021 AutoAux def2-TZVPPD DEFGRID3 NoCOSX "
        "EnGrad VERYTIGHTSCF",
        "method": "REVDSD-PBEP86-D4_2021",
        "basis": "def2-TZVPPD",
    }
    assert timing_model.flag_tokens(record, flags) == {
        "DEFGRID3",
        "NOCOSX",
        "VERYTIGHTSCF",
    }
    record = {
        "keywords": "B3LYP def2-QZVPPD def2/J RIJCOSX TightOpt keepdensity",
        "method": "B3LYP",
        "basis": "def2-QZVPPD",
    }
    assert timing_model.flag_tokens(record, flags) == {"RIJCOSX", "KEEPDENSITY"}
