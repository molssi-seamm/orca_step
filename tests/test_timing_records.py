# -*- coding: utf-8 -*-

"""The timing records ORCA runs write (seamm_exec campaign 2026-10-05)."""

from pathlib import Path

import seamm_exec
from seamm_exec.timing import read_timings

from orca_step.orca_base import method_class, task_kind, timing_descriptors

OUT = Path(__file__).parent / "data" / "run_path" / "orca.out"

INPUT = """! B3LYP def2-SVP D3BJ RIJCOSX TightSCF EnGrad
%pal nprocs 4 end
%maxcore 3000
* xyz 0 1
O 0 0 0
H 0 0 1
H 0 1 0
*
"""


def test_task_kind():
    assert task_kind("B3LYP def2-SVP") == "energy"
    assert task_kind("B3LYP def2-SVP EnGrad") == "gradient"
    assert task_kind("B3LYP def2-SVP NumGrad") == "numgrad"
    assert task_kind("Opt B3LYP def2-SVP") == "opt"
    assert task_kind("TightOpt B3LYP def2-SVP EnGrad") == "opt"
    assert task_kind("Opt Freq B3LYP def2-SVP") == "freq"
    assert task_kind("NumFreq B3LYP def2-SVP") == "numfreq"


def test_method_class_from_name():
    assert method_class("B3LYP") == "global hybrid"
    assert method_class("PBE") == "GGA"
    assert method_class("CAM-B3LYP") == "range-separated hybrid"
    assert method_class("wb97m-d4") == "range-separated hybrid"
    assert method_class("B2PLYP") == "global double-hybrid"
    assert method_class("HF") == "HF"
    assert method_class("MP2") == "MP2"
    assert method_class("RI-MP2") == "MP2"
    assert method_class("CCSD(T)") == "CC"
    assert method_class("DLPNO-CCSD(T)") == "DLPNO-CC"


def test_method_class_from_keywords():
    assert method_class(None, "B3LYP def2-SVP D3BJ") == "global hybrid"
    assert method_class("", "RI-MP2 def2-TZVP def2-TZVP/C") == "MP2"
    assert method_class("", "DLPNO-CCSD(T) cc-pVTZ cc-pVTZ/C") == "DLPNO-CC"
    assert method_class("", "CCSD(T) cc-pVTZ") == "CC"
    assert method_class("", "HF def2-SVP") == "HF"
    assert method_class("", "XTB2") == "semiempirical"
    assert method_class("", "nonsense") == ""


def test_timing_descriptors():
    d = timing_descriptors(
        INPUT,
        OUT.read_text(),
        model="DFT@B3LYP/def2-SVP",
        n_atoms=3,
        n_heavy=1,
        charge=0,
        multiplicity=1,
    )
    assert d["task"] == "gradient"
    assert d["method_class"] == "global hybrid"
    assert d["method"] == "B3LYP" and d["basis"] == "def2-SVP"
    assert d["nprocs"] == 4 and d["maxcore_mb"] == 3000
    assert d["ri"] == "RIJCOSX" and d["dispersion"] == "D3BJ"
    assert d["n_atoms"] == 3 and d["n_heavy"] == 1 and d["n_ghosts"] == 0
    # Parsed from the output
    assert d["n_electrons"] == 10
    assert d["nbf"] == 24
    assert d["hf_type"] == "RHF"
    assert d["scf_runs"] == 1 and d["scf_cycles"] == 11
    assert d["geometry_steps"] == 0
    assert abs(d["code_seconds"] - 1.24) < 1e-9
    assert d["terminated_normally"] is True


def test_descriptors_without_output():
    d = timing_descriptors(INPUT, None, model="HF/def2-SVP")
    assert d["method_class"] == "HF" and d["task"] == "gradient"
    assert "nbf" not in d


def test_record_through_seamm_exec(tmp_path, monkeypatch):
    """The descriptors land in orca.csv after the common columns."""
    monkeypatch.setattr(seamm_exec.timing, "DEFAULT_DIRECTORY", tmp_path)
    task = seamm_exec.Task(
        key="k",
        program="orca",
        cmd=["x"],
        files={"orca.inp": INPUT},
        resources=seamm_exec.Resources(ntasks=4),
        estimated_seconds=2.0,
    )
    result = seamm_exec.TaskResult(
        key="k",
        state="finished",
        returncode=0,
        history=[{"started": 10.0, "finished": 11.5}],
    )
    d = timing_descriptors(INPUT, OUT.read_text(), model="DFT@B3LYP/def2-SVP")
    seamm_exec.record_task_timing(task, result, d)
    (row,) = read_timings("orca", directory=tmp_path)
    assert row["program"] == "orca" and row["wall"] == "1.500"
    assert row["nbf"] == "24" and row["scf_cycles"] == "11"
    assert row["task"] == "gradient" and row["method_class"] == "global hybrid"
    assert row["terminated_normally"] == "1"


def test_estimated_basis_functions():
    from orca_step.orca_base import estimated_basis_functions

    # Water at def2-SVP: 14 + 2 x 5 = 24, as the ORCA output says
    assert estimated_basis_functions("def2-SVP", [8, 1, 1]) == 24
    assert estimated_basis_functions("bse:def2-SVP", [8, 1, 1]) == 24
    # A ghost centre carries its basis functions too
    assert estimated_basis_functions("def2-SVP", [8, 1, 1, 8]) == 38
    # An unknown basis falls back to a double-zeta count
    assert estimated_basis_functions("no-such-basis", [6, 1, 1, 1, 1]) == 14 + 4 * 5
    assert estimated_basis_functions("", [26]) == 32


def test_predicted_seconds_uses_the_model_or_falls_back(monkeypatch):
    from seamm_exec import timing_model

    from orca_step.orca_base import estimated_seconds, predicted_seconds

    line = "B3LYP def2-SVP TightSCF EnGrad"
    seen = {}

    def fake_predict(program, descriptors, ntasks=1, quantile=0.95, **kw):
        seen.update(program=program, descriptors=descriptors, ntasks=ntasks, q=quantile)
        return {"seconds": 42.0}

    monkeypatch.setattr(timing_model, "predict", fake_predict)
    t = predicted_seconds(
        line, [8, 1, 1], model="DFT@B3LYP/def2-SVP", ntasks=4, charge=0, multiplicity=1
    )
    assert t == 42.0
    d = seen["descriptors"]
    assert seen["program"] == "orca" and seen["ntasks"] == 4 and seen["q"] == 0.5
    assert d["task"] == "gradient" and d["method_class"] == "global hybrid"
    assert d["method"] == "B3LYP" and d["basis"] == "def2-SVP"
    assert d["nbf"] == 24 and d["n_electrons"] == 10 and d["n_atoms"] == 3

    # No model: the hand estimate
    monkeypatch.setattr(timing_model, "predict", lambda *a, **k: None)
    assert predicted_seconds(line, [8, 1, 1], model="DFT@B3LYP/def2-SVP") == (
        estimated_seconds(line, 3)
    )

    # A failure inside the model never stops the step
    def boom(*a, **k):
        raise RuntimeError("no")

    monkeypatch.setattr(timing_model, "predict", boom)
    assert predicted_seconds(line, [8, 1, 1]) == estimated_seconds(line, 3)


def test_timing_spec_is_passed_when_recording(monkeypatch):
    from orca_step import orca_base

    assert orca_base.TIMING_SPEC["size"][0] == "nbf"
    assert orca_base.TIMING_SPEC["units"] == "scf_runs"
    assert orca_base._record_kwargs() == {"spec": orca_base.TIMING_SPEC}
    monkeypatch.delattr(seamm_exec, "TimingSpec", raising=False)
    assert orca_base._record_kwargs() == {}  # an older seamm-exec: no spec
