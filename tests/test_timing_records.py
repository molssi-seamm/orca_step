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
