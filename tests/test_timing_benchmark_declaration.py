# -*- coding: utf-8 -*-
"""The step's timing benchmark declaration (seamm_exec.timing_benchmark)."""

import pytest

import orca_step


def test_declaration_shape():
    d = orca_step.TIMING_BENCHMARK
    assert d["program"] == "orca" and d["step"] == "ORCA" and d["parallel"]
    sizes = [s["size"] for s in d["systems"]]
    assert sizes == sorted(sizes) and sizes[0] == 3 and sizes[-1] == 302
    assert all(s["steps"][0]["FromSMILESStep"]["smiles string"] for s in d["systems"])
    assert "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPPD" in d["chemistries"]
    for limits in (*d["chemistries"].values(), *d["tasks"].values()):
        assert limits["quick"] <= limits["full"]
    assert {"results": {"gradients": {}}} in d["variants"]["Energy"]


def test_spec_builds():
    tb = pytest.importorskip("seamm_exec.timing_benchmark")
    if not hasattr(tb, "declarations"):
        pytest.skip("seamm_exec without benchmark discovery")
    assert tb.declarations(refresh=True).get("orca") is orca_step.TIMING_BENCHMARK
    text = tb.build_spec(("orca",), "quick")
    assert (
        text.count("REVDSD-PBEP86-D4_2021") == 8
    )  # water..toluene x (energy, gradient) + 2 opts
    assert "Optimization" in text
