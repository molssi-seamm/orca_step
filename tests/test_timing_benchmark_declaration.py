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
    # The basis-set ladder: SVP to QZVPPD for a hybrid and the double hybrid
    for method in ("B3LYP", "REVDSD-PBEP86-D4_2021"):
        for basis in ("def2-SVP", "def2-TZVP", "def2-TZVPPD", "def2-QZVPPD"):
            assert f"ORCA:DFT@{method}/{basis}" in d["chemistries"]
    assert d["chemistries"]["ORCA:DFT@B3LYP/def2-QZVPPD"]["tasks"] == ["Energy"]


def test_spec_builds():
    tb = pytest.importorskip("seamm_exec.timing_benchmark")
    if not hasattr(tb, "declarations"):
        pytest.skip("seamm_exec without benchmark discovery")
    assert tb.declarations(refresh=True).get("orca") is orca_step.TIMING_BENCHMARK
    text = tb.build_spec(("orca",), "quick")
    # TZVPPD: water..toluene x (energy, gradient) + 2 optimizations = 8; the
    # ladder's SVP, TZVP and QZVPPD: 3 x 3 molecules x (energy, gradient) = 18
    # (seamm-exec before 2026.10.7.2 also optimizes at the ladder's levels)
    assert text.count("REVDSD-PBEP86-D4_2021") >= 26
    assert text.count("def2-QZVPPD") >= 12
    assert "Optimization" in text


def test_records_carry_neighbours_and_the_whole_basis():
    """The timing record gets the structure's neighbour count, and a functional
    whose ORCA keyword has a "/" (REVDSD-PBEP86-D4/2021) no longer spills into
    the basis."""
    import seamm_exec

    from orca_step.orca_base import _neighbours, timing_descriptors

    if not hasattr(seamm_exec, "neighbour_count"):
        pytest.skip("seamm_exec without neighbour_count")
    g = seamm_exec.Geometry([6] * 40, [(1.5 * i, 0.0, 0.0) for i in range(40)])
    assert 9.0 < _neighbours(g)["neighbours"] < 10.0
    assert _neighbours(g, atom_indices=[0, 1])["neighbours"] == 1.0
    d = timing_descriptors(
        "! REVDSD-PBEP86-D4/2021 def2-TZVPPD EnGrad",
        None,
        model="DFT@REVDSD-PBEP86-D4/2021/def2-TZVPPD",
    )
    assert d["basis"] == "def2-TZVPPD" and d["method"] == "REVDSD-PBEP86-D4/2021"
    assert "neighbours" in orca_step.TIMING_SPEC["size"]
