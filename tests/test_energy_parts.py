# -*- coding: utf-8 -*-

"""The parts of ORCA 6.1.1 energies (energy_parts), from real outputs: water with
cc-pVTZ(-F12) at each kind of run, and the OH radical (UHF, QRO)."""

import gzip
from pathlib import Path

import pytest

from orca_step.energy_parts import SUM_TOLERANCE, energy_parts, property_blocks

DATA = Path(__file__).parent / "data" / "energy_parts"


def _read(name):
    out = gzip.open(DATA / f"{name}.out.gz", "rt").read()
    prop = gzip.open(DATA / f"{name}.property.txt.gz", "rt").read()
    return out, prop


def _final(prop):
    return property_blocks(prop)["Single_Point_Data"]["FinalEnergy"]


def _summed(result):
    return sum(result["parts"][k] for k in result["total"])


@pytest.mark.parametrize(
    "name, expected",
    [
        ("water_mp2", {"HF", "MP2 correlation"}),
        ("water_rimp2", {"HF", "RI-MP2 correlation"}),
        ("water_ccsdt", {"HF", "CCSD correlation", "(T)", "MP2 correlation"}),
        (
            "water_f12",
            {
                "HF",
                "CABS singles",
                "CCSD correlation",
                "F12 correction",
                "(T)",
                "MP2 correlation",
            },
        ),
        (
            "water_dlpno",
            {"HF", "DLPNO-MP2 correlation", "DLPNO-CCSD correlation", "DLPNO-(T0)"},
        ),
        (
            "water_dlpnof12",
            {
                "HF",
                "CABS singles",
                "DLPNO-CCSD correlation",
                "F12 correction",
                "DLPNO-(T0)",
                "DLPNO-MP2 correlation",
            },
        ),
    ],
)
def test_parts_of_each_run(name, expected):
    out, prop = _read(name)
    result = energy_parts(out, prop)
    assert result is not None
    assert set(result["parts"]) == expected
    assert result["reference"] == "RHF"
    assert result["units"] == "E_h"
    assert result["frozen core orbitals"] == 1  # O 1s
    assert result["relativistic"] == "none"
    # The parts of the total (not the MP2 by-product of a CC run) are the energy
    assert abs(_summed(result) - _final(prop)) <= SUM_TOLERANCE


def test_f12_parts():
    out, prop = _read("water_f12")
    r = energy_parts(out, prop)
    p = r["parts"]
    # Values as ORCA prints them for this run
    assert p["CABS singles"] == pytest.approx(-0.001072393, abs=1e-12)
    assert p["F12 correction"] == pytest.approx(-0.016427302, abs=1e-12)
    assert p["(T)"] == pytest.approx(-0.0088895434833934240, abs=1e-15)
    # ORCA's printed f, and its definition
    f = r["factors"]["(T) scale f"]
    assert f == pytest.approx(1.059246828, abs=1e-12)
    defined = (p["CCSD correlation"] + p["F12 correction"]) / p["CCSD correlation"]
    assert f == pytest.approx(defined, abs=1e-8)
    assert r["sources"]["(T) scale f"] == "output"
    # HF excludes the CABS singles: SCF energy
    scf = property_blocks(prop)["SCF_Energy"]["totalEnergy"]
    assert p["HF"] == pytest.approx(scf, abs=1e-9)
    assert r["diagnostics"]["T1"] > 0


def test_dlpno_f12_factor_is_computed():
    out, prop = _read("water_dlpnof12")
    r = energy_parts(out, prop)
    p = r["parts"]
    assert r["sources"]["(T) scale f"] == "computed"
    assert r["factors"]["(T) scale f"] == pytest.approx(
        (p["DLPNO-CCSD correlation"] + p["F12 correction"])
        / p["DLPNO-CCSD correlation"]
    )


def test_open_shell_uhf_and_qro():
    out, prop = _read("oh_cc_uhf")
    uhf = energy_parts(out, prop)
    assert uhf["reference"] == "UHF"
    assert "MP2 correlation" not in uhf["parts"]  # not a standalone method here
    assert uhf["diagnostics"]["<S**2>"] == pytest.approx(0.756049)
    assert abs(_summed(uhf) - _final(prop)) <= SUM_TOLERANCE

    out, prop = _read("oh_cc_qro")
    qro = energy_parts(out, prop)
    assert qro["reference"] == "QRO"
    # HF is the QRO determinant's energy, not the UHF SCF energy
    assert qro["parts"]["HF"] == pytest.approx(-75.414436838104464, abs=1e-12)
    assert qro["parts"]["HF"] > uhf["parts"]["HF"]
    assert abs(_summed(qro) - _final(prop)) <= SUM_TOLERANCE

    out, prop = _read("oh_f12_uhf")
    f12 = energy_parts(out, prop)
    assert f12["reference"] == "UHF"
    assert {"CABS singles", "F12 correction", "(T)"} <= set(f12["parts"])
    assert abs(_summed(f12) - _final(prop)) <= SUM_TOLERANCE


def test_inconsistent_parts_are_not_given():
    """Parts that do not sum to the final energy (here, the final energy moved
    by 1e-8 E_h) are not given."""
    import re

    out, prop = _read("water_ccsdt")
    final = _final(prop)
    tampered = re.sub(
        r"(&FinalEnergy \[&Type \"Double\"\]\s+)(\S+)",
        lambda m: m.group(1) + repr(final + 1e-8),
        prop,
    )
    assert _final(tampered) == pytest.approx(final + 1e-8, abs=1e-12)
    assert energy_parts(out, tampered) is None


def test_no_parts_without_the_property_file():
    out, _ = _read("water_ccsdt")
    assert energy_parts(out, "") is None


def test_analyze_task_reports_the_parts(tmp_path):
    """A finished F12 task's results carry its energy parts."""
    from orca_step.batch import analyze_task
    from seamm_exec import Geometry, TaskResult

    out, prop = _read("water_f12")
    (tmp_path / "orca.out").write_text(out)
    (tmp_path / "orca.property.txt").write_text(prop)
    result = TaskResult(
        key="w",
        state="finished",
        directory=tmp_path,
        files={"orca.out": "orca.out", "orca.property.txt": "orca.property.txt"},
    )
    mc = {
        "level": "ORCA:QC@CCSD(T)-F12D_RI/cc-pVTZ-F12",
        "type": "QC",
        "method": "CCSD(T)-F12D_RI",
        "basis": "cc-pVTZ-F12",
        "step": "ORCA",
        "options": {},
    }
    water = Geometry(
        [8, 1, 1], [[0, 0, 0.117], [0, 0.757, -0.468], [0, -0.757, -0.468]]
    )
    data = analyze_task(result, mc, water, properties=("energy",))
    assert data["energy parts"]["reference"] == "RHF"
    assert "CABS singles" in data["energy parts"]["parts"]
    # Without the property file there are no parts, and nothing else changes
    (tmp_path / "orca.property.txt").unlink()
    result = TaskResult(
        key="w", state="finished", directory=tmp_path, files={"orca.out": "orca.out"}
    )
    data = analyze_task(result, mc, water, properties=("energy",))
    assert "energy parts" not in data and "energy" in data
