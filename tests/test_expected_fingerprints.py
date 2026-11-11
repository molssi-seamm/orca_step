# -*- coding: utf-8 -*-

"""``batch.expected_fingerprints``: the settings fingerprints a calculation of a
structure would record, predicted without running ORCA.

The golden data are real records from SEAMM job 5684 (ORCA 6.1.1, see
tests/data/fingerprints_job5684.json): methane and borane on the same r2SCAN-3c
geometry key (sha256:6b46…, COSX moot for r2SCAN-3c), and their revDSD single
points, which differ because boron gets exact exchange (NoCOSX, orca_step#44):
sha256:1e16… for methane, sha256:7ca1… for borane.
"""

import json
from pathlib import Path

import pytest

from seamm_exec import Geometry
from seamm_util import element_data

from orca_step.batch import expected_fingerprints, get_task
from orca_step.fingerprint import orca_calculation_record

DATA = Path(__file__).parent / "data" / "fingerprints_job5684.json"

#: The Model Chemistry step's wrappers for the presets' levels (2026-11-10)
R2SCAN = {
    "level": "ORCA:CMP@R2SCAN-3C",
    "owner": "ORCA",
    "type": "CMP",
    "method": "R2SCAN-3C",
    "basis": None,
    "cutoff": None,
    "step": "ORCA",
    "options": {
        "mdi_capable": True,
        "mdi_method_arg": "R2SCAN-3C",
        "mdi_basis_arg": None,
        "prefers_batch": True,
        "analytic_hessian": True,
    },
}
REVDSD = {
    "level": "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP",
    "owner": "ORCA",
    "type": "DFT",
    "method": "REVDSD-PBEP86-D4_2021",
    "basis": "def2-TZVPP",
    "cutoff": None,
    "step": "ORCA",
    "options": {
        "mdi_capable": True,
        "mdi_method_arg": "REVDSD-PBEP86-D4/2021",
        "mdi_basis_arg": "def2-TZVPP",
        "prefers_batch": True,
        "analytic_hessian": False,
    },
}
B3LYP = {
    "level": "ORCA:DFT@B3LYP/def2-SVP",
    "owner": "ORCA",
    "type": "DFT",
    "method": "B3LYP",
    "basis": "def2-SVP",
    "step": "ORCA",
    "options": {
        "mdi_capable": True,
        "mdi_method_arg": "B3LYP",
        "mdi_basis_arg": "def2-SVP",
        "prefers_batch": True,
    },
}


def _species():
    out = {}
    for s in json.loads(DATA.read_text())["species"]:
        numbers = [element_data[e]["atomic number"] for e in s["symbols"]]
        geometry = Geometry(
            numbers,
            s["coordinates"],
            charge=s["charge"],
            multiplicity=s["multiplicity"],
        )
        out[s["species"]] = (geometry, s)
    return out


SPECIES = _species()


@pytest.mark.parametrize("species", ["C", "B"])
def test_revdsd_energy_fingerprints_of_job_5684(species):
    geometry, recorded = SPECIES[species]
    got = expected_fingerprints(
        geometry, REVDSD, kind="energy", version=recorded["energy_version"]
    )
    assert recorded["energy_fingerprint"] in got
    # The as-written fingerprint comes first: revDSD has exact exchange, so COSX
    # is not moot and the stored one is the first
    assert got[0] == recorded["energy_fingerprint"]


def test_boron_and_methane_differ_at_revdsd():
    (c, _), (b, _) = SPECIES["C"], SPECIES["B"]
    fc = expected_fingerprints(c, REVDSD, version="6.1.1")
    fb = expected_fingerprints(b, REVDSD, version="6.1.1")
    assert len(fc) == 1  # no NoCOSX: moot or not, the same input
    assert len(fb) == 2 and fc[0] not in fb[:1]
    # With COSX moot the NoCOSX of boron drops out, and the two agree
    assert fb[1] == fc[0]


@pytest.mark.parametrize("species", ["C", "B"])
def test_r2scan_hessian_fingerprint_of_job_5684(species):
    """The geometry key: the MDI engine's analytic Hessian at r2SCAN-3c. COSX is
    moot for r2SCAN-3c, so methane and borane share it."""
    geometry, recorded = SPECIES[species]
    got = expected_fingerprints(
        geometry, R2SCAN, kind="hessian", version=recorded["geometry_version"]
    )
    assert recorded["frequency_fingerprint"] in got


def test_one_center_and_sodium_get_exact_exchange():
    water = Geometry(
        [8, 1, 1], [[0, 0, 0.117], [0, 0.757, -0.467], [0, -0.757, -0.467]]
    )
    plain = expected_fingerprints(water, B3LYP, version="6.1.1")
    for structure in (
        Geometry([8], [[0, 0, 0]], multiplicity=3),
        Geometry([11, 17], [[0, 0, 0], [2.4, 0, 0]]),
        Geometry([15, 1, 1, 1], [[0, 0, 0], [1.4, 0, 0], [0, 1.4, 0], [0, 0, 1.4]]),
    ):
        got = expected_fingerprints(structure, B3LYP, version="6.1.1")
        assert len(got) == 2 and got[0] not in plain and got[1] == plain[0]
    # The element set can be changed, as for the run
    nacl = Geometry([11, 17], [[0, 0, 0], [2.4, 0, 0]])
    assert expected_fingerprints(
        nacl, B3LYP, version="6.1.1", options={"exact_exchange_elements": "none"}
    ) == (plain)


def test_the_prediction_is_the_runs_record():
    """The energy prediction hashes exactly the files the run's task writes, as
    analyze_task's record does with ORCA's output."""
    nacl = Geometry([11, 17], [[0, 0, 0], [2.4, 0, 0]])
    task = get_task(nacl, B3LYP, key="x", properties=("energy",))
    hybrid_out = (
        "Program Version 6.1.1\n Hybrid DFT is turned on\n"
        "   Fraction HF Exchange ScalHFX         ....  0.200000\n"
    )
    record = orca_calculation_record(task.files, "6.1.1", out_text=hybrid_out)
    assert expected_fingerprints(nacl, B3LYP, version="6.1.1")[0] == (
        record["fingerprint"]
    )


def test_the_version_is_part_of_the_fingerprint():
    """Major.minor is hashed (seamm_util.fingerprint.hashed_version)."""
    geometry, _ = SPECIES["C"]
    ours = expected_fingerprints(geometry, REVDSD, version="6.1.1")
    assert ours == expected_fingerprints(geometry, REVDSD, version="6.1.0")
    assert ours != expected_fingerprints(geometry, REVDSD, version="6.0.1")


def test_unknown_kind():
    geometry, _ = SPECIES["C"]
    with pytest.raises(ValueError):
        expected_fingerprints(geometry, REVDSD, kind="gradient", version="6.1.1")
