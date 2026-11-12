# -*- coding: utf-8 -*-

"""Frozen-core counts and scalar-relativistic Hamiltonians for model-chemistry
calculations (composite phase 2 §3)."""

import gzip
from pathlib import Path

import pytest

from orca_step.frozen_core import (
    ORCA_DEFAULT_CORE_ELECTRONS,
    check_frozen_core,
    frozen_core_input,
    relativistic_input,
)

WATER = [8, 1, 1]


def test_orca_defaults_as_printed():
    # From ORCA 6.1.1's "Setting NCore(X) from N to ..." for every element
    for z, n in {1: 0, 3: 0, 4: 0, 5: 2, 8: 2, 11: 2, 12: 2, 13: 10, 17: 10,
                 19: 10, 20: 10, 31: 18, 35: 18, 72: 46, 86: 68}.items():  # fmt: skip
        assert ORCA_DEFAULT_CORE_ELECTRONS[z] == n


def test_default_and_none():
    assert frozen_core_input(None, WATER) == ([], "", None)
    assert frozen_core_input("default", WATER) == ([], "", None)
    assert frozen_core_input("none", WATER) == (["NoFrozenCore"], "", 0)


def test_orcas_own_counts_write_nothing():
    """Asking for ORCA's default (W4's valence core for H-F) changes nothing in the
    input, so the settings fingerprint is that of a calculation that does not say."""
    keywords, block, expected = frozen_core_input({"H": 0, "O": 1}, WATER)
    assert (keywords, block, expected) == ([], "", 1)


def test_other_counts_are_written_by_atomic_number():
    _, block, expected = frozen_core_input({"H": 0, "O": 0, "Hf": 23}, WATER)
    assert "NewNCore 8 0 end" in block and "CheckFrozenCore true" in block
    assert "Hf" not in block  # only the elements present; never by symbol
    assert expected == 0
    # Hafnium by atomic number (ORCA 6.1.1 reads 'NewNCore Hf' as hydrogen)
    _, block, expected = frozen_core_input({"Hf": 30, "H": 0}, [72, 1, 1])
    assert "NewNCore 72 60 end" in block
    assert expected == 30
    # Na, Mg: a valence core of 5 orbitals is not ORCA's default of 1
    _, block, _ = frozen_core_input({"Na": 5, "H": 0}, [11, 1])
    assert "NewNCore 11 10 end" in block


def test_an_element_without_a_count_is_refused():
    with pytest.raises(ValueError, match="atomic number 8"):
        frozen_core_input({"H": 0}, WATER)


def test_relativistic_input():
    assert relativistic_input(None) == ([], "")
    assert relativistic_input("DKH2", "point") == (["DKH2"], "%rel FiniteNuc false end")
    assert relativistic_input("X2C", "finite") == (["X2C"], "%rel FiniteNuc true end")
    with pytest.raises(ValueError, match="nucleus model"):
        relativistic_input("DKH2")
    with pytest.raises(ValueError, match="needs a relativistic"):
        relativistic_input(None, "point")


def test_check_frozen_core():
    assert check_frozen_core(None, None) is None
    assert check_frozen_core(1, {"frozen core orbitals": 1}) is None
    assert "froze 1" in check_frozen_core(0, {"frozen core orbitals": 1})
    assert "no count" in check_frozen_core(0, None)


# -- through the task route ------------------------------------------------------


def _mc(method="CCSD(T)", basis="cc-pVTZ", **options):
    return {
        "level": f"ORCA:QC@{method}/{basis}",
        "type": "QC",
        "method": method,
        "basis": basis,
        "step": "ORCA",
        "options": options,
    }


def _water():
    from seamm_exec import Geometry

    return Geometry(WATER, [[0, 0, 0.117], [0, 0.757, -0.468], [0, -0.757, -0.468]])


def test_task_input_and_fingerprint():
    from orca_step.batch import get_task

    plain = get_task(_water(), _mc(), key="w", properties=("energy",))
    same = get_task(
        _water(), _mc(frozen_core={"H": 0, "O": 1}), key="w", properties=("energy",)
    )
    assert same.files["orca.inp"] == plain.files["orca.inp"]
    inner = get_task(
        _water(), _mc(frozen_core={"H": 0, "O": 0}), key="w", properties=("energy",)
    )
    assert "NewNCore 8 0 end" in inner.files["orca.inp"]
    assert inner.files["orca.inp"] != plain.files["orca.inp"]
    rel = get_task(
        _water(),
        _mc(basis="cc-pVTZ-DK", relativistic="DKH2", nucleus="point"),
        key="w",
        properties=("energy",),
    )
    first = rel.files["orca.inp"].splitlines()[0]
    assert "DKH2" in first.split()
    assert "%rel FiniteNuc false end" in rel.files["orca.inp"]


def test_a_frozen_core_mismatch_fails_the_task(tmp_path):
    """The water CCSD(T) run froze 1 orbital (O 1s); asking for all electrons
    there is caught after the run."""
    from seamm_exec import AnalysisError, TaskResult

    from orca_step.batch import analyze_task

    data = Path(__file__).parent / "data" / "energy_parts"
    for name in ("out", "property.txt"):
        text = gzip.open(data / f"water_ccsdt.{name}.gz", "rt").read()
        (tmp_path / f"orca.{name}").write_text(text)
    result = TaskResult(
        key="w",
        state="finished",
        directory=tmp_path,
        files={"orca.out": "orca.out", "orca.property.txt": "orca.property.txt"},
    )
    ok = analyze_task(
        result, _mc(frozen_core={"H": 0, "O": 1}), _water(), properties=("energy",)
    )
    assert ok["energy parts"]["frozen core orbitals"] == 1
    with pytest.raises(AnalysisError, match="froze 1 core orbitals, not the 0"):
        analyze_task(
            result, _mc(frozen_core={"H": 0, "O": 0}), _water(), properties=("energy",)
        )
