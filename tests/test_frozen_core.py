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


def test_decontract():
    from orca_step.batch import get_task
    from orca_step.frozen_core import decontract_block

    assert decontract_block(None) == "" and decontract_block(False) == ""
    assert decontract_block(True) == "%basis Decontract true end"
    with pytest.raises(ValueError):
        decontract_block("yes")
    plain = get_task(_water(), _mc(), key="w", properties=("energy",))
    decon = get_task(_water(), _mc(decontract=True), key="w", properties=("energy",))
    assert "%basis Decontract true end" in decon.files["orca.inp"]
    assert decon.files["orca.inp"] != plain.files["orca.inp"]


# -- the reference of open-shell calculations ---------------------------------------


def _oh(multiplicity=2):
    from seamm_exec import Geometry

    return Geometry([8, 1], [[0, 0, 0], [0, 0, 0.9697]], multiplicity=multiplicity)


def test_reference_input():
    from orca_step.frozen_core import reference_input

    assert reference_input("QRO", 2) == "%mdci UseQROs true end"
    # Closed shells, and ORCA's own UHF, write nothing
    assert reference_input("QRO", 1) == ""
    assert reference_input("UHF", 2) == "" and reference_input(None, 2) == ""
    with pytest.raises(ValueError, match="must be one of"):
        reference_input("ROHF", 2)


def test_qro_changes_only_open_shell_inputs():
    from orca_step.batch import get_task

    plain = get_task(_water(), _mc(), key="w", properties=("energy",))
    qro = get_task(_water(), _mc(reference="QRO"), key="w", properties=("energy",))
    assert qro.files["orca.inp"] == plain.files["orca.inp"]
    uhf = get_task(_oh(), _mc(), key="oh", properties=("energy",))
    same = get_task(_oh(), _mc(reference="UHF"), key="oh", properties=("energy",))
    assert same.files["orca.inp"] == uhf.files["orca.inp"]
    qro = get_task(_oh(), _mc(reference="QRO"), key="oh", properties=("energy",))
    assert "%mdci UseQROs true end" in qro.files["orca.inp"]
    assert qro.files["orca.inp"] != uhf.files["orca.inp"]


def test_qro_joins_an_existing_mdci_block():
    """DLPNO-F12 has its own %mdci block (TCutPNO); QRO goes into it."""
    from orca_step.batch import get_task
    from orca_step.orca_step import mc_method_alias

    method = mc_method_alias("DLPNO-CCSD(T)-F12D(TightPNO,TCutPNO=1e-8)")
    mc = _mc(method, "cc-pVTZ-F12", reference="QRO", aux_basis="cc-pVQZ/C")
    text = get_task(_oh(), mc, key="oh", properties=("energy",)).files["orca.inp"]
    mdci = [line for line in text.splitlines() if line.lower().startswith("%mdci")]
    assert mdci == ["%mdci TCutPNO 1e-8 UseQROs true end"]


@pytest.mark.parametrize("fixture, ok", [("oh_cc_qro", True), ("oh_cc_uhf", False)])
def test_the_reference_is_checked_after_the_run(tmp_path, fixture, ok):
    from seamm_exec import AnalysisError, TaskResult

    from orca_step.batch import analyze_task

    data = Path(__file__).parent / "data" / "energy_parts"
    for name in ("out", "property.txt"):
        text = gzip.open(data / f"{fixture}.{name}.gz", "rt").read()
        (tmp_path / f"orca.{name}").write_text(text)
    result = TaskResult(
        key="oh",
        state="finished",
        directory=tmp_path,
        files={"orca.out": "orca.out", "orca.property.txt": "orca.property.txt"},
    )
    if ok:
        out = analyze_task(result, _mc(reference="QRO"), _oh(), properties=("energy",))
        assert out["energy parts"]["reference"] == "QRO"
    else:
        with pytest.raises(AnalysisError, match="quasi-restricted orbitals"):
            analyze_task(result, _mc(reference="QRO"), _oh(), properties=("energy",))
        # UHF asked for, UHF run: fine
        analyze_task(result, _mc(reference="UHF"), _oh(), properties=("energy",))


# -- one-electron systems: nothing to correlate --------------------------------------


def _h(charge=0, n=1):
    from seamm_exec import Geometry

    xyz = [[0, 0, 0], [0, 0, 1.06]][:n]
    multiplicity = 2 if (n - charge) % 2 else 1
    return Geometry([1] * n, xyz, charge=charge, multiplicity=multiplicity)


def test_which_model_chemistries_are_correlated():
    from orca_step.frozen_core import is_correlated

    assert is_correlated(_mc())
    assert is_correlated(_mc("CCSD(T)-F12D_RI", "cc-pVTZ-F12"))
    assert is_correlated(_mc("RI-MP2"))
    assert not is_correlated(_mc("HF"))
    assert not is_correlated({**_mc("B3LYP"), "type": "DFT"})


def test_electron_pairs():
    from orca_step.frozen_core import check_electron_pairs, correlation_electrons

    assert check_electron_pairs({"H": 0}, [1], 0)  # H
    assert check_electron_pairs(None, [1, 1], 1)  # H2+
    assert not check_electron_pairs({"H": 0}, [1, 1], 0)  # H2: one pair
    assert not check_electron_pairs("default", [3], 0)  # Li: ORCA freezes nothing
    assert correlation_electrons({"Li": 1}, [3], 0) == (3, 1)
    with pytest.raises(ValueError, match="leaves 1 of the 3 electrons correlated"):
        check_electron_pairs({"Li": 1}, [3], 0)


@pytest.mark.parametrize(
    "mc",
    [
        _mc("CCSD(T)", "cc-pwCVTZ", frozen_core={"H": 0}, reference="QRO"),
        _mc("CCSD(T)-F12D_RI", "cc-pVTZ-F12", reference="QRO", aux_basis="cc-pVQZ/C"),
        _mc("RI-MP2", "cc-pVTZ"),
    ],
)
def test_a_one_electron_system_runs_as_hf(mc):
    """ORCA's MP2 and MDCI modules refuse one electron, so the task is HF in the
    same basis: no correlation keywords, frozen core, CABS or QRO block."""
    from orca_step.batch import get_task

    for system in (_h(), _h(charge=1, n=2)):
        text = get_task(system, mc, key="h", properties=("energy",)).files["orca.inp"]
        line = text.splitlines()[0]
        assert line.split()[:2] == ["!", "HF"], line
        assert f" {mc['basis']} " in f"{line} "
        for word in ("CCSD", "MP2", "F12D", "CABS", "AutoAux", "UseQROs", "NewNCore"):
            assert word not in text
    # Two electrons are correlated as asked
    text = get_task(_h(n=2), mc, key="h2", properties=("energy",)).files["orca.inp"]
    assert text.splitlines()[0].split()[1] != "HF"


def test_the_hamiltonian_is_kept_for_one_electron():
    from orca_step.batch import get_task

    mc = _mc(
        "CCSD(T)", "cc-pVTZ-DK", relativistic="DKH2", nucleus="point", decontract=True
    )
    text = get_task(_h(), mc, key="h", properties=("energy",)).files["orca.inp"]
    assert text.splitlines()[0].split()[:3] == ["!", "HF", "DKH2"]
    assert "%rel FiniteNuc false end" in text
    assert "%basis Decontract true end" in text


def test_a_frozen_core_without_pairs_is_refused():
    from seamm_exec import Geometry

    from orca_step.batch import get_task

    li = Geometry([3], [[0, 0, 0]], multiplicity=2)
    with pytest.raises(ValueError, match="no electron pairs"):
        get_task(li, _mc(frozen_core={"Li": 1}), key="li", properties=("energy",))


@pytest.mark.parametrize("method", ["CCSD(T)", "CCSD(T)-F12D_RI"])
def test_one_electron_parts(tmp_path, method):
    """The H atom's HF/cc-pVTZ-DK (DKH2) run, analysed as the correlated
    calculation asked for: its energy is the HF part and the correlation is zero."""
    from seamm_exec import TaskResult

    from orca_step.batch import analyze_task

    data = Path(__file__).parent / "data" / "energy_parts"
    for name in ("out", "property.txt"):
        text = gzip.open(data / f"h_atom_hf_dkh2.{name}.gz", "rt").read()
        (tmp_path / f"orca.{name}").write_text(text)
    result = TaskResult(
        key="h",
        state="finished",
        directory=tmp_path,
        files={"orca.out": "orca.out", "orca.property.txt": "orca.property.txt"},
    )
    mc = _mc(
        method,
        "cc-pVTZ-DK",
        frozen_core={"H": 0},
        reference="QRO",
        relativistic="DKH2",
        nucleus="point",
    )
    out = analyze_task(result, mc, _h(), properties=("energy",))
    parts = out["energy parts"]
    assert parts["parts"] == {"HF": pytest.approx(-0.499816346412, abs=1e-11)}
    assert parts["total"] == ["HF"] and parts["one electron"]
    assert out["energy"] == pytest.approx(-0.499816346412 * 2625.4996394798254)
    # The reference asked for (one electron: UHF = ROHF = QRO), no frozen core
    assert parts["reference"] == "QRO"
    assert parts["frozen core orbitals"] == 0
    assert parts["relativistic"] == "DKH2"
    assert ("omitted" in parts) == ("F12" in method)
