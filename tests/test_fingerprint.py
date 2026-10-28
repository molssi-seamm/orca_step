# -*- coding: utf-8 -*-

"""The fingerprints of the settings of ORCA calculations (orca_step.fingerprint):
what is hashed, and that the ORCA step and the MDI engine give the same
fingerprint for the same calculation."""

import importlib.util
from pathlib import Path

import pytest

from orca_step.fingerprint import (
    basis_file_blocks,
    orca_calculation_record,
    orca_version,
    settings_text,
)

SOURCE = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "run_path"
ENGINE = SOURCE / "orca_step" / "data" / "orca_mdi.py"

STEP_INPUT = """\
! B3LYP D3BJ def2-TZVP AutoAux TightSCF DEFGRID3 EnGrad NoCOSX
%pal nprocs 8 end
%maxcore 4000
%scf
  Guess MORead
  SThresh 1.000e-07
end
%moinp "orca.gbw"
%output Print[ P_Hirshfeld ] 1 end
* xyz 0 1
P 0.0 0.0 0.0
F 1.6 0.0 0.0
*
"""


def fingerprint(text, version="6.1.1", **kwargs):
    return orca_calculation_record({"orca.inp": text}, version, **kwargs)["fingerprint"]


def test_what_is_hashed():
    lines, elements = settings_text(STEP_INPUT)
    assert lines == [
        "! AUTOAUX B3LYP D3BJ DEF2-TZVP DEFGRID3 NOCOSX TIGHTSCF",
        "%SCF",
        "STHRESH 1.000E-07",
        "END",
    ]
    assert elements == {}


def test_same_settings_written_differently():
    """Order and case of keywords, one-line blocks, the task, the cores and
    memory, the guess, the geometry, charge and multiplicity do not count."""
    other = """\
! b3lyp autoaux def2-tzvp D3BJ NoCOSX TIGHTSCF DEFGRID3 Opt PAL4
%SCF SThresh 1.000e-07 end
*xyz -1 2
C 0 0 0
*
"""
    assert fingerprint(STEP_INPUT) == fingerprint(other)
    # a patch release of ORCA too, but not a minor one
    assert fingerprint(STEP_INPUT) == fingerprint(other, version="6.1.0")
    assert fingerprint(STEP_INPUT) != fingerprint(other, version="6.0.1")


@pytest.mark.parametrize(
    "change",
    [
        ("DEFGRID3", "DEFGRID2"),
        ("NoCOSX", ""),
        ("TightSCF", "VeryTightSCF"),
        ("SThresh 1.000e-07", "SThresh 1.000e-06"),
        ("def2-TZVP", "def2-SVP"),
        ("B3LYP", "B3LYP UHF"),
    ],
)
def test_settings_that_count(change):
    old, new = change
    assert fingerprint(STEP_INPUT) != fingerprint(STEP_INPUT.replace(old, new, 1))


def test_element_blocks():
    text = """\
! HF
%basis
 NewGTO H "STO-3G" end
 NewECP Pb
  N_core 60
  lmax f
 end
end
* xyz 0 1
H 0 0 0
*
"""
    lines, elements = settings_text(text)
    assert lines == ["! HF"]
    assert elements == {
        "H": ["NEWGTO", '"STO-3G"'],
        "Pb": ["NEWECP", "N_CORE 60", "LMAX F"],
    }


def test_bse_basis_file():
    """A Basis Set Exchange basis is hashed per element, so ethane and methanol
    share the fingerprint and those of H and C."""
    bse = pytest.importorskip("basis_set_exchange")
    text = '! B3LYP\n%basis GTOName "basis.bas" end\n* xyz 0 1\nC 0 0 0\n*\n'
    hco = bse.get_basis("def2-SVP", elements=[1, 6, 8], fmt="orca", header=False)
    hc = bse.get_basis("def2-SVP", elements=[1, 6], fmt="orca", header=False)
    assert sorted(basis_file_blocks(hco)) == ["C", "H", "O"]
    a = orca_calculation_record({"orca.inp": text, "basis.bas": hco}, "6.1.1")
    b = orca_calculation_record({"orca.inp": text, "basis.bas": hc}, "6.1.1")
    assert a["fingerprint"] == b["fingerprint"]
    assert a["files"] == {"orca.inp": "! B3LYP"}
    for symbol in ("H", "C"):
        assert a["element fingerprints"][symbol] == b["element fingerprints"][symbol]
    assert "O" not in b["element fingerprints"]


def test_orca_version():
    assert orca_version((DATA / "orca.out").read_text()) == "6.1.1"
    assert orca_version("no version here") is None


def _engine():
    spec = importlib.util.spec_from_file_location("_orca_mdi_engine", ENGINE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_engine_exchange_scheme():
    engine = _engine()
    assert engine.exchange_method("B3LYP", ["P", "F"]) == "B3LYP NoCOSX"
    assert engine.exchange_method("B3LYP", ["C", "H"]) == "B3LYP"
    assert engine.exchange_method("B3LYP", ["P", "F"], cosx="on") == "B3LYP"
    assert engine.exchange_method("B3LYP", ["C", "H"], cosx="off") == "B3LYP NoCOSX"
    assert engine.exchange_method("B3LYP RIJK", ["C"], cosx="off") == "B3LYP RIJK"
    with pytest.raises(ValueError):
        engine.exchange_method("B3LYP", ["C"], cosx="maybe")


SPEC = """\
title: ORCA fingerprint test
steps:
- Water: {}
- ORCA:
    steps:
    - Energy:
        use model chemistry: 'no'
        method: DFT
        basis: def2-SVP
        grid: DEFGRID3
"""


def test_step_and_engine_agree(tmp_path):
    """The ORCA step stores the fingerprint that the MDI engine gives for the
    same calculation (B3LYP/def2-SVP, TightSCF, DEFGRID3, AutoAux)."""
    import molsystem
    from seamm_exec.testing import fake_program, run_spec

    code = fake_program(
        tmp_path / "bin" / "orca",
        files=[
            DATA / "orca.property.txt",
            DATA / "orca.bibtex",
        ],
        stdout=DATA / "orca.out",
    )
    job = run_spec(
        tmp_path,
        SPEC,
        inis={"orca": {"code": code}},
        args=["--ncores", "1"],
        source=SOURCE,
    )

    db = molsystem.SystemDB(filename=str(job / "seamm.db"))
    try:
        configuration = db.system.configuration
        stored = configuration.properties.get("settings fingerprint#ORCA#*")
        assert len(stored) == 1
        ((name, value),) = stored.items()
        step_fingerprint = value["value"]
        model = name.split("#", 2)[2]
        record = configuration.properties.get(f"calculation#ORCA#{model}")
        record = record[f"calculation#ORCA#{model}"]["value"]
        # The energy is stored though no result was selected
        assert f"total energy#ORCA#{model}" in configuration.properties.get(
            "total energy#ORCA#*"
        )
        symbols = configuration.atoms.symbols
        coordinates = configuration.atoms.get_coordinates(fractionals=False)
    finally:
        db.close()
    assert record["version"] == "6.1.1"
    assert record["fingerprint"] == step_fingerprint

    engine = _engine()
    text = engine.orca_input(
        engine.exchange_method("B3LYP AutoAux", symbols),
        "def2-SVP",
        0,
        1,
        symbols,
        coordinates,
        ncores=4,
    )
    engine_record = orca_calculation_record({"orca.inp": text}, "6.1.1")
    assert engine_record["fingerprint"] == step_fingerprint, (
        record["files"],
        engine_record["files"],
    )
