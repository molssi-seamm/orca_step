# -*- coding: utf-8 -*-

"""The F12 model chemistries (P4b): canonical CCSD(T)-F12D/RI and thermo's
DLPNO-CCSD(T)-F12D(TightPNO,TCutPNO=1e-8), with cc-pVTZ-F12, its CABS and
cc-pVQZ/C, run through the task route exactly as the Energy sub-step ran them.

The golden fingerprints are those of the fragment store's canonical F12 rows,
harvested from jobs 5608-5620 (Energy sub-step, ORCA 6.1.1): their hashed input is
'! CC-PVQZ/C CC-PVTZ-F12 CC-PVTZ-F12-CABS CCSD(T)-F12D/RI TIGHTSCF', sha256:119baa11…,
and with NoCOSX (boron) sha256:36760b4f…. The fingerprint is of the settings only,
so any structure without (or with) NoCOSX elements gives the same one.
"""

import pytest

from seamm_exec import Geometry

from orca_step.batch import expected_fingerprints, get_task
from orca_step.orca_step import ORCAStep, mc_method_alias, mc_method_unalias

CANONICAL = "ORCA:QC@CCSD(T)-F12D_RI/cc-pVTZ-F12"
DLPNO = "ORCA:QC@DLPNO-CCSD(T)-F12D(TightPNO,TCutPNO=1e-8)/cc-pVTZ-F12"

ETHANE = Geometry(
    [6, 6, 1, 1, 1, 1, 1, 1],
    [
        [0.0, 0.0, 0.765],
        [0.0, 0.0, -0.765],
        [1.02, 0.0, 1.16],
        [-0.51, 0.88, 1.16],
        [-0.51, -0.88, 1.16],
        [-1.02, 0.0, -1.16],
        [0.51, -0.88, -1.16],
        [0.51, 0.88, -1.16],
    ],
)
BF3 = Geometry(
    [5, 9, 9, 9],
    [[0, 0, 0], [1.31, 0, 0], [-0.655, 1.134, 0], [-0.655, -1.134, 0]],
)
METHYL = Geometry(
    [6, 1, 1, 1],
    [[0, 0, 0], [1.08, 0, 0], [-0.54, 0.935, 0], [-0.54, -0.935, 0]],
    multiplicity=2,
)


def _wrapper(level):
    options = ORCAStep.get_model_chemistry_options()[level]
    return {
        "level": level,
        "owner": "ORCA",
        "type": options["type"],
        "method": options["method"],
        "basis": options["basis"],
        "cutoff": None,
        "step": "ORCA",
        "options": options,
    }


def _first_line(task):
    return task.files["orca.inp"].splitlines()[0]


def test_both_are_advertised_and_not_mdi():
    options = ORCAStep.get_model_chemistry_options()
    for level in (CANONICAL, DLPNO):
        assert level in options
        assert options[level]["aux_basis"] == "cc-pVQZ/C"
        assert options[level]["mdi_capable"] is False
    mdi = ORCAStep.get_model_chemistry_options(mdi_only=True)
    assert CANONICAL not in mdi and DLPNO not in mdi
    # NormalPNO DLPNO-F12 and other F12 bases are not offered (not validated)
    assert not [k for k in options if "F12" in k and k not in (CANONICAL, DLPNO)]


def test_the_alias_round_trips():
    assert mc_method_alias("CCSD(T)-F12D/RI") == "CCSD(T)-F12D_RI"
    assert mc_method_unalias("CCSD(T)-F12D_RI") == "CCSD(T)-F12D/RI"
    assert mc_method_unalias("REVDSD-PBEP86-D4_2021") == "REVDSD-PBEP86-D4/2021"
    assert mc_method_unalias("B3LYP") == "B3LYP"


def test_canonical_input_is_the_energy_substeps():
    task = get_task(ETHANE, _wrapper(CANONICAL), key="k", properties=("energy",))
    words = set(_first_line(task)[1:].split())
    assert words == {
        "CCSD(T)-F12D/RI",
        "cc-pVTZ-F12",
        "cc-pVTZ-F12-CABS",
        "cc-pVQZ/C",
        "TIGHTSCF",
    }


@pytest.mark.parametrize(
    "structure, fingerprint",
    [
        (
            ETHANE,
            "sha256:119baa115802c193899287ebef85607cf38e6ee54f34de4dfa3a68c46809c6f0",
        ),
        (BF3, "sha256:36760b4f"),
    ],
)
def test_canonical_fingerprints_are_the_stores(structure, fingerprint):
    got = expected_fingerprints(
        structure, _wrapper(CANONICAL), kind="energy", version="6.1.1"
    )
    assert got[0].startswith(fingerprint)


def test_dlpno_input_is_the_validations():
    """As in part A's validation runs (5516, 5536, 5560): TightPNO, cc-pVQZ/C, the
    CABS, TightSCF and %mdci TCutPNO 1e-8 end; no grid."""
    task = get_task(ETHANE, _wrapper(DLPNO), key="k", properties=("energy",))
    words = set(_first_line(task)[1:].split())
    assert words == {
        "DLPNO-CCSD(T)-F12D",
        "TightPNO",
        "cc-pVTZ-F12",
        "cc-pVTZ-F12-CABS",
        "cc-pVQZ/C",
        "TIGHTSCF",
    }
    assert "%mdci TCutPNO 1e-8 end" in task.files["orca.inp"]


def test_open_shell_f12_is_not_refused_as_a_dlpno_double_hybrid():
    for level in (CANONICAL, DLPNO):
        task = get_task(METHYL, _wrapper(level), key="k", properties=("energy",))
        assert "DEFGRID3" not in _first_line(task)
