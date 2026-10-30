"""A 3c composite (r2SCAN-3c, B97-3c, PBEh-3c, B3LYP-3c, wB97X-3c) runs only with
its own basis: ORCA uses a basis on the '!' line instead of the composite's, and
the result is then not the composite (thermo, TC 8315864 / 8316022, 2026-10-30)."""

import pytest

import orca_step
from orca_step.orca_base import composite_basis
from orca_step.orca_step import ORCAStep

from .test_orca_step import _load_orca_mdi

OWN = {
    "R2SCAN-3C": "def2-mTZVPP",
    "B97-3C": "def2-mTZVP",
    "PBEH-3C": "def2-mSVP",
    "B3LYP-3C": "def2-mSVP",
    "WB97X-3C": "vDZP",
}


@pytest.mark.parametrize("method, own", OWN.items())
def test_own_basis_filled_in_or_kept(method, own):
    assert composite_basis(method, None) == own
    assert composite_basis(method, "") == own
    assert composite_basis(method, own) == own
    assert composite_basis(method.lower(), own.upper()) == own
    with pytest.raises(ValueError, match="its own basis"):
        composite_basis(method, "def2-TZVP")


def test_other_methods_keep_their_basis():
    assert composite_basis("R2SCAN", "def2-TZVP") == "def2-TZVP"
    assert composite_basis("B3LYP D3BJ", None) is None


def test_engine_table_matches_metadata():
    mdi = _load_orca_mdi()
    for name, rec in orca_step.metadata["functionals"].items():
        if rec.get("basis"):
            assert mdi.COMPOSITE_BASES[name] == rec["basis"]
    assert set(mdi.COMPOSITE_BASES) == set(OWN)
    assert mdi.engine_basis("R2SCAN-3C", None) == "def2-mTZVPP"
    assert mdi.engine_basis("B3LYP", None) == "def2-SVP"
    with pytest.raises(SystemExit):
        mdi.engine_basis("R2SCAN-3C", "def2-TZVP")


def test_offered_only_with_own_basis():
    options = ORCAStep.get_model_chemistry_options()
    for method, own in OWN.items():
        keys = [k for k in options if k.startswith(f"ORCA:DFT@{method}/")]
        assert keys == [f"ORCA:DFT@{method}/{own}"], keys


def test_keyword_line_uses_own_basis():
    node = orca_step.Energy()
    P = {
        "use model chemistry": "no",
        "method": "DFT",
        "functional": "R2SCAN-3C",
        "basis": "def2-mTZVPP",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "extra keywords": "",
    }
    assert node.keyword_line(P).startswith("R2SCAN-3C def2-mTZVPP")
    P["basis"] = "def2-TZVP"
    with pytest.raises(ValueError, match="its own basis"):
        node.keyword_line(P)


def test_parameters_narrow_and_imply_the_basis():
    P = orca_step.EnergyParameters()
    values = P.current_values()
    values.update(
        {"use model chemistry": "no"},
        method="DFT",
        functional="PBEH-3C",
        basis="def2-TZVP",
    )
    assert P.choices("basis", values) == ("def2-mSVP",)
    implied = P.implied(values)
    assert implied["basis"] == "def2-mSVP"
    assert implied["basis source"] == "ORCA internal"
    assert not P.applies("basis set extrapolation", values)
    assert not P.applies("basis source", values)
