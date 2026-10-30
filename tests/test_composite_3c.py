"""The "3c" composites (r2SCAN-3c, B97-3c, PBEh-3c, B3LYP-3c, wB97X-3c) are methods
of their own, ORCA:CMP@<name>, run by their keyword alone: ORCA uses a basis on the
'!' line instead of the composite's, and the result is then not the composite
(thermo, TC 8315864 / 8316022, 2026-10-30)."""

import pytest

import orca_step
from orca_step.orca_base import composite_basis, composite_method
from orca_step.orca_step import ORCAStep

from .test_orca_step import _load_orca_mdi

OWN = {
    "R2SCAN-3C": "def2-mTZVPP",
    "B97-3C": "def2-mTZVP",
    "PBEH-3C": "def2-mSVP",
    "B3LYP-3C": "def2-mSVP",
    "WB97X-3C": "vDZP",
}


def test_metadata():
    assert orca_step.metadata["composite methods"] == OWN
    for name in OWN:
        assert orca_step.metadata["methods"][name]["type"] == "CMP"


@pytest.mark.parametrize("method, own", OWN.items())
def test_no_basis_words_and_another_is_refused(method, own):
    assert composite_method(method.lower()) == method
    assert composite_basis(method, None) == ""
    assert composite_basis(method, own) == ""  # its own reproduces the keyword
    with pytest.raises(ValueError, match="its own basis"):
        composite_basis(method, "def2-TZVP")


def test_other_methods_keep_their_basis():
    assert composite_method("R2SCAN") is None
    assert composite_basis("R2SCAN", "def2-TZVP") == "def2-TZVP"


def test_engine():
    mdi = _load_orca_mdi()
    assert mdi.COMPOSITE_BASES == OWN
    assert mdi.engine_basis("R2SCAN-3C", None) == ""
    assert mdi.engine_basis("B3LYP", None) == "def2-SVP"
    with pytest.raises(SystemExit):
        mdi.engine_basis("R2SCAN-3C", "def2-TZVP")


def test_offered_as_cmp_without_a_basis():
    options = ORCAStep.get_model_chemistry_options()
    for method in OWN:
        assert f"ORCA:CMP@{method}" in options
        option = options[f"ORCA:CMP@{method}"]
        assert option["basis"] is None and option["mdi_method_arg"] == method
        assert not [k for k in options if k.startswith(f"ORCA:DFT@{method}")]


def test_engine_command_has_no_basis(monkeypatch):
    monkeypatch.setattr(
        ORCAStep,
        "get_executor_config",
        classmethod(lambda cls, e, o: {"mdi_script": "orca_mdi.py", "code": "orca"}),
    )
    argv = ORCAStep.get_mdi_engine_command(None, {}, method="R2SCAN-3C", port=1)
    assert "--basis" not in argv
    assert argv[argv.index("--auxiliary-basis") + 1] == "none"
    argv = ORCAStep.get_mdi_engine_command(None, {}, method="B3LYP", port=1)
    assert argv[argv.index("--basis") + 1] == "def2-SVP"


def test_batch_method_words():
    from orca_step.batch import method_and_blocks

    mc = {
        "level": "ORCA:CMP@R2SCAN-3C",
        "type": "CMP",
        "method": "R2SCAN-3C",
        "basis": None,
        "options": {"mdi_method_arg": "R2SCAN-3C", "mdi_basis_arg": None},
    }
    words, basis, blocks = method_and_blocks(mc)
    assert words == "R2SCAN-3C" and basis == "" and blocks == []


def test_keyword_line_is_the_keyword_alone():
    node = orca_step.Energy()
    P = {
        "use model chemistry": "no",
        "method": "R2SCAN-3C",
        "basis": "def2-TZVP",  # does not apply to a composite: ignored
        "basis source": "ORCA internal",
        "auxiliary basis": "AutoAux",
        "basis set extrapolation": "none",
        "extra keywords": "",
    }
    line = node.keyword_line(P)
    assert line.split()[0] == "R2SCAN-3C"
    assert "def2-TZVP" not in line and "AutoAux" not in line


def test_parameters_rules():
    P = orca_step.EnergyParameters()
    values = P.current_values()
    values.update({"use model chemistry": "no", "method": "PBEH-3C"})
    for key in ("basis", "basis source", "auxiliary basis", "basis set extrapolation"):
        assert not P.applies(key, values), key
    assert "R2SCAN-3C" in P["method"].enumeration
    assert "R2SCAN-3C" not in P["functional"].enumeration


def _model_chemistry_node(monkeypatch, basis=None):
    """An Energy node seeing the model chemistry ORCA:CMP@R2SCAN-3C[/basis]."""
    level = "ORCA:CMP@R2SCAN-3C" + (f"/{basis}" if basis else "")
    mc = {
        "level": level,
        "owner": "ORCA",
        "type": "CMP",
        "method": "R2SCAN-3C",
        "basis": basis,
    }
    node = orca_step.Energy()
    monkeypatch.setattr(
        node, "variable_exists", lambda name: name == "_model_chemistry"
    )
    monkeypatch.setattr(node, "get_variable", lambda name: mc)
    return node


_MC_P = {
    "use model chemistry": "yes",
    "method": "B3LYP",
    "basis": "def2-TZVP",  # the step's default: must not be used for a composite
    "basis source": "ORCA internal",
    "auxiliary basis": "AutoAux",
    "basis set extrapolation": "none",
    "extra keywords": "",
}


def test_model_chemistry_without_basis_ignores_the_step_basis(monkeypatch):
    """ORCA:CMP@R2SCAN-3C through 'use model chemistry' once took the step's default
    basis and then refused it (ChemAI 5502, Frequencies, 2026-10-30)."""
    node = _model_chemistry_node(monkeypatch)
    assert node._resolve_method_basis(_MC_P) == ("R2SCAN-3C", "")
    line = node.keyword_line(_MC_P)
    assert line.split()[0] == "R2SCAN-3C"
    assert "def2-TZVP" not in line and "AutoAux" not in line
    assert node._model_string(_MC_P) == "CMP@R2SCAN-3C"


def test_model_chemistry_with_its_own_basis(monkeypatch):
    node = _model_chemistry_node(monkeypatch, "def2-mTZVPP")
    assert node.keyword_line(_MC_P).split()[0] == "R2SCAN-3C"
    node = _model_chemistry_node(monkeypatch, "def2-TZVP")
    with pytest.raises(ValueError, match="its own basis"):
        node.keyword_line(_MC_P)


def test_explicit_composite_level_strings():
    node = orca_step.Energy()
    P = dict(_MC_P, **{"use model chemistry": "no", "method": "R2SCAN-3C"})
    assert node._resolve_method_basis(P) == ("R2SCAN-3C", "")
    assert node._model_string(P) == "CMP@R2SCAN-3C"
