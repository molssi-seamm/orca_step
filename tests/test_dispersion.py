"""The D4-corrected functionals offered as model chemistries ("R2SCAN-D4",
written "! R2SCAN D4"), and the accuracy keywords on every batch/MDI input."""

from seamm_exec import Geometry

import orca_step
from orca_step import ORCAStep
from orca_step.batch import get_task

WATER = Geometry([8, 1, 1], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.96], [0.9, 0.0, -0.3]])


class _FakeExecutor:
    def __init__(self, name):
        self.name = name


def _offering(level):
    options = ORCAStep.get_model_chemistry_options()
    return options[level]


def test_d4_functionals_are_orca_functionals():
    functionals = orca_step.metadata["functionals"]
    assert orca_step.D4_FUNCTIONALS <= set(functionals)
    for name in orca_step.D4_FUNCTIONALS:
        assert "double-hybrid" not in functionals[name]["category"]


def test_d4_variants_are_offered():
    offering = _offering("ORCA:DFT@R2SCAN-D4/def2-TZVP")
    assert offering["method"] == "R2SCAN-D4"
    assert offering["mdi_capable"]
    assert offering["mdi_method_arg"] == "R2SCAN D4"
    assert offering["analytic_hessian"]
    options = ORCAStep.get_model_chemistry_options()
    # ORCA's own WB97X-D4 is offered once, as itself
    assert options["ORCA:DFT@WB97X-D4/def2-SVP"]["mdi_method_arg"] == "WB97X-D4"
    # Functionals with dispersion of their own get no second one
    for level in (
        "ORCA:DFT@B97M-D4-D4/def2-SVP",
        "ORCA:DFT@WB97X-D3-D4/def2-SVP",
        "ORCA:DFT@R2SCAN-3C-D4/def2-SVP",
        "ORCA:DFT@REVDSD-PBEP86-D4_2021-D4/def2-SVP",
        "ORCA:DFT@B2PLYP-D4/def2-SVP",
    ):
        assert level not in options


def test_split_dispersion():
    assert orca_step.split_dispersion("r2SCAN-D4") == ("R2SCAN", "D4")
    assert orca_step.split_dispersion("PBE0-D4") == ("PBE0", "D4")
    assert orca_step.split_dispersion("B97M-D4") == ("B97M-D4", "")
    assert orca_step.split_dispersion("REVDSD-PBEP86-D4/2021") == (
        "REVDSD-PBEP86-D4/2021",
        "",
    )
    assert orca_step.orca_method_keyword("r2SCAN-D4") == "R2SCAN D4"
    assert orca_step.orca_method_keyword("WB97X-D4") == "WB97X-D4"


def test_batch_input_for_a_d4_variant():
    """Matched (options from the offering) or typed (no options), the input is
    the functional plus D4, with the accuracy keywords."""
    offering = _offering("ORCA:DFT@R2SCAN-D4/def2-TZVP")
    matched = {
        "level": "ORCA:DFT@R2SCAN-D4/def2-TZVPPD",
        "type": "DFT",
        "method": "R2SCAN-D4",
        "basis": "def2-TZVPPD",
        "step": "ORCA",
        "options": {**offering, "mdi_basis_arg": "def2-TZVPPD"},
    }
    typed = {**matched, "method": "r2SCAN-D4", "options": {}}
    for mc in (matched, typed):
        line = get_task(WATER, mc, key="w").files["orca.inp"].splitlines()[0]
        assert line == "! R2SCAN D4 AutoAux def2-TZVPPD TIGHTSCF DEFGRID3 EnGrad"


def test_high_level_input_has_the_accuracy_keywords():
    mc = {
        "level": "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPPD",
        "type": "DFT",
        "method": "REVDSD-PBEP86-D4_2021",
        "basis": "def2-TZVPPD",
        "step": "ORCA",
        "options": {},
    }
    line = get_task(WATER, mc, key="w").files["orca.inp"].splitlines()[0]
    # A double-hybrid gradient runs with VERYTIGHTSCF (orca_base.double_hybrid_scf)
    assert line == (
        "! REVDSD-PBEP86-D4/2021 AutoAux def2-TZVPPD DEFGRID3 EnGrad VERYTIGHTSCF"
    )


def test_mdi_engine_command_for_a_d4_variant(tmp_path):
    (tmp_path / "orca.ini").write_text("[local]\ncode = /opt/orca/orca\n")
    argv = ORCAStep.get_mdi_engine_command(
        _FakeExecutor("local"),
        {"root": str(tmp_path)},
        method="R2SCAN D4",
        basis="def2-TZVP",
        port=8021,
    )
    assert argv[argv.index("--method") + 1] == "R2SCAN D4"
    assert argv[argv.index("--hessian") + 1] == "yes"


def test_energy_keyword_line_for_a_d4_variant():
    node = orca_step.Energy()
    line = node.keyword_line(
        {
            "use model chemistry": "no",
            "method": "R2SCAN-D4",
            "basis": "def2-SVP",
            "basis source": "ORCA internal",
            "auxiliary basis": "AutoAux",
            "extra keywords": "",
        }
    )
    assert line == "R2SCAN D4 def2-SVP AutoAux"


def test_accuracy_keywords_not_doubled():
    """A method that already names an SCF preset or a grid keeps its own; ORCA
    refuses two of a kind."""
    from orca_step.batch import engine_helpers

    helpers = engine_helpers()
    assert helpers.accuracy_keywords("B3LYP") == "TIGHTSCF DEFGRID3"
    assert helpers.accuracy_keywords("B3LYP VeryTightSCF") == "DEFGRID3"
    assert helpers.accuracy_keywords("B3LYP DefGrid2") == "TIGHTSCF"
    assert helpers.accuracy_keywords("B3LYP TIGHTSCF DEFGRID3") == ""
    text = helpers.orca_input(
        "B3LYP VeryTightSCF", "def2-SVP", 0, 1, ["H"], [[0, 0, 0]]
    )
    assert text.splitlines()[0] == "! B3LYP VeryTightSCF def2-SVP DEFGRID3 EnGrad"
