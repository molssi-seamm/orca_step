#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the `orca_step` package."""

import pytest  # noqa: F401

import orca_step


def test_main_factory():
    """The main stevedore helper class."""
    step = orca_step.ORCAStep()
    assert step.description()["name"] == "ORCA"
    assert step.description()["group"] == "Simulations"


def test_substep_factories():
    """The Energy and Optimization sub-step helpers."""
    assert orca_step.EnergyStep().description()["name"] == "Energy"
    assert orca_step.OptimizationStep().description()["name"] == "Optimization"


def test_main_node_has_subflowchart():
    """The main ORCA node owns a sub-flowchart."""
    node = orca_step.ORCA()
    assert node.subflowchart is not None
    assert str(type(node)) == "<class 'orca_step.orca.ORCA'>"


def test_energy_defaults():
    """Energy defaults: model chemistry on, DLPNO-CCSD(T)/def2-TZVP, AutoAux."""
    P = orca_step.EnergyParameters()
    assert P["use model chemistry"].value == "yes"
    assert P["method"].value == "DLPNO-CCSD(T)"
    assert P["basis"].value == "def2-TZVP"
    assert P["auxiliary basis"].value == "AutoAux"


def test_optimization_extends_energy():
    """Optimization is an Energy with the convergence control."""
    assert issubclass(orca_step.Optimization, orca_step.Energy)
    P = orca_step.OptimizationParameters()
    assert P["optimization convergence"].value == "NormalOpt"
    # Inherits the energy parameters
    assert P["method"].value == "DLPNO-CCSD(T)"


def test_keyword_line_explicit():
    """The explicit keyword line assembles method/basis/aux/extra."""
    node = orca_step.Energy()
    line = node.keyword_line(
        {
            "use model chemistry": "no",
            "method": "B3LYP",
            "basis": "def2-SVP",
            "basis source": "ORCA internal",
            "auxiliary basis": "AutoAux",
            "extra keywords": "TightSCF",
        }
    )
    assert line == "B3LYP def2-SVP AutoAux TightSCF"
    # 'none' aux is omitted
    line2 = node.keyword_line(
        {
            "use model chemistry": "no",
            "method": "HF",
            "basis": "6-31G*",
            "basis source": "ORCA internal",
            "auxiliary basis": "none",
            "extra keywords": "",
        }
    )
    assert line2 == "HF 6-31G*"


def test_keyword_line_bse_omits_basis():
    """With the Basis Set Exchange source, the basis name is left off the '!'
    line (it is supplied as an embedded file instead)."""
    node = orca_step.Energy()
    line = node.keyword_line(
        {
            "use model chemistry": "no",
            "method": "HF",
            "basis": "cc-pVDZ",
            "basis source": "Basis Set Exchange",
            "auxiliary basis": "AutoAux",
            "extra keywords": "TightSCF",
        }
    )
    assert line == "HF AutoAux TightSCF"


def test_keyword_line_gradients():
    """Requesting the 'gradients' result appends EnGrad to the keyword line."""
    node = orca_step.Energy()
    base = {
        "use model chemistry": "no",
        "method": "B3LYP",
        "basis": "def2-SVP",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "extra keywords": "",
    }
    assert "EnGrad" not in node.keyword_line(base)
    assert node.keyword_line({**base, "results": {"gradients": {}}}) == (
        "B3LYP def2-SVP EnGrad"
    )


def test_keyword_line_keepdensity():
    """Requesting the wavefunction adds 'keepdensity' (orca_2aim then makes wfx)."""
    node = orca_step.Energy()
    base = {
        "use model chemistry": "no",
        "method": "PBE0",
        "basis": "def2-SVP",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "extra keywords": "",
    }
    assert "keepdensity" not in node.keyword_line(base)
    line = node.keyword_line({**base, "save wavefunction": "yes"})
    assert "keepdensity" in line


def test_parse_gradients_engrad(tmp_path):
    """The gradient is read from orca.engrad as an [n_atoms, 3] array (E_h/bohr)."""
    (tmp_path / "orca.engrad").write_text(
        "#\n# Number of atoms\n#\n 2\n"
        "#\n# The current total energy in Eh\n#\n   -1.5\n"
        "#\n# The current gradient in Eh/bohr\n#\n"
        "   0.1\n   0.2\n   0.3\n  -0.1\n  -0.2\n  -0.3\n"
        "#\n# The atomic numbers and current coordinates in Bohr\n#\n"
        "  1   0.0 0.0 0.0\n  1   0.0 0.0 1.4\n"
    )
    grad = orca_step.Energy._parse_gradients(orca_step.Energy(), tmp_path)
    assert grad == [[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]]


def test_bse_shorthand_forces_bse():
    """A 'bse:NAME' basis forces the Basis Set Exchange even when the source
    toggle is 'ORCA internal': the name is dropped from the '!' line (embedded
    as a file) and the shorthand strips to the bare basis name."""
    node = orca_step.Energy()
    P = {
        "use model chemistry": "no",
        "method": "HF",
        "basis": "bse:cc-pVDZ",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "extra keywords": "",
    }
    assert node._using_bse(P)
    assert node._strip_bse("bse:cc-pVDZ") == "cc-pVDZ"
    assert node._strip_bse("def2-SVP") == "def2-SVP"
    assert node.keyword_line(P) == "HF"


def test_bse_basis_file():
    """The BSE basis file is ORCA-format ($DATA/$END) for the requested
    elements, with no '!'-comment header."""
    content = orca_step.Energy._bse_basis_file("cc-pVDZ", [1, 8])
    assert "$DATA" in content
    assert not content.lstrip().startswith("!")  # header omitted
    assert "OXYGEN" in content.upper()


def test_get_model_chemistry_options():
    """ORCA advertises ORCA:<type>@<method>/<basis> and nothing periodic/MDI."""
    opts = orca_step.ORCAStep().get_model_chemistry_options()
    key = "ORCA:QC@DLPNO-CCSD(T)/def2-TZVP"
    assert key in opts
    assert opts[key]["method"] == "DLPNO-CCSD(T)"
    assert opts[key]["basis"] == "def2-TZVP"
    assert opts[key]["type"] == "QC"
    assert orca_step.ORCAStep().get_model_chemistry_options(periodic_only=True) == {}
    assert orca_step.ORCAStep().get_model_chemistry_options(mdi_only=True) == {}


def _stub_model_chemistry(node, mc):
    """Make the node read `mc` as the global model chemistry, without needing
    the flowchart-variable store that only exists during a real run.
    """
    node.variable_exists = lambda v: True
    node.get_variable = lambda v: mc


def test_consume_model_chemistry():
    """An ORCA-owned model chemistry yields the right keyword line."""
    node = orca_step.Energy()
    _stub_model_chemistry(
        node,
        {
            "level": "ORCA:DFT@B3LYP/def2-TZVP",
            "owner": "ORCA",
            "type": "DFT",
            "method": "B3LYP",
            "basis": "def2-TZVP",
        },
    )
    line = node.keyword_line(
        {
            "use model chemistry": "yes",
            "method": "IGNORED",
            "basis": "IGNORED",
            "basis source": "ORCA internal",
            "auxiliary basis": "AutoAux",
            "extra keywords": "TightSCF",
        }
    )
    assert line == "B3LYP def2-TZVP AutoAux TightSCF"


def test_consume_rejects_foreign_owner():
    """A model chemistry owned by another program is rejected clearly."""
    node = orca_step.Energy()
    _stub_model_chemistry(
        node,
        {
            "level": "MOPAC:SQM@PM6",
            "owner": "MOPAC",
            "type": "SQM",
            "method": "PM6",
            "basis": "",
        },
    )
    with pytest.raises(RuntimeError, match="MOPAC"):
        node.keyword_line(
            {
                "use model chemistry": "yes",
                "method": "X",
                "basis": "Y",
                "auxiliary basis": "none",
                "extra keywords": "",
            }
        )


def test_bibliography_loads_libraries():
    """data/references.bib parses and the support-library + BSE entries load."""
    node = orca_step.Energy()
    assert "libint2" in node._bibliography
    assert "libxc" in node._bibliography
    assert "bse" in node._bibliography  # cited only when BSE supplies the basis


def test_dft_functional_citations():
    """DFT functionals map to citation keys that exist in the bibliography
    (bibtex in references.bib, keys in metadata -- no duplication)."""
    node = orca_step.Energy()
    funcs = orca_step.metadata["dft functionals"]
    assert {"B3LYP", "PBE0", "M062X"} <= set(funcs)
    for fkey in ("B3LYP", "PBE0"):
        assert funcs[fkey]  # has at least one citation key
        for key in funcs[fkey]:
            assert key in node._bibliography


def test_metadata_methods_and_bases():
    """Metadata exposes methods (incl. DLPNO-CCSD(T)) and the basis families."""
    assert "DLPNO-CCSD(T)" in orca_step.metadata["methods"]
    bases = orca_step.metadata["basis sets"]
    assert "def2-TZVP" in bases and "cc-pVTZ" in bases and "6-311++G**" in bases
