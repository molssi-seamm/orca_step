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
            "auxiliary basis": "none",
            "extra keywords": "",
        }
    )
    assert line2 == "HF 6-31G*"


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


def test_metadata_methods_and_bases():
    """Metadata exposes methods (incl. DLPNO-CCSD(T)) and the basis families."""
    assert "DLPNO-CCSD(T)" in orca_step.metadata["methods"]
    bases = orca_step.metadata["basis sets"]
    assert "def2-TZVP" in bases and "cc-pVTZ" in bases and "6-311++G**" in bases
