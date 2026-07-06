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


def test_consume_model_chemistry_aliased_functional():
    """A '/'-aliased DFT functional from the model chemistry is translated back
    to the real ORCA keyword on the '!' line."""
    node = orca_step.Energy()
    _stub_model_chemistry(
        node,
        {
            "level": "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVP",
            "owner": "ORCA",
            "type": "DFT",
            "method": "REVDSD-PBEP86-D4_2021",  # aliased ('/' -> '_')
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
    assert line == "REVDSD-PBEP86-D4/2021 def2-TZVP AutoAux TightSCF"


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


def test_keyword_line_dft_functional():
    """With method 'DFT' the '!' keyword is the chosen functional, not 'DFT'."""
    node = orca_step.Energy()
    line = node.keyword_line(
        {
            "use model chemistry": "no",
            "method": "DFT",
            "functional type": "global double-hybrid",
            "functional": "REVDSD-PBEP86-D4/2021",
            "basis": "def2-TZVP",
            "basis source": "ORCA internal",
            "auxiliary basis": "AutoAux",
            "extra keywords": "TightSCF",
        }
    )
    assert line == "REVDSD-PBEP86-D4/2021 def2-TZVP AutoAux TightSCF"


def test_gradient_keyword_analytic_vs_numeric():
    """Analytic-gradient functionals get EnGrad; wB97M(2) and (DLPNO-)CCSD(T)
    have no analytic gradient, so they get NumGrad."""
    node = orca_step.Energy()
    base = {
        "use model chemistry": "no",
        "basis": "def2-SVP",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "extra keywords": "",
        "results": {"gradients": {}},
    }
    # Analytic double hybrid -> EnGrad.
    analytic = node.keyword_line(
        {
            **base,
            "method": "DFT",
            "functional type": "global double-hybrid",
            "functional": "B2PLYP",
        }
    )
    assert "EnGrad" in analytic and "NumGrad" not in analytic
    # Non-self-consistent wB97M(2) -> NumGrad.
    numeric = node.keyword_line(
        {
            **base,
            "method": "DFT",
            "functional type": "range-separated double-hybrid",
            "functional": "WB97M(2)",
        }
    )
    assert "NumGrad" in numeric and "EnGrad" not in numeric
    # DLPNO-CCSD(T): no analytic (T) gradient -> NumGrad.
    cc = node.keyword_line({**base, "method": "DLPNO-CCSD(T)"})
    assert "NumGrad" in cc and "EnGrad" not in cc


def test_metadata_functionals_catalog():
    """The functional catalog carries category, gradient availability, and
    citation keys, and drives the DFT method (no per-functional method entries)."""
    md = orca_step.metadata
    funcs = md["functionals"]
    # The old JSON is gone; the catalog lives in metadata now.
    assert len(funcs) > 100
    assert {"B3LYP", "REVDSD-PBEP86-D4/2021", "WB97M(2)"} <= set(funcs)
    # Every functional's category is one of the declared categories.
    cats = set(md["functional categories"])
    assert all(rec["category"] in cats for rec in funcs.values())
    # Only the non-self-consistent double hybrids are numeric-gradient.
    numeric = {n for n, r in funcs.items() if r["gradients"] == "numeric"}
    assert numeric == {"WB97M(2)", "WB97X-2"}
    # DFT is now a single method; the functional is chosen separately.
    assert "DFT" in md["methods"]
    assert "B3LYP" not in md["methods"] and "B3LYP" in funcs
    # The citation alias still resolves keys that exist in the bibliography.
    assert funcs["B3LYP"]["citations"] == md["dft functionals"]["B3LYP"]


def test_results_have_registered_properties():
    """Every DB-storable result declares a 'property' whose template is
    registered in data/properties.csv (so store_results can create it)."""
    import csv as _csv
    import importlib.resources

    md = orca_step.metadata
    csv_path = importlib.resources.files("orca_step") / "data" / "properties.csv"
    registered = set()
    with csv_path.open(encoding="utf-8") as fh:
        for record in _csv.DictReader(fh):
            registered.add(record["Property"])

    # Results that must now be storable to the database.
    for key in ("gradients", "dipole moment", "mulliken charges", "energy"):
        assert "property" in md["results"][key], f"{key} is missing a property"

    # Every declared property template must be registered in the CSV.
    for key, entry in md["results"].items():
        if "property" in entry:
            assert entry["property"] in registered, f"{key}: {entry['property']}"


def test_dehumanize_bytes():
    """The memory-string parser handles SI, binary, bare-number, and bad input."""
    from orca_step.orca_base import _dehumanize_bytes

    assert _dehumanize_bytes("1000") == 1000
    assert _dehumanize_bytes("2 GB") == 2_000_000_000
    assert _dehumanize_bytes("512MB") == 512_000_000
    assert _dehumanize_bytes("1Gi") == 1024**3
    import pytest

    with pytest.raises(ValueError):
        _dehumanize_bytes("lots of memory")


def test_model_chemistry_advertises_functionals():
    """DFT functionals are advertised as model chemistries; keywords containing
    '/' are aliased with '_', and the alias round-trips to the real keyword."""
    opts = orca_step.ORCAStep.get_model_chemistry_options()

    # A plain functional is advertised under its own name.
    assert any(k.startswith("ORCA:DFT@B3LYP/") for k in opts)
    # A '/'-containing functional is aliased (no reserved '/' inside the method).
    assert any("REVDSD-PBEP86-D4_2021" in k for k in opts)
    assert not any("@REVDSD-PBEP86-D4/2021/" in k for k in opts)
    # Non-DFT methods are still advertised.
    assert any(k.startswith("ORCA:QC@DLPNO-CCSD(T)/") for k in opts)

    # The alias round-trips; non-aliased methods pass through unchanged.
    assert (
        orca_step.mc_method_unalias("REVDSD-PBEP86-D4_2021") == "REVDSD-PBEP86-D4/2021"
    )
    assert orca_step.mc_method_unalias("B3LYP") == "B3LYP"
    assert orca_step.mc_method_unalias("DLPNO-CCSD(T)") == "DLPNO-CCSD(T)"


def test_library_path_vars():
    """The MPI library path is exported for parallel runs, prepended to any
    existing value, and skipped for serial runs or when unset."""
    from orca_step.orca_base import _library_path_vars

    # Serial or no path -> nothing to set.
    assert _library_path_vars(1, "/opt/openmpi/lib") == []
    assert _library_path_vars(4, "") == []

    # Parallel with a path -> both loader variables, no prior value.
    pairs = dict(_library_path_vars(4, "/opt/openmpi/lib", environ={}))
    assert pairs["DYLD_LIBRARY_PATH"] == "/opt/openmpi/lib"
    assert pairs["LD_LIBRARY_PATH"] == "/opt/openmpi/lib"

    # An existing value is preserved after the new directory.
    pairs = dict(
        _library_path_vars(4, "/opt/openmpi/lib", environ={"DYLD_LIBRARY_PATH": "/x"})
    )
    assert pairs["DYLD_LIBRARY_PATH"] == "/opt/openmpi/lib:/x"
