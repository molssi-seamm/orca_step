#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the `orca_step` package."""

import importlib.resources
import importlib.util
import sys

import pytest  # noqa: F401

import orca_step


def _load_orca_mdi():
    """Load the standalone data/orca_mdi.py engine as a module (it is not part
    of an importable package; `mdi` is imported lazily inside main, so importing
    the module here does not require it)."""
    path = importlib.resources.files("orca_step") / "data" / "orca_mdi.py"
    spec = importlib.util.spec_from_file_location("orca_mdi", str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _FakeExecutor:
    def __init__(self, name):
        self.name = name


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


def test_keyword_line_cbs_extrapolation():
    """CBS extrapolation replaces the fixed basis with ORCA's Extrapolate(...)."""
    node = orca_step.Energy()
    line = node.keyword_line(
        {
            "use model chemistry": "no",
            "method": "DLPNO-CCSD(T)",
            "basis": "def2-TZVP",  # ignored when extrapolating
            "basis source": "ORCA internal",
            "auxiliary basis": "AutoAux",
            "extra keywords": "",
            "basis set extrapolation": "2/3",
            "extrapolation family": "cc",
        }
    )
    assert line == "DLPNO-CCSD(T) Extrapolate(2/3,cc) AutoAux"


def test_cbs_extrapolation_blocks_gradients():
    """Requesting gradients with CBS extrapolation is a clear error -- ORCA has
    no gradient for an extrapolated energy."""
    node = orca_step.Energy()
    P = {
        "use model chemistry": "no",
        "method": "HF",
        "basis": "def2-SVP",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "extra keywords": "",
        "basis set extrapolation": "3/4",
        "extrapolation family": "def2",
        "results": {"gradients": {}},
    }
    with pytest.raises(RuntimeError, match="extrapolat"):
        node.keyword_line(P)


def test_optimization_ignores_cbs():
    """The Optimization step never extrapolates: CBS has no gradient, so even a
    hand-set extrapolation is ignored (the control is hidden in its GUI)."""
    node = orca_step.Optimization()
    line = node.keyword_line(
        {
            "use model chemistry": "no",
            "method": "HF",
            "basis": "def2-SVP",
            "basis source": "ORCA internal",
            "auxiliary basis": "none",
            "extra keywords": "",
            "basis set extrapolation": "2/3",  # would extrapolate on an Energy step
            "extrapolation family": "cc",
        }
    )
    assert "Extrapolate" not in line
    assert line == "HF def2-SVP"


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


def test_keyword_line_grid():
    """The integration-grid preset is appended; 'default' emits nothing."""
    node = orca_step.Energy()
    base = {
        "use model chemistry": "no",
        "method": "B3LYP",
        "basis": "def2-SVP",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "extra keywords": "",
    }
    assert node.keyword_line({**base, "grid": "default"}) == "B3LYP def2-SVP"
    assert node.keyword_line({**base, "grid": "DEFGRID3"}) == "B3LYP def2-SVP DEFGRID3"


def test_keyword_line_scf_convergence():
    """The SCF convergence preset is appended; 'default' emits nothing."""
    node = orca_step.Energy()
    base = {
        "use model chemistry": "no",
        "method": "HF",
        "basis": "def2-SVP",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "extra keywords": "",
    }
    assert node.keyword_line({**base, "scf convergence": "default"}) == "HF def2-SVP"
    assert (
        node.keyword_line({**base, "scf convergence": "TIGHTSCF"})
        == "HF def2-SVP TIGHTSCF"
    )


def test_scf_and_extra_keyword_defaults():
    """The SCF tolerance now has its own control (default TIGHTSCF); the extra-
    keywords default is empty so they do not double up."""
    P = orca_step.EnergyParameters()
    assert P["scf convergence"].value == "TIGHTSCF"
    assert P["extra keywords"].value == ""


def test_basis_name_forms():
    """The basis name is extracted whether the value arrives as the dict (run
    pass), its string repr (the pre-run description pass), or a plain name."""
    E = orca_step.Energy
    assert E._basis_name({"name": "$basis", "elements": []}) == "$basis"
    assert E._basis_name("{'name': '$basis', 'elements': []}") == "$basis"
    assert E._basis_name("def2-SVP") == "def2-SVP"


def test_basis_variable_expansion(monkeypatch):
    """A $variable typed into the basis field is expanded against the flowchart
    variables (the basis is a special dict parameter, so SEAMM does not expand
    it the way it does ordinary string parameters)."""
    import seamm
    from seamm.variables import Variables

    monkeypatch.setattr(seamm, "flowchart_variables", Variables(basis="def2-SVP"))
    node = orca_step.Energy()
    assert node._expand_variables("$basis") == "def2-SVP"
    assert node._expand_variables("def2-TZVP") == "def2-TZVP"  # plain name unchanged
    line = node.keyword_line(
        {
            "use model chemistry": "no",
            "method": "HF",
            "basis": {"name": "$basis", "elements": []},
            "basis source": "ORCA internal",
            "auxiliary basis": "none",
            "extra keywords": "",
        }
    )
    assert line == "HF def2-SVP"


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
    """ORCA advertises ORCA:<type>@<method>/<basis>; nothing for periodic. The
    full list includes methods without an analytic gradient (e.g.
    DLPNO-CCSD(T)); mdi_only keeps only the MDI-drivable subset."""
    opts = orca_step.ORCAStep().get_model_chemistry_options()
    key = "ORCA:QC@DLPNO-CCSD(T)/def2-TZVP"
    assert key in opts
    assert opts[key]["method"] == "DLPNO-CCSD(T)"
    assert opts[key]["basis"] == "def2-TZVP"
    assert opts[key]["type"] == "QC"
    assert orca_step.ORCAStep().get_model_chemistry_options(periodic_only=True) == {}
    # mdi_only now returns the analytic-gradient subset (non-empty), and it is a
    # strict subset of the full list.
    mdi = orca_step.ORCAStep().get_model_chemistry_options(mdi_only=True)
    assert 0 < len(mdi) < len(opts)


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


# --- MDI engine wrapper (data/orca_mdi.py) ---------------------------------


def test_orca_mdi_input():
    """The engine input is a single-point EnGrad job; %pal appears only in
    parallel; the geometry is written in Angstrom."""
    mod = _load_orca_mdi()
    text = mod.orca_input(
        "B3LYP AutoAux",
        "def2-SVP",
        0,
        1,
        ["O", "H", "H"],
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.96], [0.9, 0.0, -0.3]],
        ncores=1,
    )
    assert text.startswith("! B3LYP AutoAux def2-SVP EnGrad")
    assert "* xyz 0 1" in text and text.rstrip().endswith("*")
    assert "%pal" not in text
    assert "%pal nprocs 4 end" in mod.orca_input(
        "HF", "def2-SVP", 0, 1, ["H", "H"], [[0, 0, 0], [0, 0, 0.74]], ncores=4
    )


def test_orca_mdi_parse_energy_and_gradient():
    """Parse the final energy (last match) and the engrad gradient array."""
    mod = _load_orca_mdi()
    out = (
        "FINAL SINGLE POINT ENERGY   -76.100000\n"
        "...\nFINAL SINGLE POINT ENERGY   -76.400000\n"
    )
    assert mod.parse_energy(out) == -76.4
    assert mod.parse_energy("nothing here") is None

    engrad = (
        "#\n# Number of atoms\n#\n 2\n"
        "#\n# The current total energy in Eh\n#\n  -1.500000\n"
        "#\n# The current gradient in Eh/bohr\n#\n"
        " 0.1\n 0.2\n 0.3\n -0.1\n -0.2\n -0.3\n"
        "#\n# coordinates in Bohr\n#\n 8 0.0 0.0 0.0\n 1 0.0 0.0 1.4\n"
    )
    grad = mod.parse_engrad(engrad, 2)
    assert grad.tolist() == [[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]]


def test_orca_mdi_hessian_input():
    """The MDI Hessian input requests an analytic Hessian (AnFreq)."""
    mod = _load_orca_mdi()
    text = mod.orca_hessian_input(
        "HF", "def2-SVP", 0, 1, ["H", "H"], [[0, 0, 0], [0, 0, 0.74]], ncores=1
    )
    assert text.startswith("! HF def2-SVP AnFreq")
    assert "* xyz 0 1" in text


def test_orca_mdi_parse_hessian():
    """Parse the ORCA .hess $hessian block, including its 5-column blocking."""
    mod = _load_orca_mdi()
    hess = (
        "$orca_hessian_file\n\n$hessian\n6\n"
        "        0      1      2      3      4\n"
        "  0   1.0    0.1    0.0    0.0    0.0\n"
        "  1   0.1    2.0    0.0    0.0    0.0\n"
        "  2   0.0    0.0    3.0    0.0    0.0\n"
        "  3   0.0    0.0    0.0    4.0    0.0\n"
        "  4   0.0    0.0    0.0    0.0    5.0\n"
        "  5   0.0    0.0    0.0    0.0    0.0\n"
        "        5\n"
        "  0   0.0\n  1   0.0\n  2   0.0\n  3   0.0\n  4   0.0\n  5   6.0\n"
        "$end\n"
    )
    H = mod.parse_hessian(hess, 2)  # 2 atoms -> 6x6
    assert H.shape == (6, 6)
    assert H[0, 0] == 1.0 and H[5, 5] == 6.0
    assert H[0, 1] == 0.1 and H[1, 0] == 0.1  # off-diagonal, both blocks


def test_method_has_analytic_hessian():
    """Analytic Hessian: yes for HF/MP2/ordinary DFT; no for double hybrids and
    (DLPNO-)CCSD(T) (which have an analytic gradient but no analytic Hessian)."""
    from orca_step.orca_step import method_has_analytic_hessian as has

    assert has("HF") is True
    assert has("MP2") is True
    assert has("B3LYP") is True
    assert has("REVDSD-PBEP86-D4/2021") is False  # double hybrid
    assert has("DLPNO-CCSD(T)") is False
    assert has("CCSD(T)-F12D/RI") is False


def test_mdi_engine_command_advertises_hessian(tmp_path):
    """The launcher passes --hessian yes only for analytic-Hessian methods."""
    (tmp_path / "orca.ini").write_text("[local]\ncode = /opt/orca/orca\n")

    def argv(method):
        return orca_step.ORCAStep.get_mdi_engine_command(
            _FakeExecutor("local"),
            {"root": str(tmp_path)},
            method=method,
            basis="def2-TZVP",
            port=8021,
        )

    hf = argv("HF")
    assert hf[hf.index("--hessian") + 1] == "yes"
    rev = argv("REVDSD-PBEP86-D4/2021")
    assert rev[rev.index("--hessian") + 1] == "no"


def test_get_mdi_engine_command(tmp_path):
    """The engine command runs orca_mdi.py with the orca binary and the
    method/basis flags, over TCP with the driver's port."""
    (tmp_path / "orca.ini").write_text("[local]\ncode = /opt/orca/orca\n")
    argv = orca_step.ORCAStep.get_mdi_engine_command(
        _FakeExecutor("local"),
        {"root": str(tmp_path)},
        method="B3LYP",
        basis="def2-TZVP",
        port=8021,
        charge=0,
        multiplicity=1,
        ncores=2,
    )
    assert argv[0] == sys.executable
    assert any("orca_mdi.py" in a for a in argv)
    assert argv[argv.index("--orca") + 1] == "/opt/orca/orca"
    assert argv[argv.index("--method") + 1] == "B3LYP"
    assert argv[argv.index("--basis") + 1] == "def2-TZVP"
    assert argv[argv.index("--ncores") + 1] == "2"
    init = argv[argv.index("-mdi") + 1]
    assert "-role ENGINE" in init and "-port 8021" in init


def test_model_chemistry_mdi_only():
    """mdi_only keeps analytic-gradient methods (with the real keyword) and
    drops the ones with no analytic gradient."""
    opts = orca_step.ORCAStep.get_model_chemistry_options(mdi_only=True)
    assert all(o["mdi_capable"] for o in opts.values())
    # An analytic-gradient functional is present.
    assert any(k.startswith("ORCA:DFT@B3LYP/") for k in opts)
    # DLPNO-CCSD(T) has no analytic gradient -> not MDI-capable.
    assert not any("DLPNO-CCSD(T)" in k for k in opts)
    # An aliased '/' functional carries the real keyword for the engine.
    revs = [o for k, o in opts.items() if "REVDSD-PBEP86-D4_2021" in k]
    assert revs and revs[0]["mdi_method_arg"] == "REVDSD-PBEP86-D4/2021"


# --------------------------------------------------------------------------
# BSSE (counterpoise) sub-step
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Explicitly-correlated (F12) methods
# --------------------------------------------------------------------------
def _f12_P(method, basis, extra=""):
    return {
        "use model chemistry": "no",
        "method": method,
        "basis": basis,
        "basis source": "ORCA internal",
        "auxiliary basis": "AutoAux",
        "extra keywords": extra,
    }


def test_f12_methods_and_bases_in_metadata():
    """The F12 methods and their orbital bases are offered, with numeric
    gradients (no analytic gradient in ORCA)."""
    m = orca_step.metadata["methods"]
    # Canonical F12 needs the '/RI' suffix; DLPNO-F12 must not have it.
    assert "CCSD(T)-F12D/RI" in m and "DLPNO-CCSD(T)-F12D" in m
    assert "CCSD(T)-F12D" not in m
    assert m["DLPNO-CCSD(T)-F12D"]["gradients"] == "numeric"
    bases = orca_step.metadata["basis sets"]
    assert {"cc-pVDZ-F12", "cc-pVTZ-F12", "cc-pVQZ-F12"} <= set(bases)


def test_f12_methods_excluded_from_model_chemistry():
    """F12 methods need F12-specific bases, so they are not advertised as
    generic model chemistries (which pair with non-F12 bases)."""
    opts = orca_step.ORCAStep.get_model_chemistry_options()
    assert not any("F12" in k.upper() for k in opts)


def test_cabs_keyword_derivation():
    """The CABS is derived from the F12 basis (tracking DZ/TZ/QZ), skipped for
    non-F12 methods/bases, and not duplicated if the user supplied one."""
    node = orca_step.Energy()
    assert (
        node._cabs_keyword(_f12_P("DLPNO-CCSD(T)-F12D", "cc-pVTZ-F12"))
        == "cc-pVTZ-F12-CABS"
    )
    assert (
        node._cabs_keyword(_f12_P("CCSD(T)-F12D/RI", "cc-pVDZ-F12"))
        == "cc-pVDZ-F12-CABS"
    )
    assert node._cabs_keyword(_f12_P("DLPNO-CCSD(T)", "cc-pVTZ")) == ""  # not F12
    assert (
        node._cabs_keyword(_f12_P("CCSD(T)-F12D/RI", "cc-pVTZ")) == ""
    )  # non-F12 basis
    # Already supplied by the user -> not added twice.
    assert (
        node._cabs_keyword(_f12_P("CCSD(T)-F12D/RI", "cc-pVTZ-F12", "cc-pVTZ-F12-CABS"))
        == ""
    )


def test_keyword_line_f12_adds_cabs():
    """keyword_line puts both the F12 basis and its CABS on the '!' line."""
    line = orca_step.Energy().keyword_line(_f12_P("DLPNO-CCSD(T)-F12D", "cc-pVTZ-F12"))
    assert "cc-pVTZ-F12" in line and "cc-pVTZ-F12-CABS" in line


def test_check_f12_rejects_non_f12_basis():
    """An F12 method with a non-F12 basis (and no manual CABS) fails early."""
    node = orca_step.Energy()
    with pytest.raises(RuntimeError, match="F12"):
        node._check_f12(_f12_P("CCSD(T)-F12D/RI", "cc-pVTZ"))
    # OK with an F12 basis, or with a hand-supplied CABS.
    node._check_f12(_f12_P("CCSD(T)-F12D/RI", "cc-pVTZ-F12"))
    node._check_f12(_f12_P("CCSD(T)-F12D/RI", "cc-pVQZ", "cc-pVQZ-F12-CABS"))


def test_bsse_compound_input_adds_f12_cabs():
    """The BSSE compound input also injects the CABS for an F12 method."""
    node = orca_step.BSSE()
    P = {
        "use model chemistry": "no",
        "method": "DLPNO-CCSD(T)-F12D",
        "basis": "cc-pVTZ-F12",
        "basis source": "ORCA internal",
        "auxiliary basis": "AutoAux",
        "grid": "default",
        "scf convergence": "default",
        "extra keywords": "TightPNO",
        "optimize monomers": "no",
    }
    block, _, _ = node._compound_input(P, "bsse.xyz", "bssenergy.cmp")
    assert "cc-pVTZ-F12-CABS" in block


def test_max_am_from_bse():
    """Angular momentum is read from the Basis Set Exchange (offline)."""
    assert orca_step.Energy._max_am_from_bse("cc-pVTZ", [8]) == 3  # f
    assert orca_step.Energy._max_am_from_bse("cc-pV5Z", [8]) == 5  # h


def test_auto_grid_threshold():
    """cc-pV5Z (h) reaches the auto-DEFGRID3 threshold; cc-pVTZ (f) does not."""
    from orca_step.energy import _HIGH_ANGULAR_MOMENTUM

    assert orca_step.Energy._max_am_from_bse("cc-pV5Z", [8]) >= _HIGH_ANGULAR_MOMENTUM
    assert orca_step.Energy._max_am_from_bse("cc-pVTZ", [8]) < _HIGH_ANGULAR_MOMENTUM


def test_bsse_compound_input_wavefunction():
    """'save wavefunction' toggles ProduceWavefunction in the Compound block."""
    node = orca_step.BSSE()
    base = {
        "use model chemistry": "no",
        "method": "B3LYP",
        "basis": "def2-SVP",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "grid": "default",
        "scf convergence": "default",
        "extra keywords": "",
        "optimize monomers": "no",
    }
    off, _, _ = node._compound_input({**base, "save wavefunction": "no"}, "b.xyz", "s")
    on, _, _ = node._compound_input({**base, "save wavefunction": "yes"}, "b.xyz", "s")
    assert "ProduceWavefunction = false;" in off
    assert "ProduceWavefunction = true;" in on


def test_bsse_cmp_scripts_support_wavefunction():
    """Both shipped Compound scripts keep the dimer density on request."""
    import importlib.resources

    for name in ("bssenergy.cmp", "bssegradient.cmp"):
        text = (importlib.resources.files("orca_step") / "data" / name).read_text()
        assert "ProduceWavefunction" in text
        assert "KeepDensity" in text


# --------------------------------------------------------------------------
# Frequencies (Hessian / vibrational) sub-step
# --------------------------------------------------------------------------
_FREQ_OUT = """
VIBRATIONAL FREQUENCIES
-----------------------

Scaling factor for frequencies =  1.000000000  (already applied!)

     0:       0.00 cm**-1
     5:       0.00 cm**-1
     6:     -50.00 cm**-1
     7:    1790.72 cm**-1
     8:    4062.03 cm**-1

IR SPECTRUM
-----------

 Mode   freq       eps      Int      T**2
       cm**-1   L/(mol*cm) km/mol    a.u.
------------------------------------------------------------
  7:   1790.72   0.015824   79.97  0.002758  ( 0.0  0.0  0.05)
  8:   4062.03   0.013551   68.48  0.001041  ( 0.0  0.03 0.0)

Zero point energy                ...      0.02238377 Eh      14.05 kcal/mol
Total Enthalpy                    ...    -75.93482201 Eh
Final Gibbs free energy         ...    -75.95623266 Eh
"""


def test_frequencies_factory():
    """The Frequencies sub-step helper."""
    assert orca_step.FrequenciesStep().description()["name"] == "Frequencies"


def test_frequencies_extends_energy():
    """Frequencies is an Energy with the second-derivative + temperature controls."""
    assert issubclass(orca_step.Frequencies, orca_step.Energy)
    node = orca_step.Frequencies()
    assert node._calculation == "frequencies"
    P = orca_step.FrequenciesParameters()
    assert P["second derivatives"].value == "analytic"
    assert P["temperature"].value == "298.15"
    assert P["method"].value == "DLPNO-CCSD(T)"  # inherits the energy params


def test_frequencies_extra_input_temperature():
    """extra_input adds the '%freq Temp' thermochemistry block."""
    node = orca_step.Frequencies()
    P = {
        "use model chemistry": "no",
        "method": "HF",
        "basis": "def2-SVP",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "grid": "default",
        "scf convergence": "default",
        "extra keywords": "",
        "Hirshfeld charges": "no",
        "polarizability": "no",
        "temperature": 350.0,
    }
    blocks, _ = node.extra_input(P)
    assert "%freq Temp 350.0000 end" in blocks


def test_frequencies_parsers():
    """Frequencies, IR intensities, and thermochemistry parse from ORCA output."""
    node = orca_step.Frequencies()
    # _parse_frequencies returns all 3N modes in order (classification is separate).
    assert node._parse_frequencies(_FREQ_OUT) == [0.0, 0.0, -50.0, 1790.72, 4062.03]
    assert node._parse_ir_intensities(_FREQ_OUT) == [79.97, 68.48]
    # Thermochemistry is converted from E_h to kJ/mol (SEAMM's SI default).
    thermo = node._parse_thermochemistry(_FREQ_OUT)
    assert thermo == pytest.approx(
        {
            "zero point energy": 58.76858,
            "enthalpy": -199366.84781,
            "gibbs energy": -199423.06147,
        }
    )


def test_frequencies_classify():
    """The 6 nominally-zero modes (non-linear) are split off; the imaginary mode
    stays, and the largest zero-mode magnitude is reported."""

    class _Atoms:
        def get_coordinates(self, fractionals=False):
            # A bent (non-linear) 3-atom geometry -> 6 zero modes expected.
            return [[0.0, 0.4, 0.0], [-0.8, -0.2, 0.0], [0.8, -0.2, 0.0]]

    class _Config:
        n_atoms = 3
        atoms = _Atoms()

    node = orca_step.Frequencies()
    # 9 modes: an imaginary (-50), 6 near-zero (incl. a 7.3 residual), 2 real.
    all_freqs = [-50.0, -7.3, -0.1, 0.0, 0.0, 0.2, 5.0, 1790.72, 4062.03]
    vibrational, max_zero = node._classify_frequencies(all_freqs, _Config())
    # The 6 smallest-|f| are the trans/rot; -50 (imaginary) and the two real
    # stretches remain as vibrational.
    assert vibrational == [-50.0, 1790.72, 4062.03]
    assert max_zero == pytest.approx(7.3)


def test_frequencies_is_linear():
    """A diatomic / collinear geometry is detected as linear (5 zero modes)."""

    class _Atoms:
        def __init__(self, xyz):
            self._xyz = xyz

        def get_coordinates(self, fractionals=False):
            return self._xyz

    class _Config:
        def __init__(self, xyz):
            self.atoms = _Atoms(xyz)

    linear = _Config([[0.0, 0.0, 0.0], [0.0, 0.0, 1.1], [0.0, 0.0, 2.2]])
    bent = _Config([[0.0, 0.4, 0.0], [-0.8, -0.2, 0.0], [0.8, -0.2, 0.0]])
    assert orca_step.Frequencies._is_linear(linear) is True
    assert orca_step.Frequencies._is_linear(bent) is False


def test_frequencies_results_in_metadata():
    """The frequency/thermochemistry results are 'frequencies'-only."""
    results = orca_step.metadata["results"]
    for key in ("frequencies", "IR intensities", "zero point energy", "gibbs energy"):
        assert results[key]["calculation"] == ["frequencies"]


def test_bsse_factory():
    """The BSSE sub-step helper."""
    assert orca_step.BSSEStep().description()["name"] == "BSSE"


def test_bsse_extends_energy():
    """BSSE is an Energy with the fragment controls and the 'bsse' calculation."""
    assert issubclass(orca_step.BSSE, orca_step.Energy)
    node = orca_step.BSSE()
    assert node._calculation == "bsse"
    P = orca_step.BSSEParameters()
    assert P["fragments"].value == "auto (2 molecules)"
    assert P["optimize monomers"].value == "no"
    # Inherits the energy parameters.
    assert P["method"].value == "DLPNO-CCSD(T)"


def test_bsse_never_extrapolates():
    """CBS extrapolation is off for BSSE (an extrapolated energy has no
    gradient)."""
    node = orca_step.BSSE()
    assert not node._extrapolating(
        {"use model chemistry": "no", "basis set extrapolation": "3/4"}
    )


def test_bsse_gradients_result_available():
    """The gradients result is offered for the 'bsse' calculation."""
    assert "bsse" in orca_step.metadata["results"]["gradients"]["calculation"]


def test_bsse_extra_results_registered():
    """The uncorrected energy and BSSE correction are 'bsse'-only results with
    registered property templates."""
    results = orca_step.metadata["results"]
    for key in ("uncorrected energy", "bsse correction"):
        assert results[key]["calculation"] == ["bsse"]
    import importlib.resources

    csv = (
        importlib.resources.files("orca_step") / "data" / "properties.csv"
    ).read_text()
    assert "uncorrected energy#ORCA#{model}" in csv
    assert "BSSE correction#ORCA#{model}" in csv


def test_bsse_parse_compound_energies(tmp_path):
    """The five sub-calculation totals come from each COMPOUND JOB block's LAST
    FINAL SINGLE POINT ENERGY (so an optimized monomer's converged value wins)."""
    (tmp_path / "orca.out").write_text(
        "COMPOUND JOB MainJOB\n"
        "COMPOUND JOB  1\nFINAL SINGLE POINT ENERGY      -76.100000\n"  # fragA(AB)
        "COMPOUND JOB  2\n"  # monomer A -- optimized, two cycles; last wins
        "FINAL SINGLE POINT ENERGY      -76.040000\n"
        "FINAL SINGLE POINT ENERGY      -76.050000\n"
        "COMPOUND JOB  3\nFINAL SINGLE POINT ENERGY      -76.110000\n"  # fragB(AB)
        "COMPOUND JOB  4\nFINAL SINGLE POINT ENERGY      -76.060000\n"  # monomer B
        "COMPOUND JOB  5\nFINAL SINGLE POINT ENERGY     -152.220000\n"  # dimer
        "COMPOUND: System command to be executed: rm ...\n"
    )
    node = orca_step.BSSE()
    e = node._parse_compound_energies(tmp_path)
    assert e == [-76.10, -76.05, -76.11, -76.06, -152.22]
    # E_CP = tot - (fragA - monA) - (fragB - monB); the in-dimer-basis fragments
    # are lower, so the correction raises the energy (BSSE removed).
    corrected = e[4] - (e[0] - e[1]) - (e[2] - e[3])
    assert corrected == pytest.approx(-152.22 + 0.05 + 0.05)
    # Fewer than five jobs -> None (cannot trust the energy).
    (tmp_path / "short.out").write_text("COMPOUND JOB  1\n")
    assert node._parse_compound_energies(tmp_path.parent / "missing") is None


def test_bsse_parse_indices():
    """1-based index/range lists parse to unique 0-based indices; out-of-range
    is an error."""
    assert orca_step.BSSE._parse_indices("1-3, 5 7", 8) == [0, 1, 2, 4, 6]
    assert orca_step.BSSE._parse_indices("2 2 3", 5) == [1, 2]  # de-duplicated
    with pytest.raises(RuntimeError):
        orca_step.BSSE._parse_indices("9", 8)


def test_bsse_compound_input():
    """The %Compound block injects the level of theory and options via `with`."""
    node = orca_step.BSSE()
    P = {
        "use model chemistry": "no",
        "method": "B3LYP",
        "basis": "def2-SVP",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "grid": "DEFGRID3",
        "scf convergence": "TIGHTSCF",
        "extra keywords": "D3BJ",
        "optimize monomers": "no",
    }
    block, method, basis = node._compound_input(P, "bsse.xyz", "bssegradient.cmp")
    assert method == "B3LYP" and basis == "def2-SVP"
    assert '%Compound "bssegradient.cmp"' in block
    assert "with" in block and block.rstrip().endswith("end")
    assert 'molecule       = "bsse.xyz";' in block
    assert 'method         = "B3LYP";' in block
    assert 'basis          = "def2-SVP";' in block
    # aux 'none' is omitted; grid, scf, and extra flow into restOfInput.
    assert 'restOfInput    = "DEFGRID3 TIGHTSCF D3BJ";' in block
    assert "DoOptimization = false;" in block


def test_bsse_compound_input_optimize():
    """'optimize monomers' maps to DoOptimization = true."""
    node = orca_step.BSSE()
    P = {
        "use model chemistry": "no",
        "method": "HF",
        "basis": "6-31G*",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "grid": "default",
        "scf convergence": "default",
        "extra keywords": "",
        "optimize monomers": "yes",
    }
    block, _, _ = node._compound_input(P, "bsse.xyz", "bssenergy.cmp")
    assert '%Compound "bssenergy.cmp"' in block
    assert "DoOptimization = true;" in block
    assert 'restOfInput    = "";' in block


def test_bsse_energy_only_script_shipped():
    """The gradient-free Compound script is shipped and requests no EnGrad."""
    import importlib.resources

    text = (
        importlib.resources.files("orca_step") / "data" / "bssenergy.cmp"
    ).read_text()
    # Look only at the active directives, not the explanatory comments.
    active = "\n".join(
        ln for ln in text.splitlines() if not ln.lstrip().startswith("#")
    )
    assert "CreateBSSE" in active  # still splits the fragments
    assert "EnGrad" not in active  # but never asks for a gradient
    assert "Nuclear_Gradient" not in active


def _neutral_singlet():
    """A minimal stand-in configuration for _check_supported."""

    class _Cfg:
        charge = 0
        spin_multiplicity = 1

    return _Cfg()


def test_bsse_rejects_numeric_gradient_method():
    """With the gradient requested, (DLPNO-)CCSD(T) is refused (numerical
    gradient only). Double hybrids and MP2 (analytic) are NOT refused."""
    node = orca_step.BSSE()
    P = {
        "use model chemistry": "no",
        "method": "DLPNO-CCSD(T)",
        "basis": "def2-SVP",
        "basis source": "ORCA internal",
        "auxiliary basis": "AutoAux",
        "extra keywords": "",
        "compute gradient": "yes",
    }
    with pytest.raises(RuntimeError, match="analytic gradient"):
        node._check_supported(P, _neutral_singlet())


def test_bsse_energy_only_allows_ccsdt():
    """Energy-only lifts the analytic-gradient requirement, so CCSD(T) is
    accepted (for gold-standard counterpoise interaction energies)."""
    node = orca_step.BSSE()
    P = {
        "use model chemistry": "no",
        "method": "DLPNO-CCSD(T)",
        "basis": "def2-SVP",
        "basis source": "ORCA internal",
        "auxiliary basis": "AutoAux",
        "extra keywords": "",
        "compute gradient": "no",
    }
    # Should not raise.
    node._check_supported(P, _neutral_singlet())


def test_bsse_rejects_bse_basis():
    """BSSE (Phase 1) refuses a Basis Set Exchange basis."""
    node = orca_step.BSSE()
    P = {
        "use model chemistry": "no",
        "method": "HF",
        "basis": "bse:cc-pVDZ",
        "basis source": "ORCA internal",
        "auxiliary basis": "none",
        "extra keywords": "",
    }
    with pytest.raises(RuntimeError, match="Basis Set Exchange"):
        node._check_supported(P, None)


def test_read_engrad(tmp_path):
    """The Compound EnGrad file yields (energy, gradient); blank-line separated
    like the script writes it."""
    (tmp_path / "result.engrad").write_text(
        "\n\n\n 2\n\n\n\n   -1.5\n\n\n\n"
        "   0.1\n   0.2\n   0.3\n  -0.1\n  -0.2\n  -0.3\n"
    )
    energy, grad = orca_step.BSSE._read_engrad(tmp_path / "result.engrad")
    assert energy == -1.5
    assert grad == [[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]]
    # Missing file -> (None, None).
    assert orca_step.BSSE._read_engrad(tmp_path / "missing.engrad") == (None, None)
