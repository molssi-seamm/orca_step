# -*- coding: utf-8 -*-

"""Metadata describing ORCA's methods, basis sets, and results.

This drives the explicit method/basis GUI (like the Gaussian step) and will also
back the model-chemistry protocol (``get_model_chemistry_options``) used by the
Model Chemistry step. It is kept deliberately small for now and will be filled
out as more of ORCA's capabilities are exposed.
"""

metadata = {}

# Methods offered in the explicit GUI. Each is an ORCA "!" keyword. "type"
# groups them for the model-chemistry protocol (QC = wavefunction, DFT, etc.).
metadata["methods"] = {
    "HF": {"type": "HF", "description": "Hartree-Fock"},
    "MP2": {"type": "MP2", "description": "second-order Moller-Plesset (MP2)"},
    "RI-MP2": {"type": "MP2", "description": "RI-MP2 (resolution of identity)"},
    "CCSD(T)": {"type": "QC", "description": "coupled cluster CCSD(T)"},
    "DLPNO-CCSD(T)": {
        "type": "QC",
        "description": "near-linear-scaling DLPNO-CCSD(T)",
        "needs_aux": True,
    },
    "B3LYP": {"type": "DFT", "description": "B3LYP hybrid functional"},
    "PBE": {"type": "DFT", "description": "PBE GGA functional"},
    "PBE0": {"type": "DFT", "description": "PBE0 hybrid functional"},
    "wB97X-D3": {
        "type": "DFT",
        "description": "range-separated wB97X-D3 functional",
    },
    "M062X": {"type": "DFT", "description": "M06-2X meta-hybrid functional"},
}

# Curated basis-set families (ORCA's built-in names). Free text is also allowed
# in the GUI; these are the guided choices.
metadata["basis sets"] = [
    # Pople
    "6-31G",
    "6-31G*",
    "6-31G**",
    "6-31+G*",
    "6-311G**",
    "6-311++G**",
    # Dunning correlation-consistent
    "cc-pVDZ",
    "cc-pVTZ",
    "cc-pVQZ",
    "aug-cc-pVDZ",
    "aug-cc-pVTZ",
    # Karlsruhe def2
    "def2-SVP",
    "def2-TZVP",
    "def2-TZVPP",
    "def2-QZVP",
]

# Auxiliary / fitting basis choices. AutoAux generates a fitting basis for any
# method/basis and is the robust default for correlated (DLPNO/MP2) methods.
metadata["auxiliary basis sets"] = [
    "AutoAux",
    "none",
    "def2/J",
    "def2/JK",
]

"""Results ORCA can produce. Same recognized fields as the other QM steps.

`methods` (when present) gates the result to certain levels of theory.
"""
metadata["results"] = {
    "energy": {
        "description": "The total (final single point) energy",
        "dimensionality": "scalar",
        "property": "total energy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
    },
    "scf energy": {
        "description": "The SCF (HF/DFT reference) energy",
        "dimensionality": "scalar",
        "type": "float",
        "units": "E_h",
    },
    "mp2 energy": {
        "description": "The MP2 total energy",
        "dimensionality": "scalar",
        "type": "float",
        "units": "E_h",
        "methods": ["MP2", "RI-MP2"],
    },
    "ccsd energy": {
        "description": "The CCSD total energy",
        "dimensionality": "scalar",
        "type": "float",
        "units": "E_h",
        "methods": ["CCSD", "CCSD(T)", "DLPNO-CCSD", "DLPNO-CCSD(T)"],
    },
    "ccsd(t) energy": {
        "description": "The CCSD(T) total energy",
        "dimensionality": "scalar",
        "type": "float",
        "units": "E_h",
        "methods": ["CCSD(T)", "DLPNO-CCSD(T)"],
    },
    "HOMO energy": {
        "description": "Energy of the highest occupied molecular orbital",
        "dimensionality": "scalar",
        "type": "float",
        "units": "E_h",
    },
    "LUMO energy": {
        "description": "Energy of the lowest unoccupied molecular orbital",
        "dimensionality": "scalar",
        "type": "float",
        "units": "E_h",
    },
    "HOMO-LUMO gap": {
        "description": "The HOMO-LUMO gap",
        "dimensionality": "scalar",
        "type": "float",
        "units": "eV",
    },
    "dipole moment": {
        "description": "The dipole moment vector",
        "dimensionality": [3],
        "type": "float",
        "units": "debye",
    },
    "dipole magnitude": {
        "description": "The magnitude of the dipole moment",
        "dimensionality": "scalar",
        "type": "float",
        "units": "debye",
    },
    "rotational constants": {
        "description": "The rotational constants A, B, C",
        "dimensionality": [3],
        "type": "float",
        "units": "GHz",
    },
    "S^2": {
        "description": "Expectation value of the total spin operator <S^2>",
        "dimensionality": "scalar",
        "type": "float",
    },
    "isotropic polarizability": {
        "description": "The isotropic dipole polarizability",
        "dimensionality": "scalar",
        "type": "float",
        "units": "a.u.",
        "methods": ["HF", "B3LYP", "PBE", "PBE0", "wB97X-D3", "M062X"],
    },
    "mulliken charges": {
        "description": "The Mulliken atomic charges",
        "dimensionality": ["n_atoms"],
        "type": "float",
        "units": "e",
    },
    "löwdin charges": {
        "description": "The Löwdin atomic charges",
        "dimensionality": ["n_atoms"],
        "type": "float",
        "units": "e",
    },
    "hirshfeld charges": {
        "description": "The Hirshfeld atomic charges",
        "dimensionality": ["n_atoms"],
        "type": "float",
        "units": "e",
    },
    "mayer valences": {
        "description": "The Mayer total valence of each atom",
        "dimensionality": ["n_atoms"],
        "type": "float",
    },
    "gradients": {
        "description": "The gradient of the energy on each atom",
        "dimensionality": ["n_atoms", 3],
        "calculation": ["energy", "optimization"],
        "type": "float",
        "units": "E_h/bohr",
    },
}

# Placeholder for the model-chemistry protocol; populated in the Model Chemistry
# integration phase via get_model_chemistry_options() on ORCAStep.
metadata["computational models"] = {}

# Map of DFT functional keyword -> list of citation keys (mined from the ORCA
# manual). The bibtex itself lives in data/references.bib (loaded into
# self._bibliography); this holds only the keys -- no bibtex duplication.
import importlib.resources  # noqa: E402
import json  # noqa: E402

_dft_file = importlib.resources.files("orca_step") / "data" / "dft_functionals.json"
try:
    metadata["dft functionals"] = json.loads(_dft_file.read_text())
except Exception:  # pragma: no cover - missing/!json data is non-fatal
    metadata["dft functionals"] = {}
