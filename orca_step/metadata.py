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
    "M06-2X": {"type": "DFT", "description": "M06-2X meta-hybrid functional"},
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

"""Results ORCA can produce. Same recognized fields as the other QM steps."""
metadata["results"] = {
    "energy": {
        "description": "The total energy",
        "dimensionality": "scalar",
        "property": "total energy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
    },
}

# Placeholder for the model-chemistry protocol; populated in the Model Chemistry
# integration phase via get_model_chemistry_options() on ORCAStep.
metadata["computational models"] = {}
