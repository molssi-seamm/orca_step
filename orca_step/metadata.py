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
# The top-level method. For DFT the specific functional is chosen separately
# (see metadata["functionals"]); "DFT" here is just the Kohn-Sham method. The
# "gradients" key records whether ORCA has an ANALYTIC nuclear gradient for the
# method ("analytic") or whether a numerical gradient is required ("numeric").
# For "DFT" the gradient availability depends on the functional, so it is looked
# up per functional instead.
metadata["methods"] = {
    "HF": {"type": "HF", "description": "Hartree-Fock", "gradients": "analytic"},
    "MP2": {
        "type": "MP2",
        "description": "second-order Moller-Plesset (MP2)",
        "gradients": "analytic",
    },
    "RI-MP2": {
        "type": "MP2",
        "description": "RI-MP2 (resolution of identity)",
        "gradients": "analytic",
    },
    "CCSD(T)": {
        "type": "QC",
        "description": "coupled cluster CCSD(T)",
        # No analytic gradient for the (T) triples -> ORCA needs NUMGRAD.
        "gradients": "numeric",
    },
    "DLPNO-CCSD(T)": {
        "type": "QC",
        "description": "near-linear-scaling DLPNO-CCSD(T)",
        "needs_aux": True,
        # No analytic gradient for the (T) triples -> ORCA needs NUMGRAD.
        "gradients": "numeric",
    },
    "CCSD(T)-F12D/RI": {
        "type": "QC",
        "description": "explicitly-correlated CCSD(T)-F12D (needs an F12 basis)",
        "needs_aux": True,
        # F12 + (T) has no analytic gradient in ORCA. The '/RI' is required --
        # canonical F12 uses the RI approximation for the F12 integrals; ORCA
        # rejects the bare 'CCSD(T)-F12D'. (DLPNO-CCSD(T)-F12D implies RI, and
        # ORCA conversely rejects a '/RI' on it.)
        "gradients": "numeric",
    },
    "DLPNO-CCSD(T)-F12D": {
        "type": "QC",
        "description": "DLPNO-CCSD(T)-F12D (explicitly correlated; F12 basis)",
        "needs_aux": True,
        "gradients": "numeric",
    },
    "DFT": {
        "type": "DFT",
        "description": "Kohn-Sham density functional theory (choose a functional)",
    },
}

# Curated basis-set families (ORCA's built-in names). Free text is also allowed
# in the GUI; these are the guided choices. Ordering matters: it is what the user
# sees in the dropdown. They are grouped by family, and within a family split into
# separate "ladders" -- valence, then polarization-augmented, then
# diffuse-augmented -- each ordered by zeta rung (DZ -> TZ -> QZ -> 5Z), so a
# sensible progression is a single ladder read top to bottom. Keep this order.
metadata["basis sets"] = [
    # --- Minimal / small ---
    "STO-3G",
    "3-21G",
    # --- Pople (6-31G = DZ, 6-311G = TZ) ---
    # valence
    "6-31G",
    "6-311G",
    # + polarization
    "6-31G*",
    "6-31G**",
    "6-311G*",
    "6-311G**",
    # + diffuse (and polarization)
    "6-31+G*",
    "6-31+G**",
    "6-31++G**",
    "6-311+G*",
    "6-311++G**",
    "6-311++G(2d,2p)",
    "6-311++G(3df,3pd)",
    # --- Dunning correlation-consistent ---
    # valence
    "cc-pVDZ",
    "cc-pVTZ",
    "cc-pVQZ",
    "cc-pV5Z",
    # + diffuse (augmented)
    "aug-cc-pVDZ",
    "aug-cc-pVTZ",
    "aug-cc-pVQZ",
    "aug-cc-pV5Z",
    # core-valence
    "cc-pCVDZ",
    "cc-pCVTZ",
    "cc-pCVQZ",
    # weighted core-valence
    "cc-pwCVDZ",
    "cc-pwCVTZ",
    "cc-pwCVQZ",
    # explicitly-correlated (F12): use ONLY with an F12 method; the step adds the
    # matching '<basis>-CABS' complementary basis automatically.
    "cc-pVDZ-F12",
    "cc-pVTZ-F12",
    "cc-pVQZ-F12",
    # --- Karlsruhe def2 (SVP = DZ, TZVP = TZ, QZVP = QZ) ---
    # valence (polarization is built in; increasing within each zeta)
    "def2-SV(P)",
    "def2-SVP",
    "def2-TZVP(-f)",
    "def2-TZVP",
    "def2-TZVPP",
    "def2-QZVP",
    "def2-QZVPP",
    # + diffuse (the 'D' sets; property/anion work)
    "def2-SVPD",
    "def2-TZVPD",
    "def2-TZVPPD",
    "def2-QZVPD",
    "def2-QZVPPD",
    # + diffuse, minimally augmented (Truhlar 'ma-')
    "ma-def2-SVP",
    "ma-def2-TZVP",
    "ma-def2-TZVPP",
    "ma-def2-QZVPP",
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
        "property": "SCF energy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
    },
    "mp2 energy": {
        "description": "The MP2 total energy",
        "dimensionality": "scalar",
        "property": "MP2 energy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
        "methods": ["MP2", "RI-MP2"],
    },
    "ccsd energy": {
        "description": "The CCSD total energy",
        "dimensionality": "scalar",
        "property": "CCSD energy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
        "methods": ["CCSD", "CCSD(T)", "DLPNO-CCSD", "DLPNO-CCSD(T)"],
    },
    "ccsd(t) energy": {
        "description": "The CCSD(T) total energy",
        "dimensionality": "scalar",
        "property": "CCSD(T) energy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
        "methods": ["CCSD(T)", "DLPNO-CCSD(T)"],
    },
    "HOMO energy": {
        "description": "Energy of the highest occupied molecular orbital",
        "dimensionality": "scalar",
        "property": "HOMO energy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
    },
    "LUMO energy": {
        "description": "Energy of the lowest unoccupied molecular orbital",
        "dimensionality": "scalar",
        "property": "LUMO energy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
    },
    "HOMO-LUMO gap": {
        "description": "The HOMO-LUMO gap",
        "dimensionality": "scalar",
        "property": "HOMO-LUMO gap#ORCA#{model}",
        "type": "float",
        "units": "eV",
    },
    "dipole moment": {
        "description": "The dipole moment vector",
        "dimensionality": [3],
        "property": "dipole moment vector#ORCA#{model}",
        "type": "float",
        "units": "debye",
    },
    "dipole magnitude": {
        "description": "The magnitude of the dipole moment",
        "dimensionality": "scalar",
        "property": "dipole magnitude#ORCA#{model}",
        "type": "float",
        "units": "debye",
    },
    "rotational constants": {
        "description": "The rotational constants A, B, C",
        "dimensionality": [3],
        "property": "rotational constants#ORCA#{model}",
        "type": "float",
        "units": "GHz",
    },
    "S^2": {
        "description": "Expectation value of the total spin operator <S^2>",
        "dimensionality": "scalar",
        "property": "S^2#ORCA#{model}",
        "type": "float",
    },
    "isotropic polarizability": {
        "description": "The isotropic dipole polarizability",
        "dimensionality": "scalar",
        "property": "isotropic polarizability#ORCA#{model}",
        "type": "float",
        "units": "a.u.",
        "methods": ["HF", "DFT"],
    },
    "mulliken charges": {
        "description": "The Mulliken atomic charges",
        "dimensionality": ["n_atoms"],
        "property": "Mulliken charges#ORCA#{model}",
        "type": "float",
        "units": "e",
    },
    "löwdin charges": {
        "description": "The Löwdin atomic charges",
        "dimensionality": ["n_atoms"],
        "property": "Löwdin charges#ORCA#{model}",
        "type": "float",
        "units": "e",
    },
    "hirshfeld charges": {
        "description": "The Hirshfeld atomic charges",
        "dimensionality": ["n_atoms"],
        "property": "Hirshfeld charges#ORCA#{model}",
        "type": "float",
        "units": "e",
    },
    "mayer valences": {
        "description": "The Mayer total valence of each atom",
        "dimensionality": ["n_atoms"],
        "property": "Mayer valences#ORCA#{model}",
        "type": "float",
    },
    "gradients": {
        "description": "The gradient of the energy on each atom",
        "dimensionality": ["n_atoms", 3],
        "calculation": ["energy", "optimization", "bsse"],
        "property": "gradients#ORCA#{model}",
        "type": "float",
        "units": "E_h/bohr",
    },
    "uncorrected energy": {
        "description": "The uncorrected (raw) total energy of the complex",
        "dimensionality": "scalar",
        "calculation": ["bsse"],
        "property": "uncorrected energy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
    },
    "bsse correction": {
        "description": "The counterpoise (BSSE) correction to the energy "
        "(corrected minus uncorrected)",
        "dimensionality": "scalar",
        "calculation": ["bsse"],
        "property": "BSSE correction#ORCA#{model}",
        "type": "float",
        "units": "E_h",
    },
    "frequencies": {
        "description": "The harmonic vibrational frequencies",
        "dimensionality": ["n_dof"],
        "calculation": ["frequencies"],
        "type": "float",
        "units": "1/cm",
    },
    "IR intensities": {
        "description": "The IR intensities of the vibrational modes",
        "dimensionality": ["n_dof"],
        "calculation": ["frequencies"],
        "type": "float",
        "units": "km/mol",
    },
    "n imaginary frequencies": {
        "description": "The number of imaginary vibrational frequencies",
        "dimensionality": "scalar",
        "calculation": ["frequencies"],
        "property": "number of imaginary frequencies#ORCA#{model}",
        "type": "integer",
    },
    "zero point energy": {
        "description": "The zero-point vibrational energy",
        "dimensionality": "scalar",
        "calculation": ["frequencies"],
        "property": "zero point energy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
    },
    "enthalpy": {
        "description": "The total enthalpy (H) from the thermochemistry",
        "dimensionality": "scalar",
        "calculation": ["frequencies"],
        "property": "enthalpy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
    },
    "gibbs energy": {
        "description": "The Gibbs free energy (G) from the thermochemistry",
        "dimensionality": "scalar",
        "calculation": ["frequencies"],
        "property": "Gibbs energy#ORCA#{model}",
        "type": "float",
        "units": "E_h",
    },
}

# Placeholder for the model-chemistry protocol; populated in the Model Chemistry
# integration phase via get_model_chemistry_options() on ORCAStep.
metadata["computational models"] = {}

# The DFT functionals ORCA supports, keyed by the exact ORCA '!' keyword. Each
# record carries:
#   category   -- ORCA's own classification (see metadata["functional categories"])
#   gradients  -- "analytic" if ORCA has an analytic nuclear gradient, or
#                 "numeric" if only a numerical gradient (NUMGRAD) is possible
#                 (e.g. wB97M(2)/wB97X-2, which are evaluated non-self-consistently)
#   citations  -- reference keys mined from the ORCA manual; the bibtex itself
#                 lives in data/references.bib (loaded into self._bibliography)
#   note       -- optional short annotation (composite "3c" method, %HFX, ...)
# Categorization and gradient availability follow the ORCA 6.1 manual's DFT
# model-chemistries page.
metadata["functional categories"] = [
    "local",
    "GGA",
    "meta-GGA",
    "global hybrid",
    "range-separated hybrid",
    "global double-hybrid",
    "range-separated double-hybrid",
]

metadata["functionals"] = {
    # --- local ---
    "PWLDA": {
        "category": "local",
        "gradients": "analytic",
        "citations": ["orca_dft_id724"],
    },
    "VWN": {
        "category": "local",
        "gradients": "analytic",
        "citations": ["orca_dft_id924"],
        "note": "defaults to VWN5",
    },
    "VWN3": {
        "category": "local",
        "gradients": "analytic",
        "citations": ["orca_dft_id924"],
    },
    "VWN5": {
        "category": "local",
        "gradients": "analytic",
        "citations": ["orca_dft_id924"],
    },
    # --- GGA ---
    "B97-3C": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id348"],
        "note": "composite 3c method",
    },
    "BLYP": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id109"],
    },
    "BP86": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id87", "orca_dft_id718"],
    },
    "GLYP": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id317"],
    },
    "MPWLYP": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id77"],
    },
    "MPWPW": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id77"],
    },
    "OLYP": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id389"],
    },
    "PBE": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id725"],
    },
    "PW91": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id721"],
    },
    "PWP": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id721"],
    },
    "REVPBE": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id964"],
    },
    "RPBE": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id382"],
    },
    "RPW86PBE": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id794"],
    },
    "XLYP": {
        "category": "GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id322"],
    },
    # --- meta-GGA ---
    "B97M-D3BJ": {
        "category": "meta-GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id327"],
    },
    "B97M-D4": {
        "category": "meta-GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id650"],
    },
    "B97M-V": {
        "category": "meta-GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id409"],
    },
    "M06L": {
        "category": "meta-GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id586"],
    },
    "R2SCAN": {
        "category": "meta-GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id297"],
    },
    "R2SCAN-3C": {
        "category": "meta-GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id354"],
        "note": "composite 3c method",
    },
    "REVTPSS": {
        "category": "meta-GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id775", "orca_dft_id776"],
    },
    "RSCAN": {
        "category": "meta-GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id79"],
    },
    "SCANFUNC": {
        "category": "meta-GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id728"],
        "note": "SCAN",
    },
    "TPSS": {
        "category": "meta-GGA",
        "gradients": "analytic",
        "citations": ["orca_dft_id899"],
    },
    # --- global hybrid ---
    "B1LYP": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id7"],
    },
    "B1P86": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id718"],
    },
    "B3LYP": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id65"],
    },
    "B3LYP-3C": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id355"],
        "note": "composite 3c method",
    },
    "B3LYP/G": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id65"],
        "note": "Gaussian B3LYP variant",
    },
    "B3P86": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id88"],
    },
    "B3PW91": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id88"],
    },
    "BHANDHLYP": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id88"],
        "note": "50% HFX",
    },
    "M06": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id585"],
        "note": "meta-GGA hybrid",
    },
    "M062X": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id585"],
        "note": "meta-GGA hybrid, 54% HFX",
    },
    "MPW1LYP": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id77"],
    },
    "MPW1PW": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id77"],
    },
    "O3LYP": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id388"],
    },
    "PBE0": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id8"],
    },
    "PBEH-3C": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id362"],
        "note": "composite 3c method, 42% HFX",
    },
    "PW1PW": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id8"],
    },
    "PW6B95": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id965"],
        "note": "meta-GGA hybrid",
    },
    "R2SCAN0": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id123"],
        "note": "25% HFX",
    },
    "R2SCAN50": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id123"],
        "note": "50% HFX",
    },
    "R2SCANH": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id123"],
        "note": "10% HFX",
    },
    "REVPBE0": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id342"],
    },
    "REVPBE38": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id342"],
        "note": "37.5% HFX",
    },
    "TPSS0": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id338"],
        "note": "25% HFX",
    },
    "TPSSH": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id899"],
        "note": "10% HFX",
    },
    "X3LYP": {
        "category": "global hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id322"],
    },
    # --- range-separated hybrid ---
    "CAM-B3LYP": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id390"],
    },
    "LC-BLYP": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id438"],
    },
    "LC-PBE": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id437"],
    },
    "WB97": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id407"],
    },
    "WB97M-D3BJ": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id327"],
        "note": "meta-GGA range-separated hybrid",
    },
    "WB97M-D4": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id650"],
        "note": "meta-GGA range-separated hybrid",
    },
    "WB97M-D4REV": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id358"],
        "note": "meta-GGA range-separated hybrid",
    },
    "WB97M-V": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id410"],
        "note": "meta-GGA range-separated hybrid",
    },
    "WB97X": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id407"],
    },
    "WB97X-3C": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id358"],
        "note": "composite 3c method",
    },
    "WB97X-D3": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id152"],
    },
    "WB97X-D3BJ": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id327"],
    },
    "WB97X-D4": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id650"],
    },
    "WB97X-D4REV": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id358"],
    },
    "WB97X-V": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id408"],
    },
    "WR2SCAN": {
        "category": "range-separated hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id972"],
        "note": "omega-r2SCAN",
    },
    # --- global double-hybrid ---
    "B2GP-PLYP": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id62"],
    },
    "B2K-PLYP": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id63"],
    },
    "B2NC-PLYP": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id64"],
    },
    "B2PLYP": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id339"],
    },
    "B2T-PLYP": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id63"],
    },
    "DSD-BLYP": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id600"],
    },
    "DSD-BLYP/2013": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id602"],
    },
    "DSD-PBEB95": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id602"],
    },
    "DSD-PBEP86": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id601"],
    },
    "DSD-PBEP86/2013": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id602"],
    },
    "KPR2SCAN": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id972"],
        "note": "kappa-Pr2SCAN50 (manual keyword KPR2SCAN50)",
    },
    "mPW2PLYP": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id832"],
    },
    "PBE-QIDH": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id712"],
    },
    "PBE0-DH": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id711"],
    },
    "PR2SCAN50": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id972"],
        "note": "Pr2SCAN50",
    },
    "PR2SCAN69": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id972"],
        "note": "Pr2SCAN69",
    },
    "PWPB95": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id325"],
    },
    "R2SCAN-CIDH": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id972"],
    },
    "R2SCAN-QIDH": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id972"],
    },
    "R2SCAN0-2": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id972"],
    },
    "R2SCAN0-DH": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id972"],
    },
    "REVDOD-PBEP86-D4/2021": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id603"],
        "note": "opposite-spin DOD variant",
    },
    "REVDOD-PBEP86/2021": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id603"],
        "note": "opposite-spin DOD variant",
    },
    "REVDSD-PBEP86-D4/2021": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id603"],
    },
    "REVDSD-PBEP86/2021": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id603"],
    },
    "SCS-B2GP-PLYP21": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SCS-PBE-QIDH": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SCS/SOS-B2PLYP21": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SOS-B2GP-PLYP21": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SOS-PBE-QIDH": {
        "category": "global double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    # --- range-separated double-hybrid ---
    "RSX-0DH": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id795"],
    },
    "RSX-QIDH": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id796"],
    },
    "SCS-RSX-QIDH": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SCS-WB2GP-PLYP": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SCS-WB88PP86": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SCS-WPBEPP86": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SCS/SOS-WB2PLYP": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SOS-RSX-QIDH": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SOS-WB2GP-PLYP": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SOS-WB88PP86": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "SOS-WPBEPP86": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id835"],
        "note": "excited-state optimized",
    },
    "WB2GP-PLYP": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id328"],
    },
    "WB2PLYP": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id328"],
    },
    "WB88PP86": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id328"],
    },
    "WB97M(2)": {
        "category": "range-separated double-hybrid",
        "gradients": "numeric",
        "citations": ["orca_dft_id594"],
        "note": "non-self-consistent; no analytic gradient/density, NUMGRAD only",
    },
    "WB97X-2": {
        "category": "range-separated double-hybrid",
        "gradients": "numeric",
        "citations": ["orca_dft_id151"],
        "note": "non-self-consistent; no analytic gradient, NUMGRAD only",
    },
    "WPBEPP86": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id328"],
    },
    "WPR2SCAN50": {
        "category": "range-separated double-hybrid",
        "gradients": "analytic",
        "citations": ["orca_dft_id972"],
        "note": "omega-Pr2SCAN50",
    },
}

# Backward-compatible view used by the citation code: functional keyword -> list
# of reference keys. Derived from metadata["functionals"] so there is a single
# source of truth (the old data/dft_functionals.json is no longer needed).
metadata["dft functionals"] = {
    name: rec["citations"] for name, rec in metadata["functionals"].items()
}
