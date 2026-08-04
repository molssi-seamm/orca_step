#!/usr/bin/env python
"""M3 sanity check: charged-fragment counterpoise on Na+/Cl-/H2O two-body
pairs, real ORCA, compared to approximate literature binding energies.

This is a magnitude-level sanity check of the per-fragment-charge wiring
(seamm_bsse.Fragment with independent charge, validated in M2's unit tests
only with stubbed ORCA) on real charged systems -- NOT a rigorous benchmark.
Na+...H2O and Cl-...H2O geometries are literature-informed (typical M-O /
Cl...H distances and orientations) but NOT optimized at this level of
theory, so a few kcal/mol mismatch vs. literature is expected. Na+...Cl-
points come from the already-built 40-point log-spaced R-scan
(~/Dropbox/GM/TrainingData/Molecular/Water/Na-Cl.sdf).

Run inside the seamm-dev environment:
    python validate_bsse_m3.py [workdir]
"""
import math
import re
import sys
from pathlib import Path

import molsystem
import orca_step
import seamm
import seamm_exec
from seamm.variables import Variables

seamm.flowchart_variables = Variables()

_HARTREE_TO_KCAL = 627.509474

WORKDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/bsse_m3_validate")
WORKDIR.mkdir(parents=True, exist_ok=True)

executor = seamm_exec.get_executor("local")

# A fast, dispersion-corrected hybrid -- good enough for a magnitude sanity
# check, not the production revDSD-PBEP86-D4/def2-TZVPPD (too slow to run
# ad hoc here). CP throughout: this is exactly what the correction is for.
P = {
    "use model chemistry": "no",
    "method": "DFT",
    "functional": "B3LYP",
    "basis": "def2-TZVP",
    "basis source": "ORCA internal",
    "auxiliary basis": "AutoAux",
    "grid": "default",
    "scf convergence": "default",
    "extra keywords": "D3BJ",
    "basis set extrapolation": "none",
    "save wavefunction": "no",
    "compute gradient": "no",
    "optimize monomers": "no",
    # "specified", not "auto (molecules)": these configurations are built by
    # hand (atoms.append only, no bond table), so find_molecules() would see
    # every atom as isolated. Each call below supplies the right groups.
    "fragments": "specified",
    "fragment atoms": "",
    "fragment charges": "",
}


class _Parent:
    options = {"ncores": "1", "memory": "available", "max-atoms-to-print": 25}
    global_options = {"root": str(Path.home() / "SEAMM"), "ncores": "1"}


class _Flowchart:
    def __init__(self, root_directory, executor):
        self.root_directory = root_directory
        self.executor = executor


def new_configuration(db_name):
    db = molsystem.SystemDB(filename=f"file:{db_name}?mode=memory&cache=shared")
    system = db.create_system(name=db_name)
    return system.create_configuration(name=db_name)


def run_bsse(label, configuration, fragment_atoms, fragment_charges):
    node = orca_step.BSSE()
    node._id = (label,)
    node.flowchart = _Flowchart(str(WORKDIR), executor)
    node.parent = _Parent()
    node.get_system_configuration = lambda arg: (None, configuration)
    node._cite_references = lambda P: None
    node._cite_bsse = lambda: None
    node.next = lambda: None

    this_P = dict(
        P, **{"fragment atoms": fragment_atoms, "fragment charges": fragment_charges}
    )
    for key, value in this_P.items():
        if key in node.parameters:
            node.parameters[key].value = value

    captured = {}
    node.analyze = lambda **kwargs: captured.update(kwargs)
    node.run()
    return captured["data"]


# ---------------------------------------------------------------------
# Na+...Cl-: points from the already-built 40-point log-spaced R-scan.
# ---------------------------------------------------------------------
print("=" * 70)
print("Na+ ... Cl-  (points from Na-Cl.sdf)")
print("=" * 70)

nacl_sdf = Path("/Users/psaxe/Dropbox/GM/TrainingData/Molecular/Water/Na-Cl.sdf")
text = nacl_sdf.read_text()
records = text.split("$$$$")
r_values = []
for rec in records:
    m = re.search(r"dimer separation.*?\n([0-9.]+)", rec, re.S)
    if m:
        r_values.append(float(m.group(1)))

# Bracket the literature R_e (~2.36 A) plus a couple of longer-range points
# to see the interaction decay sensibly.
indices = [7, 8, 9, 10, 16, 24]
print(f"{'R (A)':>8s} {'CP interaction (kcal/mol)':>28s}")
for i in indices:
    r = r_values[i]
    configuration = new_configuration(f"nacl_{i}")
    configuration.atoms.append(symbol=["Na", "Cl"], x=[0.0, 0.0], y=[0.0, 0.0], z=[0.0, r])
    configuration.charge = 0
    configuration.spin_multiplicity = 1
    data = run_bsse(f"nacl_{i}", configuration, "1; 2", "1, -1")
    kcal = data["interaction energy"] / 4.184  # kJ/mol -> kcal/mol
    print(f"{r:8.3f} {kcal:28.3f}")

print(
    "\nApprox. literature reference (Born-Haber from D0(NaCl)=97.5, "
    "IP(Na)=118.5, EA(Cl)=83.3 kcal/mol): De(Na+ + Cl- -> NaCl) ~ -133 kcal/mol "
    "at R_e ~ 2.36 A."
)

# ---------------------------------------------------------------------
# Na+...H2O and Cl-...H2O: literature-informed (NOT optimized) geometries.
# ---------------------------------------------------------------------
r_oh = 0.9572
half_angle = math.radians(104.5 / 2.0)
h1 = (r_oh * math.sin(half_angle), 0.0, r_oh * math.cos(half_angle))
h2 = (-r_oh * math.sin(half_angle), 0.0, r_oh * math.cos(half_angle))
o = (0.0, 0.0, 0.0)

print()
print("=" * 70)
print("Na+ ... H2O  (Na+ on the O lone-pair side, Na-O = 2.30 A)")
print("=" * 70)
na_o_distance = 2.30
na = (0.0, 0.0, -na_o_distance)
configuration = new_configuration("na_water")
configuration.atoms.append(
    symbol=["O", "H", "H", "Na"],
    x=[o[0], h1[0], h2[0], na[0]],
    y=[o[1], h1[1], h2[1], na[1]],
    z=[o[2], h1[2], h2[2], na[2]],
)
configuration.charge = 1
configuration.spin_multiplicity = 1
data = run_bsse("na_water", configuration, "1-3; 4", "0, 1")
kcal = data["interaction energy"] / 4.184
print(f"CP interaction energy = {kcal:.3f} kcal/mol")
print("Approx. literature reference (Dzidic & Kebarle 1970 et al.): ~ -24 kcal/mol")

print()
print("=" * 70)
print("Cl- ... H2O  (Cl- along the O-H1 bond extension, H...Cl = 2.20 A)")
print("=" * 70)
oh1_unit = (h1[0] / r_oh, h1[1] / r_oh, h1[2] / r_oh)
h_cl_distance = 2.20
cl_distance_from_o = r_oh + h_cl_distance
cl = (
    oh1_unit[0] * cl_distance_from_o,
    oh1_unit[1] * cl_distance_from_o,
    oh1_unit[2] * cl_distance_from_o,
)
configuration = new_configuration("cl_water")
configuration.atoms.append(
    symbol=["O", "H", "H", "Cl"],
    x=[o[0], h1[0], h2[0], cl[0]],
    y=[o[1], h1[1], h2[1], cl[1]],
    z=[o[2], h1[2], h2[2], cl[2]],
)
configuration.charge = -1
configuration.spin_multiplicity = 1
data = run_bsse("cl_water", configuration, "1-3; 4", "0, -1")
kcal = data["interaction energy"] / 4.184
print(f"CP interaction energy = {kcal:.3f} kcal/mol")
print("Approx. literature reference (Hiraoka et al. ion-clustering thermochemistry): ~ -13 kcal/mol")
