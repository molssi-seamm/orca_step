#!/usr/bin/env python
"""Optimize the Na+...H2O and Cl-...H2O geometries (B3LYP-D3BJ/def2-TZVP,
matching the M3 level of theory) starting from M3's literature-informed,
unoptimized geometries, then re-run the ORCA BSSE sub-step at the optimized
geometry to sharpen the CP-corrected interaction energy against literature.

Run inside the seamm-dev environment:
    python optimize_ion_water.py [workdir]
"""
import math
import sys
from pathlib import Path

import molsystem
import orca_step
import seamm
import seamm_exec
from seamm.variables import Variables

seamm.flowchart_variables = Variables()

WORKDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/ion_water_opt")
WORKDIR.mkdir(parents=True, exist_ok=True)

executor = seamm_exec.get_executor("local")


class _Parent:
    options = {"ncores": "1", "memory": "available", "max-atoms-to-print": 25}
    global_options = {"root": str(Path.home() / "SEAMM"), "ncores": "1"}


class _Flowchart:
    def __init__(self, root_directory, executor):
        self.root_directory = root_directory
        self.executor = executor


P_LEVEL = {
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
}


def new_configuration(db_name):
    db = molsystem.SystemDB(filename=f"file:{db_name}?mode=memory&cache=shared")
    system = db.create_system(name=db_name)
    return system.create_configuration(name=db_name)


def make_node(cls, label, configuration):
    node = cls()
    node._id = (label,)
    node.flowchart = _Flowchart(str(WORKDIR), executor)
    node.parent = _Parent()
    node.get_system_configuration = lambda *a, **k: (configuration.system, configuration)
    node._cite_references = lambda P: None
    node.next = lambda: None
    for key, value in P_LEVEL.items():
        if key in node.parameters:
            node.parameters[key].value = value
    return node


def optimize(label, configuration):
    node = make_node(orca_step.Optimization, label, configuration)
    node.parameters["optimization convergence"].value = "TightOpt"
    node.run()
    print(f"[{label}] optimization finished.")


def bsse(label, configuration, fragment_atoms, fragment_charges):
    node = make_node(orca_step.BSSE, label, configuration)
    node.parameters["fragments"].value = "specified"
    node.parameters["fragment atoms"].value = fragment_atoms
    node.parameters["fragment charges"].value = fragment_charges
    node.parameters["compute gradient"].value = "no"
    captured = {}
    node.analyze_orig = node.analyze
    node.analyze = lambda **kwargs: (captured.update(kwargs), node.analyze_orig(**kwargs))
    node.run()
    return captured["data"]


def report(label, coords, symbols):
    print(f"\n[{label}] optimized geometry (A):")
    for symbol, (x, y, z) in zip(symbols, coords):
        print(f"  {symbol:2s} {x:12.6f} {y:12.6f} {z:12.6f}")


# ---------------------------------------------------------------------
# Na+...H2O -- same starting geometry as validate_bsse_m3.py.
# ---------------------------------------------------------------------
r_oh = 0.9572
half_angle = math.radians(104.5 / 2.0)
h1 = (r_oh * math.sin(half_angle), 0.0, r_oh * math.cos(half_angle))
h2 = (-r_oh * math.sin(half_angle), 0.0, r_oh * math.cos(half_angle))
o = (0.0, 0.0, 0.0)

print("=" * 70)
print("Na+ ... H2O  (optimizing from the M3 starting geometry)")
print("=" * 70)
na = (0.0, 0.0, -2.30)
configuration = new_configuration("na_water_opt")
configuration.atoms.append(
    symbol=["O", "H", "H", "Na"],
    x=[o[0], h1[0], h2[0], na[0]],
    y=[o[1], h1[1], h2[1], na[1]],
    z=[o[2], h1[2], h2[2], na[2]],
)
configuration.charge = 1
configuration.spin_multiplicity = 1
optimize("na_water_opt", configuration)
report(
    "Na+...H2O",
    configuration.atoms.get_coordinates(fractionals=False, in_cell=True),
    configuration.atoms.symbols,
)
data = bsse("na_water_bsse", configuration, "1-3; 4", "0, 1")
kcal = data["interaction energy"] / 4.184
print(f"\nCP interaction energy (optimized) = {kcal:.3f} kcal/mol")
print("M3 (unoptimized geometry): -26.110 kcal/mol")
print("Approx. literature reference: ~ -24 kcal/mol")

# ---------------------------------------------------------------------
# Cl-...H2O -- same starting geometry as validate_bsse_m3.py.
# ---------------------------------------------------------------------
print()
print("=" * 70)
print("Cl- ... H2O  (optimizing from the M3 starting geometry)")
print("=" * 70)
oh1_unit = (h1[0] / r_oh, h1[1] / r_oh, h1[2] / r_oh)
cl_distance_from_o = r_oh + 2.20
cl = tuple(oh1_unit[i] * cl_distance_from_o for i in range(3))
configuration = new_configuration("cl_water_opt")
configuration.atoms.append(
    symbol=["O", "H", "H", "Cl"],
    x=[o[0], h1[0], h2[0], cl[0]],
    y=[o[1], h1[1], h2[1], cl[1]],
    z=[o[2], h1[2], h2[2], cl[2]],
)
configuration.charge = -1
configuration.spin_multiplicity = 1
optimize("cl_water_opt", configuration)
report(
    "Cl-...H2O",
    configuration.atoms.get_coordinates(fractionals=False, in_cell=True),
    configuration.atoms.symbols,
)
data = bsse("cl_water_bsse", configuration, "1-3; 4", "0, -1")
kcal = data["interaction energy"] / 4.184
print(f"\nCP interaction energy (optimized) = {kcal:.3f} kcal/mol")
print("M3 (unoptimized geometry): -15.638 kcal/mol")
print("Approx. literature reference: ~ -13 kcal/mol")
