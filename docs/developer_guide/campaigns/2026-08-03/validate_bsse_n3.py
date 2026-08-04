#!/usr/bin/env python
"""N=3 sanity check: the ORCA BSSE path on a real Na+/Cl-/H2O trimer.

First real (real ORCA, not stubbed) exercise of the 3-fragment case: 2N+1=7
jobs (cluster + 3 fragment-in-cluster-basis + 3 fragment-alone), mixed
charges (+1/-1/0), one polyatomic fragment (H2O) alongside two monatomic
ions. Confirms the wiring generalizes past N=2 for real, and that the
assembled numbers are physically sane -- not a rigorous trimer PES study
(no optimization, one hand-built geometry: a contact Na+...Cl- pair with a
water coordinating Na+ from a different direction).

Run inside the seamm-dev environment:
    python validate_bsse_n3.py [workdir]
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

WORKDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/bsse_n3_validate")
WORKDIR.mkdir(parents=True, exist_ok=True)

executor = seamm_exec.get_executor("local")

# Same fast, dispersion-corrected level as the M3 two-body checks.
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
    # Hand-built configuration (atoms.append only, no bond table) -> "auto
    # (molecules)"/find_molecules() would see every atom as its own
    # fragment; "specified" with explicit groups avoids that.
    "fragments": "specified",
    "fragment atoms": "1; 2; 3-5",
    "fragment charges": "1, -1, 0",
}


class _Parent:
    options = {"ncores": "1", "memory": "available", "max-atoms-to-print": 25}
    global_options = {"root": str(Path.home() / "SEAMM"), "ncores": "1"}


class _Flowchart:
    def __init__(self, root_directory, executor):
        self.root_directory = root_directory
        self.executor = executor


# --- Geometry: a contact Na+...Cl- pair (R = 2.44 A, the M3-validated
# near-equilibrium point) with a water coordinating Na+ from a different
# direction (130 deg off the Na-Cl axis) through its O lone pair. ---
na = (0.0, 0.0, 0.0)
cl = (0.0, 0.0, 2.44)

approach_angle = math.radians(130.0)
approach_dir = (math.sin(approach_angle), 0.0, math.cos(approach_angle))  # from Na
na_o_distance = 2.35
o = tuple(na[i] + na_o_distance * approach_dir[i] for i in range(3))

r_oh = 0.9572
half_angle = math.radians(104.5 / 2.0)
local_x = (0.0, 1.0, 0.0)  # orthogonal to approach_dir (which lies in the xz-plane)
local_z = approach_dir  # H's point away from Na, O's lone pair faces Na


def offset(local_hx, local_hz):
    return tuple(o[i] + local_hx * local_x[i] + local_hz * local_z[i] for i in range(3))


h1 = offset(r_oh * math.sin(half_angle), r_oh * math.cos(half_angle))
h2 = offset(-r_oh * math.sin(half_angle), r_oh * math.cos(half_angle))

db = molsystem.SystemDB(filename="file:bsse_n3_validate?mode=memory&cache=shared")
system = db.create_system(name="na-cl-water")
configuration = system.create_configuration(name="na-cl-water")
configuration.atoms.append(
    symbol=["Na", "Cl", "O", "H", "H"],
    x=[na[0], cl[0], o[0], h1[0], h2[0]],
    y=[na[1], cl[1], o[1], h1[1], h2[1]],
    z=[na[2], cl[2], o[2], h1[2], h2[2]],
)
configuration.charge = 0
configuration.spin_multiplicity = 1
print(f"Loaded {configuration.n_atoms} atoms: {list(configuration.atoms.symbols)}")
print(f"Na-Cl distance   = {cl[2] - na[2]:.3f} A")
print(f"Na-O distance    = {na_o_distance:.3f} A (approach {math.degrees(approach_angle):.0f} deg off the Na-Cl axis)")

node = orca_step.BSSE()
node._id = ("na-cl-water",)
node.flowchart = _Flowchart(str(WORKDIR), executor)
node.parent = _Parent()
node.get_system_configuration = lambda arg: (None, configuration)
node._cite_references = lambda P: None
node._cite_bsse = lambda: None
node.next = lambda: None
for key, value in P.items():
    if key in node.parameters:
        node.parameters[key].value = value

captured = {}
node.analyze = lambda **kwargs: captured.update(kwargs)
node.run()
data = captured["data"]

print()
print("=" * 70)
print("Na+ ... Cl- ... H2O  (N=3, 7 ORCA jobs)")
print("=" * 70)
print(f"uncorrected energy       = {data['uncorrected energy']:.8f} Eh")
print(f"corrected energy         = {data['energy']:.8f} Eh")
print(f"bsse correction          = {data['bsse correction'] * 627.509474:.3f} kcal/mol")
print(f"uncorrected interaction  = {data['uncorrected interaction energy'] / 4.184:.3f} kcal/mol")
print(f"corrected interaction    = {data['interaction energy'] / 4.184:.3f} kcal/mol")
print()
print(
    "Sanity range: the Na+...Cl- pair alone binds ~-136 kcal/mol (M3); the "
    "3-fragment interaction here should be comparably large and negative "
    "(dominated by that term, plus the weaker Na+...water and Cl-...water "
    "contributions and any 3-body non-additivity) -- not a quantitative "
    "target, just a magnitude/sign check."
)
