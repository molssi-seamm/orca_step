#!/usr/bin/env python
"""Does adding diffuse functions (def2-TZVPPD, the campaign's actual
production basis) close the Cl-...H2O gap to literature that def2-TZVP left
open, at the SAME (already def2-TZVP-optimized) geometry? Isolates the
basis effect from any geometry-relaxation effect.

Run inside the seamm-dev environment:
    python check_diffuse_basis.py [workdir]
"""
import sys
from pathlib import Path

import molsystem
import orca_step
import seamm
import seamm_exec
from seamm.variables import Variables

seamm.flowchart_variables = Variables()

WORKDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/diffuse_basis_check")
WORKDIR.mkdir(parents=True, exist_ok=True)

executor = seamm_exec.get_executor("local")


class _Parent:
    options = {"ncores": "1", "memory": "available", "max-atoms-to-print": 25}
    global_options = {"root": str(Path.home() / "SEAMM"), "ncores": "1"}


class _Flowchart:
    def __init__(self, root_directory, executor):
        self.root_directory = root_directory
        self.executor = executor


def make_node(cls, label, configuration, basis):
    node = cls()
    node._id = (label,)
    node.flowchart = _Flowchart(str(WORKDIR), executor)
    node.parent = _Parent()
    node.get_system_configuration = lambda *a, **k: (configuration.system, configuration)
    node._cite_references = lambda P: None
    node.next = lambda: None
    P = {
        "use model chemistry": "no",
        "method": "DFT",
        "functional": "B3LYP",
        "basis": basis,
        "basis source": "ORCA internal",
        "auxiliary basis": "AutoAux",
        "grid": "default",
        "scf convergence": "default",
        "extra keywords": "D3BJ",
        "basis set extrapolation": "none",
        "save wavefunction": "no",
    }
    for key, value in P.items():
        if key in node.parameters:
            node.parameters[key].value = value
    return node


def bsse(label, configuration, fragment_atoms, fragment_charges, basis):
    node = make_node(orca_step.BSSE, label, configuration, basis)
    node.parameters["fragments"].value = "specified"
    node.parameters["fragment atoms"].value = fragment_atoms
    node.parameters["fragment charges"].value = fragment_charges
    node.parameters["compute gradient"].value = "no"
    captured = {}
    orig = node.analyze
    node.analyze = lambda **kwargs: (captured.update(kwargs), orig(**kwargs))
    node.run()
    return captured["data"]


def new_configuration(db_name):
    db = molsystem.SystemDB(filename=f"file:{db_name}?mode=memory&cache=shared")
    system = db.create_system(name=db_name)
    return system.create_configuration(name=db_name)


# Optimized geometries from optimize_ion_water.py.
print("=" * 70)
print("Cl- ... H2O (def2-TZVP-optimized geometry), basis comparison")
print("=" * 70)
configuration = new_configuration("cl_water_diffuse")
configuration.atoms.append(
    symbol=["O", "H", "H", "Cl"],
    x=[-0.038552, 0.831586, -0.651766, 2.355097],
    y=[0.0, -0.0, 0.0, 0.0],
    z=[-0.022974, 0.459218, 0.718902, 1.949776],
)
configuration.charge = -1
configuration.spin_multiplicity = 1

for basis in ("def2-TZVP", "def2-TZVPPD"):
    data = bsse(f"cl_water_{basis}", configuration, "1-3; 4", "0, -1", basis)
    kcal = data["interaction energy"] / 4.184
    print(f"{basis:12s} CP interaction energy = {kcal:.3f} kcal/mol")
print("Approx. literature reference: ~ -13 kcal/mol")

print()
print("=" * 70)
print("Na+ ... H2O (def2-TZVP-optimized geometry), basis comparison")
print("=" * 70)
configuration = new_configuration("na_water_diffuse")
configuration.atoms.append(
    symbol=["O", "H", "H", "Na"],
    x=[-0.0, 0.767862, -0.767862, 0.0],
    y=[-0.0, 0.0, 0.0, 0.0],
    z=[-0.022620, 0.564173, 0.564173, -2.233698],
)
configuration.charge = 1
configuration.spin_multiplicity = 1

for basis in ("def2-TZVP", "def2-TZVPPD"):
    data = bsse(f"na_water_{basis}", configuration, "1-3; 4", "0, 1", basis)
    kcal = data["interaction energy"] / 4.184
    print(f"{basis:12s} CP interaction energy = {kcal:.3f} kcal/mol")
print("Approx. literature reference: ~ -24 kcal/mol")
