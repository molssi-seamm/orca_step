#!/usr/bin/env python
"""M2 acceptance-gate check: does the new seamm_bsse-driven BSSE path
reproduce the old ORCA Compound-script path, on real ORCA, for a given
dimer?

Run inside the seamm-dev environment:
    python validate_bsse_m2.py <sdf-path> [workdir] [label]

The first structure record of `sdf-path` is used (so a multi-configuration
production SDF works too -- only its first entry is read). Examples used for
the M2 gate:
    water: /Users/psaxe/structures/dimers/H2O-H2O.sdf
    FEC:   /Users/psaxe/structures/dimers/FEC-FEC dimer opt PM6-ORG.sdf
    EC:    /Users/psaxe/Dropbox/GM/TrainingData/Molecular/EC/EC_dimers_run_1.sdf
"""
import sys
from pathlib import Path

import molsystem
import orca_step
import seamm
import seamm_exec
from seamm.variables import Variables

seamm.flowchart_variables = Variables()

_HARTREE_TO_KJMOL = 2625.499639

SDF = Path(sys.argv[1])
WORKDIR = (
    Path(sys.argv[2])
    if len(sys.argv) > 2
    else Path("/tmp/bsse_m2_validate") / SDF.stem
)
LABEL = sys.argv[3] if len(sys.argv) > 3 else SDF.stem
WORKDIR.mkdir(parents=True, exist_ok=True)

# --- A real molsystem Configuration for the dimer (first record of SDF) ---
db = molsystem.SystemDB(filename=f"file:bsse_m2_validate_{abs(hash(LABEL))}?mode=memory&cache=shared")
system = db.create_system(name=LABEL)
configuration = system.create_configuration(name=LABEL)
configuration.from_sdf(SDF)
configuration.charge = 0
configuration.spin_multiplicity = 1
print(f"[{LABEL}] Loaded {configuration.n_atoms} atoms: {list(configuration.atoms.symbols)}")
molecules = configuration.find_molecules(as_indices=True)
print(f"[{LABEL}] find_molecules -> {len(molecules)} molecule(s), sizes {[len(m) for m in molecules]}")
assert len(molecules) == 2, f"expected two fragments, got {len(molecules)}"


class _Parent:
    options = {"ncores": "1", "memory": "available", "max-atoms-to-print": 25}
    global_options = {"root": str(Path.home() / "SEAMM"), "ncores": "1"}


class _Flowchart:
    def __init__(self, root_directory, executor):
        self.root_directory = root_directory
        self.executor = executor


executor = seamm_exec.get_executor("local")

# Deliberately cheap/fast (HF/def2-SVP, no RI): this run is about verifying
# the WIRING (job charges/ghosts/combine), not chemistry accuracy.
P = {
    "use model chemistry": "no",
    "method": "HF",
    "basis": "def2-SVP",
    "basis source": "ORCA internal",
    "auxiliary basis": "none",
    "grid": "default",
    "scf convergence": "default",
    "extra keywords": "",
    "basis set extrapolation": "none",
    "save wavefunction": "no",
    "compute gradient": "yes",
    "optimize monomers": "no",
    "fragments": "auto (molecules)",
    "fragment atoms": "",
    "fragment charges": "",
}


def make_node(subdir):
    node = orca_step.BSSE()
    node._id = (subdir,)
    node.flowchart = _Flowchart(str(WORKDIR), executor)
    node.parent = _Parent()
    node.get_system_configuration = lambda arg: (None, configuration)
    node._cite_references = lambda P: None
    node._cite_bsse = lambda: None
    node.next = lambda: None
    return node


# --- OLD path: the ORCA Compound script (bssegradient.cmp), the trusted
# Phase-1 reference. Calls bsse.py's still-present helper methods directly
# (they are no longer wired into run(), but are unchanged). ---
old_node = make_node("old")
fragA, fragB = sorted(molecules[0]), sorted(molecules[1])
xyz_filename = "bsse.xyz"
xyz_text = old_node._ghost_xyz(configuration, ghost_atoms=fragB)
compound_block, _, _ = old_node._compound_input(P, xyz_filename, "bssegradient.cmp")
extra_files = {
    "bssegradient.cmp": old_node._compound_script("bssegradient.cmp"),
    xyz_filename: xyz_text,
}
_, old_gradient = old_node.run_orca_compound(compound_block, extra_files=extra_files)
old_energies = old_node._parse_compound_energies(old_node.directory)
if old_energies is None:
    raise RuntimeError("OLD path: could not parse the 5 Compound sub-energies")
e_fragA, e_monA, e_fragB, e_monB, e_total = old_energies
old_corrected = e_total - (e_fragA - e_monA) - (e_fragB - e_monB)
old_uncorrected_interaction = e_total - e_monA - e_monB
old_corrected_interaction = e_total - e_fragA - e_fragB

print("\n=== OLD (Compound script) ===")
print(f"e_total (raw dimer)      = {e_total:.10f} Eh")
print(f"corrected energy         = {old_corrected:.10f} Eh")
print(f"bsse correction          = {old_corrected - e_total:.10f} Eh")
print(f"uncorrected interaction  = {old_uncorrected_interaction:.10f} Eh")
print(f"corrected interaction    = {old_corrected_interaction:.10f} Eh")
print(f"gradient[0]              = {old_gradient[0] if old_gradient else None}")

# --- NEW path: the seamm_bsse-driven run(). ---
new_node = make_node("new")
for key, value in P.items():
    if key in new_node.parameters:
        new_node.parameters[key].value = value
    else:
        print(f"(note: '{key}' is not a BSSEParameters key, skipped)")
captured = {}
new_node.analyze = lambda **kwargs: captured.update(kwargs)
new_node.run()
new_data = captured["data"]

new_uncorrected_interaction_eh = new_data["uncorrected interaction energy"] / _HARTREE_TO_KJMOL
new_corrected_interaction_eh = new_data["interaction energy"] / _HARTREE_TO_KJMOL

print("\n=== NEW (seamm_bsse) ===")
print(f"e_total (uncorrected)    = {new_data['uncorrected energy']:.10f} Eh")
print(f"corrected energy         = {new_data['energy']:.10f} Eh")
print(f"bsse correction          = {new_data['bsse correction']:.10f} Eh")
print(f"uncorrected interaction  = {new_uncorrected_interaction_eh:.10f} Eh")
print(f"corrected interaction    = {new_corrected_interaction_eh:.10f} Eh")
print(f"gradient[0]              = {new_data['gradients'][0]}")

print("\n=== DIFF (new - old) ===")
print(f"uncorrected energy       = {new_data['uncorrected energy'] - e_total:.3e} Eh")
print(f"corrected energy         = {new_data['energy'] - old_corrected:.3e} Eh")
print(
    "uncorrected interaction  = "
    f"{new_uncorrected_interaction_eh - old_uncorrected_interaction:.3e} Eh"
)
print(
    "corrected interaction    = "
    f"{new_corrected_interaction_eh - old_corrected_interaction:.3e} Eh"
)
if old_gradient and new_data.get("gradients"):
    diffs = [
        abs(a - b)
        for row_o, row_n in zip(old_gradient, new_data["gradients"])
        for a, b in zip(row_o, row_n)
    ]
    print(f"max |gradient diff|      = {max(diffs):.3e} Eh/bohr")

print(f"\nWorking directories: {WORKDIR}/old/ and {WORKDIR}/new/")
