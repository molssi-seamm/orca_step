# -*- coding: utf-8 -*-

"""ORCA's side of the Model Chemistry batch contract.

``get_task`` writes the same ORCA input the MDI engine writes for a structure
(the helpers are those of ``data/orca_mdi.py``: method, AutoAux, basis --
``bse:`` included -- ``EnGrad``, the DLPNO block), and ``analyze_task`` reads it
back with the engine's parsers, so the batch path and the MDI path give the same
numbers. ``options`` adds what fragments need: ``atom_indices``,
``ghost_atoms``, ``charge`` and ``multiplicity``.
"""

import importlib.util
from pathlib import Path
import tempfile

import numpy as np

import seamm_exec
from seamm_exec.evaluator import AnalysisError, check_properties, structure_data
from seamm_util import Q_

from .orca_base import _fingerprint, estimated_seconds

_ENGINE = Path(__file__).parent / "data" / "orca_mdi.py"
_engine = None

#: The defaults for one calculation: ranks and memory per rank (MB)
DEFAULT_NTASKS = 1
DEFAULT_MAXCORE = 2000


def engine_helpers():
    """The MDI engine's input builder and parsers (``data/orca_mdi.py``)."""
    global _engine
    if _engine is None:
        spec = importlib.util.spec_from_file_location("_orca_mdi_helpers", _ENGINE)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _engine = module
    return _engine


def method_and_blocks(model_chemistry):
    """The method words (keyword plus AutoAux), the basis, and any '%' block for
    a model chemistry -- exactly as ORCA's MDI engine is launched for it."""
    from .orca_step import dlpno_parent, mc_method_unalias, orca_method_keyword

    options = model_chemistry.get("options") or {}
    # The real ORCA keyword: the engine's, or (for methods with no MDI engine,
    # e.g. CCSD(T), which the batch path can still run) the un-aliased method.
    method = options.get("mdi_method_arg") or mc_method_unalias(
        model_chemistry.get("method")
    )
    basis = options.get("mdi_basis_arg") or model_chemistry.get("basis")
    blocks = []
    parent = dlpno_parent(method)
    keyword = orca_method_keyword(method)
    if parent is not None:
        blocks.append("%mp2 DLPNO true end")
    return f"{keyword} AutoAux", basis or "def2-SVP", blocks


def get_task(
    configuration,
    model_chemistry,
    *,
    key,
    properties=("energy", "gradients"),
    options=None,
    resources=None,
):
    """A :class:`seamm_exec.Task` computing ORCA's energy (and gradient) for one
    structure. See the module docstring."""
    options = dict(options or {})
    helpers = engine_helpers()
    data = structure_data(configuration)

    atom_indices = options.get("atom_indices")
    if atom_indices is None:
        atom_indices = list(range(len(data["atomic_numbers"])))
    ghosts = set(options.get("ghost_atoms") or ())
    symbols = []
    numbers = []
    for i in atom_indices:
        symbol = data["symbols"][i]
        symbols.append(f"{symbol}:" if i in ghosts else symbol)
        numbers.append(data["atomic_numbers"][i])
    coordinates = data["coordinates"][list(atom_indices)]
    charge = int(options.get("charge", data["charge"]))
    multiplicity = int(options.get("multiplicity", data["multiplicity"]))

    method, basis, blocks = method_and_blocks(model_chemistry)
    if multiplicity != 1 and blocks:
        raise ValueError(
            f"DLPNO {model_chemistry.get('method')} needs a closed-shell system "
            f"(ORCA has DLPNO-MP2 gradients only for RHF), not multiplicity "
            f"{multiplicity}."
        )
    files = {}
    with tempfile.TemporaryDirectory() as tmp:
        basis, basis_block = helpers.basis_keyword_and_block(basis, numbers, tmp)
        if basis_block:
            files["basis.bas"] = (Path(tmp) / "basis.bas").read_text()
            blocks.insert(0, basis_block)

    ntasks = DEFAULT_NTASKS
    maxcore = DEFAULT_MAXCORE
    if resources is not None:
        if resources.ntasks:
            ntasks = int(resources.ntasks)
        if resources.mem_per_cpu:
            maxcore = max(256, int(resources.mem_per_cpu / 1_000_000))
    blocks.append(f"%maxcore {maxcore}")

    gradients = "gradients" in properties
    text = helpers.orca_input(
        method,
        basis,
        charge,
        multiplicity,
        symbols,
        coordinates,
        ntasks,
        blocks="\n".join(blocks),
    )
    if not gradients:
        text = text.replace(" EnGrad\n", "\n", 1)
    files["orca.inp"] = text

    if resources is None:
        resources = seamm_exec.Resources(
            ntasks=ntasks, mem_per_cpu=int(maxcore * 1_000_000 / 0.8)
        )
    return seamm_exec.Task(
        key=key,
        program="orca",
        cmd=["{code}", "orca.inp", ">", "orca.out", "2>", "orca.err"],
        shell=True,
        files=files,
        return_files=["orca.out", "orca.err", "orca.engrad"],
        resources=resources,
        estimated_seconds=estimated_seconds(text.splitlines()[0], len(symbols)),
        fingerprint=_fingerprint(files, False),
        success_text={"orca.out": "ORCA TERMINATED NORMALLY"},
    )


def analyze_task(
    result,
    model_chemistry,
    configuration,
    *,
    properties=("energy", "gradients"),
    options=None,
):
    """The energy (kJ/mol) and gradients ((n, 3) kJ/mol/Å) of a finished task,
    converted exactly as the MDI path converts the engine's atomic units."""
    helpers = engine_helpers()
    options = dict(options or {})
    out = _text(result, "orca.out")
    data = {}
    energy = helpers.parse_energy(out) if out else None
    if energy is not None:
        data["energy"] = float(Q_(energy, "hartree").m_as("kJ/mol"))
    if "gradients" in properties:
        engrad = _text(result, "orca.engrad")
        if engrad:
            atom_indices = options.get("atom_indices")
            n = (
                len(atom_indices)
                if atom_indices is not None
                else len(configuration.atoms.atomic_numbers)
            )
            gradient = helpers.parse_engrad(engrad, n)
            data["gradients"] = Q_(
                np.asarray(gradient, dtype=float), "hartree/bohr"
            ).m_as("kJ/mol/Å")
    check_properties(data, properties, f"The ORCA calculation '{result.key}'")
    return data


def _text(result, name):
    data = result.files.get(name)
    if data is None and result.directory is not None:
        path = Path(result.directory) / name
        if path.exists():
            data = path.read_text(errors="replace")
    if isinstance(data, bytes):
        data = data.decode(errors="replace")
    return data


__all__ = ["get_task", "analyze_task", "AnalysisError"]
