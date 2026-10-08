# -*- coding: utf-8 -*-

"""ORCA's side of the Model Chemistry batch contract.

``get_task`` writes the same ORCA input the MDI engine writes for a structure
(the helpers are those of ``data/orca_mdi.py``: method, AutoAux, basis --
``bse:`` included -- ``EnGrad``, the DLPNO block), and ``analyze_task`` reads it
back with the engine's parsers, so the batch path and the MDI path give the same
numbers. ``options`` adds what fragments need: ``atom_indices``,
``ghost_atoms``, ``charge`` and ``multiplicity``.

A job with one center, or with any atom of Na, Mg, Zn, B or P (ghosts
included), uses exact exchange (``NoCOSX``), as on every ORCA path
(orca_step#44). ``options["exact_exchange_elements"]``, or the same key in the
model chemistry's options, changes the element set ("none" for one-center jobs
only).
"""

import importlib.util
import logging
from pathlib import Path
import tempfile

import numpy as np

import seamm_exec
from seamm_exec.evaluator import AnalysisError, check_properties, structure_data
from seamm_util import Q_

from .orca_base import _fingerprint, gradient_batching_problem, predicted_seconds

logger = logging.getLogger(__name__)

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

    from seamm_exec.evaluator import mdi_method_and_basis

    # The basis rule lives once, in seamm_exec. The real ORCA keyword: the
    # engine's, or (for methods with no MDI engine, e.g. CCSD(T), which the
    # batch path can still run) the un-aliased method.
    method, basis = mdi_method_and_basis(model_chemistry)
    method = mc_method_unalias(method)
    blocks = []
    parent = dlpno_parent(method)
    keyword = orca_method_keyword(method)
    if parent is not None:
        blocks.append("%mp2 DLPNO true end")
    return f"{keyword} AutoAux", basis or "def2-SVP", blocks


def can_run_task(configuration, model_chemistry, *, options=None):
    """Whether :func:`get_task` can run this structure: molecules only (ORCA has
    no periodic mode here), and no open-shell DLPNO double hybrid."""
    options = dict(options or {})
    data = structure_data(configuration)
    if data["periodicity"] != 0:
        return False
    multiplicity = int(options.get("multiplicity", data["multiplicity"]))
    _, _, blocks = method_and_blocks(model_chemistry)
    return not (blocks and multiplicity != 1)


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
    if data["periodicity"] != 0:
        raise ValueError(
            "ORCA calculations are molecular: a periodic structure cannot run as "
            "an ORCA task."
        )
    if options.get("guess") is not None:
        raise ValueError("ORCA tasks do not take an initial guess yet.")

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

    # The exact-exchange guard's elements: the task's option, else the model
    # chemistry's, else the default (orca_step#44)
    elements = options.get(
        "exact_exchange_elements",
        (model_chemistry.get("options") or {}).get("exact_exchange_elements"),
    )

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
        helpers.exact_exchange_method(method, symbols, elements),
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
        estimated_seconds=predicted_seconds(
            text.splitlines()[0].lstrip("! "),
            [z for z, sym in zip(numbers, symbols) if not sym.endswith(":")],
            ghost_numbers=[z for z, sym in zip(numbers, symbols) if sym.endswith(":")],
            model=f"{method.split()[0]}/{basis}",
            charge=charge,
            multiplicity=multiplicity,
            ntasks=ntasks,
        ),
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
    task=None,
):
    """The energy (kJ/mol) and gradients ((n, 3) kJ/mol/Å) of a finished task,
    converted exactly as the MDI path converts the engine's atomic units.

    Given the ``task`` that produced the result, its timing record is appended
    too (``seamm_exec.record_task_timing``), as the ORCA step's own runs do."""
    helpers = engine_helpers()
    options = dict(options or {})
    out = _text(result, "orca.out")
    if task is not None:
        _record_timing(task, result, out, model_chemistry, configuration, options)
    if "gradients" in properties:
        keyword_line = ""
        if task is not None:
            inp = task.files.get("orca.inp") or ""
            keyword_line = inp.splitlines()[0] if inp else ""
        problem = gradient_batching_problem(out, keyword_line)
        if problem is not None:
            raise AnalysisError(f"The ORCA calculation '{result.key}': {problem}")
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


def _record_timing(task, result, out, model_chemistry, configuration, options):
    """The timing record of a model-chemistry task; never raises."""
    try:
        from seamm_exec.evaluator import mdi_method_and_basis

        from .orca_base import (
            _heavy_atoms,
            _neighbours,
            _record_kwargs,
            timing_descriptors,
        )
        from .orca_step import mc_method_unalias

        method, basis = mdi_method_and_basis(model_chemistry)
        method = mc_method_unalias(method)
        atom_indices = options.get("atom_indices")
        ghosts = options.get("ghost_atoms")
        n_atoms, n_heavy, n_ghosts = _heavy_atoms(configuration, atom_indices, ghosts)
        descriptors = timing_descriptors(
            task.files.get("orca.inp") or _text(result, "orca.inp"),
            out,
            model=f"{method}/{basis}" if basis else method,
            n_ghosts=n_ghosts,
            n_atoms=n_atoms,
            n_heavy=n_heavy,
            charge=options.get("charge", configuration.charge),
            multiplicity=options.get("multiplicity", configuration.spin_multiplicity),
        )
        descriptors.update(_neighbours(configuration, atom_indices))
        seamm_exec.record_task_timing(task, result, descriptors, **_record_kwargs())
    except Exception as e:  # pragma: no cover
        logger.warning(f"Could not record the timing of ORCA task {task.key}: {e}")


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
