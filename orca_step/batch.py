# -*- coding: utf-8 -*-

"""ORCA's side of the Model Chemistry batch contract.

``get_task`` writes the same ORCA input the MDI engine writes for a structure
(the helpers are those of ``data/orca_mdi.py``: method, AutoAux, basis --
``bse:`` included -- ``EnGrad``, the DLPNO block), and ``analyze_task`` reads it
back with the engine's parsers, so the batch path and the MDI path give the same
numbers. ``options`` adds what fragments need: ``atom_indices``,
``ghost_atoms``, ``charge`` and ``multiplicity``.

``options["aux_basis"]`` (or the model chemistry's) replaces AutoAux with explicit
auxiliary basis sets, e.g. "def2/J,def2-TZVPP/C" for a double hybrid -- commas,
since MBE's task options separate settings with spaces -- or "explicit" for the
sets that suit the method and basis, for structures where AutoAux makes the
auxiliary basis linearly dependent.

A job with one center, or with any atom of Na, Mg, Zn, B or P (ghosts
included), uses exact exchange (``NoCOSX``), as on every ORCA path
(orca_step#44). ``options["exact_exchange_elements"]``, or the same key in the
model chemistry's options, changes the element set ("none" for one-center jobs
only). ``options["precision"] = "tight"``, for finite differences of gradients,
adds VeryTightSCF and exact exchange (``NoCOSX``) for every structure.

An SCF that fails with the basis set nearly linearly dependent is run once more
with SThresh 1e-6 (see ``orca_base.sthresh_retry_command``); the task keeps its
fingerprint, and ``orca.sthresh-retry`` in its directory records the retry.
``options["sthresh_retry"]`` (or the model chemistry's) sets that SThresh, or
"no" turns the retry off; ``options["sthresh"]`` sets the SThresh of every run.
"""

import importlib.util
import logging
from pathlib import Path
import tempfile

import numpy as np

import seamm_exec
from seamm_exec.evaluator import AnalysisError, check_properties, structure_data
from seamm_util import Q_

from .orca_base import (
    DEFAULT_RETRY_STHRESH,
    STHRESH_FAILED_OUTPUT,
    STHRESH_RETRY_MARKER,
    _RETRY_INPUT,
    _fingerprint,
    _neighbours,
    _threshold,
    auxiliary_basis_words,
    composite_basis,
    composite_method,
    double_hybrid_scf,
    engrad_gradient,
    gradient_batching_problem,
    net_force_problem,
    predicted_seconds,
    sthresh_retry_command,
    sthresh_retry_input,
    sthresh_retry_note,
    sthresh_retry_value,
)

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


def method_and_blocks(model_chemistry, aux_basis=None):
    """The method words (keyword plus the auxiliary basis, by default AutoAux), the
    basis, and any '%' block for a model chemistry -- exactly as ORCA's MDI engine
    is launched for it. ``aux_basis`` is as :func:`auxiliary_basis_words` takes."""
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
    if composite_method(keyword) is not None:
        # A 3c composite is its keyword alone: no basis, no fitting basis (see
        # composite_basis, which refuses a basis that is not its own)
        return keyword, composite_basis(keyword, basis), blocks
    basis = basis or "def2-SVP"
    aux = auxiliary_basis_words(aux_basis, keyword, basis)
    return f"{keyword} {aux}".strip(), basis, blocks


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
    from .orca_step import mc_method_alias

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

    def option(name, default=None):
        """A task option, else the model chemistry's, else ``default``."""
        if name in options:
            return options[name]
        return (model_chemistry.get("options") or {}).get(name, default)

    method, basis, blocks = method_and_blocks(model_chemistry, option("aux_basis"))
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
            # ORCA's %maxcore is a target that it overshoots, and the leader rank
            # needs more besides, so keep 20 % of the memory per core in reserve,
            # as for the default resources below and the ORCA step's own runs:
            # at 2 ranks x 2 GB, %maxcore 2000 was killed for exceeding its 4 GB.
            maxcore = max(256, int(0.8 * resources.mem_per_cpu / 1_000_000))
    blocks.append(f"%maxcore {maxcore}")
    sthresh = _threshold(option("sthresh"))
    if sthresh is not None:
        blocks.append(f"%scf SThresh {sthresh:.3e} end")
    retry = sthresh_retry_value(option("sthresh_retry", DEFAULT_RETRY_STHRESH), sthresh)

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
    # A double-hybrid gradient runs with VeryTightSCF (see double_hybrid_scf)
    first, rest = text.split("\n", 1)
    line, _ = double_hybrid_scf(first[1:].strip())
    # precision="tight" (finite differences of gradients, e.g. VPT2): VeryTightSCF
    # and exact exchange (NoCOSX), unless an exchange scheme is already chosen
    if option("precision") == "tight":
        from .orca_step import tight_exchange_words

        words = line.split()
        if "VERYTIGHTSCF" not in (w.upper() for w in words):
            words = [w for w in words if w.upper() != "TIGHTSCF"] + ["VERYTIGHTSCF"]
        line = tight_exchange_words(" ".join(words))
    text = f"! {line}\n{rest}"
    files["orca.inp"] = text
    # The fingerprint is the calculation asked for, without the retry's input,
    # so that adding the retry reruns nothing already finished.
    fingerprint = _fingerprint(files, False)
    cmd = ["{code}", "orca.inp", ">", "orca.out", "2>", "orca.err"]
    return_files = ["orca.out", "orca.err", "orca.engrad"]
    if retry is not None:
        files[_RETRY_INPUT] = sthresh_retry_input(text, retry)
        cmd += sthresh_retry_command(retry)
        return_files += [STHRESH_RETRY_MARKER, STHRESH_FAILED_OUTPUT]

    if resources is None:
        resources = seamm_exec.Resources(
            ntasks=ntasks, mem_per_cpu=int(maxcore * 1_000_000 / 0.8)
        )
    estimate, estimate_source = predicted_seconds(
        text.splitlines()[0].lstrip("! "),
        [z for z, sym in zip(numbers, symbols) if not sym.endswith(":")],
        ghost_numbers=[z for z, sym in zip(numbers, symbols) if sym.endswith(":")],
        # A model-chemistry string, so the '/' in a keyword such as
        # REVDSD-PBEP86-D4/2021 is spelled '_' as in the model chemistry
        model=f"{mc_method_alias(method.split()[0])}/{basis}",
        basis=basis,
        charge=charge,
        multiplicity=multiplicity,
        ntasks=ntasks,
        neighbours=_neighbours(configuration, atom_indices).get("neighbours"),
        with_source=True,
    )
    return seamm_exec.Task(
        key=key,
        program="orca",
        cmd=cmd,
        shell=True,
        files=files,
        return_files=return_files,
        resources=resources,
        estimated_seconds=estimate,
        estimate_source=estimate_source,
        fingerprint=fingerprint,
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
    note = sthresh_retry_note(_text(result, STHRESH_RETRY_MARKER))
    if note is not None:
        logger.warning(f"The ORCA calculation '{result.key}': {note}")
    if task is not None and note is None:
        # Not a retried run: its time is that of two runs
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
            # Every centre, ghosts included: only their sum must vanish.
            problem = net_force_problem(engrad_gradient(engrad), out)
            if problem is not None:
                raise AnalysisError(f"The ORCA calculation '{result.key}': {problem}")
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
    if task is not None:
        record = calculation_record(task, configuration, options, out)
        if record is not None:
            data["calculation"] = record
    return data


def calculation_record(task, configuration, options, out):
    """The record of the task's settings (``orca_step.fingerprint``): the settings
    fingerprint, the element fingerprints, and the fingerprint of the geometry the
    calculation ran on, so that a consumer can store the energy with its level of
    theory and check that it pairs with the right geometry. None if it cannot be
    made; never raises."""
    from .fingerprint import orca_calculation_record, orca_version

    try:
        files = {
            name: text
            for name, text in task.files.items()
            if name != _RETRY_INPUT and isinstance(text, str)
        }
        data = structure_data(configuration)
        indices = (options or {}).get("atom_indices")
        if indices is None:
            indices = range(len(data["atomic_numbers"]))
        indices = list(indices)
        return orca_calculation_record(
            files,
            orca_version(out),
            [data["atomic_numbers"][i] for i in indices],
            data["coordinates"][indices],
            out_text=out,
        )
    except Exception as e:  # pragma: no cover - never fail the result
        logger.warning(f"No calculation record for '{getattr(task, 'key', '?')}': {e}")
        return None


def expected_fingerprints(
    configuration,
    model_chemistry,
    *,
    kind="energy",
    properties=("energy",),
    options=None,
    version=None,
):
    """The settings fingerprints a calculation of this structure would record,
    without running it.

    ORCA's settings depend on the structure through its elements: a one-center
    job, or one with Na, Mg, Zn, B or P, uses exact exchange (``NoCOSX``,
    orca_step#44). So the species of one model chemistry do not all share a
    fingerprint, and a consumer that looks up stored results (e.g. the fragment
    store of the thermochemistry step) must ask for each species' own.

    The input is built by the same code as the run: :func:`get_task` for an
    energy (the batch path), and the MDI engine's input builder for a Hessian
    (``orca_hessian_input``, as the engine's ``AnFreq`` run writes it). Whether
    COSX is moot (a functional without exact exchange, a basis without an
    auxiliary Coulomb basis) is known only from ORCA's output, so both
    fingerprints are returned: the stored one is always one of them, and for a
    structure without ``NoCOSX`` they are the same.

    Parameters
    ----------
    configuration : molsystem configuration or seamm_exec.Geometry
        The structure.
    model_chemistry : dict
        The model chemistry, as for :func:`get_task`.
    kind : str
        "energy" (the energy or energy + gradient of :func:`get_task`) or
        "hessian" (the MDI engine's analytic Hessian).
    properties : iterable of str
        For "energy": the properties asked for, as for :func:`get_task`.
    options : dict, optional
        The task options, as for :func:`get_task` (e.g.
        ``exact_exchange_elements``, ``precision``).
    version : str
        The ORCA version, e.g. "6.1.1", which is part of the fingerprint: take
        it from a record of the same run (one ORCA per installation).

    Returns
    -------
    [str]
        The fingerprints, without duplicates: as written first, then with COSX
        moot.
    """
    from .fingerprint import orca_calculation_record

    options = dict(options or {})
    if kind == "energy":
        task = get_task(
            configuration,
            model_chemistry,
            key="expected-fingerprint",
            properties=tuple(properties),
            options=options,
        )
        files = {
            name: text
            for name, text in task.files.items()
            if name != _RETRY_INPUT and isinstance(text, str)
        }
    elif kind == "hessian":
        files = _engine_hessian_files(configuration, model_chemistry, options)
    else:
        raise ValueError(f"kind must be 'energy' or 'hessian', not '{kind}'")
    out = []
    for moot in (False, True):
        record = orca_calculation_record(files, version, cosx_moot=moot)
        if record["fingerprint"] not in out:
            out.append(record["fingerprint"])
    return out


def _engine_hessian_files(configuration, model_chemistry, options):
    """The input files of the MDI engine's analytic Hessian (``AnFreq``) for a
    structure, as ``data/orca_mdi.py`` writes them when launched by
    ``ORCAStep.mdi_engine_command`` for this model chemistry."""
    from seamm_exec.evaluator import mdi_method_and_basis

    from .orca_base import composite_basis, composite_method
    from .orca_step import _engine_method_words, mc_method_unalias

    helpers = engine_helpers()
    data = structure_data(configuration)
    symbols = list(data["symbols"])
    numbers = list(data["atomic_numbers"])
    method, basis = mdi_method_and_basis(model_chemistry)
    method = mc_method_unalias(method)

    def option(name, default=None):
        if name in options:
            return options[name]
        return (model_chemistry.get("options") or {}).get(name, default)

    words = _engine_method_words(method, option("precision"))
    # As the engine: a composite has no basis or fitting-basis words
    # (orca_mdi.engine_basis; launched with --auxiliary-basis none)
    if composite_method(method) is not None:
        composite_basis(method, basis)  # refuses a basis that is not its own
        basis = helpers.engine_basis(words, None)
    else:
        basis = helpers.engine_basis(words, basis)
        words = f"{words} AutoAux"
    words = helpers.exchange_method(
        words, symbols, option("cosx", "by-element"), option("exact_exchange_elements")
    )
    files = {}
    with tempfile.TemporaryDirectory() as tmp:
        basis, block = helpers.basis_keyword_and_block(basis, numbers, tmp)
        if block:
            files["basis.bas"] = (Path(tmp) / "basis.bas").read_text()
    files["orca.inp"] = helpers.orca_hessian_input(
        words,
        basis,
        int(options.get("charge", data["charge"])),
        int(options.get("multiplicity", data["multiplicity"])),
        symbols,
        data["coordinates"],
        1,
        blocks=block,
    )
    return files


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
    """A returned file's text, read now: the result holds where it is, not its
    contents (from disk, or from the task archive)."""
    data = result.read(name)
    if data is None and result.directory is not None:
        path = Path(result.directory) / name
        if path.exists():
            data = path.read_text(errors="replace")
    return data


__all__ = ["get_task", "analyze_task", "AnalysisError"]
