# -*- coding: utf-8 -*-

"""Shared base class for the ORCA step and its sub-steps.

Holds the command-line/seamm.ini parser, the geometry block, and the routine
that writes the ORCA input, runs ORCA through the flowchart executor, and parses
the output. The main ``ORCA`` node and the ``Energy`` sub-step both inherit from
this (``Optimization`` inherits from ``Energy``), mirroring ``MOPACBase`` in the
MOPAC step.
"""

import configparser
import hashlib
import importlib
import logging
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

import seamm
import seamm_exec
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

from .orca_step import full_orca_path

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("ORCA")

# Families of '!' keywords of which ORCA takes only one: the SCF convergence and the
# integration grid presets.
_KEYWORD_FAMILIES = (
    (
        "SCF convergence",
        {
            "SLOPPYSCF",
            "LOOSESCF",
            "NORMALSCF",
            "STRONGSCF",
            "TIGHTSCF",
            "VERYTIGHTSCF",
            "EXTREMESCF",
        },
    ),
    ("integration grid", {"DEFGRID1", "DEFGRID2", "DEFGRID3"}),
)


def tidy_keyword_line(keyword_line):
    """Remove repeated and conflicting keywords from an ORCA '!' line.

    ORCA refuses a keyword given twice, ignoring case (e.g. the 'TIGHTSCF' from
    the SCF convergence setting plus a 'TightSCF' in the extra keywords). An exact
    repeat is dropped. Of several presets from one family (SCF convergence,
    integration grid) only the last is kept: the extra keywords come last, so a
    preset typed there overrides the setting.

    Returns
    -------
    (str, [str])
        The tidied line, and a note for each preset that was overridden.
    """
    words = keyword_line.split()
    keep = [True] * len(words)
    seen = set()
    for i, word in enumerate(words):
        if word.upper() in seen:
            keep[i] = False
        seen.add(word.upper())
    notes = []
    for family, members in _KEYWORD_FAMILIES:
        found = [i for i, w in enumerate(words) if keep[i] and w.upper() in members]
        for i in found[:-1]:
            keep[i] = False
            notes.append(
                f"The {family} preset {words[i]} is overridden by "
                f"{words[found[-1]]} in the keywords."
            )
    return " ".join(w for w, k in zip(words, keep) if k), notes


# Byte multipliers for parsing a memory string (SI 'GB' and binary 'GiB').
_MEMORY_UNITS = {
    "": 1,
    "b": 1,
    "k": 1000,
    "kb": 1000,
    "ki": 1024,
    "kib": 1024,
    "m": 1000**2,
    "mb": 1000**2,
    "mi": 1024**2,
    "mib": 1024**2,
    "g": 1000**3,
    "gb": 1000**3,
    "gi": 1024**3,
    "gib": 1024**3,
    "t": 1000**4,
    "tb": 1000**4,
    "ti": 1024**4,
    "tib": 1024**4,
}


def _dehumanize_bytes(text):
    """Parse a memory string ('3 GB', '512MB', '2Gi', or a bare number of bytes)
    into an integer number of bytes. Raises ValueError on anything unparseable."""
    m = re.fullmatch(r"\s*([0-9]*\.?[0-9]+)\s*([a-zA-Z]*)\s*", str(text))
    if not m:
        raise ValueError(f"cannot parse memory '{text}'")
    unit = m.group(2).lower()
    if unit not in _MEMORY_UNITS:
        raise ValueError(f"unknown memory unit '{m.group(2)}'")
    return int(float(m.group(1)) * _MEMORY_UNITS[unit])


def _library_path_vars(n_cores, library_path, environ=None):
    """The dynamic-library search variables to set for a parallel ORCA run.

    Returns a list of ``(name, value)`` pairs prepending ``library_path`` to
    ``DYLD_LIBRARY_PATH`` and ``LD_LIBRARY_PATH`` (preserving any existing
    value). Empty for a serial run or when no path is given. ORCA's OpenMPI
    launcher needs these to find e.g. ``libmpi.40.dylib``.
    """
    if n_cores <= 1 or not library_path:
        return []
    environ = os.environ if environ is None else environ
    pairs = []
    for var in ("DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
        existing = environ.get(var, "")
        pairs.append((var, library_path + (os.pathsep + existing if existing else "")))
    return pairs


# Hartree -> the energy unit ORCA reports (E_h); kept explicit for clarity.


# The rough cost of one ORCA calculation, relative to a hybrid DFT one, by the
# words in the "!" line (the first match wins) and by the basis set family.
_METHOD_COST = (
    (("ccsd(t)", "ccsd", "qcisd"), 30.0),
    (("dlpno",), 8.0),
    (("mp2", "dsd", "b2plyp", "pwpb95", "b2gp", "dsd-"), 4.0),
    (("xtb", "am1", "pm3", "mndo", "zindo"), 0.01),
    (("hf-3c", "pbeh-3c", "r2scan-3c", "b97-3c"), 0.3),
)
_BASIS_COST = (
    (("qzv", "cc-pvqz", "pcseg-3"), 15.0),
    (("tzv", "cc-pvtz", "pcseg-2"), 4.0),
)


def estimated_seconds(keyword_line, n_atoms):
    """A rough estimate of an ORCA run's time, for the inline rule.

    It only needs to tell a sub-minute job from a long one: about a second for
    water at hybrid DFT/def2-SVP, scaling as the atom count to the 2.5 power, more
    for correlated methods, larger and diffuse basis sets, and gradients.
    """
    words = keyword_line.lower()
    method = 1.0
    for keys, cost in _METHOD_COST:
        if any(k in words for k in keys):
            method = cost
            break
    basis = 1.0
    for keys, cost in _BASIS_COST:
        if any(k in words for k in keys):
            basis = cost
            break
    if "aug-" in words or re.search(r"def2-\w+d\b", words):
        basis *= 2.0
    extra = 3.0 if ("engrad" in words or "opt" in words or "freq" in words) else 1.0
    n = max(1, int(n_atoms))
    return 1.0 + 0.03 * n**2.5 * method * basis * extra


# ----------------------------------------------------------------------
# Predicted time (seamm_exec.timing_model; campaign seamm_exec 2026-10-05, Phase 3)
# ----------------------------------------------------------------------
_nbf_cache = {}

#: Basis functions per atom when the basis is not known to the Basis Set
#: Exchange: a double-zeta-plus-polarization count (def2-SVP: 5 for H, 14 for
#: C-Ne, 18 for Na-Ar, more beyond)
_NBF_FALLBACK = {1: 5, 2: 5}


def _fallback_nbf(z):
    if z in _NBF_FALLBACK:
        return _NBF_FALLBACK[z]
    if z <= 10:
        return 14
    if z <= 18:
        return 18
    if z <= 36:
        return 32
    return 40


def estimated_basis_functions(basis, atomic_numbers):
    """The number of (spherical) basis functions ``basis`` gives these atoms,
    counted from the Basis Set Exchange's definition when it has the basis,
    else from a double-zeta-plus-polarization estimate. Ghost centres count
    like atoms. Cached per basis and element."""
    name = (basis or "").strip()
    if name.lower().startswith("bse:"):
        name = name[4:]
    total = 0
    for z in atomic_numbers:
        key = (name.lower(), int(z))
        if key not in _nbf_cache:
            count = None
            if name:
                try:
                    import basis_set_exchange as bse

                    data = bse.get_basis(name, elements=[int(z)])
                    shells = data["elements"][str(int(z))]["electron_shells"]
                    count = sum(
                        len(shell["coefficients"]) * (2 * am + 1)
                        for shell in shells
                        for am in shell["angular_momentum"]
                    )
                except Exception:
                    count = None
            _nbf_cache[key] = count if count else _fallback_nbf(int(z))
        total += _nbf_cache[key]
    return total


def _atomic_numbers(configuration):
    """The configuration's atomic numbers, from its symbols if it has no
    numbers (a test's stand-in), else carbon for every atom: this only feeds
    a time estimate."""
    try:
        numbers = configuration.atoms.atomic_numbers
        if numbers is not None:
            return [int(z) for z in numbers]
    except Exception:
        pass
    try:
        from molsystem.elements import to_atnos

        return [int(z) for z in to_atnos(list(configuration.atoms.symbols))]
    except Exception:
        pass
    try:
        return [6] * int(configuration.n_atoms)
    except Exception:
        return []


def predicted_seconds(
    keyword_line,
    atomic_numbers,
    *,
    ghost_numbers=(),
    model=None,
    basis=None,
    charge=0,
    multiplicity=1,
    ntasks=1,
    quantile=0.5,
):
    """The expected time of an ORCA run, from the fitted cost model when there
    is one (``seamm_exec.timing_model.predict``), else :func:`estimated_seconds`.

    The task layer adds its own margin to a task's ``estimated_seconds`` (the
    bundle walltime is twice the estimate plus ten minutes), so the median is
    returned by default; ask for ``quantile=0.95`` for a figure to submit as is.
    The descriptors are those the run's record will carry, with the basis
    functions estimated from the basis and the atoms (ghosts included).

    Parameters
    ----------
    keyword_line : str
        The '!' line.
    atomic_numbers : sequence of int
        The real atoms.
    ghost_numbers : sequence of int
        The ghost centres' atomic numbers (basis functions, no electrons).
    model : str, optional
        The ``type@method/basis`` string; the method and basis come from it
        unless ``basis`` is given.
    """
    numbers = [int(z) for z in atomic_numbers]
    n_atoms = len(numbers)
    try:
        from seamm_exec import timing_model

        method = ""
        if model:
            level = str(model).split("@", 1)[-1]
            method, _, model_basis = level.partition("/")
            basis = basis or model_basis
        descriptors = {
            "task": task_kind(keyword_line),
            "method_class": method_class(method, keyword_line),
            "method": method,
            "basis": basis or "",
            "n_atoms": n_atoms,
            "n_heavy": sum(1 for z in numbers if z > 1),
            "n_ghosts": len(ghost_numbers),
            "charge": charge,
            "multiplicity": multiplicity,
            "n_electrons": sum(numbers) - int(charge or 0),
            "nbf": estimated_basis_functions(basis, [*numbers, *ghost_numbers]),
        }
        result = timing_model.predict(
            "orca", descriptors, ntasks=ntasks, quantile=quantile
        )
        if result is not None:
            return float(result["seconds"])
    except Exception as e:  # the hand estimate is the fallback
        logger.debug(f"No model prediction for the ORCA run: {e}")
    return estimated_seconds(keyword_line, n_atoms)


# ----------------------------------------------------------------------
# ORCA 6.1.1: a correlated gradient with exact exchange is wrong when batched
# ----------------------------------------------------------------------
_BATCHES = re.compile(r"Number of batches necessary\s*\.*\s*(\d+)")
_MOS_PER_BATCH = re.compile(r"Number of MOs treated per batch\s*\.*\s*(\d+)")
_MEMORY_NEEDED = re.compile(r"Memory needed for all in one shot\s*\.*\s*(\d+)\s*MB")
_MEMORY_GIVEN = re.compile(r"Memory devoted for MP2\s*\.*\s*(\d+)\s*MB")


def gradient_batching_problem(output_text, keyword_line=""):
    """Why an ORCA gradient in ``output_text`` cannot be trusted, or None.

    ORCA 6.1.1 computes the RI-MP2 / double-hybrid gradient wrongly -- the
    energy is unchanged to the last digit -- whenever ``%maxcore`` makes it
    split the gradient's MO loop into more than one batch ("Number of batches
    necessary ... N" with N > 1 in the output). Seen with exact exchange
    (NoCOSX, which this step chooses for Na, Mg, Zn, B and P): 417 meV/Å rms on
    Na+(H2O)6 at revDSD/def2-TZVPPD with 4 batches, forces ~2 eV/Å wrong on
    35-atom clusters (2026-10-06). A RIJCOSX run does not print the section.
    The remedy is more memory per process: fewer ranks, a larger %maxcore.
    """
    if not output_text:
        return None
    batches = [int(m) for m in _BATCHES.findall(output_text)]
    if not batches or max(batches) <= 1:
        return None
    mos = _MOS_PER_BATCH.findall(output_text)
    needed = _MEMORY_NEEDED.findall(output_text)
    given = _MEMORY_GIVEN.findall(output_text)
    exchange = "NoCOSX" if "nocosx" in (keyword_line or "").lower() else "this exchange"
    memory = ""
    if needed and given:
        memory = (
            f" It needed {needed[-1]} MB per process for one batch and had "
            f"{given[-1]} MB (%maxcore)."
        )
    return (
        f"ORCA split the correlated gradient into {max(batches)} batches of MOs"
        + (f" ({mos[0]} per batch)" if mos else "")
        + ", and ORCA 6.1.1 then computes the RI-MP2/double-hybrid gradient "
        f"wrongly with {exchange} while the energy stays correct.{memory} Give ORCA "
        "more memory per process so one batch holds every MO: fewer ranks and a "
        "larger 'memory' per core in orca.ini or the [orca-step] options (e.g. 4 "
        "ranks x 8 GB instead of 16 x 2 GB). The gradient of this run is not used."
    )


# ----------------------------------------------------------------------
# Timing records (seamm_exec.timing; campaign seamm_exec 2026-10-05)
# ----------------------------------------------------------------------
#: What ORCA's cost model is made of (seamm_exec.timing_model.Spec as plain
#: data): the size variables, the method class, the task and the unit count.
#: Written beside the records when a run is recorded.
TIMING_SPEC = {
    "size": ["nbf", "n_electrons", "n_atoms"],
    "klass": ["method_class"],
    "task": "task",
    "units": "scf_runs",
    "multiplier": None,
    "default_alpha": 0.8,
}


def _record_kwargs():
    """``spec=`` for seamm-exec releases that take it (2026.10.6.1 on)."""
    return {"spec": TIMING_SPEC} if hasattr(seamm_exec, "TimingSpec") else {}


#: ORCA keywords that mark a correlated or semiempirical method, for
#: :func:`method_class` when the method name is not known
_CC_PREFIXES = ("ccsd", "qcisd", "cepa", "ri-ccsd", "cisd", "ncisd")
_SEMIEMPIRICAL = {
    "xtb",
    "xtb0",
    "xtb1",
    "xtb2",
    "gfn-xtb",
    "gfn2-xtb",
    "am1",
    "pm3",
    "mndo",
}
_HF = {"hf", "rhf", "uhf", "rohf", "hf-3c", "ri-hf", "rijk-hf"}


def task_kind(keyword_line):
    """What kind of calculation an ORCA '!' line asks for, for the timing
    records: ``numfreq``, ``freq``, ``opt``, ``numgrad``, ``gradient`` or
    ``energy`` (the most expensive that applies)."""
    words = set(keyword_line.lower().split())
    if "numfreq" in words:
        return "numfreq"
    if "freq" in words or "anfreq" in words:
        return "freq"
    if any(w.endswith("opt") or w.startswith("opt") for w in words):
        return "opt"
    if "numgrad" in words:
        return "numgrad"
    if "engrad" in words:
        return "gradient"
    return "energy"


def method_class(method=None, keyword_line=""):
    """The class of a method, for the timing records: ``HF``, a DFT category
    from the step's metadata (``local``, ``GGA``, ``meta-GGA``, ``global
    hybrid``, ``range-separated hybrid``, ``global double-hybrid``,
    ``range-separated double-hybrid``), ``MP2``, ``DLPNO-CC``, ``CC``,
    ``semiempirical``, or ``""`` if it cannot be told.

    ``method`` is the step's method name (``B3LYP``, ``DLPNO-CCSD(T)``, ...);
    when it is not known the '!' line is scanned instead.
    """
    import orca_step

    functionals = orca_step.metadata["functionals"]
    methods = orca_step.metadata["methods"]
    name = (method or "").strip()
    if name:
        if name not in functionals and name not in methods:
            # The metadata's spelling, whatever the case given (wB97X-D3)
            lower = name.lower()
            name = next(
                (k for k in (*functionals, *methods) if k.lower() == lower), name
            )
        if name in functionals:
            return functionals[name]["category"]
        if name in methods:
            kind = methods[name]["type"]
            if kind == "QC":
                return "DLPNO-CC" if "dlpno" in name.lower() else "CC"
            if kind == "DFT":
                return "DFT"
            return kind
    words = keyword_line.lower().split()
    for word in words:
        if word.startswith("dlpno-cc"):
            return "DLPNO-CC"
        if word.startswith(_CC_PREFIXES):
            return "CC"
    by_lower = {k.lower(): v["category"] for k, v in functionals.items()}
    for word in words:
        if word in by_lower:
            return by_lower[word]
    for word in words:
        if "mp2" in word:
            return "MP2"
        if word in _SEMIEMPIRICAL:
            return "semiempirical"
        if word in _HF:
            return "HF"
    return ""


def timing_descriptors(
    input_text,
    output_text,
    *,
    model=None,
    n_atoms=None,
    n_heavy=None,
    n_ghosts=0,
    charge=None,
    multiplicity=None,
):
    """The descriptors of an ORCA run for its timing record: the numbers a
    cost model is fitted to (see the seamm_exec campaign of 2026-10-05).

    Parameters
    ----------
    input_text : str or None
        ``orca.inp``: the '!' line, ``%pal`` and ``%maxcore``.
    output_text : str or None
        ``orca.out``: electrons, basis functions, SCF runs and cycles, ORCA's
        own run time.
    model : str, optional
        The step's ``type@method/basis`` string, for the method and basis.
    n_atoms, n_heavy, n_ghosts, charge, multiplicity : optional
        Of the system as run (ghosts not counted in the atoms).
    """
    d = {}
    keyword_line = ""
    if input_text:
        for line in input_text.splitlines():
            if line.lstrip().startswith("!"):
                keyword_line = (keyword_line + " " + line.lstrip()[1:].strip()).strip()
        m = re.search(r"%pal\s+nprocs\s+(\d+)", input_text, re.IGNORECASE)
        d["nprocs"] = int(m.group(1)) if m else 1
        m = re.search(r"%maxcore\s+(\d+)", input_text, re.IGNORECASE)
        if m:
            d["maxcore_mb"] = int(m.group(1))

    method = basis = ""
    if model:
        level = model.split("@", 1)[-1]
        method, _, basis = level.partition("/")
    d["task"] = task_kind(keyword_line)
    d["method_class"] = method_class(method, keyword_line)
    d["method"] = method
    d["basis"] = basis
    d["model"] = model or ""
    d["keywords"] = keyword_line

    words = set(keyword_line.lower().split())
    if "nocosx" in words:
        d["ri"] = "NoCOSX"
    elif "rijcosx" in words:
        d["ri"] = "RIJCOSX"
    elif "rijk" in words or "ri-jk" in words:
        d["ri"] = "RIJK"
    elif "nori" in words:
        d["ri"] = "NoRI"
    else:
        d["ri"] = ""
    d["dispersion"] = next(
        (
            w
            for w in keyword_line.split()
            if w.lower() in ("d4", "d3bj", "d3", "d3zero")
        ),
        "",
    )

    d["n_atoms"] = n_atoms
    d["n_heavy"] = n_heavy
    d["n_ghosts"] = n_ghosts
    d["charge"] = charge
    d["multiplicity"] = multiplicity

    if output_text:
        text = output_text
        m = re.search(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", text)
        d["n_electrons"] = int(m.group(1)) if m else None
        m = re.search(r"Number of basis functions\s+\.+\s+(\d+)", text)
        d["nbf"] = int(m.group(1)) if m else None
        m = re.search(r"Number of auxiliary basis functions\s*\.+\s+(\d+)", text)
        d["n_aux"] = int(m.group(1)) if m else None
        m = re.search(r"Hartree-Fock type\s+HFTyp\s+\.+\s+(\w+)", text)
        d["hf_type"] = m.group(1) if m else ""
        cycles = [
            int(n)
            for n in re.findall(
                r"SCF (?:CONVERGED|NOT CONVERGED) AFTER\s+(\d+)\s+CYCLES", text
            )
        ]
        d["scf_runs"] = len(cycles)
        d["scf_cycles"] = sum(cycles)
        d["geometry_steps"] = len(
            re.findall(r"GEOMETRY OPTIMIZATION CYCLE\s+\d+", text)
        )
        m = re.search(
            r"TOTAL RUN TIME:\s+(\d+) days (\d+) hours (\d+) minutes (\d+) seconds"
            r" (\d+) msec",
            text,
        )
        if m:
            days, hours, minutes, seconds, msec = (int(x) for x in m.groups())
            d["code_seconds"] = (
                ((days * 24 + hours) * 60 + minutes) * 60 + seconds + msec / 1000.0
            )
        else:
            d["code_seconds"] = None
        d["terminated_normally"] = "ORCA TERMINATED NORMALLY" in text
    return d


def _heavy_atoms(configuration, atom_indices=None, ghost_atoms=None):
    """(n_atoms, n_heavy, n_ghosts) of the atoms ORCA was given."""
    numbers = list(configuration.atoms.atomic_numbers)
    indices = range(len(numbers)) if atom_indices is None else atom_indices
    ghosts = set(ghost_atoms or ())
    real = [i for i in indices if i not in ghosts]
    n_heavy = sum(1 for i in real if numbers[i] > 1)
    return len(real), n_heavy, len(ghosts)


def _orca_2aim(config):
    """The command for orca_2aim, which lives beside the orca binary.

    ``{code_dir}`` when ``code`` is a path, so the command names no absolute
    path; a bare ``orca`` (e.g. from conda) means orca_2aim is on the PATH.
    """
    if Path(config["code"]).expanduser().parent != Path("."):
        return "{code_dir}/orca_2aim"
    return "orca_2aim"


def _fingerprint(files, make_wfx):
    """Identify an ORCA run's inputs for restart, ignoring ``%pal``/``%maxcore``.

    Those lines come from the cores and memory of the machine the job happens
    to run on, which must not force a finished calculation to be redone.
    """
    h = hashlib.sha256()
    h.update(f"orca wfx={bool(make_wfx)}".encode())
    for name in sorted(files):
        data = files[name]
        if name == "orca.inp" and isinstance(data, str):
            data = "\n".join(
                line
                for line in data.splitlines()
                if not line.lstrip().lower().startswith(("%pal", "%maxcore"))
            )
        h.update(b"\0" + name.encode() + b"\0")
        h.update(data if isinstance(data, bytes) else str(data).encode())
    return "orca:" + h.hexdigest()


def mpi_env(n_cores, config):
    """Return ``(env, lib_prefix)`` for launching a (possibly parallel) ORCA.

    Parallel ORCA needs a matching OpenMPI runtime. Two things must line up
    and they are handled differently:

      1. mpirun (PATH): ORCA launches its workers with whatever ``mpirun`` it
         finds on PATH; it MUST be the same OpenMPI as the libraries the
         workers link, or the ranks corrupt each other's data (ORCA aborts
         with a "BLAS-ERROR: incompatible matrices"). We prepend the OpenMPI
         bin (the sibling of the lib dir given by library-path) so a
         different mpirun on PATH (e.g. a newer Homebrew OpenMPI) is not used.
         PATH is an ordinary variable, so it reaches ORCA's children.

      2. libmpi (dynamic-loader path): on Linux, LD_LIBRARY_PATH is honored
         and inherited, so exporting it is enough. On macOS, ORCA does NOT
         pass DYLD_* to the MPI sub-processes it spawns (and SIP strips DYLD_*
         through /bin/sh anyway), so the OpenMPI libraries must instead be on
         dyld's default search path -- e.g. symlink libmpi.*.dylib into
         /usr/local/lib. That is a one-time machine setup, documented in the
         User Guide; nothing here can substitute for it. We still export the
         loader variables for Linux.

      3. Core binding: ORCA calls its own ``mpirun`` internally once it
         parses ``%pal``, so we never see that command line -- but its
         environment is ours to set. On an unmanaged, interactive host
         (no SLURM etc. coordinating which cores belong to which job),
         every independent ``mpirun`` invocation applies OpenMPI's
         default binding policy (bind each rank to a core, starting from
         a low core number) with no knowledge of other concurrently-
         running ORCA jobs, so several jobs launched at once all pile
         onto the same one or two cores instead of spreading across the
         machine. Disabling binding (``OMPI_MCA_hwloc_base_binding_policy
         =none``, the environment-variable spelling of ``mpirun
         --bind-to none``, same on Linux and macOS -- no ``DYLD_*``/SIP
         wrinkle here) lets the OS scheduler load-balance ranks from
         concurrent jobs across all cores instead. Under a real
         scheduler this is unnecessary and can be counter-productive: a
         SLURM allocation already restricts the job to specific cores
         (cgroups), so we leave OpenMPI's binding on in that case, the
         same ``SLURM_JOB_ID`` check seamm_exec's ``in_situ`` auto-
         detection and ``computational_environment()`` use.

    The OpenMPI library directory is part of *how to run ORCA*, so it comes
    from the executor config (~/SEAMM/orca.ini), not the user-facing
    [orca-step] options. configparser lower-cases keys.
    """
    env = {}
    lib_prefix = []
    library_path = config.get("library-path", "") or ""
    quoted_library_path = shlex.quote(library_path)
    for var, value in _library_path_vars(n_cores, library_path):
        env[var] = value
        # Prepend at shell *runtime*, not with this Python-computed
        # snapshot of os.environ: an installation=modules code
        # (seamm_exec.Local.exec()) runs `module load ...` earlier in
        # this same generated script, which Python's os.environ can't
        # see yet at this point -- a static value here would silently
        # discard whatever that module load just added (confirmed for
        # real: ORCA's own liborca_tools_*.so disappeared from a
        # parallel run's LD_LIBRARY_PATH this way).
        #
        # Deliberately brace-free (an if/then/else, not
        # ``${VAR:+...}`` parameter expansion): this text is later run
        # through seamm_exec.local.Local.exec()'s own
        # ``command.format(**config, **ce)`` (for ITS "{code}"-style
        # placeholders), called in a loop *until the output stops
        # changing* -- a single-braced ``${VAR:+...}`` here survives
        # exactly one pass (mistaken for one of its own fields on the
        # next), and even doubling the braces only buys one extra
        # pass before the same KeyError (confirmed for real, both
        # ways). A construct with no braces at all is immune
        # regardless of how many passes run.
        lib_prefix.append(
            f'if [ -n "${var}" ]; then '
            f"export {var}={quoted_library_path}:${var}; "
            f"else export {var}={quoted_library_path}; fi;"
        )
    if n_cores > 1 and library_path:
        bindir = Path(library_path).expanduser().parent / "bin"
        if bindir.is_dir():
            lib_prefix.insert(0, f"export PATH={shlex.quote(str(bindir))}:$PATH;")
    if n_cores > 1 and "SLURM_JOB_ID" not in os.environ:
        env["OMPI_MCA_hwloc_base_binding_policy"] = "none"
    return env, lib_prefix


def _mpirun_version(mpirun):
    """``(implementation, major)`` of an ``mpirun``, e.g. ("Open MPI", 5), or
    ``(None, None)`` if it cannot be told."""
    try:
        result = subprocess.run(
            [str(mpirun), "--version"], capture_output=True, text=True, timeout=20
        )
    except (OSError, subprocess.SubprocessError):
        return None, None
    match = re.search(r"\(Open MPI\)\s+(\d+)\.", result.stdout + result.stderr)
    if match is None:
        return None, None
    return "Open MPI", int(match.group(1))


def check_mpirun(n_cores, config, path=None):
    """Stop a parallel run that ORCA would start with an unusable ``mpirun``.

    ORCA 6 is built with OpenMPI 4.1 and starts its parallel workers with the
    first ``mpirun`` on the PATH; with OpenMPI 5 -- Homebrew's, on a Mac -- the
    workers abort in ORCA's start-up with an obscure error. ``mpi_env`` puts
    the OpenMPI beside ``library-path`` first on the PATH, so this checks that
    one if ``orca.ini`` gives it, else the first on the PATH, and raises a
    clear error instead.

    Not checked: one core; ``installation = modules``, whose module, loaded in
    the job's script, supplies the ``mpirun`` Python cannot see beforehand; and
    an MPI that cannot be identified (not Open MPI).

    Parameters
    ----------
    n_cores : int
        The processes ORCA will be asked for (``%pal``).
    config : dict
        This machine's ``orca.ini`` section.
    path : str, optional
        The PATH to search; default the environment's.
    """
    if n_cores <= 1 or (config.get("installation") or "").strip() == "modules":
        return
    library_path = (config.get("library-path") or "").strip()
    mpirun = None
    if library_path:
        candidate = Path(library_path).expanduser().parent / "bin" / "mpirun"
        if candidate.is_file():
            mpirun = candidate
    if mpirun is None:
        mpirun = shutil.which("mpirun", path=path)
    where = "the [local] section of orca.ini"
    if mpirun is None:
        raise RuntimeError(
            f"ORCA was asked for {n_cores} processes, but there is no 'mpirun' on "
            "the PATH. Parallel ORCA needs OpenMPI 4.1: set 'library-path' in "
            f"{where} to the 'lib' directory of an OpenMPI 4.1 installation, or "
            "run ORCA on one core."
        )
    implementation, major = _mpirun_version(mpirun)
    if implementation == "Open MPI" and major is not None and major >= 5:
        raise RuntimeError(
            f"ORCA was asked for {n_cores} processes, but the 'mpirun' it would "
            f"use, {mpirun}, is Open MPI {major}. ORCA 6 needs OpenMPI 4.1; with "
            "a newer one its workers abort in start-up. Set 'library-path' in "
            f"{where} to the 'lib' directory of an OpenMPI 4.1 installation (its "
            "'bin/mpirun' is then used), or run ORCA on one core."
        )


class ORCABase(seamm.Node):
    """Common functionality for ORCA nodes."""

    def create_parser(self):
        """Set up the command-line / seamm.ini parser for the ORCA step.

        All ORCA nodes share the ``[orca-step]`` section (step_type), so the
        options are registered once here.
        """
        parser_name = self.step_type  # 'orca-step'
        parser = self.flowchart.parser

        parser_exists = parser.exists(parser_name)

        result = super().create_parser(name=parser_name)

        if parser_exists:
            return result

        # User-configurable run options (the [orca-step] section of the main
        # seamm.ini). How to find/launch ORCA -- its executable path and the
        # OpenMPI library directory -- lives in ~/SEAMM/orca.ini instead (see
        # data/orca.ini and _orca_config).
        parser.add_argument(
            parser_name,
            "--ncores",
            default="available",
            help=(
                "How many cores/processes ORCA may use (via %%pal). 'available' "
                "(the default) uses all cores the job/machine provides; give an "
                "integer to cap it, or '1' to force serial. Parallel ORCA needs a "
                "matching OpenMPI runtime, whose location is set with "
                "'library-path' in orca.ini."
            ),
        )
        parser.add_argument(
            parser_name,
            "--memory",
            default="available",
            help=(
                "Memory ORCA may use per process (its %%maxcore). 'available' (the "
                "default) scales to the memory-per-core of the machine; 'all' uses "
                "the whole node divided among the processes; or give an explicit "
                "amount such as '3 GB' (per process)."
            ),
        )
        parser.add_argument(
            parser_name,
            "--max-atoms-to-print",
            default=25,
            help="Maximum number of atoms for which to print detailed results.",
        )

        return result

    @property
    def is_runable(self):
        """Whether this node actually runs (vs. only contributing input)."""
        return True

    @property
    def version(self):
        """The semantic version of this module."""
        import orca_step

        return orca_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        import orca_step

        return orca_step.__git_revision__

    # ------------------------------------------------------------------
    # Input generation
    # ------------------------------------------------------------------
    def geometry_block(
        self, configuration, charge, multiplicity, atom_indices=None, ghost_atoms=None
    ):
        """Return the ORCA coordinate block for `configuration`.

        Parameters
        ----------
        atom_indices : Sequence[int] | None
            0-based atom indices to include, in this order. ``None`` (the
            default) is every atom, in the configuration's own order -- a BSSE
            sub-job passes a fragment's (sub-)list instead.
        ghost_atoms : Container[int] | None
            The subset of `atom_indices` to write as ORCA ghost centres (the
            element symbol with a trailing ``:``, e.g. ``O:`` -- basis
            functions only, no nucleus/electrons). ``None``/empty writes every
            atom as real.
        """
        lines = [f"* xyz {charge} {multiplicity}"]
        symbols = configuration.atoms.symbols
        xyzs = configuration.atoms.get_coordinates(fractionals=False, in_cell=True)
        indices = range(len(symbols)) if atom_indices is None else atom_indices
        ghosts = set(ghost_atoms) if ghost_atoms else set()
        for i in indices:
            symbol = symbols[i]
            x, y, z = xyzs[i]
            label = f"{symbol}:" if i in ghosts else symbol
            lines.append(f"{label:4s} {x:15.8f} {y:15.8f} {z:15.8f}")
        lines.append("*")
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------
    def run_orca(self, keyword_line, extra_blocks="", extra_files=None, make_wfx=False):
        """Write the ORCA input for the node's own configuration, run ORCA,
        and return the parsed results.

        A thin wrapper around :meth:`run_orca_job` for the common case: the
        whole of the node's own system, at its own charge/multiplicity, in
        the node's own directory. See :meth:`run_orca_job` for a sub-job
        against an explicit geometry/charge/multiplicity/directory (e.g. one
        of a BSSE correction's counterpoise sub-jobs).

        Parameters
        ----------
        keyword_line : str
            The contents of the ORCA "!" simple-input line (without the "!").
        extra_blocks : str
            Any additional ``%`` blocks to place before the geometry.
        extra_files : dict | None
            Extra input files to write into the run directory (e.g. an external
            basis file referenced by ``%basis GTOName ... end``).
        make_wfx : bool
            After ORCA finishes, run ``orca_2aim`` (shipped alongside ORCA) to
            convert the density (retained by the ``keepdensity`` keyword) into an
            AIMPAC ``orca.wfx`` wavefunction file for a following Atomic Charges
            step. Requires the ORCA input to include ``keepdensity``.

        Returns
        -------
        dict
            Parsed results, at least ``{"energy": <E_h>, "success": bool}``.
        """
        _, configuration = self.get_system_configuration(None)
        return self.run_orca_job(
            keyword_line,
            configuration,
            configuration.charge,
            configuration.spin_multiplicity,
            extra_blocks=extra_blocks,
            extra_files=extra_files,
            make_wfx=make_wfx,
        )

    def run_orca_job(
        self,
        keyword_line,
        configuration,
        charge,
        multiplicity,
        atom_indices=None,
        ghost_atoms=None,
        directory=None,
        extra_blocks="",
        extra_files=None,
        make_wfx=False,
    ):
        """Write an ORCA input for an explicit geometry/charge/multiplicity,
        run ORCA, and return the parsed results.

        The primitive ``run_orca`` and a BSSE counterpoise correction's
        per-fragment sub-jobs both drive: geometry/charge/multiplicity come
        from arguments instead of ``self.get_system_configuration()``/
        ``self.directory``, since a BSSE sub-job is a different atom subset
        (and, for a ghost-augmented fragment, a different charge) than the
        node's own system, run in its own sub-directory so the 2N + 1 jobs of
        one BSSE correction do not collide.

        Parameters
        ----------
        keyword_line : str
            The contents of the ORCA "!" simple-input line (without the "!").
        configuration : molsystem.Configuration
            The system to take atoms/coordinates from.
        charge, multiplicity : int
            The charge/multiplicity for *this job* -- not necessarily
            ``configuration.charge``/``configuration.spin_multiplicity`` (a
            BSSE fragment sub-job has its own).
        atom_indices : Sequence[int] | None
            0-based atom indices to include, in this order. ``None`` (the
            default) is every atom, in the configuration's own order.
        ghost_atoms : Container[int] | None
            The subset of `atom_indices` to write as ORCA ghost centres. See
            :meth:`geometry_block`.
        directory : str | Path | None
            Where to run this job. ``None`` (the default) is ``self.directory``
            -- the node's own job directory, as ``run_orca`` uses.
        extra_blocks : str
            Any additional ``%`` blocks to place before the geometry.
        extra_files : dict | None
            Extra input files to write into the run directory.
        make_wfx : bool
            After ORCA finishes, run ``orca_2aim`` to convert the retained
            density into an AIMPAC ``orca.wfx`` in *this job's* directory.

        Returns
        -------
        dict
            Parsed results, at least ``{"energy": <E_h>, "success": bool}``.
        """
        task = self.orca_job_task(
            keyword_line,
            configuration,
            charge,
            multiplicity,
            atom_indices=atom_indices,
            ghost_atoms=ghost_atoms,
            directory=directory,
            extra_blocks=extra_blocks,
            extra_files=extra_files,
            make_wfx=make_wfx,
        )
        run_directory = Path(task.directory)
        result = seamm_exec.run_task(task, node=self, directory=run_directory)
        self._check_task_result(result, run_directory)
        self._report_run_location(result, run_directory)
        self.record_timing(
            task,
            result,
            run_directory,
            configuration,
            charge,
            multiplicity,
            atom_indices=atom_indices,
            ghost_atoms=ghost_atoms,
        )

        data = self._parse_output(run_directory / "orca.out")
        problem = gradient_batching_problem(
            (
                (run_directory / "orca.out").read_text(errors="replace")
                if (run_directory / "orca.out").exists()
                else ""
            ),
            keyword_line,
        )
        if problem is not None and task_kind(keyword_line) in (
            "gradient",
            "opt",
            "freq",
            "numgrad",
            "numfreq",
        ):
            raise RuntimeError(f"ORCA in {run_directory}: {problem}")
        return data

    def record_timing(
        self,
        task,
        result,
        run_directory,
        configuration,
        charge,
        multiplicity,
        atom_indices=None,
        ghost_atoms=None,
    ):
        """Append this ORCA run's timing record (``~/.seamm.d/timing/orca.csv``)
        -- the common columns from the task layer, the descriptors from
        :func:`timing_descriptors`. Never raises; a restored result is not
        recorded. See seamm_exec's campaign of 2026-10-05.
        """
        try:
            if result.restored:
                return
            n_atoms, n_heavy, n_ghosts = _heavy_atoms(
                configuration, atom_indices, ghost_atoms
            )
            out = Path(run_directory) / "orca.out"
            descriptors = timing_descriptors(
                task.files.get("orca.inp"),
                out.read_text(errors="replace") if out.exists() else None,
                model=getattr(self, "model", None),
                n_atoms=n_atoms,
                n_heavy=n_heavy,
                n_ghosts=n_ghosts,
                charge=charge,
                multiplicity=multiplicity,
            )
            seamm_exec.record_task_timing(task, result, descriptors, **_record_kwargs())
        except Exception as e:  # pragma: no cover - must never stop the step
            logger.warning(f"Could not record the timing of the ORCA run: {e}")

    def orca_job_task(
        self,
        keyword_line,
        configuration,
        charge,
        multiplicity,
        atom_indices=None,
        ghost_atoms=None,
        directory=None,
        extra_blocks="",
        extra_files=None,
        make_wfx=False,
        key="orca",
    ):
        """The :class:`seamm_exec.Task` that :meth:`run_orca_job` runs, built
        without running it, so that a step can run several at once (a BSSE
        correction's sub-jobs) in one ``TaskSet``. ``key`` names it in that
        ``TaskSet``; the other arguments are those of :meth:`run_orca_job`.
        """
        run_directory = (
            Path(directory) if directory is not None else Path(self.directory)
        )
        run_directory.mkdir(parents=True, exist_ok=True)

        # Resources (cores + per-process memory), shared with run_orca_compound.
        n_cores, memory_mb = self._resources()

        keyword_line, notes = tidy_keyword_line(keyword_line)
        for note in notes:
            printer.normal(__(f"Note: {note}", indent=self.indent + 4 * " "))
        lines = [f"! {keyword_line}"]
        if n_cores > 1:
            lines.append(f"%pal nprocs {n_cores} end")
        lines.append(f"%maxcore {memory_mb}")
        if extra_blocks.strip() != "":
            lines.append(extra_blocks.rstrip())
        lines.append(
            self.geometry_block(
                configuration,
                charge,
                multiplicity,
                atom_indices=atom_indices,
                ghost_atoms=ghost_atoms,
            )
        )
        lines.append("")
        input_text = "\n".join(lines)

        files = {"orca.inp": input_text}
        if extra_files:
            files.update(extra_files)
        logger.debug(f"orca.inp ({run_directory}):\n" + input_text)

        # The task names only the program; where it runs, the 'orca' resolver
        # (orca_step.resolver) supplies ORCA's full path -- ORCA must be invoked
        # by it to find its sub-programs -- and, for a parallel run, the
        # OpenMPI paths, from that machine's orca.ini. When a wavefunction file
        # is wanted, chain orca_2aim (which lives next to the orca binary) in
        # the same shell so it runs in this directory right after ORCA, reading
        # the just-written orca.gbw/orca.densities.
        cmd = ["{code}", "orca.inp", ">", "orca.out", "2>", "orca.err"]
        return_files = [
            "orca.out",
            "orca.err",
            "orca.gbw",
            "orca.xyz",
            "orca.hess",
            "*.bibtex",
            "*.txt",
            "*.engrad",
        ]
        if make_wfx:
            # {orca_2aim} comes from the resolver: beside the orca binary.
            cmd += ["&&", "{orca_2aim}", "orca", ">", "orca_2aim.out", "2>&1"]
            return_files += ["orca.wfx", "orca_2aim.out"]

        all_numbers = _atomic_numbers(configuration)
        indices = range(len(all_numbers)) if atom_indices is None else atom_indices
        ghosts = set(ghost_atoms or ())
        real_numbers = [all_numbers[i] for i in indices if i not in ghosts]
        ghost_numbers = [all_numbers[i] for i in indices if i in ghosts]
        estimate = predicted_seconds(
            keyword_line,
            real_numbers,
            ghost_numbers=ghost_numbers,
            model=getattr(self, "model", None),
            charge=charge,
            multiplicity=multiplicity,
            ntasks=n_cores,
        )
        task = seamm_exec.Task(
            key=key,
            program="orca",
            cmd=cmd,
            directory=run_directory,
            files=files,
            # Wildcards: ORCA writes orca.bibtex (suggested citations) and
            # orca.property.txt (and other *.txt depending on options); the
            # executor discards files not requested. orca.engrad appears only
            # when gradients are requested (! EnGrad).
            return_files=return_files,
            # None: run in node-local scratch under a scheduler (SLURM's job
            # dir is NFS, unsafe for ORCA's MPI scratch I/O), in place
            # otherwise. See molssi-seamm/orca_step#20.
            in_situ=None,
            shell=True,
            resources=seamm_exec.Resources(ntasks=n_cores),
            estimated_seconds=estimate,
            fingerprint=_fingerprint(files, make_wfx),
            # ORCA exits 0 even after an error termination.
            success_text={"orca.out": "ORCA TERMINATED NORMALLY"},
        )
        return task

    def _check_task_result(self, result, directory):
        """Report a failed ORCA task, and stop if its output can't be used.

        A task that could not be run, or that was not run again because it has
        used up its attempts (its output is stale), raises. Otherwise a failure
        is noted and the output is parsed as always, so ORCA's own messages
        decide what happens next.
        """
        if result.ok:
            return
        reason = result.reason or "unknown"
        if result.returncode is None or reason.startswith("attempts exhausted"):
            raise RuntimeError(
                f"ORCA in {directory} did not run: {reason}\n{result.stderr}"
            )
        printer.normal(
            __(
                f"ORCA failed in {directory}: {reason}.",
                indent=self.indent + 4 * " ",
            )
        )

    def _report_run_location(self, result, directory):
        """Note in step.out where ORCA actually ran -- the job directory, or
        (under a scheduler) node-local scratch, with only results copied back.
        See molssi-seamm/orca_step#20: this distinction matters for diagnosing
        NFS-related failures and for knowing where any leftover scratch
        (e.g. after a crash) actually landed.

        ``result`` is the :class:`seamm_exec.TaskResult` of the run.
        """
        if result.restored:
            printer.normal(
                __(
                    "ORCA had already finished this calculation, with the same "
                    f"input, in {directory}; using those results.",
                    indent=self.indent + 4 * " ",
                )
            )
        elif result.in_situ in (None, True):
            printer.normal(
                __(
                    f"Ran ORCA directly in the job directory, {directory}.",
                    indent=self.indent + 4 * " ",
                )
            )
        else:
            printer.normal(
                __(
                    "Ran ORCA in node-local scratch "
                    f"({result.run_directory}), not the job directory, "
                    "because this job is running under a batch scheduler and "
                    "ORCA's MPI scratch I/O is not safe on NFS; only the "
                    "requested result files were copied back here.",
                    indent=self.indent + 4 * " ",
                )
            )

    def _resources(self):
        """Resolve (n_cores, memory_mb) for an ORCA run from the executor's
        computational environment and the [orca-step] user options.

        run_orca (and run_orca_compound) are called by a sub-step, whose options
        live on its parent (the main ORCA node), mirroring the MOPAC step.
        """
        ce = seamm_exec.computational_environment()
        options = self.parent.options
        seamm_options = self.parent.global_options

        # --- Cores/processes (ORCA %pal) ---
        available_cores = max(1, int(ce.get("NTASKS", 1) or 1))
        ncores_opt = str(options.get("ncores", "available")).strip().lower()
        if ncores_opt in ("available", "default", "all", ""):
            n_cores = available_cores
        else:
            try:
                n_cores = min(available_cores, int(ncores_opt))
            except ValueError:
                n_cores = available_cores
        # Respect a global cap from the SEAMM options, if one is set.
        global_ncores = str(seamm_options.get("ncores", "available")).strip().lower()
        if global_ncores not in ("available", "default", "all", ""):
            try:
                n_cores = min(n_cores, int(global_ncores))
            except ValueError:
                pass
        n_cores = max(1, n_cores)

        # --- Memory per process (ORCA %maxcore, in MB) ---
        mem_per_node = int(ce.get("MEM_PER_NODE", 0) or 0)  # bytes
        mem_per_cpu = int(ce.get("MEM_PER_CPU", 0) or 0)  # bytes
        memory_opt = str(options.get("memory", "available")).strip().lower()
        if memory_opt in ("available", "default", ""):
            # ~80% of the per-core memory, leaving headroom for the OS/driver.
            per_proc_bytes = int(0.8 * mem_per_cpu) if mem_per_cpu else 0
        elif memory_opt == "all":
            per_proc_bytes = int(mem_per_node / n_cores) if mem_per_node else 0
        else:
            try:
                per_proc_bytes = _dehumanize_bytes(options["memory"])
            except ValueError:
                per_proc_bytes = 0
        # ORCA's %maxcore is per-process MB; fall back to a conservative default
        # if the machine's memory could not be determined.
        memory_mb = int(per_proc_bytes / 1_000_000) if per_proc_bytes else 2000
        memory_mb = max(256, memory_mb)
        return n_cores, memory_mb

    def _mpi_env(self, n_cores, config):
        """``(env, lib_prefix)`` for a (possibly parallel) ORCA: see
        :func:`mpi_env`."""
        return mpi_env(n_cores, config)

    def run_orca_compound(
        self,
        compound_block,
        extra_files=None,
        engrad="result.engrad",
        make_wfx=False,
        wfx_step="orca_Compound_5",
    ):
        """Write and run an ORCA *Compound* job, returning ``(energy, gradient)``.

        Unlike ``run_orca`` (a single ``! keywords`` + inline geometry job), a
        Compound job runs several calculations from a script and writes its own
        result. The ``%pal``/``%maxcore`` preamble is emitted here (applying to
        every sub-calculation); the caller supplies the ``%Compound ... end``
        block and any files it references (the ``.cmp`` script, the geometry).

        Parameters
        ----------
        compound_block : str
            The ``%Compound "..." ... end`` block for ``orca.inp``.
        extra_files : dict | None
            Files the Compound job reads (e.g. ``{"bssegradient.cmp": ...,
            "bsse.xyz": ...}``), written into the run directory.
        engrad : str | None
            The EnGrad file the Compound writes with the final energy and
            gradient. Pass ``None`` for an energy-only job that writes no EnGrad
            (e.g. a method with no analytic gradient); then ``(None, None)`` is
            returned and the caller gets the energy elsewhere.

        Returns
        -------
        tuple(float | None, list | None)
            The energy (E_h) and gradient (``[n_atoms][3]``, E_h/bohr) from the
            Compound's EnGrad file, or ``(None, None)`` when ``engrad`` is None.
        """
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        n_cores, memory_mb = self._resources()

        lines = []
        if n_cores > 1:
            lines.append(f"%pal nprocs {n_cores} end")
        lines.append(f"%maxcore {memory_mb}")
        lines.append(compound_block.rstrip())
        lines.append("")
        input_text = "\n".join(lines)

        files = {"orca.inp": input_text}
        if extra_files:
            files.update(extra_files)
        logger.debug("orca.inp (Compound):\n" + input_text)

        return_files = ["orca.out", "orca.err", "*.bibtex", "*.txt", "*.engrad"]
        if engrad:
            return_files.append(engrad)
        # ORCA's full path and OpenMPI paths come from the resolver (run_orca_job)
        cmd = ["{code}", "orca.inp", ">", "orca.out", "2>", "orca.err"]
        # Convert the dimer step's retained density into a .wfx (as run_orca does
        # for the Energy substep), for a following Atomic Charges (DDEC6) step.
        # wfx_step is the Compound sub-job that kept its density (the dimer).
        if make_wfx:
            cmd += [
                "&&",
                "{orca_2aim}",
                wfx_step,
                ">",
                "orca_2aim.out",
                "2>&1",
                "&&",
                "cp",
                f"{wfx_step}.wfx",
                "orca.wfx",
            ]
            return_files += ["orca.wfx", "orca_2aim.out"]
        task = seamm_exec.Task(
            key="orca",
            program="orca",
            cmd=cmd,
            directory=directory,
            files=files,
            return_files=return_files,
            # See the comment in run_orca -- molssi-seamm/orca_step#20.
            in_situ=None,
            shell=True,
            resources=seamm_exec.Resources(ntasks=n_cores),
            fingerprint=_fingerprint(files, make_wfx),
            # ORCA exits 0 even after an error termination.
            success_text={"orca.out": "ORCA TERMINATED NORMALLY"},
        )
        result = seamm_exec.run_task(task, node=self, directory=directory)
        self._check_task_result(result, directory)
        self._report_run_location(result, directory)

        # Read the corrected gradient from the EnGrad file the Compound wrote,
        # when one was requested (energy-only jobs pass engrad=None).
        if engrad:
            return self._read_engrad(directory / engrad)
        return None, None

    @staticmethod
    def _read_engrad(path):
        """Read ``(energy, gradient)`` from an ORCA EnGrad file.

        Format (comment lines start with ``#``; blank lines are ignored): the
        atom count, then the energy (E_h), then ``3*n`` gradient components
        (E_h/bohr). Returns ``(None, None)`` if it cannot be parsed.
        """
        path = Path(path)
        if not path.exists():
            return None, None
        lines = [
            ln.strip()
            for ln in path.read_text().splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        try:
            n = int(lines[0])
            energy = float(lines[1])
            vals = [float(x) for x in lines[2 : 2 + 3 * n]]  # noqa: E203
        except (ValueError, IndexError):
            return None, None
        if len(vals) != 3 * n:
            return energy, None
        gradient = [vals[3 * i : 3 * i + 3] for i in range(n)]  # noqa: E203
        return energy, gradient

    def _orca_config(self):
        """Resolve the ORCA executable configuration.

        Reads the executor's section of ``<root>/orca.ini`` (the ``code`` = full
        path to ORCA, and the optional ``library-path`` for parallel runs);
        falls back to locating ``orca`` on the PATH. ORCA must be invoked by its
        full path so it can find its sub-programs.
        """
        executor_type = self.flowchart.executor.name
        seamm_options = self.parent.global_options
        ini_dir = Path(seamm_options["root"]).expanduser()
        ini_path = ini_dir / "orca.ini"

        full_config = configparser.ConfigParser()
        if ini_path.exists():
            full_config.read(ini_path)

        if (
            executor_type in full_config
            and full_config[executor_type].get("code", "") != ""
        ):
            config = dict(full_config.items(executor_type))
            # A bare name such as 'orca' must become the full path.
            config["code"] = full_orca_path(config["code"])
            return config

        # Fall back to finding ORCA on the PATH.
        code = shutil.which("orca")
        if code is None:
            raise RuntimeError(
                "Could not find the 'orca' executable. Put it on your PATH, or "
                "set 'code' (the full path) in the relevant section of "
                f"{ini_path}."
            )
        return {"code": code, "installation": "local"}

    # ------------------------------------------------------------------
    # Output parsing
    # ------------------------------------------------------------------
    def _parse_output(self, path):
        """Parse ORCA output: final single-point energy and success."""
        data = {"success": False, "energy": None}
        if not path.exists():
            raise RuntimeError(f"ORCA produced no output ({path}).")
        text = path.read_text()
        for line in text.splitlines():
            if "FINAL SINGLE POINT ENERGY" in line:
                try:
                    data["energy"] = float(line.split()[-1])
                except (ValueError, IndexError):
                    pass
        data["success"] = "ORCA TERMINATED NORMALLY" in text
        if not data["success"]:
            raise RuntimeError(
                f"ORCA did not terminate normally; see {path} and orca.err."
            )
        return data

    @staticmethod
    def _add_properties():
        """Register this package's properties with molsystem, if present."""
        import molsystem

        path = importlib.resources.files("orca_step") / "data"
        csv_file = path / "properties.csv"
        if path.exists() and csv_file.exists():
            molsystem.add_properties_from_file(csv_file)


# Run once at import to register any custom properties.
try:
    ORCABase._add_properties()
except Exception:  # pragma: no cover - non-fatal
    pass

# Re-export for convenience
__all__ = [
    "ORCABase",
    "tidy_keyword_line",
    "task_kind",
    "method_class",
    "timing_descriptors",
    "predicted_seconds",
    "estimated_basis_functions",
    "gradient_batching_problem",
    "printer",
    "job",
    "os",
    "__",
]
