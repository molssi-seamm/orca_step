# -*- coding: utf-8 -*-

"""Frozen-core control and scalar-relativistic Hamiltonians for model-chemistry
calculations (composite phase 2 §3).

A recipe gives the frozen core as a count of ORBITALS per element (or ``none``, or
``default`` for ORCA's own choice). ORCA counts core ELECTRONS
(``%method NewNCore Z n end``). Only the elements whose count differs from ORCA's
default are written, so a calculation that asks for ORCA's default has the same
input, and the same settings fingerprint, as one that does not say. When anything
is written, ``CheckFrozenCore`` is pinned on. The elements are written by atomic
number: ORCA 6.1.1 reads ``NewNCore Hf n`` as hydrogen.

After the run the number of frozen orbitals ORCA reports is checked against the
count asked for (:func:`check_frozen_core`).
"""

#: ORCA 6.1.1's default frozen core, in electrons, by atomic number, as ORCA
#: prints it ("Setting NCore(X) from N to ..."), 2026-11-12. Li and Be are not
#: frozen; Na and Mg freeze only 1s; K-Zn freeze 1s-2p.
ORCA_DEFAULT_CORE_ELECTRONS = {}
for _first, _last, _n in (
    (1, 4, 0),  # H-Be
    (5, 12, 2),  # B-Mg
    (13, 30, 10),  # Al-Zn
    (31, 38, 18),  # Ga-Sr
    (39, 48, 28),  # Y-Cd
    (49, 70, 36),  # In-Yb
    (71, 80, 46),  # Lu-Hg
    (81, 86, 68),  # Tl-Rn
):
    for _z in range(_first, _last + 1):
        ORCA_DEFAULT_CORE_ELECTRONS[_z] = _n

#: The scalar-relativistic Hamiltonians and their ORCA keywords
RELATIVISTIC = {"none": None, "DKH2": "DKH2", "X2C": "X2C"}
#: The nucleus models, written explicitly with any relativistic Hamiltonian
NUCLEUS = {"point": "false", "finite": "true"}


_SYMBOLS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu "
    "Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs "
    "Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl "
    "Pb Bi Po At Rn"
).split()
ATOMIC_NUMBER = {s: i + 1 for i, s in enumerate(_SYMBOLS)}


def _atomic_number(element):
    if isinstance(element, int):
        return element
    if element not in ATOMIC_NUMBER:
        raise ValueError(f"Unknown element '{element}' in the frozen core")
    return ATOMIC_NUMBER[element]


def _orbitals_by_number(frozen_core):
    """{Z: orbitals} from {symbol or Z: orbitals}."""
    out = {}
    for element, n in frozen_core.items():
        if isinstance(element, str) and element.isdigit():
            element = int(element)
        out[_atomic_number(element)] = int(n)
    return out


def frozen_core_input(frozen_core, atomic_numbers):
    """The ORCA input for a frozen core, and the frozen orbitals expected.

    Parameters
    ----------
    frozen_core : None, "default", "none" or {element: orbitals}
        None and "default": ORCA's default; "none": all electrons correlated.
    atomic_numbers : [int]
        The real atoms of the calculation (not ghosts).

    Returns
    -------
    (keywords, block, expected)
        Keywords for the '!' line ([] or ["NoFrozenCore"]), a ``%method`` block
        (or ""), and the number of frozen orbitals expected (None for ORCA's
        default, which is not checked).
    """
    if frozen_core in (None, "default"):
        return [], "", None
    if frozen_core == "none":
        return ["NoFrozenCore"], "", 0
    if not isinstance(frozen_core, dict):
        raise ValueError(
            f"The frozen core must be 'default', 'none' or a count per element, not "
            f"{frozen_core!r}"
        )
    counts = _orbitals_by_number(frozen_core)
    present = sorted(set(int(z) for z in atomic_numbers))
    missing = [z for z in present if z not in counts]
    if missing:
        raise ValueError(
            "The frozen core gives no count for the element(s) with atomic number "
            + ", ".join(str(z) for z in missing)
        )
    lines = []
    for z in present:
        electrons = 2 * counts[z]
        if electrons != ORCA_DEFAULT_CORE_ELECTRONS.get(z):
            lines.append(f"  NewNCore {z} {electrons} end")
    expected = sum(counts[int(z)] for z in atomic_numbers)
    if not lines:
        return [], "", expected
    block = "%method\n  CheckFrozenCore true\n" + "\n".join(lines) + "\nend"
    return [], block, expected


def relativistic_input(relativistic, nucleus=None):
    """The '!' keywords and the ``%rel`` block of a scalar-relativistic
    Hamiltonian. The nucleus model is always written when there is one, so the
    input, and the settings fingerprint, say which was used."""
    if relativistic in (None, "", "none"):
        if nucleus not in (None, ""):
            raise ValueError("A nucleus model needs a relativistic Hamiltonian")
        return [], ""
    if relativistic not in RELATIVISTIC:
        raise ValueError(
            f"Unknown relativistic Hamiltonian '{relativistic}'; ORCA's model "
            f"chemistries offer {', '.join(RELATIVISTIC)}"
        )
    if nucleus not in NUCLEUS:
        raise ValueError(
            f"The relativistic Hamiltonian {relativistic} needs the nucleus model "
            f"given explicitly: {' or '.join(NUCLEUS)}, not {nucleus!r}"
        )
    return [RELATIVISTIC[relativistic]], f"%rel FiniteNuc {NUCLEUS[nucleus]} end"


#: The references an open-shell correlated calculation can be built on. "default"
#: and "UHF" are ORCA's own choice (UHF); "QRO" takes quasi-restricted orbitals
#: from the UHF solution. ROHF is not offered: ORCA's canonical coupled cluster
#: refuses it, and ROHF-NoIter silently gave OH the wrong state (2026-11-12).
REFERENCES = ("default", "UHF", "QRO")


def reference_input(reference, multiplicity):
    """The block that makes an open-shell correlated calculation use quasi-
    restricted orbitals, or "".

    Closed-shell calculations, and open-shell ones on ORCA's default (UHF), get
    nothing, so asking for a reference changes neither their input nor their
    settings fingerprint.
    """
    if reference in (None, ""):
        reference = "default"
    if reference not in REFERENCES:
        raise ValueError(
            f"The reference of an open-shell calculation must be one of "
            f"{', '.join(REFERENCES)}, not {reference!r}"
        )
    if reference == "QRO" and int(multiplicity) > 1:
        return "%mdci UseQROs true end"
    return ""


def check_reference(reference, multiplicity, parts):
    """The reason the run's reference is not the one asked for, or None."""
    if reference != "QRO" or int(multiplicity) == 1:
        return None
    got = (parts or {}).get("reference")
    if got != "QRO":
        return f"quasi-restricted orbitals were asked for, but the run used {got}"
    return None


def check_frozen_core(expected, parts):
    """The reason the run's frozen core is not the one asked for, or None."""
    if expected is None:
        return None
    if not parts or parts.get("frozen core orbitals") is None:
        return (
            f"{expected} frozen-core orbitals were asked for, but the run reports "
            "no count to check them against"
        )
    got = parts["frozen core orbitals"]
    if got != expected:
        return f"ORCA froze {got} core orbitals, not the {expected} asked for"
    return None


def model_chemistry_capabilities(mtype):
    """What a model-chemistry calculation of this type can be told (composite
    recipes check their components against it): the frozen core, a
    scalar-relativistic Hamiltonian and its nucleus model, and whether the energy
    comes with its parts (correlated methods), and the reference of an open-shell
    correlated calculation."""
    return {
        "frozen_core": ["default", "none", "per-element"],
        "reference": list(REFERENCES),
        "relativistic": ["none", "DKH2", "X2C"],
        "nucleus": ["point", "finite"],
        "decontract": True,
        "energy parts": mtype == "QC",
    }


def decontract_block(decontract):
    """The ``%basis`` block that decontracts the orbital basis, or ""."""
    if decontract in (None, False):
        return ""
    if decontract is not True:
        raise ValueError(f"decontract must be true or false, not {decontract!r}")
    return "%basis Decontract true end"


# -- Systems with no electron pairs to correlate --------------------------------------

#: The method keywords of uncorrelated (SCF) model chemistries
SCF_METHODS = {"HF", "RHF", "UHF", "ROHF"}


def is_correlated(model_chemistry):
    """Whether a model chemistry is a correlated wavefunction method (MP2, CCSD(T),
    their F12 and DLPNO forms), whose energy ORCA's MP2 and MDCI modules compute."""
    if str(model_chemistry.get("type") or "").upper() != "QC":
        return False
    method = str(model_chemistry.get("method") or "").split("/")[0].strip().upper()
    return method not in SCF_METHODS


def frozen_core_electrons(frozen_core, atomic_numbers):
    """The number of core electrons frozen: ORCA's default, none, or the count per
    element asked for (see :func:`frozen_core_input`)."""
    if frozen_core in (None, "default"):
        return sum(ORCA_DEFAULT_CORE_ELECTRONS.get(int(z), 0) for z in atomic_numbers)
    if frozen_core == "none":
        return 0
    counts = _orbitals_by_number(frozen_core)
    return sum(2 * counts.get(int(z), 0) for z in atomic_numbers)


def correlation_electrons(frozen_core, atomic_numbers, charge):
    """(electrons, correlated electrons) of a calculation."""
    n = sum(int(z) for z in atomic_numbers) - int(charge)
    return n, n - frozen_core_electrons(frozen_core, atomic_numbers)


def correlated_pairs(frozen_core, atomic_numbers, charge):
    """The smallest number of correlated electron pairs ORCA's MDCI module can
    count for a calculation: m(m + 1)/2 for the m correlated electron pairs
    (pairs of occupied orbitals i <= j). Exact for a closed shell; an open
    shell has at least as many. ORCA refuses to run MDCI on more processes than
    pairs ("Number of processes (2) in parallel calculation exceeds number of
    pairs (1)" for H2)."""
    _, correlated = correlation_electrons(frozen_core, atomic_numbers, charge)
    m = max(0, correlated) // 2
    return m * (m + 1) // 2


def check_electron_pairs(frozen_core, atomic_numbers, charge):
    """Whether a correlated calculation is of a one-electron system (H, He+,
    H2+, ...), which runs as HF: there is nothing to correlate, so the correlated
    energy is the HF energy, exact within the basis. ORCA's MP2 and MDCI modules
    refuse such a system ("no electrons to correlate", "number of processes
    exceeds number of pairs (0)").

    Raises
    ------
    ValueError
        If the frozen core leaves fewer than two correlated electrons in a system
        with more than one: no pairs to correlate, but the HF energy would not be
        the correlated one either.
    """
    n, correlated = correlation_electrons(frozen_core, atomic_numbers, charge)
    if n <= 1:
        return True
    if correlated < 2:
        raise ValueError(
            f"The frozen core leaves {correlated} of the {n} electrons correlated, "
            "so there are no electron pairs to correlate. Freeze fewer core "
            "orbitals."
        )
    return False


def one_electron_parts(energy, reference, frozen_core_orbitals, relativistic, f12):
    """The ``energy parts`` of a one-electron system's correlated energy, run as HF:
    the HF energy is the total and the correlation is zero.

    The reference is reported as the one asked for: with one electron the UHF,
    ROHF and quasi-restricted determinants are the same. An F12 method's CABS
    singles correction is not included: ORCA computes it only in its F12
    correlation modules, which refuse a one-electron system (for the H atom with
    cc-pVTZ-F12 it would be about -5e-5 E_h, the basis's HF error).
    """
    parts = {
        "parts version": 1,
        "units": "E_h",
        "parts": {"HF": energy},
        "total": ["HF"],
        "factors": {},
        "reference": reference,
        "frozen core orbitals": frozen_core_orbitals,
        "relativistic": relativistic,
        "diagnostics": {},
        "sources": {"HF": "output"},
        "one electron": True,
    }
    if f12:
        parts["omitted"] = ["CABS singles"]
    return parts
