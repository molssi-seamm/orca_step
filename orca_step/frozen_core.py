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
