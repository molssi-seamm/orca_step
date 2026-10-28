# -*- coding: utf-8 -*-

"""Fingerprints of the settings of ORCA calculations.

What ORCA's input holds that does not belong to the settings is removed here,
and :func:`seamm_util.fingerprint.calculation_record` hashes the rest. The
Energy, Optimization and Frequencies sub-steps and the MDI engine
(``data/orca_mdi.py``) all use :func:`orca_calculation_record`, so a step and an
engine running the same calculation give the same fingerprint.

Removed:

* the cores and memory (``%pal``, ``PALn``, ``%maxcore``), file names (``%base``,
  ``%moinp``) and output (``%output``, ``%elprop``, ``%plots``, print keywords);
* the starting guess (``MORead``, ``AutoStart``, ``Guess`` in ``%scf``), which is
  not a setting;
* the task (``SP``, ``EnGrad``, ``Opt``, ``Freq``, ...) and the optimizer
  (``%geom``), which do not change the energy at a geometry;
* the geometry, charge and multiplicity, which belong to the species;
* the thermochemistry temperature (``%freq Temp``), which the names of the
  thermochemical properties carry.

Kept: everything else on the '!' line, in any order (ORCA's keywords are
case-insensitive and unordered), and the other blocks. NoCOSX, which the step and
engine add for some elements, is kept: a calculation with exact exchange is not
the same level as one with COSX. RHF, UHF and similar words are kept too; ORCA
chooses from the multiplicity itself, so they are on the '!' line only when the
user asked for them.

Basis sets written per element -- a Basis Set Exchange basis in ``basis.bas``,
or ``NewGTO``, ``NewAuxJGTO``, ``NewECP``, ... in ``%basis`` -- are hashed per
element, so that species with different elements can be compared.
"""

import re

from seamm_util.fingerprint import calculation_record, normalize_lines

#: '!' words that do not change the energy at a given geometry
_DROPPED_WORDS = {
    # tasks
    "SP",
    "ENERGY",
    "ENGRAD",
    "GRADIENT",
    "NUMGRAD",
    "OPT",
    "COPT",
    "ZOPT",
    "GDIIS-OPT",
    "LOOSEOPT",
    "NORMALOPT",
    "TIGHTOPT",
    "VERYTIGHTOPT",
    "OPTTS",
    "FREQ",
    "ANFREQ",
    "NUMFREQ",
    # starting guess
    "MOREAD",
    "AUTOSTART",
    "NOAUTOSTART",
    # output
    "MINIPRINT",
    "SMALLPRINT",
    "NORMALPRINT",
    "LARGEPRINT",
    "PRINTBASIS",
    "PRINTMOS",
    "PRINTGAP",
    "NOPOP",
    "ALLPOP",
    "KEEPDENSITY",
    "KEEPINTS",
    "XYZFILE",
    "PDBFILE",
    "HIRSHFELD",
    "CHELPG",
    "MAYER",
    "NBO",
}
_PAL = re.compile(r"PAL\d+")

#: Blocks dropped whole
_DROPPED_BLOCKS = {"PAL", "OUTPUT", "ELPROP", "PLOTS", "GEOM"}
#: One-line directives without an 'end'
_DROPPED_DIRECTIVES = {"MAXCORE", "BASE", "MOINP"}
#: Lines dropped within a kept block, by block
_DROPPED_LINES = {
    "SCF": ("GUESS", "MOINP", "AUTOSTART"),
    "FREQ": ("TEMP",),
    "BASIS": ("GTONAME",),
}
#: Words that open a sub-block, closed by its own 'end'
_SUB_BLOCKS = {"CONSTRAINTS", "SCAN", "MODIFY_INTERNAL", "COORDS", "HESS"}
_ELEMENT_SUB_BLOCK = re.compile(r"NEW\w*(GTO|ECP)\s+(\w+)")

#: Element names in a Basis Set Exchange ORCA-format file
_ELEMENT_NAMES = None


def orca_version(out_text):
    """The ORCA version from its output, e.g. "6.1.1", or None."""
    m = re.search(r"Program Version\s+(\S+)", out_text or "")
    return m.group(1) if m else None


def _element_symbol(name):
    """The symbol for a name in a BSE basis file ("HYDROGEN" -> "H")."""
    global _ELEMENT_NAMES
    if _ELEMENT_NAMES is None:
        from seamm_util import element_data

        _ELEMENT_NAMES = {}
        for symbol, data in element_data.items():
            if data.get("element"):
                _ELEMENT_NAMES[str(data["element"]).upper()] = symbol
            _ELEMENT_NAMES[symbol.upper()] = symbol
    return _ELEMENT_NAMES.get(name.upper(), name.capitalize())


def basis_file_blocks(text):
    """The per-element sections of an ORCA-format basis file from the Basis Set
    Exchange (``$DATA``, then a line with each element's name and its shells,
    ``$END``), as {symbol: [normalized lines]}."""
    blocks = {}
    current = None
    for line in normalize_lines(text, comment="!", upper=True):
        if line in ("$DATA", "$END"):
            current = None
            continue
        if re.fullmatch(r"[A-Z]+", line) and len(line) > 2:
            current = _element_symbol(line)
            blocks.setdefault(current, [])
            continue
        if current is not None:
            blocks[current].append(line)
    return blocks


def settings_text(input_text):
    """The lines of an ORCA input that make up its settings, and the per-element
    basis sets in it.

    Parameters
    ----------
    input_text : str
        The contents of orca.inp.

    Returns
    -------
    ([str], {str: [str]})
        The normalized lines, and the element blocks by symbol.
    """
    words = set()
    lines = []
    elements = {}
    in_geometry = False
    block = None  # the name of the current kept or dropped block
    dropping = False
    depth = 0  # sub-blocks open in the current block
    element = None  # the element sub-block of %basis being read

    for line in normalize_lines(input_text, comment="#", upper=True):
        if in_geometry:
            if line == "*":
                in_geometry = False
            continue
        if block is None:
            if line.startswith("!"):
                for word in line[1:].split():
                    if word in _DROPPED_WORDS or _PAL.fullmatch(word):
                        continue
                    words.add(word)
                continue
            if line.startswith("*"):
                # The geometry: '* xyz 0 1' ... '*', or '*xyzfile 0 1 file'
                if "FILE" not in line.split()[0] and not (
                    len(line.split()) > 1 and "FILE" in line.split()[1]
                ):
                    in_geometry = True
                continue
            if line.startswith("%"):
                tokens = line[1:].split()
                name = tokens[0] if tokens else ""
                if name in _DROPPED_DIRECTIVES:
                    continue
                single = len(tokens) > 1 and tokens[-1] == "END"
                if name in _DROPPED_BLOCKS:
                    if not single:
                        block, dropping, depth = name, True, 0
                    continue
                if single:
                    # The same as the block over several lines
                    inner = " ".join(tokens[1:-1])
                    lines.append("%" + name)
                    if inner and not any(
                        inner.startswith(x) for x in _DROPPED_LINES.get(name, ())
                    ):
                        lines.append(inner)
                    lines.append("END")
                else:
                    block, dropping, depth = name, False, 0
                    lines.append("%" + name)
                    rest = " ".join(tokens[1:])
                    if rest and not any(
                        rest.startswith(x) for x in _DROPPED_LINES.get(name, ())
                    ):
                        lines.append(rest)
                continue
            lines.append(line)  # anything else, kept as is
            continue

        # Inside a block
        first = line.split()[0]
        if element is not None:
            if line == "END":
                element = None
            else:
                elements.setdefault(element, []).append(line)
            continue
        if line == "END":
            if depth > 0:
                depth -= 1
                if not dropping:
                    lines.append(line)
            else:
                if not dropping:
                    lines.append("END")
                block = None
            continue
        if dropping:
            if first in _SUB_BLOCKS:
                depth += 1
            continue
        if block == "BASIS":
            m = _ELEMENT_SUB_BLOCK.match(line)
            if m is not None:
                element = _element_symbol(m.group(2))
                elements.setdefault(element, []).append(line.split()[0])
                rest = line[m.end() :].strip()
                closed = rest == "END" or rest.endswith(" END")
                if closed:
                    rest = rest[: -len("END")].strip()
                if rest:
                    elements[element].append(rest)
                if closed:
                    element = None
                continue
        if any(line.startswith(x) for x in _DROPPED_LINES.get(block, ())):
            continue
        if first in _SUB_BLOCKS:
            depth += 1
        lines.append(line)

    keywords = ["! " + " ".join(sorted(words))] if words else []
    # A block left with nothing in it is no setting
    cleaned = []
    for line in keywords + lines:
        if line == "END" and cleaned and cleaned[-1].startswith("%"):
            cleaned.pop()
            continue
        cleaned.append(line)
    return cleaned, elements


def orca_calculation_record(files, version, atomic_numbers=None, coordinates=None):
    """The record of an ORCA calculation's settings, with its fingerprints.

    Parameters
    ----------
    files : {str: str or bytes}
        The input files as run: orca.inp and any others, such as basis.bas.
        Binary files, e.g. an orca.gbw for the guess, are not settings.
    version : str
        The ORCA version, e.g. "6.1.1" (see :func:`orca_version`).
    atomic_numbers, coordinates : optional
        The geometry, for its fingerprint.

    Returns
    -------
    dict
        See :func:`seamm_util.fingerprint.calculation_record`.
    """
    lines, elements = settings_text(files.get("orca.inp", ""))
    if isinstance(files.get("basis.bas"), str):
        for symbol, block in basis_file_blocks(files["basis.bas"]).items():
            elements.setdefault(symbol, []).extend(block)
    return calculation_record(
        "ORCA",
        version,
        {"orca.inp": lines},
        element_blocks=elements,
        atomic_numbers=atomic_numbers,
        coordinates=coordinates,
    )
