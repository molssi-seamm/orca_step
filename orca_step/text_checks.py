# -*- coding: utf-8 -*-

"""Checks of the free text written to the ORCA input.

The 'extra keywords' and 'extra blocks' are written to orca.inp as typed, after
SEAMM has replaced each ``${expression[:format]}`` in them by its value, as in any
text parameter (e.g. ``EField 0, 0, ${Fz}`` or ``${Fz.m_as("hartree/(e*bohr)")}``).
``$${`` writes a literal ``${``.

Two mistakes would still reach ORCA unnoticed, so the text is checked as typed,
before it is expanded:

1. A bare ``$name`` is not expanded (only ``${...}`` is), so it would be written as
   is and ORCA would fail on it. It is an error that suggests ``${name}``. ORCA's
   own ``$new_job`` is allowed, and so is a value that is only ``$name``: that
   whole-value form gives the variable's value, e.g. keywords held in a variable.
2. ``${name}`` of a variable with units would be written with the units' text
   (``0.001 hartree / bohr / elementary_charge``). It is an error that suggests
   ``${name.m_as("<units>")}``, which writes the number in the given units.
"""

import re

from seamm import templates

#: ORCA's own tokens starting with '$', which are written as they are
RESERVED = ("new_job",)

# '$name' not preceded by '$' and not followed by '{' or '$'
_BARE = re.compile(r"(?<!\$)\$(?![${])([A-Za-z_][A-Za-z0-9_]*)")
# The inside of '${...}' that is just a name, possibly with !conversion or :format
_NAME = re.compile(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*(?:[!:].*)?", re.DOTALL)


class TextCheckError(RuntimeError):
    """The free text of the ORCA input has a variable that cannot be written."""


def _has_units(value):
    """Whether a value is a pint quantity with real units."""
    if not (hasattr(value, "magnitude") and hasattr(value, "units")):
        return False
    return not getattr(value, "dimensionless", False)


def _units_hint(name, value):
    return f'${{{name}.m_as("{value.units:D}")}}'


def check_text(text, variables, field="text", where="ORCA"):
    """Check the free text of the ORCA input before it is expanded.

    Parameters
    ----------
    text : str
        The text as typed, e.g. the 'extra blocks'.
    variables : dict
        The flowchart variables, name -> value.
    field : str
        The parameter the text comes from, for error messages.
    where : str
        The step, for error messages, e.g. "Step 3.2: Energy".

    Raises
    ------
    TextCheckError
        For a bare '$name', or '${name}' of a variable with units.
    """
    if not isinstance(text, str) or "$" not in text:
        return

    def fail(message):
        raise TextCheckError(f"{where}: {message} in the {field}.")

    # A whole value that is just '$name' gives the variable's value
    whole = re.fullmatch(r"\s*\$([A-Za-z_][A-Za-z0-9_]*)\s*", text)
    if whole is not None:
        name = whole.group(1)
        if _has_units(variables.get(name)):
            fail(
                f"the variable '{name}' has units ({variables[name].units:D}); "
                f"write e.g. '{_units_hint(name, variables[name])}' for its number"
            )
        return

    for match in _BARE.finditer(text):
        name = match.group(1)
        if name.lower() in RESERVED:
            continue
        hint = f"${{{name}}}"
        if _has_units(variables.get(name)):
            hint = _units_hint(name, variables[name])
        fail(
            f"'${name}' is not replaced by the variable's value: only '${{...}}' is. "
            f"Write '{hint}'"
        )

    try:
        pieces = templates.parse(text)
    except templates.TemplateError as e:
        fail(str(e))
    for inside, is_field in pieces:
        if not is_field:
            continue
        match = _NAME.fullmatch(inside)
        if match is None:
            continue
        name = match.group(1)
        if _has_units(variables.get(name)):
            fail(
                f"the variable '{name}' has units ({variables[name].units:D}), which "
                f"'${{{inside}}}' would write into the input. Write "
                f"'{_units_hint(name, variables[name])}' for its number in the units "
                "ORCA expects"
            )
