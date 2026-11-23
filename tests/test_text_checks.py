# -*- coding: utf-8 -*-

"""Flowchart variables in the ORCA step's free text ('extra keywords', 'extra
blocks'). ChemAI 5939 (2026-11-21) wrote '$Fz' from a Foreach loop into %scf
EField verbatim and ORCA failed in every iteration.

Option B (Paul, 2026-11-23): both fields expand SEAMM's ${expression[:format]}
templates (seamm 2026.11.23 lets the multi-line extra blocks opt in), and
text_checks catches what would still reach ORCA: a bare '$name' and '${name}' of a
variable with units.

* Unit tests of text_checks.check_text.
* Both fields expand through SEAMM's parameters; an undefined variable is an error.
* Every sub-step checks the text before it expands it.
* Run-path tests (seamm_exec.testing, fake ORCA): a Foreach over Fz writes each
  value into its iteration's orca.inp; a bare '$Fz' stops the step before ORCA runs.
"""

import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

from seamm import templates
from seamm_util import Q_

from orca_step.text_checks import check_text, RESERVED, TextCheckError

SOURCE = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "run_path"

FZ = Q_(0.001, "hartree/(e*bohr)")


# ------------------------------------------------------------------ the checks
@pytest.mark.parametrize(
    "text",
    [
        "",
        "%scf\n  MaxIter 400\nend",
        "EField 0, 0, ${Fz.m_as('hartree/(e*bohr)')}",
        'EField 0, 0, ${Fz.m_as("hartree/(e*bohr)"):.6f}',
        "EField 0, 0, ${x}",
        "EField 0, 0, ${x:.4f} ${2 * x}",
        "$new_job\n! HF def2-SVP",
        "$NEW_JOB",
        "literal $${x} stays",
        "$kw",  # a whole value: the variable's own value
        "  $kw  ",
        "$${",
    ],
)
def test_allowed(text):
    check_text(text, {"Fz": FZ, "x": 0.5, "kw": "TightSCF"})


@pytest.mark.parametrize(
    "text, suggestion",
    [
        ("EField 0, 0, $x", "${x}"),
        ("$x and more", "${x}"),
        ("SlowConv $kw", "${kw}"),
        ("EField 0, 0, $Fz", '${Fz.m_as("hartree / bohr / elementary_charge")}'),
        ("%scf $undefined end", "${undefined}"),
    ],
)
def test_bare_dollar_names_are_errors(text, suggestion):
    variables = {"Fz": FZ, "x": 0.5, "kw": "TightSCF"}
    with pytest.raises(TextCheckError, match="only '\\$\\{...\\}' is") as e:
        check_text(text, variables, field="'extra blocks'", where="Step 2.1: Energy")
    message = str(e.value)
    assert suggestion in message
    assert message.startswith("Step 2.1: Energy:")
    assert "'extra blocks'" in message


@pytest.mark.parametrize("text", ["${Fz}", "EField 0, 0, ${Fz}", "${ Fz :.4f}", "$Fz"])
def test_variables_with_units_need_m_as(text):
    with pytest.raises(TextCheckError, match="has units") as e:
        check_text(text, {"Fz": FZ})
    assert 'Fz.m_as("hartree / bohr / elementary_charge")' in str(e.value)


def test_dimensionless_quantity_is_fine():
    check_text("${r}", {"r": Q_(2.0, "")})


def test_template_syntax_error_is_reported():
    with pytest.raises(TextCheckError, match="not closed"):
        check_text("EField ${x", {"x": 1})


def test_reserved_list():
    assert "new_job" in RESERVED


# ---------------------------------------------------- through SEAMM's parameters
def test_both_fields_expand_through_seamm():
    from orca_step.energy_parameters import EnergyParameters

    P = EnergyParameters()
    assert P["extra keywords"].expands
    assert P["extra blocks"].expands
    P["extra keywords"].value = "SlowConv ${grid}"
    P["extra blocks"].value = (
        "%scf\n  EField 0, 0, ${Fz.m_as('hartree/(e*bohr)'):.4f}\nend\n$${literal}"
    )
    values = P.current_values_to_dict(
        context={"Fz": Q_(0.1, "V/Å"), "grid": "DefGrid3"}
    )
    assert values["extra keywords"] == "SlowConv DefGrid3"
    assert values["extra blocks"] == (
        f"%scf\n  EField 0, 0, {0.1 / 51.4220674763:.4f}\nend\n${{literal}}"
    )


def test_undefined_variable_names_the_field():
    from orca_step.energy_parameters import EnergyParameters

    P = EnergyParameters()
    P["extra blocks"].value = "%scf\n  EField 0, 0, ${Fy}\nend"
    with pytest.raises(templates.TemplateError, match="extra blocks"):
        P.current_values_to_dict(context={"Fz": 0.001})


# ------------------------------------------------------------- the sub-steps
@pytest.mark.parametrize(
    "module, name",
    [
        ("energy", "Energy"),
        ("optimization", "Optimization"),
        ("frequencies", "Frequencies"),
        ("bsse", "BSSE"),
    ],
)
def test_every_step_checks_before_expanding(module, name):
    """Optimization and Frequencies expand their parameters before Energy.run, and
    BSSE has its own run(): each must check the text first."""
    import importlib

    cls = getattr(importlib.import_module(f"orca_step.{module}"), name)
    source = inspect.getsource(cls.run)
    assert "self.check_text_fields()" in source
    assert source.index("self.check_text_fields()") < source.index(
        "current_values_to_dict"
    )


def test_check_text_fields(monkeypatch):
    import seamm
    from orca_step.energy import Energy
    from orca_step.energy_parameters import EnergyParameters

    step = object.__new__(type("S", (Energy,), {"header": "Step 2.1: Energy"}))
    step.parameters = EnergyParameters()
    monkeypatch.setattr(seamm, "flowchart_variables", SimpleNamespace(_data={"x": 1}))
    step.parameters["extra blocks"].value = "EField 0, 0, ${x}"
    step.check_text_fields()
    step.parameters["extra keywords"].value = "SlowConv $x"
    with pytest.raises(TextCheckError, match="Step 2.1: Energy.*'extra keywords'"):
        step.check_text_fields()


# ------------------------------------------------------------------ run path
LOOP_SPEC = """\
title: ORCA extra blocks with a loop variable
steps:
- Water: {{}}
- Loop:
    type: Foreach
    variable: Fz
    values: 0.0 0.001 0.002
    body:
    - ORCA:
        steps:
        - Energy:
            use model chemistry: 'no'
            method: DFT
            basis: def2-SVP
            extra blocks: |-
              %scf
                EField 0.0, 0.0, {token}
              end
"""


def run(tmp_path, token):
    from seamm_exec.testing import fake_program, run_spec

    code = fake_program(
        tmp_path / "bin" / "orca",
        files=[DATA / "orca.property.txt", DATA / "orca.bibtex"],
        stdout=DATA / "orca.out",
    )
    return run_spec(
        tmp_path,
        LOOP_SPEC.format(token=token),
        inis={"orca": {"code": code}},
        args=["--ncores", "1"],
        source=SOURCE,
    )


def test_run_path_loop_variable_in_extra_blocks(tmp_path):
    pytest.importorskip("seamm_exec.testing")
    job = run(tmp_path, "${Fz}")
    inputs = sorted(job.rglob("orca.inp"))
    fields = sorted(
        line.split()[-1]
        for path in inputs
        for line in path.read_text().splitlines()
        if "EField" in line
    )
    assert fields == ["0.0", "0.001", "0.002"]
    assert "Caught exception in loop iteration" not in (job / "job.out").read_text()


def test_run_path_formatted_expression(tmp_path):
    pytest.importorskip("seamm_exec.testing")
    job = run(tmp_path, "${float(Fz) * 2:.4f}")
    fields = sorted(
        line.split()[-1]
        for path in job.rglob("orca.inp")
        for line in path.read_text().splitlines()
        if "EField" in line
    )
    assert fields == ["0.0000", "0.0020", "0.0040"]


def test_run_path_bare_variable_stops_before_orca(tmp_path):
    pytest.importorskip("seamm_exec.testing")
    try:
        job = run(tmp_path, "$Fz")
    except AssertionError as e:  # run_spec reports a failed job
        text = str(e)
    else:
        text = (job / "job.out").read_text()
    assert "only '${...}' is" in text
    assert "${Fz}" in text
    assert not list(tmp_path.rglob("orca.out"))
