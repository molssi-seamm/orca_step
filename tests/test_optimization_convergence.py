# -*- coding: utf-8 -*-

"""Whether an ORCA optimization converged (Paul, 2026-11-23: a precondition for E1).

Before, a TightOpt that did not converge passed silently: a warning at normal print
level (not in job.out from inside a Loop) and the last geometry stored. Now:

* the properties 'optimization converged' (1/0) and 'optimization cycles' can be
  stored with the results;
* 'If not converged': 'warn and continue' (the default) warns in job.out; 'fail the
  step' stops the step before anything is stored, so a Loop counts the iteration.

The run-path tests replay a real ORCA 6.1.1 run (water, r2SCAN-3c TightOpt, 5
cycles) through a fake ORCA; "not converged" removes ORCA's convergence banner.
"""

from pathlib import Path
import sqlite3

import pytest

from orca_step.optimization import Optimization

SOURCE = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "run_path_opt"
FILES = ("orca.xyz", "orca.engrad", "orca.property.txt", "orca.bibtex")

SPEC = """\
title: ORCA optimization convergence
steps:
- Water: {{}}
- ORCA:
    steps:
    - Optimization:
        use model chemistry: 'no'
        method: R2SCAN-3C
        optimization convergence: TightOpt
        if not converged: {handling}
        results:
          energy: {{property: "total energy#ORCA#{{model}}"}}
          optimization converged: {{property: "optimization converged#ORCA#{{model}}"}}
          optimization cycles: {{property: "optimization cycles#ORCA#{{model}}"}}
"""


def test_status_from_the_output(tmp_path):
    assert Optimization._optimization_status(tmp_path) == (False, 0)
    (tmp_path / "orca.out").write_text((DATA / "orca.out").read_text())
    assert Optimization._optimization_status(tmp_path) == (True, 5)
    text = _unconverged((DATA / "orca.out").read_text())
    (tmp_path / "orca.out").write_text(text)
    assert Optimization._optimization_status(tmp_path) == (False, 5)


def _unconverged(text):
    return "\n".join(
        line
        for line in text.splitlines()
        if "THE OPTIMIZATION HAS CONVERGED" not in line and "HURRAY" not in line
    )


def _run(tmp_path, handling, converged=True):
    from seamm_exec.testing import fake_program, run_spec

    data = DATA
    if not converged:
        data = tmp_path / "unconverged"
        data.mkdir()
        (data / "orca.out").write_text(_unconverged((DATA / "orca.out").read_text()))
    code = fake_program(
        tmp_path / "bin" / "orca",
        files=[DATA / name for name in FILES],
        stdout=data / "orca.out",
    )
    return run_spec(
        tmp_path,
        SPEC.format(handling=handling),
        inis={"orca": {"code": code}},
        args=["--ncores", "1"],
        source=SOURCE,
    )


def _properties(job):
    db = sqlite3.connect(f"file:{job / 'seamm.db'}?mode=ro", uri=True)
    try:
        return {
            name: value
            for name, value in db.execute(
                "SELECT property.name, int_data.value FROM int_data, property"
                " WHERE int_data.property = property.id"
            )
        }
    finally:
        db.close()


def _value(properties, prefix):
    return next(v for k, v in properties.items() if k.startswith(prefix))


def test_run_path_converged(tmp_path):
    pytest.importorskip("seamm_exec.testing")
    job = _run(tmp_path, "fail the step")
    properties = _properties(job)
    assert _value(properties, "optimization converged#ORCA#") == 1
    assert _value(properties, "optimization cycles#ORCA#") == 5
    assert "did not converge" not in (job / "job.out").read_text()


def test_run_path_not_converged_warns(tmp_path):
    pytest.importorskip("seamm_exec.testing")
    job = _run(tmp_path, "warn and continue", converged=False)
    properties = _properties(job)
    assert _value(properties, "optimization converged#ORCA#") == 0
    assert _value(properties, "optimization cycles#ORCA#") == 5
    out = " ".join((job / "job.out").read_text().split())  # the text is wrapped
    assert "did not converge in 5 cycles" in out
    assert "The last geometry is kept" in out


def test_run_path_not_converged_fails(tmp_path):
    pytest.importorskip("seamm_exec.testing")
    try:
        job = _run(tmp_path, "fail the step", converged=False)
    except AssertionError as e:  # run_spec reports a failed job
        text = str(e)
    else:
        pytest.fail("the job should have failed:\n" + (job / "job.out").read_text())
    assert "did not converge in 5" in text
    assert "as 'If not converged' asks" in text
