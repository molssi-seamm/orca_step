# -*- coding: utf-8 -*-

"""ORCA's run path, end to end: a flowchart built from a spec and run by
``run_flowchart``, through the Energy sub-step -- the input, ``orca.ini``, the
resolver, the task, the output and its analysis (seamm_exec.testing).

* With a fake ``orca`` replaying a real run's output, so the path runs in CI.
* With the real ORCA when it is installed (skipped otherwise, as in CI).

One core: ORCA 6.1.1 aborts water in its start-up on 11 MPI processes ("the
number of points read from the grid does not match the expectation").
"""

from pathlib import Path

import pytest

from seamm_exec.testing import fake_program, find_program, run_spec, table_value

SOURCE = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "run_path"

SPEC = """\
title: ORCA run-path test
steps:
- Water: {}
- ORCA:
    steps:
    - Energy:
        use model chemistry: 'no'
        method: DFT
        basis: def2-SVP
"""

ENERGY = -76.32100225  # E_h, B3LYP/def2-SVP, ORCA 6.1.1


def check_job(job):
    run = job / "2" / "1"
    text = (run / "orca.inp").read_text()
    assert "B3LYP" in text and "def2-SVP" in text and "* xyz 0 1" in text
    assert "ORCA TERMINATED NORMALLY" in (run / "orca.out").read_text()
    assert table_value(job / "job.out", "Total energy") == pytest.approx(
        ENERGY, abs=1.0e-6
    )


def test_run_path_with_a_fake_orca(tmp_path):
    code = fake_program(
        tmp_path / "bin" / "orca",
        files=[DATA / "orca.property.txt", DATA / "orca.bibtex"],
        stdout=DATA / "orca.out",
    )
    job = run_spec(
        tmp_path,
        SPEC,
        inis={"orca": {"code": code}},
        args=["--ncores", "1"],
        source=SOURCE,
    )
    check_job(job)
    # The fake ran, not an ORCA on the PATH: the real one writes orca.gbw. (A
    # sub-step's task once lost the root, so orca.ini was not read.)
    assert not (job / "2" / "1" / "orca.gbw").exists()


@pytest.mark.skipif(find_program("orca") is None, reason="ORCA is not installed")
def test_run_path_with_the_real_orca(tmp_path):
    job = run_spec(
        tmp_path,
        SPEC,
        inis={"orca": {"code": find_program("orca")}},
        args=["--ncores", "1"],
        source=SOURCE,
    )
    check_job(job)
