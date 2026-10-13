# -*- coding: utf-8 -*-

"""The retry of a failed SCF with a larger SThresh (mbe, 2026-10-12).

The retry is a few shell commands after the ORCA run; these tests run them for
real with a stand-in for ORCA that fails the SCF, as ORCA 6.1.1 does for a
nearly linearly dependent def2-QZVPPD cluster, unless the input raises SThresh.
"""

import stat
import subprocess
import textwrap

import pytest

from seamm_exec import Geometry

from orca_step.orca_base import (
    STHRESH_FAILED_OUTPUT,
    STHRESH_RETRY_MARKER,
    _RETRY_INPUT,
    sthresh_retry_command,
    sthresh_retry_input,
    sthresh_retry_note,
    sthresh_retry_value,
)

INPUT = textwrap.dedent("""\
    ! R2SCAN D4 AutoAux def2-QZVPPD TIGHTSCF DEFGRID3 EnGrad
    %pal nprocs 4 end
    %maxcore 2000
    * xyz 0 1
    O        7.1467000000       5.2974900000       1.4502200000
    *
    """)

# What ORCA 6.1.1 printed for the failing pair (FEC32_r095_T298_f05 c1-d03_04)
FAILED = """\
Diagonalization of the overlap matrix:
Smallest eigenvalue                        ... {eigenvalue}
Threshold for overlap eigenvalues          ... 1.000e-07
Number of eigenvalues below threshold      ... 0

IMPORTANT:
Smallest eigenvalue of S below the DiffSThresh (1.000e-06) while
               *        SCF NOT CONVERGED AFTER 1305 CYCLES         *
ORCA finished by error termination in LEANSCF
"""

# The stand-in: the SCF converges only with SThresh 1e-6 in the input, or when
# it is told to ($CONVERGE). It notes whether it found an old orca.gbw, which
# ORCA would start from, and leaves one behind. Exits 0, as ORCA does.
FAKE_ORCA = """\
#!/bin/sh
if [ -e orca.gbw ]; then echo "found orca.gbw" >> calls.txt; fi
echo run >> calls.txt
echo junk > orca.gbw
if grep -q "SThresh 1.000e-06" "$1" || [ -n "$CONVERGE" ]; then
    echo "FINAL SINGLE POINT ENERGY      -883.266910000000"
    echo "****ORCA TERMINATED NORMALLY****"
else
    cat failed.txt
fi
exit 0
"""


def _run(tmp_path, eigenvalue="1.438e-07", sthresh=1.0e-06, env=None):
    """Run the task's command with the retry in ``tmp_path``; return the shell's
    return code."""
    orca = tmp_path / "orca"
    orca.write_text(FAKE_ORCA)
    orca.chmod(orca.stat().st_mode | stat.S_IEXEC)
    (tmp_path / "failed.txt").write_text(FAILED.format(eigenvalue=eigenvalue))
    (tmp_path / "orca.inp").write_text(INPUT)
    (tmp_path / _RETRY_INPUT).write_text(sthresh_retry_input(INPUT, sthresh))
    cmd = ["{code}", "orca.inp", ">", "orca.out", "2>", "orca.err"]
    cmd += sthresh_retry_command(sthresh)
    # As seamm_exec does: join, then fill {...} from the configuration until
    # nothing changes (local.py) -- the commands must survive that.
    command = " ".join(cmd)
    while True:
        filled = command.format(code=str(orca))
        if filled == command:
            break
        command = filled
    return subprocess.run(
        command, shell=True, cwd=tmp_path, env=env, executable="/bin/sh"
    ).returncode


def test_the_value():
    assert sthresh_retry_value("1.0e-06") == 1.0e-06
    assert sthresh_retry_value(1e-5) == 1e-5
    for off in ("no", "none", "", None, "default"):
        assert sthresh_retry_value(off) is None
    # It must raise the threshold, or it could remove nothing
    assert sthresh_retry_value("1e-7") is None
    assert sthresh_retry_value("1e-6", sthresh="1e-6") is None
    assert sthresh_retry_value("1e-5", sthresh="1e-6") == 1e-5
    with pytest.raises(ValueError):
        sthresh_retry_value("tight")


def test_the_input_overrides_earlier_blocks():
    text = sthresh_retry_input(
        INPUT.replace("%maxcore", "%scf SThresh 1e-7 end\n%maxcore"), 1e-6
    )
    lines = text.splitlines()
    i = lines.index("%scf SThresh 1.000e-06 end")
    assert lines[i + 1] == "* xyz 0 1"
    assert lines.index("%scf SThresh 1e-7 end") < i
    with pytest.raises(ValueError):
        sthresh_retry_input("! HF\n", 1e-6)


def test_a_failed_scf_is_retried(tmp_path):
    assert _run(tmp_path) == 0
    out = (tmp_path / "orca.out").read_text()
    assert "ORCA TERMINATED NORMALLY" in out
    assert "SCF NOT CONVERGED" in (tmp_path / STHRESH_FAILED_OUTPUT).read_text()
    assert "SThresh 1.000e-06" in (tmp_path / "orca.inp").read_text()
    assert "SThresh" not in (tmp_path / "orca.inp.scf-failed").read_text()
    # The retry started afresh, not from the failed run's orbitals
    assert (tmp_path / "calls.txt").read_text().split("\n")[:2] == ["run", "run"]
    marker = (tmp_path / STHRESH_RETRY_MARKER).read_text()
    assert marker.split()[-1] == "1.438e-07"
    note = sthresh_retry_note(marker, 1e-6)
    assert "1.44e-07" in note and "SThresh 1.0e-06" in note


def test_no_retry_if_nothing_would_be_removed(tmp_path):
    """The smallest eigenvalue is above the retry's threshold."""
    _run(tmp_path, eigenvalue="5.306e-06")
    assert "SCF NOT CONVERGED" in (tmp_path / "orca.out").read_text()
    assert not (tmp_path / STHRESH_RETRY_MARKER).exists()
    assert (tmp_path / "calls.txt").read_text().split() == ["run"]


def test_no_retry_if_the_scf_converged(tmp_path):
    _run(tmp_path, env={"CONVERGE": "1", "PATH": "/usr/bin:/bin"})
    assert "ORCA TERMINATED NORMALLY" in (tmp_path / "orca.out").read_text()
    assert not (tmp_path / STHRESH_RETRY_MARKER).exists()
    assert not (tmp_path / STHRESH_FAILED_OUTPUT).exists()


def test_the_batch_task(tmp_path):
    """On by default; the fingerprint is that of the calculation asked for; the
    options turn it off, change its SThresh, or set every run's SThresh."""
    from orca_step.batch import get_task

    water = Geometry(
        [8, 1, 1], [[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]]
    )
    mc = {
        "level": "ORCA:DFT@B3LYP/def2-SVP",
        "method": "B3LYP",
        "basis": "def2-SVP",
        "options": {"mdi_method_arg": "B3LYP", "mdi_basis_arg": "def2-SVP"},
    }
    task = get_task(water, mc, key="w")
    assert "SThresh 1.000e-06" in task.files[_RETRY_INPUT]
    assert STHRESH_RETRY_MARKER in task.return_files
    assert any("SCF NOT CONVERGED" in c for c in task.cmd)

    off = get_task(water, mc, key="w", options={"sthresh_retry": "no"})
    assert _RETRY_INPUT not in off.files
    assert off.fingerprint == task.fingerprint

    wider = get_task(water, mc, key="w", options={"sthresh_retry": 1e-5})
    assert "SThresh 1.000e-05" in wider.files[_RETRY_INPUT]
    assert wider.fingerprint == task.fingerprint

    # A campaign-wide SThresh is part of the calculation, so of its fingerprint
    mc["options"]["sthresh"] = "1e-6"
    every = get_task(water, mc, key="w")
    assert "%scf SThresh 1.000e-06 end" in every.files["orca.inp"]
    assert every.fingerprint != task.fingerprint
    assert _RETRY_INPUT not in every.files  # 1e-6 again would remove nothing


def test_the_parameters():
    from orca_step import BSSEParameters, EnergyParameters, FrequenciesParameters

    P = EnergyParameters()
    assert P.applies("retry sthresh", {"scf retry": "yes"})
    assert not P.applies("retry sthresh", {"scf retry": "no"})
    # The sub-steps that run SCFs all have it
    assert BSSEParameters().applies("scf retry", {})
    assert FrequenciesParameters().applies("scf retry", {})


def test_the_step_task(tmp_path, monkeypatch):
    """The ORCA step's own task: the retry, and the same fingerprint."""
    import orca_step.orca_base as base
    from orca_step.orca_base import ORCABase

    monkeypatch.setattr(base, "_atomic_numbers", lambda configuration: [1, 1])

    class Node:
        indent = ""
        directory = str(tmp_path)

        def _resources(self):
            return 1, 1000

        def geometry_block(self, *args, **kwargs):
            return "* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*"

    node = Node()
    plain = ORCABase.orca_job_task(node, "HF def2-SVP", None, 0, 1)
    task = ORCABase.orca_job_task(node, "HF def2-SVP", None, 0, 1, sthresh_retry=1e-6)
    assert task.fingerprint == plain.fingerprint
    assert _RETRY_INPUT in task.files and _RETRY_INPUT not in plain.files
    assert STHRESH_RETRY_MARKER in task.return_files


def test_the_note_needs_the_marker():
    assert sthresh_retry_note("") is None
    assert sthresh_retry_note(None) is None
    assert "nearly singular" in sthresh_retry_note("Smallest eigenvalue ... x")


def test_maxcore_keeps_headroom():
    """%maxcore is 80 % of the memory per core given, as for the default
    resources (mbe 2026-10-12: 2 ranks x 2 GB with %maxcore 2000 were OOM-killed)."""
    import seamm_exec
    from orca_step.batch import get_task

    water = Geometry(
        [8, 1, 1], [[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]]
    )
    mc = {
        "level": "ORCA:DFT@B3LYP/def2-SVP",
        "method": "B3LYP",
        "basis": "def2-SVP",
        "options": {"mdi_method_arg": "B3LYP", "mdi_basis_arg": "def2-SVP"},
    }
    task = get_task(
        water,
        mc,
        key="w",
        resources=seamm_exec.Resources(ntasks=2, mem_per_cpu=2_000_000_000),
    )
    assert "%maxcore 1600\n" in task.files["orca.inp"]
    assert "%pal nprocs 2 end" in task.files["orca.inp"]
    default = get_task(water, mc, key="w")
    assert "%maxcore 2000\n" in default.files["orca.inp"]
    assert default.resources.mem_per_cpu == 2_500_000_000
    # %maxcore is not part of what the task is: finished tasks are not rerun
    assert task.fingerprint == default.fingerprint
