# -*- coding: utf-8 -*-

"""check_mpirun: a parallel ORCA run with an OpenMPI 5 mpirun is stopped with a
clear error, instead of ORCA's workers aborting in its start-up."""

import os

import pytest

from orca_step.orca_base import check_mpirun


def fake_mpirun(directory, version_text):
    """An ``mpirun`` that prints ``version_text`` for ``--version``."""
    directory.mkdir(parents=True, exist_ok=True)
    mpirun = directory / "mpirun"
    mpirun.write_text(f"#!/bin/sh\necho '{version_text}'\n")
    mpirun.chmod(0o755)
    return mpirun


OMPI5 = "mpirun (Open MPI) 5.0.8"
OMPI4 = "mpirun (Open MPI) 4.1.6"


def test_openmpi_5_on_the_path_is_stopped(tmp_path):
    fake_mpirun(tmp_path / "brew", OMPI5)
    with pytest.raises(RuntimeError, match="Open MPI 5.*OpenMPI 4.1.*library-path"):
        check_mpirun(11, {"installation": "local"}, path=str(tmp_path / "brew"))


def test_openmpi_4_on_the_path_runs(tmp_path):
    fake_mpirun(tmp_path / "ompi4", OMPI4)
    check_mpirun(11, {"installation": "local"}, path=str(tmp_path / "ompi4"))


def test_library_path_mpirun_wins_over_the_path(tmp_path):
    """mpi_env puts the bin beside library-path first, so that one is checked."""
    fake_mpirun(tmp_path / "brew", OMPI5)
    fake_mpirun(tmp_path / "orca-mpi" / "bin", OMPI4)
    (tmp_path / "orca-mpi" / "lib").mkdir()
    config = {"installation": "local", "library-path": str(tmp_path / "orca-mpi/lib")}
    check_mpirun(11, config, path=str(tmp_path / "brew"))


def test_no_mpirun_is_stopped(tmp_path):
    (tmp_path / "empty").mkdir()
    with pytest.raises(RuntimeError, match="no 'mpirun'"):
        check_mpirun(4, {}, path=str(tmp_path / "empty"))


@pytest.mark.parametrize(
    "n_cores, config",
    [
        (1, {"installation": "local"}),  # one core: no MPI
        (8, {"installation": "modules"}),  # the module supplies mpirun later
    ],
)
def test_not_checked(tmp_path, n_cores, config):
    fake_mpirun(tmp_path / "brew", OMPI5)
    check_mpirun(n_cores, config, path=str(tmp_path / "brew"))


def test_an_unidentified_mpi_is_not_blocked(tmp_path):
    fake_mpirun(tmp_path / "mpich", "HYDRA build details: Version: 4.2.0")
    check_mpirun(4, {"installation": "local"}, path=str(tmp_path / "mpich"))


@pytest.mark.skipif(
    not os.path.exists("/opt/homebrew/bin/mpirun"), reason="no Homebrew mpirun"
)
def test_this_macs_homebrew_mpirun():
    """On a Mac with Homebrew's OpenMPI 5, it is the one caught."""
    from orca_step.orca_base import _mpirun_version

    implementation, major = _mpirun_version("/opt/homebrew/bin/mpirun")
    if implementation == "Open MPI" and major >= 5:
        with pytest.raises(RuntimeError, match="Open MPI"):
            check_mpirun(11, {"installation": "local"}, path="/opt/homebrew/bin")
