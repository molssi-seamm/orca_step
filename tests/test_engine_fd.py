# -*- coding: utf-8 -*-

"""The ORCA MDI engine's support for finite differences (2026-11-26, after the
warm-engine FD Hessian of 4-hydroxy-TEMPO, SEAMM job 6074): the reference
orbitals kept for every displaced point's guess, and the tight settings."""

import importlib.util
from pathlib import Path

ENGINE = Path(__file__).parents[1] / "orca_step" / "data" / "orca_mdi.py"


def _engine():
    spec = importlib.util.spec_from_file_location("_orca_mdi_engine_fd", ENGINE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_tight_method_words():
    engine = _engine()
    assert engine.tight_method_words("R2SCAN-3C") == "R2SCAN-3C VERYTIGHTSCF"
    # the engine's TIGHTSCF then gives way: ORCA refuses two SCF presets
    words = engine.tight_method_words("B3LYP")
    assert engine.accuracy_keywords(words) == "DEFGRID3"
    # an SCF preset already chosen is kept
    assert engine.tight_method_words("B3LYP EXTREMESCF") == "B3LYP EXTREMESCF"


def test_the_reference_guess_survives_the_cleanup(tmp_path):
    engine = _engine()
    for name in ("orca.gbw", "orca.out", "orca.engrad", engine.REFERENCE_GUESS):
        (tmp_path / name).write_text("x")
    removed = engine.clean_workdir(tmp_path)
    assert sorted(removed) == ["orca.engrad", "orca.out"]
    assert (tmp_path / engine.REFERENCE_GUESS).exists()
    assert (tmp_path / "orca.gbw").exists()


def test_the_commands_are_registered():
    """SAVE_GUESS, RESTORE_GUESS, TIGHT_SCF and NORMAL_SCF are advertised, so a
    driver can ask with MDI_Check_Command_Exists."""
    text = ENGINE.read_text()
    block = text.split("commands = [", 1)[1].split("]", 1)[0]
    loop = text.split("while True:", 1)[1]
    for command in ("SAVE_GUESS", "RESTORE_GUESS", "TIGHT_SCF", "NORMAL_SCF"):
        assert f'"{command}"' in block  # registered
        assert f'"{command}"' in loop  # and handled
