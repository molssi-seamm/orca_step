#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for Energy.calculate_energy_of_formation, ORCA's DfE0/DfHT/DfGT
implementation via the shared seamm_thermochemistry reference database.

A missing/broken formation-energy dependency must never break an ORCA job,
so most of this suite is about the graceful-degradation paths, not just the
happy path.
"""

import sys
from types import SimpleNamespace

import pytest

import orca_step

seamm_thermochemistry = pytest.importorskip("seamm_thermochemistry")

pytestmark = pytest.mark.skipif(
    not seamm_thermochemistry.DEFAULT_DB_PATH.exists(),
    reason=f"{seamm_thermochemistry.DEFAULT_DB_PATH} not built -- run "
    "seamm_thermochemistry/scripts/build_prototype_db.py",
)


def _fake_configuration(atomic_numbers):
    configuration = SimpleNamespace(
        atoms=SimpleNamespace(atomic_numbers=list(atomic_numbers))
    )
    configuration.PC_iupac_name = lambda fallback=None: fallback
    return configuration


WATER_ATOMIC_NUMBERS = [8, 1, 1]

# PWLDA/def2-TZVP is one of the functionals imported into the prototype
# ThermoDB for "orca" (see seamm_thermochemistry/scripts/build_prototype_db.py
# and the ORCA atom-energy-scan import), so lookups against it succeed.
DFT_P = {
    "use model chemistry": "no",
    "method": "DFT",
    "functional": "PWLDA",
    "basis": "def2-TZVP",
    "basis source": "ORCA internal",
    "basis set extrapolation": "none",
}


def test_energy_of_formation_from_energy_alone():
    # No temperature -- the case a plain single-point Energy job hits.
    node = orca_step.Energy()
    configuration = _fake_configuration(WATER_ATOMIC_NUMBERS)
    data = {"energy": -76.0}

    text = node.calculate_energy_of_formation(DFT_P, data, configuration)

    assert "E atomization" in data
    assert "DfE0" in data
    assert "DfHT" not in data  # no temperature/enthalpy available
    assert "DfGT" not in data
    assert "Thermochemistry of" in text
    assert f"{data['DfE0']:.4f}" in text


def test_energy_of_formation_adds_dfht_dfgt_when_available():
    node = orca_step.Energy()
    configuration = _fake_configuration(WATER_ATOMIC_NUMBERS)
    data = {"energy": -76.0, "enthalpy": -199990.0, "gibbs energy": -200040.0}

    node.calculate_energy_of_formation(DFT_P, data, configuration, temperature=298.15)

    assert "DfHT" in data
    assert "DfGT" in data


def test_energy_of_formation_no_temperature_skips_dfht_dfgt():
    # Even with enthalpy/gibbs energy present, no temperature -> no DfHT/DfGT
    # (calculate_energy_of_formation's temperature kwarg gates both).
    node = orca_step.Energy()
    configuration = _fake_configuration(WATER_ATOMIC_NUMBERS)
    data = {"energy": -76.0, "enthalpy": -199990.0, "gibbs energy": -200040.0}

    node.calculate_energy_of_formation(DFT_P, data, configuration)

    assert "DfHT" not in data
    assert "DfGT" not in data


def test_energy_of_formation_aliases_slash_in_functional():
    """revDSD-PBEP86-D4/2021 is looked up under its model-chemistry-safe
    ('/' -> '_') database spelling -- see orca_step.mc_method_alias."""
    node = orca_step.Energy()
    configuration = _fake_configuration(WATER_ATOMIC_NUMBERS)
    P = {
        "use model chemistry": "no",
        "method": "DFT",
        "functional": "REVDSD-PBEP86-D4/2021",
        "basis": "def2-TZVPPD",
        "basis source": "ORCA internal",
        "basis set extrapolation": "none",
    }
    data = {"energy": -76.4}

    text = node.calculate_energy_of_formation(P, data, configuration)

    assert "DfE0" in data, text


def test_unknown_basis_reports_cleanly_no_exception():
    node = orca_step.Energy()
    configuration = _fake_configuration(WATER_ATOMIC_NUMBERS)
    P = {**DFT_P, "basis": "this-basis-does-not-exist"}
    data = {"energy": -76.0}

    text = node.calculate_energy_of_formation(P, data, configuration)

    assert "DfE0" not in data
    assert "cannot calculate" in text.lower()


def test_graceful_when_seamm_thermochemistry_not_installed(monkeypatch):
    # The well-known trick: sys.modules[name] = None makes `import name`
    # raise ImportError without actually uninstalling anything.
    monkeypatch.setitem(sys.modules, "seamm_thermochemistry", None)

    node = orca_step.Energy()
    configuration = _fake_configuration(WATER_ATOMIC_NUMBERS)
    data = {"energy": -76.0}

    text = node.calculate_energy_of_formation(DFT_P, data, configuration)

    assert "DfE0" not in data
    assert "not installed" in text


def test_scalar_summary_orders_formation_quantities_first(monkeypatch):
    """The formation-referenced quantities (most chemist-relevant first) come
    before atomization energy, which comes before the raw electronic-
    structure energies -- Paul's requested reordering."""
    from orca_step import energy as energy_mod

    captured = {}
    original_tabulate = energy_mod.tabulate

    def spy(rows, **kwargs):
        captured["rows"] = rows
        return original_tabulate(rows, **kwargs)

    monkeypatch.setattr(energy_mod, "tabulate", spy)

    node = orca_step.Energy()
    node._id = ("1",)
    props = {
        "energy": -76.363,
        "scf energy": -76.217,
        "E atomization": 964.31,
        "DfE0": -285.45,
        "DfHT": -237.73,
        "DfGT": -224.41,
        "zero point energy": 56.38,
        "formation temperature": 298.15,
    }
    node._print_scalar_summary(props)

    labels = [row[0] for row in captured["rows"]]
    assert labels.index("Enthalpy of formation (298.15 K)") == 0
    assert labels.index("Gibbs energy of formation (298.15 K)") == 1
    assert labels.index("Energy of formation (0 K)") == 2
    assert labels.index("Atomization energy") == 3
    assert labels.index("Zero-point energy") == 4
    assert labels.index("Total energy") > labels.index("Zero-point energy")


def test_graceful_when_database_not_built(monkeypatch, tmp_path):
    monkeypatch.setattr(
        seamm_thermochemistry.db, "DEFAULT_DB_PATH", tmp_path / "does_not_exist.db"
    )

    node = orca_step.Energy()
    configuration = _fake_configuration(WATER_ATOMIC_NUMBERS)
    data = {"energy": -76.0}

    text = node.calculate_energy_of_formation(DFT_P, data, configuration)

    assert "DfE0" not in data
    assert "not built" in text
