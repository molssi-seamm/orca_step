"""The auxiliary (fitting) basis: AutoAux, none, or explicit sets (science,
2026-10-23: AutoAux made the auxiliary basis of a 104-atom Li+ pair linearly
dependent, which ORCA's SHARK cannot handle)."""

import pytest

from orca_step.orca_base import auxiliary_basis_words


@pytest.mark.parametrize(
    "value, words",
    [
        (None, "AutoAux"),
        ("", "AutoAux"),
        ("autoaux", "AutoAux"),
        ("none", ""),
        ("def2/J", "def2/J"),
        ("def2/J def2-TZVPP/C", "def2/J def2-TZVPP/C"),
        # MBE's task options separate settings with spaces
        ("def2/J,def2-QZVPP/C", "def2/J def2-QZVPP/C"),
        ("def2/JK+cc-pVTZ/C", "def2/JK cc-pVTZ/C"),
    ],
)
def test_words(value, words):
    assert auxiliary_basis_words(value) == words


@pytest.mark.parametrize(
    "value, message",
    [
        ("AutoAux def2/J", "combines AutoAux"),
        ("def2-TZVPP", "not the name of an auxiliary basis"),
        ("def2/J RIJCOSX", "not the name of an auxiliary basis"),
    ],
)
def test_refused(value, message):
    with pytest.raises(ValueError, match=message):
        auxiliary_basis_words(value)


def test_energy_parameters_check_the_aux_basis():
    """The flowchart tools and the check before a job see a bad name at once."""
    from orca_step import EnergyParameters

    P = EnergyParameters()
    assert not [p for p in P.problems() if "auxiliary" in p]
    P["auxiliary basis"].value = "def2/J def2-TZVPP/C"
    assert not [p for p in P.problems() if "auxiliary" in p]
    P["auxiliary basis"].value = "def2-TZVPP"
    assert any("auxiliary basis" in p for p in P.problems())
    P["auxiliary basis"].value = "$aux"
    assert not [p for p in P.problems() if "auxiliary" in p]
