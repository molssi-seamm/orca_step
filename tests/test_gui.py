# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialogs: create them and re-lay them out for every choice
that drives the layout, checking that the controls shown are those that the
parameters' rules say apply. Skipped when no display is available."""

import pytest

SUBSTEPS = ("Energy", "Optimization", "Frequencies", "BSSE")


@pytest.fixture()
def orca(request):
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw

    Pmw.initialise(root)
    yield root
    root.destroy()


def make(root, substep):
    import seamm

    flowchart = seamm.Flowchart(namespace="org.molssi.seamm.orca", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.orca.tk"
    )
    node = flowchart.create_node(substep)
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get(substep)
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    tk_node.create_dialog()
    return tk_node


LEVEL = (
    "use model chemistry",
    "method",
    "functional type",
    "functional",
    "basis set extrapolation",
    "extrapolation family",
    "basis",
    "basis source",
)


def check(tk_node):
    """The shown controls are exactly those that apply, of the ones the dialog
    lays out."""
    tk_node.reset_dialog()
    P = tk_node.node.parameters
    values = tk_node._widget_values()
    laid_out = set(LEVEL) | set(tk_node._run_detail_keys())
    if "initial guess" in laid_out:
        laid_out |= {"specified orbitals", "if wavefunction not found"}
    for key in P:
        if key == "results" or key not in tk_node:
            continue
        shown = tk_node[key].grid_info() != {}
        if shown:
            assert P.applies(key, values), f"{key} is shown but does not apply"
        elif key in laid_out:
            assert not P.applies(key, values), f"{key} applies but is not shown"
    return values


@pytest.mark.parametrize("substep", SUBSTEPS)
def test_layouts_follow_the_rules(orca, substep):
    tk_node = make(orca, substep)
    P = tk_node.node.parameters
    for use in ("yes", "no"):
        tk_node["use model chemistry"].set(use)
        check(tk_node)
        if use == "yes":
            assert tk_node["method"].grid_info() == {}
            continue
        for method in P["method"].enumeration:
            tk_node["method"].set(method)
            check(tk_node)
            # F12: only the F12 bases, from ORCA itself
            if "F12" in method.upper():
                basis = tk_node["basis"].get()
                name = basis["name"] if isinstance(basis, dict) else basis
                assert name.endswith("-F12")
                assert tk_node["basis source"].grid_info() == {}
            # DFT: the functional list is that of the functional type
            if method == "DFT":
                for ftype in P["functional type"].enumeration:
                    tk_node["functional type"].set(ftype)
                    tk_node.reset_functionals()
                    assert tk_node["functional"].get() in P.choices(
                        "functional", tk_node._widget_values()
                    )
        tk_node["method"].set("MP2")
        check(tk_node)
        can_extrapolate = P.extrapolation
        assert (tk_node["basis set extrapolation"].grid_info() != {}) == can_extrapolate
        if can_extrapolate:
            tk_node["basis set extrapolation"].set("2/3")
            check(tk_node)
            assert tk_node["basis"].grid_info() == {}
            assert tk_node["extrapolation family"].grid_info() != {}
            tk_node["basis set extrapolation"].set("none")
    if "initial guess" in P and "initial guess" not in P.unused:
        for guess in P["initial guess"].enumeration:
            tk_node["initial guess"].set(guess)
            check(tk_node)
    if substep == "BSSE":
        for fragments in P["fragments"].enumeration:
            tk_node["fragments"].set(fragments)
            check(tk_node)
