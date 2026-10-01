# -*- coding: utf-8 -*-

"""Flowcharts saved before BSSE took N fragments (2026-08-04) can be read."""

import pytest

import orca_step
from orca_step.bsse import BSSE


def test_old_auto_fragments():
    P = orca_step.BSSEParameters(
        data={
            "fragments": {"value": "auto (2 molecules)", "units": None},
            "fragment A atoms": {"value": "", "units": None},
        }
    )
    assert P["fragments"].value == "auto (molecules)"
    assert P["fragment atoms"].value == ""


def test_old_specified_fragment_a():
    P = orca_step.BSSEParameters(
        data={
            "fragments": {"value": "specified", "units": None},
            "fragment A atoms": {"value": "1-3", "units": None},
        }
    )
    assert P["fragment atoms"].value == "1-3; rest"
    assert BSSE._parse_fragment_groups(P["fragment atoms"].value, 6) == [
        [0, 1, 2],
        [3, 4, 5],
    ]


def test_rest_group():
    assert BSSE._parse_fragment_groups("1; 2-3; rest", 5) == [[0], [1, 2], [3, 4]]
    with pytest.raises(RuntimeError, match="only be the last"):
        BSSE._parse_fragment_groups("rest; 1-3", 5)
    with pytest.raises(RuntimeError, match="no atoms are left"):
        BSSE._parse_fragment_groups("1-5; rest", 5)


def test_an_old_flowchart_reads(tmp_path):
    """The whole path: a 2.0 flowchart with the old parameters converts and loads."""
    import json

    import seamm

    node = {
        "item": "object",
        "module": "orca_step.bsse",
        "class": "BSSE",
        "version": "2026.7.10",
        "extension": "BSSE",
        "attributes": {
            "_uuid": 7,
            "_title": "BSSE",
            "extension": "BSSE",
            "parameters": {
                "__class__": "BSSEParameters",
                "__module__": "orca_step.bsse_parameters",
                "fragments": {"value": "specified", "units": None},
                "fragment A atoms": {"value": "1-3", "units": None},
            },
        },
    }
    orca = {
        "item": "object",
        "module": "orca_step.orca",
        "class": "ORCA",
        "version": "2026.7.10",
        "extension": "ORCA",
        "attributes": {
            "_uuid": 5,
            "_title": "ORCA",
            "extension": "ORCA",
            "parameters": None,
        },
        "subflowchart": {
            "nodes": [
                {
                    "module": "seamm.start_node",
                    "class": "StartNode",
                    "version": "2026.7.1",
                    "extension": None,
                    "attributes": {"_uuid": 1, "parameters": None},
                },
                node,
            ],
            "edges": [
                {
                    "node1": 1,
                    "node2": 7,
                    "edge_type": "execution",
                    "edge_subtype": "next",
                    "attributes": {},
                }
            ],
        },
    }
    data = {
        "nodes": [
            {
                "module": "seamm.start_node",
                "class": "StartNode",
                "version": "2026.7.1",
                "extension": None,
                "attributes": {"_uuid": 1, "parameters": None},
            },
            orca,
        ],
        "edges": [
            {
                "node1": 1,
                "node2": 5,
                "edge_type": "execution",
                "edge_subtype": "next",
                "attributes": {},
            }
        ],
    }
    text = (
        "!MolSSI flowchart 2.0\n#metadata\n{}\n#flowchart\n"
        + json.dumps(data)
        + "\n#end\n"
    )
    flowchart = seamm.Flowchart()
    flowchart.from_text(text)
    bsse = flowchart.get_nodes()[1].subflowchart.get_node("1").next()
    assert bsse.parameters["fragment atoms"].value == "1-3; rest"


@pytest.mark.parametrize("cls", ["EnergyParameters", "BSSEParameters"])
def test_renamed_method(cls):
    """CCSD(T)-F12D was renamed CCSD(T)-F12D/RI (ORCA rejects the bare keyword)."""
    P = getattr(orca_step, cls)(
        data={"method": {"value": "CCSD(T)-F12D", "units": None}}
    )
    assert P["method"].value == "CCSD(T)-F12D/RI"
