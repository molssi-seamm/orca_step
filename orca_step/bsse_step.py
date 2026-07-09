# -*- coding: utf-8 -*-

"""Stevedore helper for the ORCA BSSE (counterpoise) sub-step."""

import orca_step


class BSSEStep(object):
    my_description = {
        "description": "Counterpoise (BSSE) corrected energy and gradient with ORCA",
        "group": "Calculation",
        "name": "BSSE",
    }

    def __init__(self, flowchart=None, gui=None):
        pass

    def description(self):
        return BSSEStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        return orca_step.BSSE(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        return orca_step.TkBSSE(canvas=canvas, **kwargs)
