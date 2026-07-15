# -*- coding: utf-8 -*-

"""Stevedore helper for the ORCA Frequencies sub-step."""

import orca_step


class FrequenciesStep(object):
    my_description = {
        "description": "Vibrational frequencies (Hessian) with ORCA",
        "group": "Calculation",
        "name": "Frequencies",
    }

    def __init__(self, flowchart=None, gui=None):
        pass

    def description(self):
        return FrequenciesStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        return orca_step.Frequencies(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        return orca_step.TkFrequencies(canvas=canvas, **kwargs)
