# -*- coding: utf-8 -*-

"""Stevedore helper for the ORCA Energy sub-step."""

import orca_step


class EnergyStep(object):
    my_description = {
        "description": "Single-point energy with ORCA",
        "group": "Calculation",
        "name": "Energy",
    }

    def __init__(self, flowchart=None, gui=None):
        pass

    def description(self):
        return EnergyStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        return orca_step.Energy(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        return orca_step.TkEnergy(canvas=canvas, **kwargs)
