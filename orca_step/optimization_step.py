# -*- coding: utf-8 -*-

"""Stevedore helper for the ORCA Optimization sub-step."""

import orca_step


class OptimizationStep(object):
    my_description = {
        "description": "Geometry optimization with ORCA",
        "group": "Calculation",
        "name": "Optimization",
    }

    def __init__(self, flowchart=None, gui=None):
        pass

    def description(self):
        return OptimizationStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        return orca_step.Optimization(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        return orca_step.TkOptimization(canvas=canvas, **kwargs)
