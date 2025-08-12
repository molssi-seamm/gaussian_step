# -*- coding: utf-8 -*-

"""Main module."""

import gaussian_step


class WavefunctionStabilityStep(object):
    my_description = {
        "description": "Wavefunction stability analysis using Gaussian",
        "group": "Calculation",
        "name": "Wavefunction Stability",
    }

    def __init__(self, flowchart=None, gui=None):
        """Initialize this helper class, which is used by
        the application via stevedore to get information about
        and create node objects for the flowchart
        """
        pass

    def description(self):
        """Return a description of what this extension does"""
        return WavefunctionStabilityStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        """Return the new node object"""
        return gaussian_step.WavefunctionStability(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        """Return the graphical Tk node object"""
        return gaussian_step.TkWavefunctionStability(canvas=canvas, **kwargs)
