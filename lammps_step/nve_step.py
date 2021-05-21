# -*- coding: utf-8 -*-

"""Main module."""

import lammps_step


class NVEStep(object):
    my_description = {
        "description": "NVE (microcanonical) dynamics",
        "group": "Calculations",
        "name": "NVE dynamics",
    }

    def __init__(self, flowchart=None, gui=None):
        """Initialize this helper class, which is used by
        the application via stevedore to get information about
        and create node objects for the flowchart
        """
        pass

    def description(self):
        """Return a description of what this extension does"""
        return NVEStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        """Return the new node object"""
        return lammps_step.NVE(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        """Return the graphical Tk node object"""
        return lammps_step.TkNVE(canvas=canvas, **kwargs)
