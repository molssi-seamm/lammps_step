# -*- coding: utf-8 -*-

"""Main module."""

import lammps_step


class CustomStep(object):
    my_description = {
        "description": (
            "A custom step for LAMMPS, where the user types in LAMMPS "
            "commmands directly"
        ),
        "group": "Customize",
        "name": "Custom",
    }

    def __init__(self, flowchart=None, gui=None):
        """Initialize this helper class, which is used by
        the application via stevedore to get information about
        and create node objects for the flowchart
        """
        pass

    def description(self):
        """Return a description of what this extension does"""
        return CustomStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        """Return the new node object"""
        return lammps_step.Custom(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        """Return the graphical Tk node object"""
        return lammps_step.TkCustom(canvas=canvas, **kwargs)
