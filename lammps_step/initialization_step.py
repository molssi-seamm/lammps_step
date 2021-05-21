# -*- coding: utf-8 -*-

"""Main module."""

import lammps_step


class InitializationStep(object):
    my_description = {
        "description": (
            "The initialization of LAMMPS, specifying e.g. nonbond "
            "methods and cutoffs"
        ),
        "group": "Calculations",
        "name": "Initialization",
    }

    def __init__(self, flowchart=None, gui=None):
        """Initialize this helper class, which is used by
        the application via stevedore to get information about
        and create node objects for the flowchart
        """
        pass

    def description(self):
        """Return a description of what this extension does"""
        return InitializationStep.my_description

    def create_node(self, flowchart=None, **kwargs):
        """Return the new node object"""
        return lammps_step.Initialization(flowchart=flowchart, **kwargs)

    def create_tk_node(self, canvas=None, **kwargs):
        """Return the graphical Tk node object"""
        return lammps_step.TkInitialization(canvas=canvas, **kwargs)
