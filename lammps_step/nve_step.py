# -*- coding: utf-8 -*-

"""Main module."""

import lammps_step


class NVEStep(object):
    my_description = {
        'description':
        'NVE (microcanonical) dynamics',
        'group': 'Calculations',
        'name': 'NVE dynamics'
    }

    def __init__(self, workflow=None, gui=None):
        """Initialize this helper class, which is used by
        the application via stevedore to get information about
        and create node objects for the workflow
        """
        pass

    def description(self):
        """Return a description of what this extension does
        """
        return NVEStep.my_description

    def factory(self, graphical=False, workflow=None, canvas=None, **kwargs):
        """Return the node object or graphical node object"""
        if graphical:
            return lammps_step.TkNVE(canvas=canvas, **kwargs)
        else:
            return lammps_step.NVE(workflow=workflow, **kwargs)
