# -*- coding: utf-8 -*-

"""Main module."""

import lammps_step


class MinimizationStep(object):
    my_description = {
        'description':
        'Minimization',
        'group': 'Calculations',
        'name': 'Minimization'
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
        return MinimizationStep.my_description

    def factory(self, graphical=False, workflow=None, canvas=None, **kwargs):
        """Return the node object or graphical node object"""
        if graphical:
            return lammps_step.TkMinimization(canvas=canvas, **kwargs)
        else:
            return lammps_step.Minimization(workflow=workflow, **kwargs)
