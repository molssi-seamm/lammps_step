# -*- coding: utf-8 -*-

"""Main module."""

import lammps_step


class NVTStep(object):
    my_description = {
        'description':
        'NVTE (canonical) dynamics',
        'group': 'Calculations',
        'name': 'NVT dynamics'
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
        return NVTStep.my_description

    def factory(self, graphical=False, workflow=None, canvas=None, **kwargs):
        """Return the node object or graphical node object"""
        if graphical:
            return lammps_step.TkNVT(canvas=canvas, **kwargs)
        else:
            return lammps_step.NVT(workflow=workflow, **kwargs)
