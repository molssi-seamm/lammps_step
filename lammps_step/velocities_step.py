# -*- coding: utf-8 -*-

"""Main module."""

import lammps_step


class VelocitiesStep(object):
    my_description = {
        'description':
        'Set the velocities in preparation for MD',
        'group': 'Calculations',
        'name': 'Velocities'
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
        return VelocitiesStep.my_description

    def factory(self, graphical=False, workflow=None, canvas=None, **kwargs):
        """Return the node object or graphical node object"""
        if graphical:
            return lammps_step.TkVelocities(canvas=canvas, **kwargs)
        else:
            return lammps_step.Velocities(workflow=workflow, **kwargs)
