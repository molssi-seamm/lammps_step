# -*- coding: utf-8 -*-
"""Set the velocities on the atoms

ToDo:
    * This does not support groups other than "all" yet. Once it does, setting
        specific velocities would also be useful.
    * Gaussian or normal distributions...
"""

import logging
import molssi_workflow
from molssi_workflow import ureg, Q_, data, units_class  # nopep8
import molssi_util.printing as printing
from molssi_util.printing import FormattedText as __
import random
from textwrap import dedent

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter('lammps')


class Velocities(molssi_workflow.Node):
    structures = {
        'current': '',
        'initial': '',
        'other': '',
    }

    def __init__(self, workflow=None, title='Velocities',
                 extension=None):
        """Initialize the node"""

        logger.debug('Creating Velocities {}'.format(self))

        super().__init__(workflow=workflow, title=title,
                         extension=extension)

        self.description = 'Set the initial velocities on the atoms'

        self.method = 'using a random distribution'
        self.seed_method = 'random'
        self.seed = 53
        self.seed_variable = ''
        self.temperature_method = 'is'
        self.temperature = Q_(298.15, ureg.K)
        self.temperature_variable = ''
        self.momentum_method = 'default'
        self.remove_linear_momentum = True
        self.remove_angular_momentum = False

    def describe(self, indent='', json_dict=None):
        """Write out information about what this node will do
        If json_dict is passed in, add information to that dictionary
        so that it can be written out by the controller as appropriate.
        """

        next_node = super().describe(indent, json_dict)

        values = {}
        if isinstance(self.temperature, units_class):
            values['temperature'] = '{:~P}'.format(self.temperature)
        else:
            values['temperature'] = self.temperature
        values['method'] = self.method

        string = dedent("""\
        Set the temperature of the system to {temperature} by setting the
        velocities {method}.""")
        if self.momentum_method == 'default':
            string += dedent("""
            By default, the linear momentum will be projected out for
            periodic systems, and both linear and angular momentum for
            molecular (non-periodic) systems.""")
        else:
            if self.remove_linear_momentum:
                if self.remove_angular_momentum:
                    string += (" Any linear or rotational momentum "
                               "will be removed.")
                string += " Any linear momentum will be removed."
            elif self.remove_angular_momentum:
                string += " Any angular momentum will be removed."

        job.job(__(string, indent=self.indent+'    ', **values))

        return next_node

    def get_input(self):
        """Get the input for setting the velocities in LAMMPS"""

        self.description = []
        self.description.append(__(self.header, indent=self.indent))

        values = {}
        values['method'] = self.get_value(self.method)

        string = dedent("""\
        Set the temperature of the system to {T} by setting the
        velocities {method}.""")

        if 'random' in self.method:
            if self.seed_method == 'random':
                seed = int(random.random() * 2**31)
                string += dedent("""
                The random number seed was picked randomly and is {seed}.""")
            elif self.seed_method == 'is':
                seed = self.get_value(self.seed)
                string += " The random number seed was set to {seed}."
            values['seed'] = seed

        lines = []

        lines.append('')
        lines.append('#     velocities')
        lines.append('')
        if self.momentum_method == 'default':
            if data.structure['periodicity'] == 3:
                remove_translations = 'yes'
                remove_rotations = 'no'
                string += dedent("""
                By default, the linear momentum will be projected out for
                this periodic systems.""")
            else:
                remove_translations = 'yes'
                remove_rotations = 'yes'
                string += dedent("""
                By default, both linear and angular momentum will be
                projected out for this molecular (non-periodic) system.""")
        elif(self.momentum_method) == 'is':
            remove_translations = 'yes' if self.remove_linear_momentum \
                                  else 'no'
            remove_rotations = 'yes' if self.remove_angular_momentum else 'no'
            if remove_translations == 'yes':
                if self.remove_rotataions == 'yes':
                    string += (" Any linear or rotational momentum "
                               "will be removed.")
                string += " Any linear momentum will be removed."
            elif remove_rotations == 'yes':
                string += " Any angular momentum will be removed."
        else:
            raise RuntimeError(
                'momentum method does not support variables yet')

        T = self.get_value(self.temperature)
        if isinstance(self.temperature, units_class):
            values['T'] = '{:~P}'.format(T)
            T = T.to('K').magnitude
        else:
            values['T'] = '{} K'.format(T)

        if 'random' in self.method:
            lines.append('velocity            '
                         'all create {} {} mom {} rot {}'.format(
                             T, seed, remove_translations, remove_rotations))
        elif 'scaling' in self.method:
            lines.append('velocity            '
                         'all scale {} mom {} rot {}'.format(
                             T, remove_translations, remove_rotations))
        else:
            raise RuntimeError(
                "Velocity method '{}' not supported yet".format(self.method))

        self.description.append(
            __(string, indent=self.indent+'    ', **values)
        )

        return lines
