# -*- coding: utf-8 -*-
"""Set the velocities on the atoms

ToDo:
    * This does not support groups other than "all" yet. Once it does, setting
        specific velocities would also be useful.
    * Gaussian or normal distributions...
"""

import molssi_workflow
from molssi_workflow import units, Q_, data  # nopep8
import logging
import random

logger = logging.getLogger(__name__)


class Velocities(molssi_workflow.Node):
    structures = {
        'current': '',
        'initial': '',
        'other': '',
    }

    def __init__(self, workflow=None, gui_object=None, title='Velocities',
                 extension=None):
        """Initialize the node"""

        logger.debug('Creating Velocities {}'.format(self))

        super().__init__(workflow=workflow, title=title, gui_object=gui_object,
                         extension=extension)

        self.description = 'Set the initial velocities on the atoms'

        self.method = 'using a random distribution'
        self.seed_method = 'random'
        self.seed = 53
        self.seed_variable = ''
        self.temperature_method = 'as given'
        self.temperature = Q_(25, units.degC)
        self.temperature_variable = ''
        self.momentum_method = 'default'
        self.remove_linear_momentum = True
        self.remove_angular_momentum = False

    def get_input(self):
        """Get the input for setting the velocities in LAMMPS"""

        lines = []

        lines.append('')
        lines.append('#     velocities')
        lines.append('')
        if self.momentum_method == 'default':
            if data.structure['periodicity'] == 3:
                remove_translations = 'yes'
                remove_rotations = 'no'
            else:
                remove_translations = 'yes'
                remove_rotations = 'yes'
        elif(self.momentum_method) == 'as given':
            remove_translations = 'yes' if self.remove_linear_momentum \
                                  else 'no'
            remove_rotations = 'yes' if self.remove_angular_momentum else 'no'
        else:
            raise RuntimeError(
                'momentum method does not support variables yet')

        if 'random' in self.method:
            if self.temperature_method == 'as given':
                T = self.temperature.to('K').magnitude
                if self.seed_method == 'random':
                    seed = int(random.random() * 2**31)
                elif self.seed_method == 'as given':
                    seed = self.seed
                else:
                    raise RuntimeError('seed does not support variables yet')
            else:
                raise RuntimeError(
                    'temperature does not support variables yet')
            lines.append('velocity            '
                         'all create {} {} mom {} rot {}'.format(
                             T, seed, remove_translations, remove_rotations))
        elif 'scaling' in self.method:
            if self.temperature_method == 'as given':
                T = self.temperature.to('K').magnitude
            else:
                raise RuntimeError(
                    'temperature does not support variables yet')
            lines.append('velocity            '
                         'all scale {} mom {} rot {}'.format(
                             T, remove_translations, remove_rotations))
        else:
            raise RuntimeError(
                "Velocity method '{}' not supported yet".format(self.method))

        return lines
