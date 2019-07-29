# -*- coding: utf-8 -*-

"""Set the velocities on the atoms

ToDo:
    * This does not support groups other than "all" yet. Once it does, setting
        specific velocities would also be useful.
    * Gaussian or normal distributions...
"""

import lammps_step
import logging
import seamm
from seamm import data
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __
import random

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter('lammps')


class Velocities(seamm.Node):

    def __init__(self, flowchart=None, title='Velocities', extension=None):
        """Initialize the node"""

        logger.debug('Creating Velocities {}'.format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.description = 'Set the initial velocities on the atoms'
        self.parameters = lammps_step.VelocitiesParameters()

    def description_text(self, P):
        """Prepare information about what this node will do
        """
        text = ('Set the velocities to give a temperature {T} ' 'by {method}.')

        if P['remove_momentum'][0] == '$':
            text += (
                ' Whether to remove translational or rotational '
                'momentum will be determined at runtime by '
                "'{remove_momentum}'"
            )
        else:
            text += ' LAMMPS will {remove_momentum}.'

        if P['method'] != 'scaling current velocities':
            if P['seed'] == 'random':
                text += (
                    ' The random number generator will be initialized '
                    'randomly.'
                )
            else:
                text += (
                    ' The random number generator will be initialized '
                    "with the seed '{seed}'."
                )

        return text

    def get_input(self):
        """Get the input for setting the velocities in LAMMPS"""
        self._long_header = ''
        self._long_header += str(__(self.header, indent=3 * ' '))
        self._long_header += '\n'

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        # Fix variables that need attention
        if 'default' in P['remove_momentum']:
            if data.structure['periodicity'] == 3:
                P['remove_momentum'] = (
                    "remove translational but not rotational momentum"
                )
            else:
                P['remove_momentum'] = (
                    "remove both translational and rotational momentum"
                )
        if P['seed'] == 'random':
            P['seed'] = int(random.random() * 2**31)

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = '{:~P}'.format(PP[key])

        self._long_header += str(
            __(self.description_text(PP), **PP, indent=7 * ' ')
        )
        self.description = [self._long_header]

        # Get the input lines
        lines = []
        lines.append('')
        lines.append('#     velocities')
        lines.append('')

        if P['remove_momentum'] == (
            "remove translational but not "
            "rotational momentum"
        ):
            remove_translations = 'yes'
            remove_rotations = 'no'
        elif P['remove_momentum'] == (
            "remove rotational but not "
            "translational momentum"
        ):
            remove_translations = 'no'
            remove_rotations = 'yes'
        elif P['remove_momentum'] == (
            "remove both translational and "
            "rotational momentum"
        ):
            remove_translations = 'yes'
            remove_rotations = 'yes'
        elif P['remove_momentum'] == (
            "remove neither translational nor "
            "rotational momentum"
        ):
            remove_translations = 'no'
            remove_rotations = 'no'
        else:
            raise RuntimeError(
                "Don't recognize 'remove_momentum' of '{}'".format(
                    P['remove_momentum']
                )
            )

        T = P['T'].to('K').magnitude

        if 'random' in P['method']:
            lines.append(
                'velocity            all create {} {} mom {} rot {}'.format(
                    T, P['seed'], remove_translations, remove_rotations
                )
            )
        elif 'scaling' in P['method']:
            lines.append(
                'velocity            all scale {} mom {} rot {}'.format(
                    T, remove_translations, remove_rotations
                )
            )
        else:
            raise RuntimeError(
                "Velocity method '{}' not supported yet".format(P['method'])
            )

        return lines
