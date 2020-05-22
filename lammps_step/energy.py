# -*- coding: utf-8 -*-

"""A single-point energy in LAMMPS"""

import logging

import lammps_step
import seamm
from seamm_util import units_class
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)


class Energy(seamm.Node):
    """Handle a singlepoint energy calculation in LAMMPS"""

    def __init__(self, flowchart=None, title='Energy', extension=None):
        """Initialize the node"""

        logger.debug('Creating Energy {}'.format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.description = 'A single point energy calculation'
        self.parameters = lammps_step.EnergyParameters()

    @property
    def header(self):
        """A printable header for this section of output"""
        return (
            'Step {}: {}'.format(
                '.'.join(str(e) for e in self._id), self.title
            )
        )

    @property
    def version(self):
        """The semantic version of this module.
        """
        return lammps_step.__version__

    @property
    def git_revision(self):
        """The git version of this module.
        """
        return lammps_step.__git_revision__

    def description_text(self, P=None):
        """Create the text description of what this step will do.
        """

        if not P:
            P = self.parameters.values_to_dict()

        text = "Single-point energy calculation."

        return self.header + '\n' + __(text, indent=4 * ' ').__str__()

    def get_input(self, extras=None):
        """Get the input for an energy calculation for LAMMPS"""

        P = self.parameters.values_to_dict()

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = '{:~P}'.format(PP[key])

        self.description = []
        self.description.append(
            __(self.description_text(PP), **PP, indent=3 * ' ')
        )

        lines = []

        lines.append('')
        lines.append('#     single-point energy')
        lines.append('')
        lines.append('run                 0')

        return lines
