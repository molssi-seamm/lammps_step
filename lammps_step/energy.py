# -*- coding: utf-8 -*-

"""A single-point energy in LAMMPS"""

import lammps_step
import seamm
import logging

logger = logging.getLogger(__name__)


class Energy(seamm.Node):
    """Handle a singlepoint energy calculation in LAMMPS"""

    def __init__(self, flowchart=None, title='Energy', extension=None):
        """Initialize the node"""

        logger.debug('Creating Energy {}'.format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.description = 'A single point energy calculation'
        self.parameters = lammps_step.EnergyParameters()

    def description_text(self):
        """Create the text description of what this step will do.
        """

        text = ("Single-point energy calculation.")

        return text

    def get_input(self):
        """Get the input for an energy calculation for LAMMPS"""

        lines = []

        lines.append('')
        lines.append('#     single-point energy')
        lines.append('')
        lines.append('run                 0')

        return lines
