# -*- coding: utf-8 -*-
"""A single-point energy in LAMMPS"""

import molssi_workflow
import logging

logger = logging.getLogger(__name__)


class Energy(molssi_workflow.Node):
    """Handle a singlepoint energy calculation in LAMMPS"""

    def __init__(self, workflow=None, title='Energy',
                 extension=None):
        """Initialize the node"""

        logger.debug('Creating Energy {}'.format(self))

        super().__init__(workflow=workflow, title=title,
                         extension=extension)

        self.description = 'A single point energy calculation'

    def get_input(self):
        """Get the input for an energy calculation for LAMMPS"""

        lines = []

        lines.append('')
        lines.append('#     single-point energy')
        lines.append('')
        lines.append('run                 0')

        return lines
