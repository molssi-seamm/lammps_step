# -*- coding: utf-8 -*-
"""A custom script for LAMMPS"""

import molssi_workflow
import logging

logger = logging.getLogger(__name__)


class Custom(molssi_workflow.Node):
    def __init__(self, workflow=None, title='Custom',
                 extension=None):
        """Initialize the node"""

        logger.debug('Creating Custom step {}'.format(self))

        super().__init__(workflow=workflow, title=title,
                         extension=extension)

        self.description = 'A custom script for LAMMPS'

        self.text = '# Custom script for LAMMPS (replace this!)\n\n'

    def get_input(self):
        """Get the custom input for LAMMPS"""

        lines = []

        lines.append('')
        lines.append('#     custom scripting')
        lines.append('')
        lines.append(self.text)

        return lines
