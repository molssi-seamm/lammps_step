# -*- coding: utf-8 -*-

"""A custom script for LAMMPS"""

import seamm
import logging

logger = logging.getLogger(__name__)


class Custom(seamm.Node):

    def __init__(self, flowchart=None, title='Custom', extension=None):
        """Initialize the node"""

        logger.debug('Creating Custom step {}'.format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.description = 'A custom script for LAMMPS'

        self.text = '# Custom script for LAMMPS (replace this!)\n\n'

    def get_input(self, extras=None):
        """Get the custom input for LAMMPS"""

        lines = []

        lines.append('')
        lines.append('#     custom scripting')
        lines.append('')
        lines.append(self.text)

        return lines
