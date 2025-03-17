# -*- coding: utf-8 -*-

"""A custom script for LAMMPS"""

import seamm
import logging

import lammps_step
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)


class Custom(seamm.Node):
    def __init__(self, flowchart=None, title="Custom", extension=None):
        """Initialize the node"""

        logger.debug("Creating Custom step {}".format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self._calculation = "energy"
        self._model = None
        self._metadata = lammps_step.metadata
        self.parameters = lammps_step.CustomParameters()
        self.description = "A custom script for LAMMPS"

    @property
    def header(self):
        """A printable header for this section of output"""
        return "Step {}: {}".format(".".join(str(e) for e in self._id), self.title)

    @property
    def version(self):
        """The semantic version of this module."""
        return lammps_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return lammps_step.__git_revision__

    def description_text(self, P=None):
        """Create the text description of what this step will do."""

        if not P:
            P = self.parameters.values_to_dict()

        text = "Custom script for LAMMPS\n\n"
        lines = P["script"].splitlines()
        if len(lines) > 5:
            text += "\n".join(lines[0:4]) + "\n...\n"
        else:
            text += P["script"]

        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def get_input(self, extras=None):
        """Get the custom input for LAMMPS"""

        P = self.parameters.values_to_dict()

        self.description = []
        self.description.append(__(self.description_text(P), **P, indent=4 * " "))

        lines = []

        lines.append("")
        lines.append("#     custom scripting")
        lines.append("")
        lines.append(P["script"])

        return {
            "script": lines,
            "postscript": None,
            "use python": False,
        }
