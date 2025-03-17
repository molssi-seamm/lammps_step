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
from seamm_util import units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __
import random

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("lammps")


class Velocities(seamm.Node):
    def __init__(self, flowchart=None, title="Velocities", extension=None):
        """Initialize the node"""

        logger.debug("Creating Velocities {}".format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.description = "Set the initial velocities on the atoms"
        self.parameters = lammps_step.VelocitiesParameters()

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
        """Return a short description of this step.

        Return a nicely formatted string describing what this step will
        do.

        Keyword arguments:
            P: a dictionary of parameter values, which may be variables
                or final values. If None, then the parameters values will
                be used as is.
        """

        if not P:
            P = self.parameters.values_to_dict()

        text = "Set the velocities to give a temperature {T} " "by {method}."

        if P["remove_momentum"][0] == "$":
            text += (
                " Whether to remove translational or rotational "
                "momentum will be determined at runtime by "
                "'{remove_momentum}'"
            )
        else:
            text += " LAMMPS will {remove_momentum}"

        if P["method"] != "scaling current velocities":
            if P["seed"] == "random":
                text += " The random number generator will be initialized randomly."
            else:
                text += (
                    " The random number generator will be initialized "
                    "with the seed '{seed}'."
                )

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def get_input(self, extras=None):
        """Get the input for setting the velocities in LAMMPS"""
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        # Fix variables that need attention
        if "default" in P["remove_momentum"]:
            system_db = self.get_variable("_system_db")
            configuration = system_db.system.configuration
            if configuration.periodicity == 3:
                P["remove_momentum"] = (
                    "remove translational but not rotational momentum"
                )
            else:
                P["remove_momentum"] = (
                    "remove both translational and rotational momentum"
                )
        if P["seed"] == "random":
            P["seed"] = int(random.random() * 2**31)

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        self.description = [str(__(self.description_text(PP), **PP, indent=4 * " "))]

        # Get the input lines
        lines = []
        lines.append("")
        lines.append(f"# {self.header}")
        lines.append("")

        if P["remove_momentum"] == "remove translational but not rotational momentum":
            remove_translations = "yes"
            remove_rotations = "no"
        elif P["remove_momentum"] == "remove rotational but not translational momentum":
            remove_translations = "no"
            remove_rotations = "yes"
        elif P["remove_momentum"] == (
            "remove both translational and rotational momentum"
        ):
            remove_translations = "yes"
            remove_rotations = "yes"
        elif P["remove_momentum"] == (
            "remove neither translational nor rotational momentum"
        ):
            remove_translations = "no"
            remove_rotations = "no"
        else:
            raise RuntimeError(
                "Don't recognize 'remove_momentum' of '{}'".format(P["remove_momentum"])
            )

        T = P["T"].to("K").magnitude

        if "random" in P["method"]:
            lines.append(
                "velocity            all create {} {} mom {} rot {}".format(
                    T, P["seed"], remove_translations, remove_rotations
                )
            )
        elif "scaling" in P["method"]:
            lines.append(
                "velocity            all scale {} mom {} rot {}".format(
                    T, remove_translations, remove_rotations
                )
            )
        else:
            raise RuntimeError(
                "Velocity method '{}' not supported yet".format(P["method"])
            )

        return {
            "script": lines,
            "postscript": None,
            "use python": False,
        }
