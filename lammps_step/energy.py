# -*- coding: utf-8 -*-

"""A single-point energy in LAMMPS"""

import logging
from pathlib import Path

from tabulate import tabulate

import lammps_step
import seamm
from seamm_util import units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
printer = printing.getPrinter("lammps")


class Energy(seamm.Node):
    """Handle a singlepoint energy calculation in LAMMPS"""

    def __init__(
        self,
        flowchart=None,
        title="Energy",
        extension=None,
        logger=logger,
    ):
        """Initialize the node"""

        logger.debug("Creating Energy {}".format(self))

        super().__init__(
            flowchart=flowchart,
            title=title,
            extension=extension,
            logger=logger,
        )

        self._calculation = "energy"
        self._model = None
        self._metadata = lammps_step.metadata
        self.parameters = lammps_step.EnergyParameters()

        self.description = "A single point energy calculation"

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

        # if not P:
        #     P = self.parameters.values_to_dict()

        text = "Single-point energy calculation."

        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def get_input(self, extras=None):
        """Get the input for an energy calculation for LAMMPS"""

        P = self.parameters.values_to_dict()

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        self.description = []
        self.description.append(__(self.description_text(PP), **PP, indent=3 * " "))

        lines = []

        filename = f"@{self._id[-1]}+forces.dump"
        lines.append("")
        lines.append(f"# {self.header}")
        lines.append("")
        lines.append(f"dump                1 all custom 1 {filename} id fx fy fz")
        lines.append("run                 0")
        lines.append("undump              1")

        return {
            "script": lines,
            "postscript": None,
            "use python": False,
        }

    def analyze(self, indent="", data={}, table=None, output=[], **kwargs):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """
        if table is not None:
            text = ""
            tmp = tabulate(
                table,
                headers="keys",
                tablefmt="simple",
                disable_numparse=True,
                colalign=(
                    "center",
                    "decimal",
                    "center",
                    "decimal",
                    "left",
                    "decimal",
                    "decimal",
                    "decimal",
                ),
            )
            length = len(tmp.splitlines()[0])
            text += "\n"
            text += "Properties".center(length)
            text += "\n"
            text += tmp
            text += "\n"

            printer.normal(__(text, indent=8 * " ", wrap=False, dedent=False))

            have_warning = False
            for value in table["Property"]:
                if len(value) > 0 and value[-1] == "*":
                    have_warning = True
                    break

            have_acf_warning = False
            for value in table["tau"]:
                if len(value) > 0 and value[-1] == "^":
                    have_warning = True
                    break

            if have_warning or have_acf_warning:
                printer.normal("\n")
            if have_warning:
                printer.normal(
                    __(
                        "          * this property has less than 100 independent "
                        "samples, so may not be accurate.",
                        wrap=False,
                        dedent=False,
                    )
                )

            if have_acf_warning:
                printer.normal(
                    __(
                        "          ^ there are not enough samples after "
                        "equilibration to plot the ACF.",
                        wrap=False,
                        dedent=False,
                    )
                )

        # Process the output to get the energy
        lines = iter(output)
        for line in lines:
            if "   Step" in line:
                tmp = line.split()
                try:
                    i = tmp.index("TotEng")
                except ValueError:
                    pass
                else:
                    line = next(lines)
                    tmp = line.split()
                    energy = float(tmp[i])
                    data["energy"] = energy
                    data["energy,units"] = "kcal/mol"

        # See if forces have been dumped
        wdir = Path(self.directory)
        path = wdir / "forces.dump"
        if path.exists():
            tmp = self.parent.get_dump(path)

            fields = tmp["fields"]
            if "fx" in fields and "fy" in fields and "fz" in fields:
                ix = fields.index("fx")
                iy = fields.index("fy")
                iz = fields.index("fz")
                fxs = tmp["data"][-1][ix]
                fys = tmp["data"][-1][iy]
                fzs = tmp["data"][-1][iz]
                gradients = [[-fx, -fy, -fz] for fx, fy, fz in zip(fxs, fys, fzs)]
                data["gradients"] = gradients
                data["gradients,units"] = "kcal/mol/angstrom"

        # Get the configuration
        _, configuration = self.get_system_configuration(None)

        # Need to set the model for properties.
        # See what type of forcefield we have and handle it
        ff = self.get_variable("_forcefield")
        if ff == "OpenKIM":
            self._model = "OpenKIM/" + self.get_variable("_OpenKIM_Potential")
        else:
            # Valence forcefield...
            self._model = ff.current_forcefield

        # Put any requested results into variables or tables
        self.store_results(
            configuration=configuration,
            data=data,
            create_tables=self.parameters["create tables"].get(),
        )
