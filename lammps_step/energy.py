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
    def model(self):
        """The name of the forcefield."""
        return self.parent.model

    @model.setter
    def model(self, value):
        self.parent.model = value

    @property
    def header(self):
        """A printable header for this section of output"""
        return "Step {}: {}".format(".".join(str(e) for e in self._id), self.title)

    @property
    def results(self):
        """The storage for the results in the main LAMMPS step."""
        return self.flowchart.parent._results

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
        self.description.append(__(self.description_text(PP), **PP, indent=4 * " "))

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
        ff = self.get_variable("_forcefield")

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

                # Check for reaxff enthalpy offset
                if ff.ff_form == "reaxff":
                    Eat = self.parent._atomic_energy_sum
                    if Eat != 0.0:
                        dHf = data["energy"] + Eat
                        data["DfH0_reax"] = dHf
                    data["energy,units"] = "kcal/mol"
            if line.startswith("Loop time of"):
                try:
                    tmp = line.split()
                    _time = round(float(tmp[3]), 2)
                    _procs = int(tmp[5])
                    _steps = int(tmp[8])
                    _natoms = int(tmp[11])
                    _type = self._calculation
                    self.results[f"t_{_type}"] = _time
                    self.results[f"np_{_type}"] = _procs
                    self.results[f"steps_{_type}"] = _steps
                    self.results[f"natoms_{_type}"] = _natoms
                except Exception as _e:
                    print(f"LAMMPS loop time: {_e}")

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

        # Add static properties such as density for e.g NVE and NVT calculations
        if "density" not in data and configuration.periodicity == 3:
            data["density"] = configuration.density
        if "a" not in data and configuration.periodicity == 3:
            data["a"] = configuration.cell.a
        if "b" not in data and configuration.periodicity == 3:
            data["b"] = configuration.cell.b
        if "c" not in data and configuration.periodicity == 3:
            data["c"] = configuration.cell.c
        if "volume" not in data and configuration.periodicity == 3:
            data["volume"] = configuration.cell.volume

        # Put any requested results into variables or tables
        self.store_results(
            configuration=configuration,
            data=data | self.results,
            create_tables=self.parameters["create tables"].get(),
        )
