# -*- coding: utf-8 -*-

"""Minimization step in LAMMPS"""

import logging
from pathlib import Path
import textwrap
import traceback

import numpy as np
from tabulate import tabulate

from molsystem import RMSD
import seamm
from seamm_util import units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __
import lammps_step

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("lammps")

minimization_style = {
    "Conjugate Gradient": "cg",
    "Steepest Descent": "sd",
    "Hessian-Free truncated Newton": "hftn",
    "QuickMin": "quickmin",
    "Fire": "fire",
}


class Minimization(lammps_step.Energy):
    def __init__(self, flowchart=None, title="Minimization", extension=None):
        """Initialize the node"""

        logger.debug("Creating Minimization {}".format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.description = "Minimization step in LAMMPS"

        self._calculation = "minimization"
        self._metadata = lammps_step.metadata
        self.parameters = lammps_step.MinimizationParameters()
        self._save = {}

    def analyze(self, data={}, properties=None, table=None, output=""):
        """Analyze the results of the simulation.

        Parameters
        ----------
        """
        ff = self.get_variable("_forcefield")

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Get the initial structure
        _, initial_configuration = self.get_system_configuration()

        # Handle the new structure as needed
        system, configuration = self.get_system_configuration(P, model=self.model)

        # Read the dump file and get the structure
        dump_file = Path(self.directory) / "minimization.dump"
        if dump_file.exists:
            try:
                xyz, fractional, cell, vxyz = self.parent.read_dump(dump_file)
            except Exception as e:
                xyz = None
                cell = None
                vxyz = None
                printer.normal("Warning: unable to read the LAMMPS dumpfile")
                logger.warning(f"The was a problem reading the LAMMPS dumpfile: {e}")
                logger.warning(traceback.format_exc())
        else:
            printer.normal("Warning: there is no 'dump' file from LAMMPS")
            xyz = None
            cell = None
            vxyz = None

        if configuration is not None:
            if cell is not None:
                configuration.cell.parameters = cell
            if xyz is not None:
                configuration.atoms.set_coordinates(xyz, fractionals=fractional)
            if vxyz is not None:
                # LAMMPS only has Cartesian velocities
                configuration.atoms.set_velocities(vxyz, fractionals=False)

            # And the name of the configuration.
            text = seamm.standard_parameters.set_names(
                system,
                configuration,
                P,
                _first=True,
                model=self.model,
            )
            printer.normal(__(text, **data, indent=self.indent + 4 * " "))
            printer.normal("")

        # Get the RMSD of the initial and final structures
        if initial_configuration.periodicity == 0 and xyz is not None:
            initial = initial_configuration.to_RDKMol()
            final = initial_configuration.to_RDKMol()
            final.GetConformer(0).SetPositions(np.array(xyz))

            result = RMSD(final, initial, symmetry=True, include_h=True, align=True)
            data["RMSD with H"] = result["RMSD"]
            data["displaced atom with H"] = result["displaced atom"]
            data["maximum displacement with H"] = result["maximum displacement"]

            # Align the structure
            if configuration is not None:
                configuration.coordinates = final.GetConformer(0).GetPositions()

            result = RMSD(final, initial, symmetry=True)
            data["RMSD"] = result["RMSD"]
            data["displaced atom"] = result["displaced atom"]
            data["maximum displacement"] = result["maximum displacement"]

        if table is None:
            table = {
                "Property": [],
                "Value": [],
                "Units": [],
            }

        first = None
        try:
            first = output.index("Minimization stats:")
        except ValueError:
            pass

        if first is not None:
            line = output[first + 1]
            if "Stopping criterion =" in line:
                criterion = " ".join(line.split()[3:])
                table["Property"].append("Converged")
                if criterion in ("energy tolerance", "force tolerance"):
                    table["Value"].append("True")
                    data["optimization is converged"] = True
                else:
                    table["Value"].append("False")
                    data["optimization is converged"] = False
                table["Units"].append("")
                table["Property"].append("Stopped because")
                table["Value"].append(criterion)
                table["Units"].append("")

            line = output[first + 2]
            if "Energy initial, next-to-last, final =" in line:
                Es = [float(e) for e in output[first + 3].split()]

                table["Property"].append("Final Energy")
                E = lammps_step.from_lammps_units(Es[2], "kcal/mol")
                table["Value"].append(f"{E.magnitude:.2f}")
                table["Units"].append("kcal/mol")
                data["energy"] = E.magnitude

                table["Property"].append("Initial Energy")
                E = lammps_step.from_lammps_units(Es[0], "kcal/mol")
                table["Value"].append(f"{E.magnitude:.2f}")
                table["Units"].append("kcal/mol")

                if ff.ff_form == "reaxff":
                    Eat = self.parent._atomic_energy_sum
                    if Eat != 0.0:
                        dHf = data["energy"] + Eat
                        data["DfH0_reax"] = dHf
                        table["Property"].append("DfH0_reax")
                        table["Value"].append(f"{dHf:.2f}")
                        table["Units"].append("kcal/mol")

                table["Property"].append("Last Energy Change")
                E = lammps_step.from_lammps_units(Es[2] - Es[1], "kcal/mol")
                table["Value"].append(f"{E.magnitude:.4f}")
                table["Units"].append("kcal/mol")
                data["energy change"] = E.magnitude

            line = output[first + 7]
            if "Iterations, force evaluations =" in line:
                F0, F1 = line.split()[4:6]

                table["Property"].append("Number of Steps")
                table["Value"].append(f"{F0}")
                table["Units"].append("")
                data["N steps optimization"] = F0

                table["Property"].append("Number of Force Evaluations")
                table["Value"].append(f"{F1}")
                table["Units"].append("")
                data["N force evaluations"] = F1

        if "RMSD" in data:
            tmp = data["RMSD"]
            table["Property"].append("RMSD in Geometry")
            table["Value"].append(f"{tmp:.2f}")
            table["Units"].append("Å")

        if "maximum displacement" in data:
            tmp = data["maximum displacement"]
            table["Property"].append("Largest Displacement")
            table["Value"].append(f"{tmp:.2f}")
            table["Units"].append("Å")

        if "displaced atom" in data:
            tmp = data["displaced atom"]
            table["Property"].append("Displaced Atom")
            table["Value"].append(f"{tmp + 1}")
            table["Units"].append("")

        if first is not None:
            line = output[first + 4]
            if "Force two-norm initial, final =" in line:
                F0, F1 = line.split()[5:7]

                table["Property"].append("Final Force Norm")
                tmp = lammps_step.from_lammps_units(float(F1), "kcal/mol/Å")
                table["Value"].append(f"{tmp.magnitude:.3f}")
                table["Units"].append("kcal/mol/Å")
                data["force norm"] = tmp.magnitude

                table["Property"].append("Initial Force Norm")
                tmp = lammps_step.from_lammps_units(float(F0), "kcal/mol/Å")
                table["Value"].append(f"{tmp.magnitude:.3f}")
                table["Units"].append("kcal/mol/Å")

            line = output[first + 5]
            if "Force max component initial, final =" in line:
                F0, F1 = line.split()[6:8]

                table["Property"].append("Final Maximum Force")
                tmp = lammps_step.from_lammps_units(float(F1), "kcal/mol/Å")
                table["Value"].append(f"{tmp.magnitude:.3f}")
                table["Units"].append("kcal/mol/Å")

                table["Property"].append("Initial Maximum Force")
                tmp = lammps_step.from_lammps_units(float(F0), "kcal/mol/Å")
                table["Value"].append(f"{tmp.magnitude:.3f}")
                table["Units"].append("kcal/mol/Å")

            line = output[first + 6]
            if "Final line search alpha, max atom move =" in line:
                F0, F1 = line.split()[8:10]

                table["Property"].append("Final Line Search Alpha")
                table["Value"].append(f"{float(F0):.3f}")
                table["Units"].append("")

                table["Property"].append("Final Maximum Displacement")
                tmp = lammps_step.from_lammps_units(float(F1), "Å")
                table["Value"].append(f"{tmp.magnitude:.3f}")
                table["Units"].append("Å")

        data["force norm threshold"] = self._save["ftol"]
        data["energy threshold"] = self._save["etol"]
        data["minimizer"] = self._save["minimizer"]

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            colalign=("center", "decimal", "left"),
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        text_lines = []
        text_lines.append("Convergence".center(length))
        text_lines.append(tmp)
        printer.normal(textwrap.indent("\n".join(text_lines), self.indent + 7 * " "))
        printer.normal("")

        if configuration is not None:
            self.store_results(configuration=configuration, data=data)

    def description_text(self, P=None):
        """Create the text description of what this step will do."""

        if P is None:
            P = self.parameters.values_to_dict()

        model = self.model

        text = "The structure will be minimized using the {minimizer} approach "
        text += "to a '{convergence}' convergence"
        if model is None:
            text += "."
        else:
            text += f" using the {model} forcefield. "
        text += "The number of steps will be limited to {nsteps} steps and no more "
        text += "than {nevaluations} energy and force evaluations."

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def get_input(self, extras=None):
        """Get the input for a minimization in LAMMPS"""

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        self.description = []
        self.description.append(__(self.description_text(PP), **PP, indent=4 * " "))

        _, configuration = self.get_system_configuration()
        nAtoms = configuration.n_atoms  # noqa: F841

        lines = []

        etol = 0.0
        ftol = None
        convergence = P["convergence"]
        if convergence == "normal":
            etol = 0.0
            ftol = 0.1
        elif convergence == "tight":
            etol = 0.0
            ftol = 0.01
        elif convergence == "loose":
            etol = 0.0
            ftol = 1
        elif convergence == "crude":
            etol = 0.0
            ftol = 10
        elif convergence == "custom":
            etol = P["etol"].magnitude
            ftol = P["ftol"].magnitude
        else:
            raise ValueError(f"Don't understand convergence '{convergence}'!")

        self._save["etol"] = etol
        self._save["ftol"] = ftol

        # maximum number of iterations
        nSteps = eval(P["nsteps"])
        nEvals = eval(P["nevaluations"])

        # Minimization style
        minimizer = P["minimizer"]
        if minimizer not in minimization_style:
            raise ValueError(f"Don't recognize minimizer '{minimizer}'")
        min_style = minimization_style[minimizer]

        self._save["minimizer"] = minimizer

        timestep = lammps_step.to_lammps_units(P["timestep"], quantity="time")

        thermo_properties = (
            "fmax fnorm press etotal ke pe ebond eangle edihed eimp evdwl etail ecoul "
            "elong"
        )
        if self.parent.have_dreiding_hbonds:
            thermo_properties += " v_N_hbond v_E_hbond"

        lines.append("")
        lines.append(f"# {self.header}")
        lines.append("")
        lines.append(f"thermo_style        custom {thermo_properties}")
        # lines.append("thermo              100")
        lines.append("thermo              1")
        lines.append(f"min_style           {min_style}")
        if min_style in ("quickmin", "fire"):
            lines.append(f"timestep            {timestep}")
        lines.append(f"minimize            {etol} {ftol} {nSteps} {nEvals}")
        lines.append("")
        filename = f"@{self._id[-1]}+minimization.dump"
        lines.append(
            f"write_dump          all custom  {filename} id xu yu zu fx fy fz"
            " modify flush yes sort id"
        )

        return {
            "script": lines,
            "postscript": None,
            "use python": False,
        }
