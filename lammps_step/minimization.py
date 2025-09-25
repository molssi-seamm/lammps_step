# -*- coding: utf-8 -*-

"""Minimization step in LAMMPS"""

import json
import logging
import textwrap
import traceback

import numpy as np
from tabulate import tabulate

from molsystem import RMSD
import seamm
from seamm_util import Q_, units_class
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
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Get the initial structure
        _, initial_configuration = self.get_system_configuration()

        # Save the initial cell parameter if periodic
        if initial_configuration.periodicity != 0:
            a0, b0, c0, alpha0, beta0, gamma0 = initial_configuration.cell.parameters

        # Handle the new structure as needed
        system, configuration = self.get_system_configuration(P, model=self.model)

        # See if there are other results in json files
        filename = self.wd / "minimization.json"
        if filename.exists():
            try:
                with filename.open("r") as fd:
                    tmp = json.load(fd)
            except Exception as e:
                printer.normal(f"Warning: error reading {filename}: {e}")
                logger.warning(f"Error with {filename}: {e}")
                pass
            else:
                data.update(tmp)

        # Read the dump file and get the structure
        dump_file = self.wd / "minimization.dump"
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

                if self.parent.ff_form() == "reaxff":
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

        # For periodic systems, the stress
        if initial_configuration.periodicity != 0:
            for key in ("Sxx", "Syy", "Szz", "Syz", "Sxz", "Sxy"):
                if key in data:
                    table["Property"].append(key)
                    value = Q_(data[key], data[key + ",units"]).m_as("GPa")
                    table["Value"].append(f"{value:.3f}")
                    table["Units"].append("GPa")

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

        # For periodic systems, the change in cell
        if (
            cell is not None
            and initial_configuration.periodicity != 0
            and P["optimize cell"]
        ):
            a, b, c, alpha, beta, gamma = cell
            ctable = {
                "": ("𝗮", "𝗯", "𝗰", "𝞪", "𝞫", "𝞬"),
                "Initial": (
                    f"{a0:.3f}",
                    f"{b0:.3f}",
                    f"{c0:.3f}",
                    f"{alpha0:.1f}",
                    f"{beta0:.1f}",
                    f"{gamma0:.1f}",
                ),
                "Final": (
                    f"{a:.3f}",
                    f"{b:.3f}",
                    f"{c:.3f}",
                    f"{alpha:.1f}",
                    f"{beta:.1f}",
                    f"{gamma:.1f}",
                ),
                "Change": (
                    f"{a - a0:.3f}",
                    f"{b - b0:.3f}",
                    f"{c - c0:.3f}",
                    f"{alpha - alpha0:.1f}",
                    f"{beta - beta0:.1f}",
                    f"{gamma - gamma0:.1f}",
                ),
                "Units": ("Å", "Å", "Å", "°", "°", "°"),
            }

            tmp = tabulate(
                ctable,
                headers="keys",
                tablefmt="rounded_outline",
                colalign=("center", "decimal", "decimal", "decimal", "center"),
                disable_numparse=True,
            )
            length = len(tmp.splitlines()[0])
            text_lines = []
            text_lines.append("Cell Parameters".center(length))
            text_lines.append(tmp)
            printer.normal(
                textwrap.indent("\n".join(text_lines), self.indent + 7 * " ")
            )
            printer.normal("")

        if configuration is not None:
            self.store_results(configuration=configuration, data=data, printer=printer)

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

        # Number of atoms may be used in the number of steps
        nAtoms = configuration.n_atoms  # noqa: F841

        # May need to force a triclinic cell
        if P["allow shear"]:
            self.parent.force_triclinic = True

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

        nfixes = 0

        lines.append("")
        lines.append(f"# {self.header}")
        lines.append("")
        lines.append(f"thermo_style        custom {thermo_properties}")
        lines.append("thermo              100")

        if configuration.periodicity != 0:
            # Attend to optimization of the cell
            if P["optimize cell"]:
                # Work out the pressure/stress part of the command
                ptext = self.get_pressure_text(P)
                if ptext != "":
                    nfixes += 1
                    lines.append(
                        f"fix                 {nfixes} all box/relax {ptext}"
                        " fixedpoint 0.0 0.0 0.0"
                    )
                    lines.append("")

        lines.append(f"min_style           {min_style}")
        if min_style in ("quickmin", "fire"):
            lines.append(f"timestep            {timestep}")
        lines.append(f"minimize            {etol} {ftol} {nSteps} {nEvals}")
        lines.append("")
        if nfixes > 0:
            for fix in range(1, nfixes + 1):
                lines.append(f"unfix               {fix}")
            lines.append("")

        if configuration.periodicity != 0:
            # Write out the stress after minimization
            filename = f"@{self._id[-1]}+minimization.json"
            units = lammps_step.lammps_units("pressure")
            lines.append(
                'print               """{\n'
                '    "Sxx": $(v_sxx:%.3f),\n'
                f'    "Sxx,units": "{units}",\n'
                '    "Syy": $(v_syy:%.3f),\n'
                f'    "Syy,units": "{units}",\n'
                '    "Szz": $(v_szz:%.3f),\n'
                f'    "Szz,units": "{units}",\n'
                '    "Syz": $(v_syz:%.3f),\n'
                f'    "Syz,units": "{units}",\n'
                '    "Sxz": $(v_sxz:%.3f),\n'
                f'    "Sxz,units": "{units}",\n'
                '    "Sxy": $(v_sxy:%.3f),\n'
                f'    "Sxy,units": "{units}"\n'
                "}"
                f'""" file {filename}'
            )
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

    def get_pressure_text(self, _P):
        """Work out and return the pressure/stress part of the
        'fix npt' or 'fix berendsen' or minimization in LAMMPS
        """
        system_type = _P["system type"]

        if system_type == "fluid":
            P = lammps_step.to_lammps_units(_P["P"], quantity="pressure")
            return f" iso {P}"

        use_stress = "pressure" not in _P["use_stress"]
        couple = _P["couple"]
        allow_shear = _P["allow shear"]

        Sxx = _P["Sxx"]
        Syy = _P["Syy"]
        Szz = _P["Szz"]
        Syz = _P["Syz"]
        Sxz = _P["Sxz"]
        Sxy = _P["Sxz"]

        if use_stress:
            if couple == "x, y and z":
                ptext = "couple xyz x {Sxx} y {Sxx} z {Sxx}"
            elif couple == "x and y":
                if Sxx == "fixed":
                    ptext = "couple xy z {Szz}"
                elif Szz == "fixed":
                    ptext = "couple xy x {Sxx} y {Sxx}"
                else:
                    ptext = "couple xy x {Sxx} y {Sxx} z {Szz}"
            elif couple == "x and z":
                if Sxx == "fixed":
                    ptext = "couple xz y {Syy}"
                elif Syy == "fixed":
                    ptext = "couple xz x {Sxx} z {Sxx}"
                else:
                    ptext = "couple xz x {Sxx} y {Syy} z {Sxx}"
            elif couple == "y and z":
                if Sxx == "fixed":
                    ptext = "couple yz y {Syy} z {Syy}"
                elif Syy == "fixed":
                    ptext = "couple yz x {Sxx}"
                else:
                    ptext = "couple yz x {Sxx} y {Syy} z {Syy}"
            else:
                ptext = "couple none"
                if Sxx != "fixed":
                    ptext += " x {Sxx}"
                if Syy != "fixed":
                    ptext += " y {Syy}"
                if Szz != "fixed":
                    ptext += " z {Szz}"

            if allow_shear:
                if Syz != "fixed":
                    ptext += " yz {Syz}"
                if Sxz != "fixed":
                    ptext += " xz {Sxz}"
                if Sxy != "fixed":
                    ptext += " xy {Sxy}"

            if ptext == "couple none":
                ptext = ""

            Tmp = {}
            for key in ("Sxx", "Syy", "Szz", "Syz", "Sxz", "Sxy"):
                if _P[key] != "fixed":
                    Tmp[key] = lammps_step.to_lammps_units(
                        -_P[key], quantity="pressure"
                    )
            ptext = ptext.format(**Tmp)
        else:
            # hydrostatic pressure applied
            if couple == "x, y and z":
                ptext = "couple xyz x {P} y {P} z {P}"
            elif couple == "x and y":
                ptext = "couple xy x {P} y {P} z {P}"
            elif couple == "x and z":
                ptext = "couple xz x {P} y {P} z {P}"
            elif couple == "y and z":
                ptext = "couple yz x {P} y {P} z {P}"
            else:
                ptext = "couple none x {P} y {P} z {P}"

            if allow_shear:
                ptext += " xy 0.0 xz 0.0 yz 0.0"

            P = lammps_step.to_lammps_units(_P["P"], quantity="pressure")
            ptext = ptext.format(P=P)

        return ptext
