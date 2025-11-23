# -*- coding: utf-8 -*-

"""A single-point energy in LAMMPS"""

from collections import Counter
import json
import logging
from math import isnan
from pathlib import Path
import pkg_resources
import textwrap
import traceback

import numpy as np
import pandas
from tabulate import tabulate

import lammps_step
from lammps_step import LAMMPS, from_lammps_units
from molsystem import elements, Cell
import seamm
from seamm_util import units_class, Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
printer = printing.getPrinter("lammps")


_subscript = {
    "0": "\N{SUBSCRIPT ZERO}",
    "1": "\N{SUBSCRIPT ONE}",
    "2": "\N{SUBSCRIPT TWO}",
    "3": "\N{SUBSCRIPT THREE}",
    "4": "\N{SUBSCRIPT FOUR}",
    "5": "\N{SUBSCRIPT FIVE}",
    "6": "\N{SUBSCRIPT SIX}",
    "7": "\N{SUBSCRIPT SEVEN}",
    "8": "\N{SUBSCRIPT EIGHT}",
    "9": "\N{SUBSCRIPT NINE}",
}


def subscript(n):
    """Return the number using Unicode subscript characters."""
    return "".join([_subscript[c] for c in str(n)])


one_half = "\N{VULGAR FRACTION ONE HALF}"
degree_sign = "\N{DEGREE SIGN}"
standard_state = {
    "H": f"{one_half}H{subscript(2)}(g)",
    "He": "He(g)",
    "Li": "Li(s)",
    "Be": "Be(s)",
    "B": "B(s)",
    "C": "C(s,gr)",
    "N": f"{one_half}N{subscript(2)}(g)",
    "O": f"{one_half}O{subscript(2)}(g)",
    "F": f"{one_half}F{subscript(2)}(g)",
    "Ne": "Ne(g)",
    "Na": "Na(s)",
    "Mg": "Mg(s)",
    "Al": "Al(s)",
    "Si": "Si(s)",
    "P": "P(s)",
    "S": "S(s)",
    "Cl": f"{one_half}Cl{subscript(2)}(g)",
    "Ar": "Ar(g)",
    "K": "K(s)",
    "Ca": "Ca(s)",
    "Sc": "Sc(s)",
    "Ti": "Ti(s)",
    "V": "V(s)",
    "Cr": "Cr(s)",
    "Mn": "Mn(s)",
    "Fe": "Fe(s)",
    "Co": "Co(s)",
    "Ni": "Ni(s)",
    "Cu": "Cu(s)",
    "Zn": "Zn(s)",
    "Ga": "Ga(s)",
    "Ge": "Ge(s)",
    "As": "As(s)",
    "Se": "Se(s)",
    "Br": f"{one_half}Br{subscript(2)}(l)",
    "Kr": "(g)",
}


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

    def analyze(self, indent="", data={}, table=None, output=[], **kwargs):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Get the initial structure
        _, initial_configuration = self.get_system_configuration()

        # Handle the new structure as needed
        system, configuration = self.get_system_configuration(P, model=self.model)

        # See if there are other results in json files
        filename = self.wd / "energy.json"
        if filename.exists():
            try:
                with filename.open("r") as fd:
                    tmp = json.load(fd)
            except Exception as e:
                printer.normal(f"Warning: error reading {filename}: {e}")
                logger.warning(f"Error with {filename}: {e}")
            else:
                data.update(tmp)

        # Get the stress vector
        stress = []
        if initial_configuration.periodicity != 0:
            for key in ("Sxx", "Syy", "Szz", "Syz", "Sxz", "Sxy"):
                if key in data:
                    value = Q_(data[key], data[key + ",units"]).m_as("atm")
                    stress.append(value)
                    data[key] = value
                    data[key + ",units"] = "atm"
        if len(stress) == 6:
            data["stress"] = stress
            data["stress,units"] = "atm"

        # Pressure
        if "P" in data:
            data["P"] = Q_(data["P"], data["P,units"]).m_as("atm")
            data["P,units"] = "atm"

        # Energies
        for key in ("Epe", "Etot"):
            if key in data:
                data[key] = Q_(data[key], data[key + ",units"]).m_as("kcal/mol")
                data[key + ",units"] = "kcal/mol"

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
                    energy = from_lammps_units(float(tmp[i]), "kcal/mol").magnitude
                    data["energy"] = energy
                    data["energy,units"] = "kcal/mol"

                # Check for reaxff enthalpy offset
                if self.parent.ff_form() == "reaxff":
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

        # Calculate the enthalpy of formation, if possible
        tmp_text = self.calculate_enthalpy_of_formation(data)
        if tmp_text != "":
            path = self.wd / "Thermochemistry.txt"
            path.write_text(tmp_text)

        # Add static properties such as density for e.g NVE and NVT calculations
        if initial_configuration.periodicity == 3:
            if "density" not in data:
                data["density"] = initial_configuration.density
            if "a" not in data:
                data["a"] = initial_configuration.cell.a
            if "b" not in data:
                data["b"] = initial_configuration.cell.b
            if "c" not in data:
                data["c"] = initial_configuration.cell.c
            if "alpha" not in data:
                data["alpha"] = initial_configuration.cell.alpha
            if "beta" not in data:
                data["beta"] = initial_configuration.cell.beta
            if "gamma" not in data:
                data["gamma"] = initial_configuration.cell.gamma
            if "volume" not in data:
                data["volume"] = initial_configuration.cell.volume

        # Print the table of results
        if table is None:
            table = {
                "Property": [],
                "Value": [],
                "Units": [],
            }
        for key in ("DfH0_reax", "DfE0", "energy"):
            if key in data:
                table["Property"].append(key)
                table["Value"].append(f"{data[key]:.2f}")
                table["Units"].append("kcal/mol")
        for key in ("a", "b", "c"):
            if key in data:
                table["Property"].append(key)
                table["Value"].append(f"{data[key]:.3f}")
                table["Units"].append("Å")
        for key, char in zip(("alpha", "beta", "gamma"), ("𝞪", "𝞫", "𝞬")):
            if key in data:
                table["Property"].append(char)
                table["Value"].append(f"{data[key]:.3f}")
                table["Units"].append("°")
        for key in ("P", "Sxx", "Syy", "Szz", "Syz", "Sxz", "Sxy"):
            if key in data:
                table["Property"].append(key)
                table["Value"].append(f"{data[key]:.2f}")
                table["Units"].append("atm")
        for key in ("density",):
            if key in data:
                table["Property"].append(key)
                table["Value"].append(f"{data[key]:.3f}")
                table["Units"].append("g/mL")
        for key in ("volume",):
            if key in data:
                table["Property"].append(key)
                table["Value"].append(f"{data[key]:.1f}")
                table["Units"].append("Å^3")

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            colalign=("center", "decimal", "left"),
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        text_lines = []
        text_lines.append("Single Point Energy".center(length))
        text_lines.append(tmp)
        printer.normal(textwrap.indent("\n".join(text_lines), self.indent + 7 * " "))
        printer.normal("")

        # Read the dump file and get the structure
        xyz = None
        cell = None
        gradients = None
        velocities = None
        dump_file = self.wd / "energy.dump"
        if dump_file.exists:
            try:
                tmp = dump_file.read_text().splitlines()
                results = self._parse_dump_frame(tmp)
                if "xyz" in results:
                    xyz = results["xyz"]
                    fractional = False
                elif "abc" in results:
                    xyz = results["abc"]
                    fractional = True
                else:
                    xyz = None
                cell = results["cell"] if "cell" in results else None
                gradients = results["gradients"] if "gradients" in results else None
                velocities = results["velocities"] if "velocities" in results else None
                if gradients is not None:
                    data["gradients"] = gradients
                if velocities is not None:
                    data["velocities"] = velocities
            except Exception as e:
                printer.normal("Warning: unable to read the LAMMPS dumpfile")
                logger.warning(f"The was a problem reading the LAMMPS dumpfile: {e}")
                logger.warning(traceback.format_exc())
        else:
            printer.normal("Warning: there is no 'dump' file from LAMMPS")

        if configuration is not None:
            if cell is not None:
                configuration.cell.parameters = cell
            if xyz is not None:
                configuration.atoms.set_coordinates(xyz, fractionals=fractional)
            if gradients is not None:
                # LAMMPS only has Cartesian gradients, in kcal/mol/Å
                factor = Q_("kcal/mol/Å").m_as("kJ/mol/Å")
                tmp = factor * np.array(gradients)
                configuration.atoms.set_gradients(tmp, fractionals=False)
            if velocities is not None:
                # LAMMPS only has Cartesian velocities
                configuration.atoms.set_velocities(velocities, fractionals=False)

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

        # Put any requested results into variables or tables
        self.store_results(
            configuration=configuration,
            data=data | self.results,
            create_tables=self.parameters["create tables"].get(),
            printer=printer,
        )

    def calculate_enthalpy_of_formation(self, data):
        """Calculate the enthalpy of formation from the results of a calculation.

        This uses tabulated values of the enthalpy of formation of the atoms for
        the elements and tabulated energies calculated for atoms with the current
        method.

        Parameters
        ----------
        data : dict
            The results of the calculation.
        """
        J_to_cal = Q_("J").m_as("cal")

        # Read the tabulated values from either user or data directory
        personal_file = Path("~/.seamm.d/data/element_energies.csv").expanduser()
        if personal_file.exists():
            personal_table = pandas.read_csv(personal_file, index_col=False)
        else:
            personal_table = None

        path = Path(pkg_resources.resource_filename(__name__, "data/"))
        csv_file = path / "element_energies.csv"
        table = pandas.read_csv(csv_file, index_col=False)

        self.logger.debug(f"self.model = {self.model}")

        # Check if have the data
        atom_formation_energy = None
        atom_energy = None
        column = self.model

        self.logger.debug(f"Looking for '{column}'")

        # atom_formation_energy is the energy of formation of the standard state,
        # per atom.
        # atom_energy is the calculated energy of the atom, which defaults to zero
        column2 = column + " atom energy"
        if personal_table is not None and column in personal_table.columns:
            atom_formation_energy = (J_to_cal * personal_table[column]).to_list()
            if column2 in personal_table.columns:
                atom_energy = (J_to_cal * personal_table[column2]).to_list()
        elif column in table.columns:
            atom_formation_energy = (J_to_cal * table[column]).to_list()
            if column2 in table.columns:
                atom_energy = (J_to_cal * table[column2]).to_list()

        if atom_formation_energy is None:
            # Not found!
            return f"There are no tabulated atom energies for {column}"

        # Assume an offset energy -- the energy of an isolated atom -- is zero if not
        # tabulated
        if atom_energy is None:
            atom_energy = [0.0] * len(atom_formation_energy)

        DfH0gas = None
        references = None
        term_symbols = None
        if personal_table is not None and "ΔfH°gas" in personal_table.columns:
            DfH0gas = personal_table["ΔfH°gas"].to_list()
            if "Reference" in personal_table.columns:
                references = personal_table["Reference"].to_list()
            if "Term Symbol" in personal_table.columns:
                term_symbols = personal_table["Term Symbols"].to_list()
        elif "ΔfH°gas" in table.columns:
            DfH0gas = table["ΔfH°gas"].to_list()
            if "Reference" in table.columns:
                references = table["Reference"].to_list()
            if "Term Symbol" in table.columns:
                term_symbols = table["Term Symbol"].to_list()

        # Get the atomic numbers and counts
        _, configuration = self.get_system_configuration(None)
        counts = Counter(configuration.atoms.atomic_numbers)

        # Get the Hill formula as a list
        symbols = sorted(elements.to_symbols(counts.keys()))
        composition = []
        if "C" in symbols:
            composition.append((6, "C", counts[6]))
            symbols.remove("C")
            if "H" in symbols:
                composition.append((1, "H", counts[1]))
                symbols.remove("H")

        for symbol in symbols:
            atno = elements.symbol_to_atno[symbol]
            composition.append((atno, symbol, counts[atno]))

        # And the reactions. First, for atomization energy
        middot = "\N{MIDDLE DOT}"
        lDelta = "\N{GREEK CAPITAL LETTER DELTA}"
        formula = ""
        tmp = []
        for atno, symbol, count in composition:
            if count == 1:
                formula += symbol
                tmp.append(f"{symbol}(g)")
            else:
                formula += f"{symbol}{subscript(count)}"
                tmp.append(f"{count}{middot}{symbol}(g)")
        gas_atoms = " + ".join(tmp)
        tmp = []
        for atno, symbol, count in composition:
            if count == 1:
                tmp.append(standard_state[symbol])
            else:
                tmp.append(f"{count}{middot}{standard_state[symbol]}")
        standard_elements = " + ".join(tmp)

        # The energy - any offsets is the negative of the atomization energy
        name = "Formula: " + formula
        try:
            name = configuration.PC_iupac_name(fallback=name)
        except Exception:
            pass

        if name is None:
            name = "Formula: " + formula

        text = f"Thermochemistry of {name} with {column}\n\n"
        text += "Atomization Energy\n"
        text += "------------------\n"
        text += textwrap.fill(
            f"The atomization energy,  {lDelta}atE{degree_sign}, is the energy to break"
            " all the bonds in the system, separating the atoms from each other."
        )
        text += f"\n\n    {formula} --> {gas_atoms}\n\n"
        text += textwrap.fill(
            "The following table shows in detail the calculation. The first line is "
            "the system and its calculated energy. The next lines are the energies "
            "of each type of atom in the system. These have been tabulated by running "
            "calculations on each atom, and are included in the SEAMM release. "
            "The line give the formation energy from atoms in kcal/mol.",
        )
        text += "\n\n"
        table = {
            "System": [],
            "Term": [],
            "Value": [],
            "Units": [],
        }

        if "Epe" in data:
            E = data["Epe"]
        elif "energy" in data:
            E = data["energy"]
        else:
            return "The energy is not in results from the calculation!"

        Eatoms = 0.0
        Ef0 = 0.0
        for atno, symbol, count in composition:
            Eatom = atom_energy[atno - 1]
            if isnan(Eatom):
                # Don't have the data for this element
                return f"Do not have tabulated atom energies for {symbol} in {column}"
            Eatoms += count * Eatom
            table["System"].append(f"{symbol}(g)")
            table["Term"].append(f"{count} * {Eatom:.2f}")
            table["Value"].append(f"{count * Eatom:.2f}")
            table["Units"].append("")

            Ef0 += count * atom_formation_energy[atno - 1]

        data["DfE0"] = E - Ef0

        table["Units"][0] = "kcal/mol"

        table["System"].append("^")
        table["Term"].append("-")
        table["Value"].append("-")
        table["Units"].append("")

        table["System"].append(formula)
        table["Term"].append(f"{-E:.2f}")
        table["Value"].append(f"{-E:.2f}")
        table["Units"].append("kcal/mol")

        data["E atomization"] = Eatoms - E

        table["System"].append("")
        table["Term"].append("")
        table["Value"].append("=")
        table["Units"].append("")

        result = f'{data["E atomization"]:.2f}'
        table["System"].append(f"{lDelta}atE")
        table["Term"].append("")
        table["Value"].append(result)
        table["Units"].append("kcal/mol")

        factor = Q_("kcal/mol").m_as("kcal/mol")
        table["System"].append("")
        table["Term"].append("")
        table["Value"].append(f'{factor * data["E atomization"]:.2f}')
        table["Units"].append("kcal/mol")

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            colalign=("center", "center", "decimal", "center"),
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        text_lines = []
        text_lines.append(f"Atomization Energy for {formula}".center(length))
        text_lines.append(tmp)
        text += textwrap.indent("\n".join(text_lines), 4 * " ")

        if "H" not in data:
            text += "\n\n"
            text += "Cannot calculate enthalpy of formation without the enthalpy"
            return text
        if DfH0gas is None:
            text += "\n\n"
            text += "Cannot calculate enthalpy of formation without the tabulated\n"
            text += "atomization enthalpies of the elements."
            return text

        # Atomization enthalpy of the elements, experimental
        table = {
            "System": [],
            "Term": [],
            "Value": [],
            "Units": [],
            "Reference": [],
        }

        E = data["energy"]

        DfH_at = 0.0
        refno = 1
        for atno, symbol, count in composition:
            DfH_atom = DfH0gas[atno - 1]
            DfH_at += count * DfH_atom
            tmp = Q_(DfH_atom, "kcal/mol").m_as("E_h")
            table["System"].append(f"{symbol}(g)")
            if count == 1:
                table["Term"].append(f"{tmp:.6f}")
            else:
                table["Term"].append(f"{count} * {tmp:.6f}")
            table["Value"].append(f"{count * tmp:.6f}")
            table["Units"].append("")
            refno += 1
            table["Reference"].append(refno)

        table["Units"][0] = "E_h"

        table["System"].append("^")
        table["Term"].append("-")
        table["Value"].append("-")
        table["Units"].append("")
        table["Reference"].append("")

        table["System"].append(standard_elements)
        table["Term"].append("")
        table["Value"].append("0.0")
        table["Units"].append("E_h")
        table["Reference"].append("")

        table["System"].append("")
        table["Term"].append("")
        table["Value"].append("=")
        table["Units"].append("")
        table["Reference"].append("")

        result = f'{Q_(DfH_at, "kcal/mol").m_as("E_h"):.6f}'
        table["System"].append(f"{lDelta}atH{degree_sign}")
        table["Term"].append("")
        table["Value"].append(result)
        table["Units"].append("E_h")
        table["Reference"].append("")

        table["System"].append("")
        table["Term"].append("")
        table["Value"].append(f"{DfH_at:.2f}")
        table["Units"].append("kcal/mol")
        table["Reference"].append("")

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            colalign=("center", "center", "decimal", "center", "center"),
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        text_lines = []
        text_lines.append(
            "Atomization enthalpy of the elements (experimental)".center(length)
        )
        text_lines.append(tmp)

        text += "\n\n"
        text += "Enthalpy of Formation\n"
        text += "---------------------\n"
        text += textwrap.fill(
            f"The enthalpy of formation, {lDelta}fHº, is the enthalpy of creating the "
            "molecule from the elements in their standard state:"
        )
        text += f"\n\n   {standard_elements} --> {formula} (1)\n\n"
        text += textwrap.fill(
            "The standard state of the element, denoted by the superscript º,"
            " is its form at 298.15 K and 1 atm pressure, e.g. graphite for carbon, "
            "H2 gas for hydrogen, etc."
        )
        text += "\n\n"
        text += textwrap.fill(
            "Since it is not easy to calculate the enthalpy of e.g. graphite we will "
            "use two sequential reactions that are equivalent. First, we will create "
            "gas phase atoms from the elements:"
        )
        text += f"\n\n    {standard_elements} --> {gas_atoms} (2)\n\n"
        text += textwrap.fill(
            "This will use the experimental values of the enthalpy of formation of the "
            "atoms in the gas phase to calculate the enthalpy of this reaction. "
            "Then we react the atoms to get the desired system:"
        )
        text += f"\n\n    {gas_atoms} --> {formula} (3)\n\n"
        text += textwrap.fill(
            "Note that this is reverse of the atomization reaction, so "
            f"{lDelta}H = -{lDelta}atH."
        )
        text += "\n\n"
        text += textwrap.fill(
            "First we calculate the enthalpy of the atomization of the elements in "
            "their standard state, using tabulated experimental values:"
        )
        text += "\n\n"
        text += textwrap.indent("\n".join(text_lines), 4 * " ")

        # And the calculated atomization enthalpy
        table = {
            "System": [],
            "Term": [],
            "Value": [],
            "Units": [],
        }

        Hatoms = 0.0
        dH = Q_(6.197, "kcal/mol").m_as("E_h")
        for atno, symbol, count in composition:
            Eatom = atom_formation_energy[atno - 1]
            # 6.197 is the H298-H0 for an atom
            Hatoms += count * (Eatom + 6.197)

            table["System"].append(f"{symbol}(g)")
            if count == 1:
                table["Term"].append(f"{-Eatom:.2f} + {dH:.2f}")
            else:
                table["Term"].append(f"{count} * ({-Eatom:.2f} + {dH:.2f})")
            table["Value"].append(f"{-count * (Eatom + dH):.2f}")
            table["Units"].append("")

        table["System"].append("^")
        table["Term"].append("-")
        table["Value"].append("-")
        table["Units"].append("")

        H = data["H"]

        table["System"].append(formula)
        table["Term"].append(f"{H:.2f}")
        table["Value"].append("")
        table["Units"].append("kcal/mol")

        data["H atomization"] = Hatoms - Q_(H, "E_h").m_as("kcal/mol")
        data["DfH0"] = DfH_at - data["H atomization"]
        table["System"].append("")
        table["Term"].append("")
        table["Value"].append("=")
        table["Units"].append("")

        table["System"].append("")
        table["Term"].append("")
        table["Value"].append(f'{data["H atomization"]:.2f}')
        table["Units"].append("kcal/mol")

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            colalign=("center", "center", "decimal", "center"),
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        text_lines = []
        text_lines.append("Atomization Enthalpy (calculated)".center(length))
        text_lines.append(tmp)
        text += "\n\n"

        text += textwrap.fill(
            "Next we calculate the atomization enthalpy of the system. We have the "
            "calculated enthalpy of the system, but need the enthalpy of gas phase "
            f"atoms at the standard state (25{degree_sign}C, 1 atm). The tabulated "
            "energies for the atoms, used above, are identical to H0 for an atom. "
            "We will add H298 - H0 to each atom, which [1] is 5/2RT = 0.002360 E_h"
        )
        text += "\n\n"
        text += textwrap.indent("\n".join(text_lines), 4 * " ")
        text += "\n\n"
        text += textwrap.fill(
            "The enthalpy change for reaction (3) is the negative of this atomization"
            " enthalpy. Putting the two reactions together with the negative for Rxn 3:"
        )
        text += "\n\n"
        text += f"{lDelta}fH{degree_sign} = {lDelta}H(rxn 2) - {lDelta}H(rxn 3)\n"
        text += f"     = {DfH_at:.2f} - {data['H atomization']:.2f}\n"
        text += f"     = {DfH_at - data['H atomization']:.2f} kcal/mol\n"

        text += "\n\n"
        text += "References\n"
        text += "----------\n"
        text += "1. https://en.wikipedia.org/wiki/Monatomic_gas\n"
        refno = 1
        for atno, symbol, count in composition:
            refno += 1
            text += f"{refno}. {lDelta}fH{degree_sign} = {DfH0gas[atno - 1]} kcal/mol"
            if term_symbols is not None:
                text += f" for {term_symbols[atno - 1]} {symbol}"
            else:
                text += f" for {symbol}"
            if references is not None:
                text += f" from {references[atno-1]}\n"

        return text

    def get_input(self, extras=None):
        """Get the input for an energy calculation for LAMMPS"""

        P = self.parameters.values_to_dict()

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        _, configuration = self.get_system_configuration()

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

        if configuration.periodicity != 0:
            # Write out the stress
            filename = f"@{self._id[-1]}+energy.json"
            units = lammps_step.lammps_units("pressure")
            eunits = lammps_step.lammps_units("energy")
            lines.append(
                'print               """{\n'
                '    "P": $(v_press:%.3f),\n'
                f'    "P,units": "{units}",\n'
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
                f'    "Sxy,units": "{units}",\n'
                '    "Etot": $(v_etotal:%.4f),\n'
                f'    "Etot,units": "{eunits}",\n'
                '    "Epe": $(v_etotal:%.4f),\n'
                f'    "Epe,units": "{eunits}"\n'
                "}"
                f'""" file {filename}'
            )
        filename = f"@{self._id[-1]}+energy.dump"
        lines.append(
            f"write_dump          all custom  {filename} id xu yu zu fx fy fz"
            " modify flush yes sort id"
        )

        return {
            "script": lines,
            "postscript": None,
            "use python": False,
        }

    def _parse_dump_frame(self, lines, line_no=0):
        """Process a single frame of a LAMMPS dump file.

        Parameters
        ----------
        lines : [str]
            The lines of the section.

        lineno : int
            The count of lines before this frame, used for error messages

        Returns
        -------
        dict(str, any)
            The results keyed by the item name.
            The units are Å, kcal/mol, and fs
        """
        results = {}
        it = iter(lines)

        def next_line(it, line_no):
            return (next(it).strip(), line_no + 1)

        line, line_no = next_line(it, line_no)

        while True:
            try:
                if not line.startswith("ITEM:"):
                    raise ValueError(
                        f"Error reading dump file: line {line_no} should be 'ITEM: ...'"
                        f"\n\t'{line}'"
                    )
                section = line[6:].strip()
                self.logger.debug("   section = " + section)

                if section == "TIMESTEP":
                    line, line_no = next_line(it, line_no)
                    results["timestep"] = int(line)
                    line, line_no = next_line(it, line_no)
                    continue
                elif section == "NUMBER OF ATOMS":
                    line, line_no = next_line(it, line_no)
                    results["n_atoms"] = n_atoms = int(line)
                    line, line_no = next_line(it, line_no)
                    continue
                elif "BOX BOUNDS" in section:
                    factor = from_lammps_units(1, "Å").magnitude
                    if len(section.split()) == 8:
                        line, line_no = next_line(it, line_no)
                        xlo_bound, xhi_bound, xy = line.split()
                        line, line_no = next_line(it, line_no)
                        ylo_bound, yhi_bound, xz = line.split()
                        line, line_no = next_line(it, line_no)
                        zlo, zhi, yz = line.split()

                        xlo_bound = float(xlo_bound) * factor
                        xhi_bound = float(xhi_bound) * factor
                        ylo_bound = float(ylo_bound) * factor
                        yhi_bound = float(yhi_bound) * factor
                        zlo = float(zlo) * factor
                        zhi = float(zhi) * factor
                        xy = float(xy) * factor
                        xz = float(xz) * factor
                        yz = float(yz) * factor

                        xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
                        xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
                        ylo = ylo_bound - min(0.0, yz)
                        yhi = yhi_bound - max(0.0, yz)
                        cell = LAMMPS.box_to_cell(
                            xhi - xlo, yhi - ylo, zhi - zlo, xy, xz, yz
                        )
                        tmp = Cell(*cell)
                        lattice = []
                        for v in tmp.vectors():
                            lattice.extend(v)
                    else:
                        line, line_no = next_line(it, line_no)
                        xlo, xhi = line.split()
                        line, line_no = next_line(it, line_no)
                        ylo, yhi = line.split()
                        line, line_no = next_line(it, line_no)
                        zlo, zhi = line.split()

                        xlo = float(xlo) * factor
                        xhi = float(xhi) * factor
                        ylo = float(ylo) * factor
                        yhi = float(yhi) * factor
                        zlo = float(zlo) * factor
                        zhi = float(zhi) * factor

                        lattice = (
                            xhi - xlo,
                            0.0,
                            0.0,
                            0.0,
                            yhi - ylo,
                            0.0,
                            0.0,
                            0.0,
                            zhi - zlo,
                        )
                        cell = (xhi - xlo, yhi - ylo, zhi - zlo, 90, 90, 90)
                    results["cell"] = cell
                    results["lattice"] = lattice
                    line, line_no = next_line(it, line_no)
                elif "ATOMS" in section:
                    xyz = []
                    f = []
                    v = []
                    keys = section.split()[1:]

                    # Coordinates
                    ixyz = None
                    if "x" in keys and "y" in keys and "z" in keys:
                        ixyz = (keys.index("x"), keys.index("y"), keys.index("z"))
                        xyz_type = "cartesian"
                        fxyz = from_lammps_units(1, "Å").magnitude
                    elif "xu" in keys and "yu" in keys and "zu" in keys:
                        ixyz = (keys.index("xu"), keys.index("yu"), keys.index("zu"))
                        xyz_type = "cartesian"
                        fxyz = from_lammps_units(1, "Å").magnitude
                    elif "xs" in keys and "ys" in keys and "zs" in keys:
                        ixyz = (keys.index("xs"), keys.index("ys"), keys.index("zs"))
                        xyz_type = "fractional"
                        fxyz = 1
                    elif "xsu" in keys and "ysu" in keys and "zsu" in keys:
                        ixyz = (keys.index("xsu"), keys.index("ysu"), keys.index("zsu"))
                        xyz_type = "fractional"
                        fxyz = 1

                    # Forces
                    fi = None
                    if "fx" in keys and "fy" in keys and "fz" in keys:
                        fi = (keys.index("fx"), keys.index("fy"), keys.index("fz"))
                        ff = -from_lammps_units(1, "kcal/mol/Å").magnitude

                    # Velocities
                    vi = None
                    if "vx" in keys and "vy" in keys and "vz" in keys:
                        vi = (keys.index("vx"), keys.index("vy"), keys.index("vz"))
                        fv = from_lammps_units(1, "Å/fs").magnitude

                    for i in range(n_atoms):
                        line, line_no = next_line(it, line_no)
                        values = line.split()
                        if ixyz is not None:
                            xyz.append(
                                [
                                    fxyz * float(values[ixyz[0]]),
                                    fxyz * float(values[ixyz[1]]),
                                    fxyz * float(values[ixyz[2]]),
                                ]
                            )
                        if fi is not None:
                            f.append(
                                [
                                    ff * float(values[fi[0]]),
                                    ff * float(values[fi[1]]),
                                    ff * float(values[fi[2]]),
                                ]
                            )
                        if vi is not None:
                            v.append(
                                [
                                    fv * float(values[vi[0]]),
                                    fv * float(values[vi[1]]),
                                    fv * float(values[vi[2]]),
                                ]
                            )
                    if ixyz is not None:
                        if xyz_type == "cartesian":
                            results["xyz"] = xyz
                        elif xyz_type == "fractional":
                            results["abc"] = xyz
                    if fi is not None:
                        results["gradients"] = f
                    if vi is not None:
                        results["velocities"] = v

                    line, line_no = next_line(it, line_no)
                else:
                    raise ValueError(
                        f"Don't recognize section {section} in the dump file"
                    )
            except StopIteration:
                break

        return results
