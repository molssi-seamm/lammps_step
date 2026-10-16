# -*- coding: utf-8 -*-

"""NVE (microcanonical) dynamics in LAMMPS"""

import bz2
import gzip
import json
import os
from pathlib import Path
import traceback

import numpy as np
from tabulate import tabulate

import lammps_step
from lammps_step import from_lammps_units
from lammps_step.md_h5md import write_md_h5md
import logging
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("lammps")


class NVE(lammps_step.Energy):
    # Molecular dynamics: the model-chemistry task for the provenance label.
    task = "MD"

    def __init__(
        self,
        flowchart=None,
        title="NVE dynamics",
        extension=None,
        logger=logger,
    ):
        """Initialize the node"""

        logger.debug("Creating NVE {}".format(self))

        super().__init__(
            flowchart=flowchart,
            title=title,
            extension=extension,
            logger=logger,
        )
        self.logger.debug("NVE.init() creating NVE_Parameters object")

        self._calculation = "nve"
        self._metadata = lammps_step.metadata
        self.parameters = lammps_step.NVE_Parameters()

    def analyze(self, data={}, properties=None, table=None, output=""):
        """Analyze the results of the simulation.

        Parameters
        ----------
        """
        # Calculate the enthalpy of formation, if possible
        tmp_text = self.calculate_enthalpy_of_formation(data)
        if tmp_text != "":
            path = self.wd / "Thermochemistry.txt"
            path.write_text(tmp_text)

        # Print the results from the analysis of the state trajectory
        for i, _property in enumerate(table["Property"]):
            if _property.startswith("Epe"):
                stderr = table["StdErr"][i]
                i += 1
                new = {k: v[0:i] for k, v in table.items()}
                for key in ("DfE0", "E atomization"):
                    if key in data:
                        new["Property"].append(key)
                        new["Value"].append(f"{data[key]:.2f}")
                        new[" "].append("±")
                        new["StdErr"].append(stderr)
                        new["Units"].append("kcal/mol")
                        new["convergence"].append("")
                        new["tau"].append("")
                        new["inefficiency"].append("")
                for key, value in table.items():
                    new[key].extend(table[key][i:])
                table = new
                break
        else:
            for key in ("DfE0", "E atomization"):
                if key in data:
                    table["Property"].append(key)
                    table["Value"].append(f"{data[key]:.2f}")
                    table[" "].append("±")
                    table["StdErr"].append("")
                    table["Units"].append("kcal/mol")
                    table["convergence"].append("")
                    table["tau"].append("")
                    table["inefficiency"].append("")

        # Print out a table of results.
        text = ""
        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
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
                have_acf_warning = True
                break

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

        if have_warning or have_acf_warning:
            printer.normal("\n")

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Get the initial structure
        _, initial_configuration = self.get_system_configuration()

        # Handle the new structure as needed
        system, configuration = self.get_system_configuration(P, model=self.model)

        # Read the dump file and get the structure
        stem = self.calculation.lower().replace(" ", "_")  # "Heat Flux" -> heat_flux
        dump_file = Path(self.directory) / f"{stem}.dump"
        if dump_file.exists():
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

        # The user's extra properties and their treatments
        if properties is not None:
            text = self._analyze_extra_properties(properties, data, configuration)
            if text != "":
                printer.normal(__(text, indent=8 * " ", wrap=False, dedent=False))
                printer.normal("")

        # Get performance data from the output
        if output != "":
            table = {
                "Metric": [],
                "Value": [],
                "Units": [],
            }
            table2 = {
                "Section": [],
                "Value": [],
                "%Variance": [],
                "%": [],
            }
            lines = iter(output)
            for line in lines:
                if line.startswith("Loop time of"):
                    tmp = line.split()
                    data["loop time"] = tmp[3]
                    data["nprocs"] = tmp[5]
                    data["nsteps"] = tmp[8]
                    data["natoms"] = tmp[11]
                    for key in ("loop time", "nprocs", "nsteps", "natoms"):
                        table["Metric"].append(key)
                        table["Value"].append(data[key])
                    table["Units"] = ["s", "", "", ""]
                elif line.startswith("Performance:"):
                    tmp = line.split()
                    data["rate per day"] = Q_(float(tmp[1]), tmp[2].rstrip(",")).m_as(
                        "ns/day"
                    )
                    data["rate timesteps/s"] = float(tmp[5])
                    if tmp[6].lower().startswith("k"):
                        data["rate timestep/s"] *= 1000
                    data["katom-steps/s"] = float(tmp[7])
                    if tmp[8].lower().startswith("m"):
                        data["katom-steps/s"] *= 1000
                    line = next(lines)
                    tmp = line.split()
                    data["%cpu"] = float(line[0].rstrip("%"))

                    table["Metric"].append("rate(per day)")
                    table["Value"].append(data["rate per day"])
                    table["Units"].append("ns/day")
                    table["Metric"].append("rate(timesteps)")
                    table["Value"].append(data["rate timesteps/s"])
                    table["Units"].append("timesteps/s")
                    table["Metric"].append("katom-steps/s")
                    table["Value"].append(data["katom-steps/s"])
                    table["Units"].append("katom*steps/s")
                    table["Metric"].append("%CPU")
                    table["Value"].append(data["%cpu"])
                    table["Units"].append("%")
                elif line.startswith("MPI task timing breakdown:"):
                    imbalance = []
                    comms = 0
                    detail = {}
                    next(lines)
                    next(lines)
                    # Pair    | 1.0275     | 1.2183     | 1.3807     |  12.1 | 24.85
                    # Other   |            | 0.04275    |            |       |  0.90
                    for line in lines:
                        if line.strip() == "":
                            break
                        (
                            section,
                            min_time,
                            avg_time,
                            max_time,
                            varavg,
                            percent,
                        ) = line.split("|")
                        section = section.strip()
                        min_time = min_time.strip()
                        avg_time = avg_time.strip()
                        max_time = max_time.strip()
                        varavg = varavg.strip()
                        percent = percent.strip()

                        detail[section] = [
                            min_time,
                            avg_time,
                            max_time,
                            varavg,
                            percent,
                        ]
                        if varavg != "" and float(varavg) > 5:
                            imbalance.append(section)
                        if section == "Comm":
                            comms = float(percent)
                        table2["Section"].append(section)
                        table2["Value"].append(avg_time)
                        table2["%Variance"].append(varavg)
                        table2["%"].append(percent)
                    break
            # Print out a table of results.
            text = ""
            tmp = tabulate(
                table,
                headers="keys",
                tablefmt="rounded_outline",
                disable_numparse=True,
                colalign=(
                    "center",
                    "decimal",
                    "center",
                ),
            )
            length = len(tmp.splitlines()[0])
            text += "\n"
            text += "Performance".center(length)
            text += "\n"
            text += tmp
            text += "\n"

            printer.normal(__(text, indent=8 * " ", wrap=False, dedent=False))

            if len(table2["Section"]) > 0:
                # Print out a table of timings per part of lammps.
                text = ""
                tmp = tabulate(
                    table2,
                    headers="keys",
                    tablefmt="rounded_outline",
                    disable_numparse=True,
                    colalign=(
                        "center",
                        "decimal",
                        "decimal",
                        "decimal",
                    ),
                )
                length = len(tmp.splitlines()[0])
                text += "\n"
                text += "Performance of Sections of LAMMPS".center(length)
                text += "\n"
                text += tmp
                text += "\n"

                printer.normal(__(text, indent=8 * " ", wrap=False, dedent=False))

                if len(imbalance) > 0 or comms > 10:
                    text = ""
                    if len(imbalance) > 0:
                        tmp = ", ".join(imbalance)
                        text += f"The variance of > 5% in the timing of {tmp}"
                        text += "between processors suggests that LAMMPS is not "
                        text += "partitioning the particles effectively. If possible "
                        text += "turn on balancing either once before the run or "
                        text += "at an interval during the run."
                        text += "\n\n"
                    if comms > 10:
                        text += f"The communication overhead of {comms}% suggests "
                        text += "that you are running on too many cores."
                        text += "\n\n"
                    printer.normal(__(text, indent=8 * " ", wrap=True, dedent=True))

        # Write the MD data as H5MD for the analysis steps
        text = write_md_h5md(
            self, P, initial_configuration, wall_time=data.get("loop time")
        )
        printer.normal(__(text, indent=self.indent + 4 * " "))
        printer.normal("")

        # Save the trajectory to a new system and its configurations
        if P["trajectory save"]:
            text = self._save_trajectory(P)
            printer.normal(__(text, indent=self.indent + 4 * " "))
            printer.normal("")

        # Import or convert the trajectory.
        if P["trajectory export"]:
            text = self._export_trajectory(P)
            printer.normal(__(text, indent=self.indent + 4 * " "))
            printer.normal("")

        if configuration is not None:
            self.store_results(configuration=configuration, data=data, printer=printer)

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

        model = self.model

        text = "{time} of microcanonical (NVE) dynamics using a "
        text += "timestep of {timestep}"
        if model is None:
            text += "."
        else:
            text += f" using the {model} forcefield. "
        text += "The trajectory will be sampled every {sampling}."

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def _export_trajectory(self, P):
        """Export the trajectory in requested formats.

        Parameters
        ----------
        P : dict(str, value)
            The control parameters for this step

        Returns
        -------
        str
            A text string for printing, containing warnings, etc.
        """
        text = ""
        dump = Path(self.directory) / "trajectory.dump_trj.gz"
        trj = Path(self.directory) / "trajectory.trj"

        if not dump.exists() and not trj.exists():
            return (
                "Neither the coordinate or energy part of the trajectory exists, so "
                "cannot export it."
            )

        if P["trajectory extxyz"]:
            text += self._convert_trajectory_to_extxyz(dump, trj, P)

        return text

    def _convert_trajectory_to_extxyz(self, dump, trj, P):
        """Convert the trajectory to the ASE extended XYZ format.

        Parameters
        ----------
        dump : Path
            The path to the dump part of the trajectory.
        trj : Path
            The path to the state (energy, pressure) part of the trajectory
        P : dict(str, value)
            The control parameters for this step

        Returns
        -------
        str
            Text for printing explaining the process
        """
        text = ""
        if not dump.exists():
            return (
                "The coordinate dump of the trajectory does not exist, so cannot "
                "create the extxyz file."
            )
        if trj.exists():
            trj_data = trj.read_text().splitlines()
            header = trj_data[0]
            if not header.startswith("!MolSSI trajectory"):
                trj_data = None
                text = (
                    "The file for the energies etc. for the trajectory exists "
                    "but is not the right type of file, so the energies, etc. "
                    "will not be in the extxyz file"
                )
            else:
                state_variables = trj_data[1].split()
                trj_data = [i.split() for i in trj_data[2:]]
        else:
            trj_data = None
            text = (
                "The energies etc. for the trajectory are missing so they will "
                "not be in the extxyz file"
            )

        # Do it!
        _, timestep = self.timestep(P["timestep"])
        _, configuration = self.get_system_configuration()
        symbols = configuration.atoms.symbols

        filename = P["trajectory extxyz filename"].strip()
        if filename.startswith("/"):
            path = self.job_path / filename[1:]
        else:
            path = self.wd / Path(filename)

        nskip = P["trajectory extxyz skip frames"]
        frame = 0
        start_line = 0
        lines = []

        append = P["trajectory extxyz append"]
        mode = "a" if append else "w"
        text += "Appended" if append else "Wrote"
        with (
            gzip.open(dump, "rt") as fdin,
            (
                gzip.open(path, mode=mode + "t")
                if path.suffix == ".gz"
                else (
                    bz2.open(path, mode=mode + "t")
                    if path.suffix == ".bz2"
                    else open(path, mode)
                )
            ) as fdout,
        ):
            for line in fdin:
                line = line.strip()
                if line.startswith("ITEM: TIMESTEP"):
                    if len(lines) > 0:
                        frame += 1
                        if frame > nskip:
                            results = self._parse_dump_frame(lines, start_line)
                            results["symbols"] = symbols
                            results.update(
                                {
                                    k: v
                                    for k, v in zip(
                                        state_variables, trj_data[frame - 1]
                                    )
                                }
                            )
                            extxyz = self._to_extxyz(results, timestep=timestep)
                            fdout.write(extxyz)
                        start_line += len(lines)
                        lines = []
                lines.append(line)
            # Process last frame
            if len(lines) > 0:
                frame += 1
                if frame > nskip:
                    results = self._parse_dump_frame(lines, start_line)
                    results["symbols"] = symbols
                    results.update(
                        {k: v for k, v in zip(state_variables, trj_data[frame - 1])}
                    )
                    extxyz = self._to_extxyz(results, timestep=timestep)
                    fdout.write(extxyz)
        if nskip > 0:
            text += (
                f" {frame - nskip} frames of the trajectory to {path}, "
                f"skipping the first {nskip} frames."
            )
        else:
            text += f" {frame - nskip} frames of the trajectory to {path}."
        return text

    def get_input(self, extras=None):
        """Get the input for an NVE dynamics run in LAMMPS"""

        # See what type of forcefield we have and handle it
        ff_form = self.parent.ff_form()

        self.description = []
        self.description.append(__(self.header, indent=4 * " "))

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        _, configuration = self.get_system_configuration()

        timestep, P["timestep"] = self.timestep(P["timestep"])
        sampling_warning = self._drop_periodic_sampling(P, configuration)

        if extras is not None and "nsteps" in extras:
            nsteps = extras["nsteps"]
        else:
            # both in LAMMPS's units: fs for "real", ps for "metal" (MLFFs)
            time = lammps_step.to_lammps_units(P["time"], quantity="time")
            nsteps = round(time / timestep)

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        self.description.append(__(self.description_text(), **PP, indent=7 * " "))
        if sampling_warning != "":
            self.description.append(__(sampling_warning, indent=7 * " "))

        form = self.parent.ff_form()
        thermo_properties = "time temp press etotal ke pe"
        # PyTorch (MACE) and MDI/QM get their energy from an external engine, so
        # LAMMPS computes no bonded/pair terms -- those columns would be all zero.
        if form not in ("PyTorch", "MDI/QM"):
            thermo_properties += " ebond eangle edihed eimp evdwl etail ecoul elong"
        if configuration.periodicity == 3:
            properties = "v_time v_temp v_press v_etotal v_ke v_pe v_emol v_epair"
            title2 = "tstep t T P Etot Eke Epe Emol Epair"
        else:
            # The pressure of a non-periodic system is meaningless
            properties = "v_time v_temp v_etotal v_ke v_pe v_emol v_epair"
            title2 = "tstep t T Etot Eke Epe Emol Epair"
        if configuration.periodicity == 3:
            properties += " v_sxx v_syy v_szz v_syz v_sxz v_sxy"
            title2 += " Sxx Syy Szz Syz Sxz Sxy"
        if self.parent.have_dreiding_hbonds:
            thermo_properties += " v_N_hbond v_E_hbond"
            properties += " v_N_hbond v_E_hbond"
            title2 += " N_hbond E_hbond"
        # The user's extra commands and properties for this sub-step
        # (not "extras": that is the argument carrying SHAKE and the steps)
        user_extras = self.extra_input(P)
        if user_extras["properties"] != "":
            properties += " " + user_extras["properties"]
            title2 += " " + user_extras["titles"]

        lines = []
        nfixes = 0
        ncomputes = 0
        ndumps = 0
        lines.append("")
        lines.append(f"# {self.header}")
        lines.append("")
        lines.append("reset_timestep      0")
        lines.append("timestep            {}".format(timestep))
        lines.append("thermo_style        custom {}".format(thermo_properties))
        lines.append("thermo              {}".format(int(nsteps / 100)))
        nfixes += 1
        lines.append("fix                 {} all nve".format(nfixes))

        lines.extend(user_extras["lines"])

        # For the heat flux, if requested, we need extra input
        if P["heat flux"] != "never":
            lines.append(self.heat_flux_input(P, ff_form))

        # The summary: 10 lines over the run, each the average over its tenth of
        # the run, sampled every 10 steps (every step of a run of fewer than 100).
        # 'off 1' keeps the time the time of the line.
        nevery = max(1, min(10, nsteps // 10))
        nfreq = int(nsteps / 10)
        nrepeat = max(1, int(nfreq / nevery))
        nfreq = nevery * nrepeat
        nfixes += 1
        filename = f"@{self._id[-1]}+nve_summary.trj"
        lines.append(
            f"fix                 {nfixes} all ave/time {nevery} {nrepeat} {nfreq} &\n"
            f"                       {properties} &\n"
            "                       off 1 &\n"
            f"                       file {filename}"
        )

        # instantaneous output written for averaging
        if P["sampling"] == "none":
            self.description.append(
                __(
                    "The run will be {nsteps:n} steps of dynamics.",
                    nsteps=nsteps,
                    indent=7 * " ",
                )
            )
        else:
            sampling = lammps_step.to_lammps_units(P["sampling"], quantity="time")
            nevery = max(1, round(sampling / timestep))
            nfreq = int(nsteps / nevery)
            nrepeat = 1
            nfreq = nevery * nrepeat
            nfixes += 1
            dt = (nevery * P["timestep"]).to_compact()
            text = json.dumps(
                {
                    "code": "LAMMPS",
                    "type": "NVE",
                    "dt": dt.magnitude,
                    "tunits": str(dt.u),
                    "nsteps": nsteps // nevery,
                },
                separators=(",", ":"),
            )
            title1 = "!MolSSI trajectory 2.0 " + text
            filename = f"@{self._id[-1]}+nve_state.trj"
            lines.append(
                f"fix                 {nfixes} all ave/time {nevery} 1 {nfreq} &\n"
                f"                       {properties} &\n"
                "                       off 1 &\n"
                f"                       title1 '{title1}' &\n"
                f"                       title2 '{title2}' &\n"
                f"                       file {filename}"
            )
            self.description.append(
                __(
                    (
                        "The run will be {nsteps:,d} steps of dynamics "
                        "sampled every {nevery:n} steps."
                    ),
                    nsteps=nsteps,
                    nevery=nevery,
                    indent=7 * " ",
                )
            )

        # Handle trajectories
        tmp, ncomputes, ndumps, nfixes = self.trajectory_input(
            P, timestep, nsteps, ncomputes, ndumps, nfixes
        )
        lines.extend(tmp)

        if extras is not None and "shake" in extras:
            nfixes += 1
            lines.append(extras["shake"].format(nfixes))

        lines.append("")
        lines.append("run                 {}".format(nsteps))
        lines.append("")
        filename = f"@{self._id[-1]}+nve.dump"
        if configuration.periodicity == 0:
            lines.append(
                f"write_dump         all custom  {filename} id xu yu zu vx vy vz"
                " modify flush yes sort id"
            )
        else:
            lines.append(
                f"write_dump         all custom  {filename} id xsu ysu zsu vx vy vz"
                " modify flush yes sort id"
            )
        lines.append("")

        for i in range(1, ncomputes + 1):
            lines.append(f"uncompute           {i}")
        for i in range(1, ndumps + 1):
            lines.append(f"undump              {i}")
        for i in range(1, nfixes + 1):
            lines.append(f"unfix               {i}")
        if P["heat flux"] != "never":
            lines.extend(self.heat_flux_cleanup(P, ff_form))
        lines.extend(user_extras["cleanup"])
        lines.append("")

        return {
            "script": lines,
            "postscript": None,
            "use python": False,
        }

    def _extras_parameters(self):
        """The extra commands and properties as given: a string and a list of
        dicts (``values_to_dict`` would give the list as text)."""
        P = self.parameters
        result = {}
        for key in ("extra commands", "extra properties"):
            value = P[key].value if key in P else None
            if isinstance(value, str) and key == "extra properties":
                value = json.loads(value) if value.strip() != "" else []
            result[key] = value
        return result

    def check_flowchart(self, preceding):
        """Check the extra commands and properties before the job runs (see
        seamm.builder.validate), so that a mistake stops the job at once."""
        text = lammps_step.extras.check(self._extras_parameters())
        if text == "":
            return []
        return [f"'{self.title}': {problem}" for problem in text.splitlines()]

    def extra_definitions(self):
        """The user's extra properties, checked, or [] if there are none or they
        are wrong (the input refuses wrong ones before the run)."""
        P = self._extras_parameters()
        try:
            return lammps_step.extras.check_properties(P["extra properties"])
        except lammps_step.extras.ExtrasError:
            return []

    def _mean_volume(self, data, configuration):
        """The mean volume of the cell, or None for a non-periodic system."""
        if "V" in data:
            return Q_(data["V"], "Å^3")
        if configuration is not None and configuration.periodicity == 3:
            return Q_(configuration.cell.volume, "Å^3")
        return None

    def _analyze_extra_properties(self, properties, data, configuration):
        """Treat, print and store the user's extra properties.

        Parameters
        ----------
        properties : dict
            The analysis of the state file, per column.
        data : dict
            The results of the step; ``extra properties`` is added.
        configuration : molsystem._Configuration
            Where to store the properties asked for.

        Returns
        -------
        str
            The table of the extra properties, or "" if there are none.
        """
        definitions = self.extra_definitions()
        if len(definitions) == 0:
            return ""
        table = {
            "Property": [],
            "Treatment": [],
            "Value": [],
            " ": [],
            "StdErr": [],
            "Units": [],
            "Note": [],
        }
        results = {}
        for d in definitions:
            name = d["name"]
            entry = {
                "expression": d["expression"],
                "treatment": d["treatment"],
                "units": d["display_units"],
            }
            results[name] = entry
            r = properties.get(d["column"])
            if r is None or "mean" not in r:
                entry["error"] = "the property did not equilibrate"
                for key, value in (
                    ("Property", name),
                    ("Treatment", d["treatment"]),
                    ("Value", ""),
                    (" ", ""),
                    ("StdErr", ""),
                    ("Units", d["display_units"]),
                    ("Note", entry["error"]),
                ):
                    table[key].append(value)
                continue
            entry["mean"] = r["mean"]
            entry["stderr"] = r["stderr"]
            value, stderr, units = r["mean"], r["stderr"], d["display_units"]
            note = ""
            dt = Q_(r["timestep"], "fs")
            # The production part of the series, kept by analyze_trajectory
            values = getattr(self, "_extra_series", {}).get(d["column"])
            try:
                if d["treatment"] == "autocorrelation":
                    ac = lammps_step.extras.autocorrelation(values, dt)
                    tau = ac["tau"].m_as("ps")
                    entry["correlation time"] = tau
                    entry["correlation time,units"] = "ps"
                    entry["acf"] = ac["acf"]
                    note = f"correlation time {tau:.3g} ps"
                    self._extra_acf_graph(d, ac)
                elif d["treatment"] == "Green-Kubo integral":
                    temperature = Q_(data["T"], "K") if "T" in data else None
                    gk = lammps_step.extras.green_kubo(
                        values,
                        dt,
                        d,
                        unit=d["display_units"],
                        volume=self._mean_volume(data, configuration),
                        temperature=temperature,
                    )
                    value = gk["value"].magnitude
                    stderr = gk["stderr"].magnitude
                    units = gk["units"]
                    entry["integral"] = value
                    entry["integral,stderr"] = stderr
                    entry["integral,units"] = units
                    entry["integral,error from"] = gk["error_from"]
                    entry["tau"] = gk["tau"]
                    entry["tau,units"] = "ps"
                    note = (
                        "Green-Kubo, error from 4 blocks"
                        if gk["error_from"] == "blocks"
                        else "Green-Kubo, error from the fit (too small: the run is "
                        "too short for blocks)"
                    )
                    self._extra_gk_graph(d, gk)
            except Exception as e:
                logger.warning(f"The {d['treatment']} of {name} failed: {e}")
                entry["error"] = f"the {d['treatment']} failed: {e}"
                note = entry["error"]

            for key, item in (
                ("Property", name),
                ("Treatment", d["treatment"]),
                ("Value", f"{value:.4g}"),
                (" ", "±"),
                ("StdErr", f"{stderr:.2g}"),
                ("Units", units),
                ("Note", note),
            ):
                table[key].append(item)

            if d["store"] and configuration is not None and "error" not in entry:
                self._store_extra_property(configuration, d, value, stderr, units)
                entry["stored as"] = f"{name}#LAMMPS#{self.model}"

        data["extra properties"] = results

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        return "\n" + "Extra Properties".center(length) + "\n" + tmp + "\n"

    def _store_extra_property(self, configuration, definition, value, stderr, units):
        """Store an extra property, and its standard error, in the configuration,
        defining them if needed."""
        properties = configuration.properties
        for suffix, item in (("", value), (",stderr", stderr)):
            name = f"{definition['name']}{suffix}#LAMMPS#{self.model}"
            if not properties.exists(name):
                what = "the standard error of " if suffix else ""
                properties.add(
                    name,
                    "float",
                    units=units,
                    description=(
                        f"{what}{definition['name']} = {definition['expression']} "
                        f"({definition['treatment']}, {self.model})"
                    ),
                    noerror=True,
                )
            properties.put(name, item)

    def _write_graph(self, figure, name):
        """Write a graph, and the other formats asked for."""
        path = Path(self.directory) / f"{name}.graph"
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.write_file(path)

    def _extra_acf_graph(self, definition, ac):
        """The graph of the autocorrelation function of an extra property."""
        dt = ac["dt"].m_as("ps")
        n = len(ac["acf"])
        ts = [i * dt for i in range(n)]
        figure = self.create_figure(
            module_path=("lammps_step", "seamm"),
            template="line.graph_template",
            title=f"Autocorrelation of {definition['name']}",
        )
        plot = figure.add_plot("acf")
        x_axis = plot.add_axis("x", label="Time (ps)")
        y_axis = plot.add_axis("y", label="ACF", anchor=x_axis)
        x_axis.anchor = y_axis
        upper = [a + 2 * e for a, e in zip(ac["acf"], ac["stderr"])]
        lower = [a - 2 * e for a, e in zip(ac["acf"], ac["stderr"])]
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name="95% confidence",
            x=ts + ts[::-1],
            xlabel="t",
            xunits="ps",
            y=upper + lower[::-1],
            ylabel="95%",
            yunits="",
            color="rgba(211,211,211,0.5)",
            fill="toself",
        )
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name="ACF",
            x0=0,
            dx=dt,
            xlabel="t",
            xunits="ps",
            y=list(ac["acf"]),
            ylabel="ACF",
            yunits="",
            color="red",
        )
        figure.grid_plots("acf")
        self._write_graph(figure, f"{definition['name']}_acf")

    def _extra_gk_graph(self, definition, gk):
        """The graph of the Green-Kubo running integral of an extra property."""
        figure = self.create_figure(
            module_path=("lammps_step", "seamm"),
            template="line.graph_template",
            title=f"Green-Kubo integral of {definition['name']}",
        )
        plot = figure.add_plot("gk")
        x_axis = plot.add_axis("x", label="Time (ps)")
        y_axis = plot.add_axis("y", label=f"Integral ({gk['units']})", anchor=x_axis)
        x_axis.anchor = y_axis
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name="integral",
            x=gk["ts"],
            xlabel="t",
            xunits="ps",
            y=gk["integral"],
            ylabel="integral",
            yunits=gk["units"],
            color="#4dbd74",
        )
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name="fit",
            x=gk["fit_ts"],
            xlabel="t",
            xunits="ps",
            y=gk["fit"],
            ylabel="fit",
            yunits=gk["units"],
            color="black",
        )
        figure.grid_plots("gk")
        self._write_graph(figure, f"{definition['name']}_gk")

    def extra_input(self, P):
        """The user's extra commands and properties for this sub-step, checked
        (see :mod:`lammps_step.extras`)."""
        try:
            return lammps_step.extras.extra_input(self._extras_parameters())
        except lammps_step.extras.ExtrasError as e:
            raise RuntimeError(
                f"{self.header}: the extra commands or properties are wrong:\n{e}"
            ) from None

    def heat_flux_input(self, P, ff_form=None):
        """The computes and variables giving the heat flux, in W/m^2, as the
        variables Jx, Jy and Jz.

        The flux is J = sum_i e_i v_i - sum_i S_i.v_i over the atoms, with e_i
        the per-atom energy and S_i the per-atom stress. The bonded terms and the
        constraint forces (SHAKE) need the centroid form of the stress, which
        compute centroid/stress/atom gives (Surblys et al., J. Appl. Phys. 2021);
        it does not handle kspace, so the pair and kspace terms are a second
        compute, and since both computes include the convective term sum_i e_i v_i,
        one copy is removed. Class2 forcefields, or "use centroid stress" off, use
        the ordinary per-atom stress of all the terms.

        Parameters
        ----------
        P : dict
            The step's parameter values ("use centroid stress").
        ff_form : str, optional
            The form of the forcefield. Default: the forcefield's.

        Returns
        -------
        str
            The LAMMPS input, to be undone by :meth:`heat_flux_cleanup`.
        """
        if ff_form is None:
            ff_form = self.parent.ff_form()
        if lammps_step.get_lammps_unit_system() == "metal":
            factor = Q_("eV/Å^2/ps")
        else:
            factor = Q_("kcal/Å^2/fs/mol") / Q_("kcal/mol") * Q_("kcal/mol").to("kJ")
        factor = factor.m_as("W/m^2")
        if ff_form == "class2" or not P["use centroid stress"]:
            # Centroid/stress/atom does not handle class2 ff ... cross-terms?
            return f"""
compute             KE all ke/atom
compute             PE all pe/atom

compute             S_p all stress/atom NULL virial
compute             flux_p all heat/flux KE PE S_p

#          Conversion from kcal/Å^2/fs/mol to W/m^2

variable            factor equal {factor}
variable            Jx equal v_factor*c_flux_p[1]/vol
variable            Jy equal v_factor*c_flux_p[2]/vol
variable            Jz equal v_factor*c_flux_p[3]/vol
"""
        return f"""
compute             KE all ke/atom
compute             PE all pe/atom

#          centroid doesn't work with kspace, so split into pair and non-pair parts.
#          The constraint forces (fix) need the centroid form.

compute             S_p all stress/atom NULL pair kspace
compute             S_b all centroid/stress/atom NULL bond angle dihedral improper fix
compute             flux_p all heat/flux KE PE S_p
compute             flux_b all heat/flux KE PE S_b

#          Conversion from kcal/Å^2/fs/mol to W/m^2

variable            factor equal {factor}
#          Both computes include the convective term, so remove one copy

variable            Jx equal v_factor*(c_flux_p[1]+c_flux_b[1]-c_flux_b[4])/vol
variable            Jy equal v_factor*(c_flux_p[2]+c_flux_b[2]-c_flux_b[5])/vol
variable            Jz equal v_factor*(c_flux_p[3]+c_flux_b[3]-c_flux_b[6])/vol
"""

    def heat_flux_cleanup(self, P, ff_form=None):
        """The lines removing what :meth:`heat_flux_input` defined."""
        if ff_form is None:
            ff_form = self.parent.ff_form()
        lines = []
        if ff_form != "class2" and P["use centroid stress"]:
            lines.append("uncompute           flux_b")
            lines.append("uncompute           S_b")
        lines.append("uncompute           flux_p")
        lines.append("uncompute           S_p")
        lines.append("uncompute           PE")
        lines.append("uncompute           KE")
        lines.append("variable            factor delete")
        lines.append("variable            Jx delete")
        lines.append("variable            Jy delete")
        lines.append("variable            Jz delete")
        return lines

    def md_provides(self):
        """What this run will write to md.h5md, from its parameters.

        Returns
        -------
        seamm_md_analysis.Provides
        """
        return lammps_step.md_provides.md_provides(self)

    def request_sampling(self, need):
        """Turn on the sampling an analysis needs, if possible.

        Parameters
        ----------
        need : seamm_md_analysis.Need

        Returns
        -------
        str or None
            A description of the change, or None if it could not be made.
        """
        return lammps_step.md_provides.request_sampling(self, need)

    def _drop_periodic_sampling(self, P, configuration):
        """Turn off sampling that needs a periodic cell, for a non-periodic system.

        The heat flux and the pressure tensor are divided by the volume of the
        cell, so they are meaningless without one. The GUI cannot hide these
        choices because it does not know the structure, so they are turned off
        here with a warning.

        Parameters
        ----------
        P : dict
            The parameters of the step, changed in place.
        configuration : molsystem._Configuration
            The configuration that will be simulated.

        Returns
        -------
        str
            A warning for the description, or "" if nothing was changed.
        """
        if configuration.periodicity != 0:
            return ""
        dropped = []
        for key in ("heat flux", "pressure tensor"):
            if P.get(key, "never") != "never":
                P[key] = "never"
                dropped.append(key)
        if len(dropped) == 0:
            return ""
        return (
            "Warning: the system is not periodic, so the "
            + " and the ".join(dropped)
            + " will not be sampled. They need the volume of a periodic cell."
        )

    def timestep(self, value):
        """Get the timestep in the correct units.

        This handles the 'normal', 'accurate' and 'coarse' values,
        which depend on the mass in an empirical fashion.

        Parameters
        ----------
        value : str or Pint quantity
            The desired timestep, which may be a Pint quantity with units
            or one of 'notmal', 'accurate but slow', or 'coarse but fast'

        Returns
        -------
        timestep : float
            The magnitude of the time step in the appropriate LAMMPS units
        """
        masses = self.parent._data["masses"]
        min_mass = min(masses)

        if self.parent.ff_form() == "reaxff":
            # ReaxFF needs a smaller timestep for the QEq part
            factor = 0.5
        else:
            # These are based on masses as a proxy for vibrational frequencies
            if min_mass < 10:
                factor = 1
            elif min_mass < 50:
                factor = 2
            else:
                factor = 4

        if value == "normal":
            timestep = 1.0 * factor
            value = Q_(timestep, ureg.fs)
        elif value == "accurate but slow":
            timestep = 0.5 * factor
            value = Q_(timestep, ureg.fs)
        elif value == "coarse but fast":
            timestep = 2.0 * factor
            value = Q_(timestep, ureg.fs)
        else:
            value = Q_(value)

        timestep = lammps_step.to_lammps_units(value, quantity="time")

        return (timestep, value)

    def trajectory_input(self, P, timestep, nsteps, ncomputes, ndumps, nfixes):
        """Create the part of the input handling the trajectories.

        Parameters
        ----------
        P : dict
            The dictionary of options
        timestep : int
            The timestep in LAMMPS units
        nsteps : int
            Total number of steps in the run
        ncomputes : int
            The counter for the computes
        ndumps : int
            The counter for the dumps
        nfixes : int
            The counter for the fixes

        Returns
        -------
        [str]
            The input lines as array
        int
            The counter for the computes
        int
            The counter for the dumps
        int
            The counter for the fixes
        """
        lines = []
        sampling = P["trajectory"]
        if sampling != "never":
            if "interval" in sampling:
                t_s = lammps_step.to_lammps_units(P["trajectory rate"], quantity="time")
                n = max(1, round(t_s / timestep))
            else:
                n = max(1, nsteps // P["trajectory number of samples"])
            ndumps += 1
            filename = f"@{self._id[-1]}+trajectory.dump_trj.gz"
            if "charges" in self.parent.eex:
                line = (
                    f"dump                {ndumps} all custom {n} {filename} id element"
                    " q xu yu zu"
                )
            else:
                line = (
                    f"dump                {ndumps} all custom {n} {filename} id element"
                    " xu yu zu"
                )
            if P["trajectory forces"]:
                line += " fx fy fz"
            if P["trajectory velocities"]:
                line += " vx vy vz"
            lines.append("\n")
            lines.append(line)
            tmp = " ".join(self.parent.eex["elements"])
            lines.append(
                f"dump_modify         {ndumps} sort id format float %.5f element {tmp}"
            )

            if False:
                # Temporarily put in dump to extxyz
                ndumps += 1
                filename = f"@{self._id[-1]}+trajectory.extxyz.gz"
                line = f"dump                {ndumps} all extxyz {n} {filename}"
                lines.append(line)
                element_string = " ".join(self.parent.eex["elements"])
                lines.append(f"dump_modify         {ndumps} element {element_string}")

            # The energies, stress, etc
            properties = "v_time v_temp v_press v_etotal v_ke v_pe v_emol v_epair"
            title2 = "tstep t T P Etot Eke Epe Emol Epair"
            if self.parent.have_dreiding_hbonds:
                properties += " v_N_hbond v_E_hbond"
                title2 += " N_hbond E_hbond"

            filename = f"@{self._id[-1]}+trajectory.trj"
            nfixes += 1
            dt = (n * P["timestep"]).to_compact()
            text = json.dumps(
                {
                    "code": "LAMMPS",
                    "type": "state for trajectory",
                    "dt": dt.magnitude,
                    "tunits": str(dt.u),
                    "nsteps": nsteps // n,
                },
                separators=(",", ":"),
            )
            title1 = "!MolSSI trajectory 2.0 " + text
            _, configuration = self.get_system_configuration()
            if configuration.periodicity == 3:
                title2 += " Sxx Syy Szz Syz Sxz Sxy"
                properties += " v_sxx v_syy v_szz v_syz v_sxz v_sxy"
            else:
                # No pressure for a non-periodic system: it is meaningless, and
                # the stress variables do not exist
                properties = properties.replace(" v_press", "")
                title2 = title2.replace(" P ", " ")
            lines.append(
                "\n"
                f"fix                 {nfixes} all ave/time {n} 1 {n} "
                f"                        {properties}    &\n"
                f"                        title1 '{title1}' &\n"
                f"                        title2 '{title2}' &\n"
                f"                        file {filename}"
            )
        sampling = P["atomic positions"]
        if sampling != "never":
            if "interval" in sampling:
                t_s = lammps_step.to_lammps_units(
                    P["atomic positions rate"], quantity="time"
                )
                n = max(1, round(t_s / timestep))
            else:
                n = max(1, nsteps // P["atomic positions number of samples"])
            ndumps += 1
            filename = f"@{self._id[-1]}+atomic_positions.dump_trj"
            lines.append(
                "\n"
                f"dump                {ndumps} all custom {n} {filename} id xu yu zu\n"
                f"dump_modify         {ndumps} sort id time yes units yes"
            )
        sampling = P["com positions"]
        if sampling != "never":
            if "interval" in sampling:
                t_s = lammps_step.to_lammps_units(
                    P["com positions rate"], quantity="time"
                )
                n = max(1, round(t_s / timestep))
            else:
                n = max(1, nsteps // P["com positions number of samples"])

            filename = f"@{self._id[-1]}+com_positions.trj"
            ncomputes += 1
            c1 = ncomputes
            ncomputes += 1
            c2 = ncomputes
            nfixes += 1
            dt = (n * P["timestep"]).to_compact()
            text = json.dumps(
                {
                    "code": "LAMMPS",
                    "type": "com positions",
                    "dt": dt.magnitude,
                    "tunits": str(dt.u),
                    "nsteps": nsteps // n,
                },
                separators=(",", ":"),
            )
            title1 = "!MolSSI vector_trajectory 2.0 " + text
            title2 = "! timestep n_molecules"
            title3 = "molecule com_x com_y com_z"
            lines.append(
                "\n"
                f"compute             {c1} all chunk/atom molecule\n"
                f"compute             {c2} all com/chunk {c1}\n"
                f"fix                 {nfixes} all ave/time {n} 1 {n} c_{c2}[*] &\n"
                f"                        title1 '{title1}' &\n"
                f"                        title2 '{title2}' &\n"
                f"                        title3 '{title3}' &\n"
                f"                        file {filename} mode vector"
            )
        sampling = P["atomic velocities"]
        if sampling != "never":
            if "interval" in sampling:
                t_s = lammps_step.to_lammps_units(
                    P["atomic velocities rate"], quantity="time"
                )
                n = max(1, round(t_s / timestep))
            else:
                n = max(1, nsteps // P["atomic velocities number of samples"])
            ndumps += 1
            filename = f"@{self._id[-1]}+atomic_velocities.dump_trj"
            lines.append(
                "\n"
                f"dump                {ndumps} all custom {n} {filename} id vx vy vz\n"
                f"dump_modify         {ndumps} sort id time yes units yes"
            )
        sampling = P["com velocities"]
        if sampling != "never":
            if "interval" in sampling:
                t_s = lammps_step.to_lammps_units(
                    P["com velocities rate"], quantity="time"
                )
                n = max(1, round(t_s / timestep))
            else:
                n = max(1, nsteps // P["com velocities number of samples"])
            filename = f"@{self._id[-1]}+com_velocities.trj"
            ncomputes += 1
            c1 = ncomputes
            ncomputes += 1
            c2 = ncomputes
            nfixes += 1
            dt = (n * P["timestep"]).to_compact()
            text = json.dumps(
                {
                    "code": "LAMMPS",
                    "type": "com velocities",
                    "dt": dt.magnitude,
                    "tunits": str(dt.u),
                    "nsteps": nsteps // n,
                },
                separators=(",", ":"),
            )
            title1 = "!MolSSI vector_trajectory 2.0 " + text
            title2 = "! timestep n_molecules"
            title3 = "molecule com_vx com_vy com_vz"
            lines.append(
                "\n"
                f"compute             {c1} all chunk/atom molecule\n"
                f"compute             {c2} all vcm/chunk {c1}\n"
                f"fix                 {nfixes} all ave/time {n} 1 {n} c_{c2}[*] &\n"
                f"                        title1 '{title1}' &\n"
                f"                        title2 '{title2}' &\n"
                f"                        title3 '{title3}' &\n"
                f"                        file {filename} mode vector"
            )
        sampling = P["heat flux"]
        if sampling != "never":
            if "interval" in sampling:
                t_s = lammps_step.to_lammps_units(P["heat flux rate"], quantity="time")
                n = max(1, round(t_s / timestep))
            else:
                n = max(1, nsteps // P["heat flux number of samples"])
            filename = f"@{self._id[-1]}+heat_flux.trj"
            nfixes += 1
            dt = (n * P["timestep"]).to_compact()
            text = json.dumps(
                {
                    "code": "LAMMPS",
                    "type": "heat flux",
                    "dt": dt.magnitude,
                    "tunits": str(dt.u),
                    "nsteps": nsteps // n,
                },
                separators=(",", ":"),
            )
            title1 = "!MolSSI trajectory 2.0 " + text
            title2 = "Jx Jy Jz"
            lines.append(
                "\n"
                f"fix                 {nfixes} all ave/time {n} 1 {n} "
                "v_Jx v_Jy v_Jz &\n"
                f"                        title1 '{title1}' &\n"
                f"                        title2 '{title2}' &\n"
                f"                        file {filename}"
            )
        sampling = P["pressure tensor"]
        if sampling != "never":
            if "interval" in sampling:
                t_s = lammps_step.to_lammps_units(
                    P["pressure tensor rate"], quantity="time"
                )
                n = max(1, round(t_s / timestep))
            else:
                n = max(1, nsteps // P["pressure tensor number of samples"])
            filename = f"@{self._id[-1]}+pressure_tensor.trj"
            nfixes += 1
            dt = (n * P["timestep"]).to_compact()
            text = json.dumps(
                {
                    "code": "LAMMPS",
                    "type": "pressure tensor",
                    "dt": dt.magnitude,
                    "tunits": str(dt.u),
                    "nsteps": nsteps // n,
                },
                separators=(",", ":"),
            )
            title1 = "!MolSSI trajectory 2.0 " + text
            title2 = "Pxx Pyy Pzz Pxy Pxz Pyz"
            lines.append(
                "\n"
                f"fix                 {nfixes} all ave/time {n} 1 {n} "
                "v_pxx v_pyy v_pzz v_pxy v_pxz v_pyz &\n"
                f"                        title1 '{title1}' &\n"
                f"                        title2 '{title2}' &\n"
                f"                        file {filename}"
            )

        return lines, ncomputes, ndumps, nfixes

    def _save_trajectory(self, P):
        """Save the trajectory to configurations.

        Parameters
        ----------
        P : dict(str, value)
            The control parameters for this step

        Returns
        -------
        str
            A text string for printing, containing warnings, etc.
        """
        text = ""
        dump = Path(self.directory) / "trajectory.dump_trj.gz"
        trj = Path(self.directory) / "trajectory.trj"

        if not dump.exists() and not trj.exists():
            return (
                "Neither the coordinate or energy part of the trajectory exists, so "
                "cannot save it."
            )
        if trj.exists():
            trj_data = trj.read_text().splitlines()
            header = trj_data[0]
            if not header.startswith("!MolSSI trajectory"):
                trj_data = None
                text = (
                    "The file for the energies etc. for the trajectory exists "
                    "but is not the right type of file, so the energies, etc. "
                    "will not be saved with the configuration. "
                )
            else:
                state_variables = trj_data[1].split()
                trj_data = [i.split() for i in trj_data[2:]]
        else:
            trj_data = None
            text = (
                "The energies etc. for the trajectory are missing so they will "
                "not be saved with the configuration. "
            )

        # Get the elements from the starting configuration
        initial_system, initial_configuration = self.get_system_configuration()
        symbols = initial_configuration.atoms.symbols

        # Make a new system
        name = P["trajectory system name"]
        system_db = initial_system.system_db
        if name == "current":
            system = system_db.system
        else:
            if system_db.system_exists(name):
                system = system_db.get_system(name)
            else:
                system = system_db.create_system(name=name)
        if P["make current"]:
            system_db.system = system

        # Read the trajectory dump file and process
        _, timestep = self.timestep(P["timestep"])
        frame = 0
        start_line = 0
        lines = []

        with gzip.open(dump, "rt") as fdin:
            for line in fdin:
                line = line.strip()
                if line.startswith("ITEM: TIMESTEP"):
                    if len(lines) > 0:
                        frame += 1
                        results = self._parse_dump_frame(lines, start_line)
                        results["symbols"] = symbols
                        results.update(
                            {k: v for k, v in zip(state_variables, trj_data[frame - 1])}
                        )
                        # Make a new configuration
                        time = timestep * results["timestep"]
                        name = f"{time:.1f~P}"
                        configuration = system.copy_configuration(
                            configuration=initial_configuration, name=name
                        )
                        self._to_configuration(configuration, results)
                    start_line += len(lines)
                    lines = []
                lines.append(line)
            # Process last frame
            if len(lines) > 0:
                frame += 1
                results = self._parse_dump_frame(lines, start_line)
                results["symbols"] = symbols
                results.update(
                    {k: v for k, v in zip(state_variables, trj_data[frame - 1])}
                )
                # Make a new configuration
                time = timestep * results["timestep"]
                name = f"{time:.1f~P}"
                configuration = system.copy_configuration(
                    configuration=initial_configuration, name=name
                )
                self._to_configuration(configuration, results)
        text += f"Created {frame} configurations of system '{system.name}'"
        if P["make current"]:
            text += ", which was made the current system,"
        text += " from the trajectory."
        return text

    def _check_property(self, properties, _property, **kwargs):
        """Check if a property exists after substitution and create if necessary.

        Parameters
        ----------
        properties : molsystem._Properties
            The properties object

        _property : str
            The name of the property

        **kwargs : {str: str}
            Optional variables for substituting in the property name.
        """
        expanded_property = _property.format_map(kwargs)
        if not properties.exists(expanded_property):
            # Get the general property's info to create the model property.
            _type, units, description = properties.metadata(_property)
            properties.add(
                expanded_property,
                _type=_type,
                units=units,
                description=description.format_map(kwargs),
            )
        return expanded_property

    def _to_configuration(self, configuration, data):
        """Update a configuration with the data from a step in a trajectory.

        Parameters
        ----------
        configuration : molsystem._Configuration
            The configuration to update
        data : dict(str, any)
            Dictionary containing coordinates and other atom properties
        """
        if configuration.periodicity == 3:
            # Set the cell
            if "cell" in data:
                configuration.cell.parameters = data["cell"]

            # Stress
            if "Sxx" in data:
                sfact = from_lammps_units(1, "atm").magnitude
                stress = [
                    sfact * float(data["Sxx"]),
                    sfact * float(data["Syy"]),
                    sfact * float(data["Szz"]),
                    sfact * float(data["Syz"]),
                    sfact * float(data["Sxz"]),
                    sfact * float(data["Sxy"]),
                ]
                _property = self._check_property(
                    configuration.properties, "stress#LAMMPS#{model}", model=self.model
                )
                configuration.properties.put(_property, stress)

        # Coordinates
        if "abc" in data:
            configuration.atoms.set_coordinates(data["abc"], fractionals=True)
        elif "xyz" in data:
            configuration.atoms.set_coordinates(data["xyz"], fractionals=False)

        if "velocities" in data:
            # LAMMPS only has Cartesian velocities. Already in Å/fs
            tmp = np.array(data["velocities"])
            configuration.atoms.set_velocities(tmp, fractionals=False)

        if "gradients" in data:
            # LAMMPS only has Cartesian forces, in kcal/mol/Å
            factor = Q_("kcal/mol/Å").m_as("kJ/mol/Å")
            tmp = factor * np.array(data["gradients"])
            configuration.atoms.set_gradients(tmp, fractionals=False)

        # Energy
        if "Epe" in data:
            _property = self._check_property(
                configuration.properties,
                "potential energy#LAMMPS#{model}",
                model=self.model,
            )
            tmp = from_lammps_units(float(data["Epe"]), "kcal/mol").magnitude
            configuration.properties.put(_property, tmp)

            # Calculate the energy of formation if we can
            tmp_data = {"Epe": tmp}
            self.calculate_enthalpy_of_formation(tmp_data)
            if "DfE0" in tmp_data:
                _property = self._check_property(
                    configuration.properties,
                    "DfE0#LAMMPS#{model}",
                    model=self.model,
                )
                configuration.properties.put(_property, tmp_data["DfE0"])

    def _to_extxyz(self, data, timestep=None):
        """Create the text of ASE extxyz format for a structure

        Parameters
        ----------
        data : dict(str, any)
            Dictionary containing symbols, coordinates, and other atom properties

        timestep : pint.units
            The timestep.

        Note
        ----
        Example:

        8
        Properties=species:S:1:pos:R:3:REF_forces:R:3 REF_energy=-9289.028102188924 pbc="F F F"
        C       -0.97142774       0.56951904       1.10325229       2.77158783      -3.19505603       0.10221193
        O       -0.27179274       1.07438707      -0.03295067      -1.40444169       0.58287428       2.02393798
        C        0.44642827       0.04566107      -0.56782067      -2.59873178       1.73475706       0.42253222
        O        0.22332226      -1.04711390       0.04433633      -1.29638003      -2.74941166       4.75917280
        C       -0.54656476      -0.74672896       1.23491836      -0.44334119       2.10839425      -2.18446253
        O        0.98560923       0.09795906      -1.57642770       3.16925789       1.02455543      -4.39760192
        H       -0.66579676      -1.59970891       1.88536036      -0.27973052       0.31498569       0.12677338
        H       -1.44834375       1.17851210       1.87373829       0.08523677       0.17712004      -0.85144738
        """  # noqa: E501
        lines = []
        lines.append(f"{data['n_atoms']}")

        # The lattice/properties line
        if "lattice" in data:
            # already in Å
            line = 'Lattice="'
            line += " ".join([f"{v:.5f}" for v in data["lattice"]])
            line += '"'
        line += " Properties=species:S:1:pos:R:3"
        if "gradients" in data:
            # currently in kcal/mol/Å
            line += ":REF_forces:R:3"
            ffact = -Q_(1, "kcal/mol/Å").m_as("eV/Å")
        if "velocities" in data:
            # currentl in Å/fs
            line += ":vel:R:3"
            vfact = Q_(1, "Å/fs").m_as("eV^0.5/amu^0.5")

        # Add in other properties, like the energy ... currently just energy and stress
        if "Epe" in data:
            # Calculate the energy of formation if we can
            tmp_data = {
                "Epe": from_lammps_units(float(data["Epe"]), "kcal/mol").magnitude
            }
            self.calculate_enthalpy_of_formation(tmp_data)
            if "DfE0" in tmp_data:
                value = Q_(tmp_data["DfE0"], "kcal/mol").m_as("eV")
                line += f" REF_energy={value:.4f}"
            else:
                # The energy is the potential, not total energy which includes kinetic
                efact = from_lammps_units(1, "eV").magnitude
                line += f" REF_energy={efact*float(data['Epe']):.4f}"
        if "Sxx" in data:
            sfact = from_lammps_units(1, "eV/Å^3").magnitude
            Sxx = f"{sfact*float(data['Sxx']):.7f}"
            Syy = f"{sfact*float(data['Syy']):.7f}"
            Szz = f"{sfact*float(data['Szz']):.7f}"
            Syz = f"{sfact*float(data['Syz']):.7f}"
            Sxz = f"{sfact*float(data['Sxz']):.7f}"
            Sxy = f"{sfact*float(data['Sxy']):.7f}"
            line += f' REF_stress="{Sxx} {Syy} {Szz} {Syz} {Sxz} {Sxy}" pbc="T T T"'
        else:
            line += ' pbc="F F F"'

        # The time
        if timestep is not None and "timestep" in data:
            time = timestep * data["timestep"]
            line += f' time="{time:.1f~P}"'

        # Add the model, jobserver, and job id if available
        line += f' model="{self.model}"'
        if "SEAMM_JOBSERVER" in os.environ:
            tmp = os.environ["SEAMM_JOBSERVER"]
            line += f' jobserver="{tmp}"'
        if "SEAMM_JOB_ID" in os.environ:
            tmp = os.environ["SEAMM_JOB_ID"]
            line += f' job_id="{tmp}"'

        lines.append(line)

        # And the coordinate lines
        for i, symbol in enumerate(data["symbols"]):
            line = f"{symbol:<2}"
            for val in data["xyz"][i]:
                line += f" {val:14.8f}"
            if "gradients" in data:
                for val in data["gradients"][i]:
                    line += f" {ffact*val:14.8f}"
            if "velocities" in data:
                for val in data["velocities"][i]:
                    line += f" {vfact*val:14.8f}"
            lines.append(line)

        # Add empty line so text ends with newline
        lines.append("")

        return "\n".join(lines)
