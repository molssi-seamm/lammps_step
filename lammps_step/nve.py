# -*- coding: utf-8 -*-

"""NVE (microcanonical) dynamics in LAMMPS"""

import gzip
import json
from pathlib import Path
import traceback

import lammps_step
from lammps_step import LAMMPS, from_lammps_units
import logging
from molsystem import Cell
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("lammps")


class NVE(lammps_step.Energy):
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
        super().analyze(data=data, properties=properties, table=table, output=output)

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Get the initial structure
        _, initial_configuration = self.get_system_configuration()

        # Handle the new structure as needed
        system, configuration = self.get_system_configuration(P, model=self.model)

        # Read the dump file and get the structure
        dump_file = Path(self.directory) / f"{self.calculation}.dump"
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

        # Import or convert the trajectory.
        if P["trajectory export"]:
            text = self._export_trajectory(P)
            printer.normal(__(text, indent=self.indent + 4 * " "))
            printer.normal("")

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
        _, configuration = self.get_system_configuration()
        symbols = configuration.atoms.symbols

        filename = P["trajectory extxyz filename"].strip()
        if filename.startswith("/"):
            path = Path(self.flowchart.root_directory) / filename[1:]
        else:
            path = Path(filename)

        nskip = P["trajectory extxyz skip frames"]
        if P["trajectory extxyz append"]:
            mode = "at"
            text += "Appended"
        else:
            mode = "wt"
            text += "Wrote"
        frame = 0
        start_line = 0
        lines = []
        with gzip.open(dump, "rt") as fdin, gzip.open(path, mode) as fdout:
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
                            extxyz = self._to_extxyz(results)
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
                    extxyz = self._to_extxyz(results)
                    fdout.write(extxyz)
        if nskip > 0:
            text += (
                f" {frame - nskip} frames to the trajectory in .extxyz format, "
                f"skipping the first {nskip} frames."
            )
        else:
            text += f" {frame - nskip} frames to the trajectory in .extxyz format."
        return text

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
                    if len(section.split()) == 8:
                        line, line_no = next_line(it, line_no)
                        xlo_bound, xhi_bound, xy = line.split()
                        line, line_no = next_line(it, line_no)
                        ylo_bound, yhi_bound, xz = line.split()
                        line, line_no = next_line(it, line_no)
                        zlo, zhi, yz = line.split()

                        xlo_bound = float(xlo_bound)
                        xhi_bound = float(xhi_bound)
                        ylo_bound = float(ylo_bound)
                        yhi_bound = float(yhi_bound)
                        zlo = float(zlo)
                        zhi = float(zhi)
                        xy = float(xy)
                        xz = float(xz)
                        yz = float(yz)

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

                        xlo = float(xlo)
                        xhi = float(xhi)
                        ylo = float(ylo)
                        yhi = float(yhi)
                        zlo = float(zlo)
                        zhi = float(zhi)

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
                        ff = from_lammps_units(1, "kcal/mol/Å").magnitude

                    # Velocities
                    vi = None
                    if "vx" in keys and "vy" in keys and "vz" in keys:
                        vi = (keys.index("vx"), keys.index("vy"), keys.index("vz"))
                        fv = from_lammps_units(1, "m/s").magnitude

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
                        results["forces"] = f
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

        if extras is not None and "nsteps" in extras:
            nsteps = extras["nsteps"]
        else:
            time = P["time"].to("fs").magnitude
            nsteps = round(time / timestep)

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        self.description.append(__(self.description_text(), **PP, indent=7 * " "))

        # time = lammps_step.to_lammps_units(P['time'], quantity='time')
        # nsteps = round(time / timestep)

        thermo_properties = (
            "time temp press etotal ke pe ebond "
            "eangle edihed eimp evdwl etail ecoul elong"
        )
        properties = "v_time v_temp v_press v_etotal v_ke v_pe v_emol v_epair"
        title2 = "tstep t T P Etot Eke Epe Emol Epair"
        if configuration.periodicity == 3:
            properties += " v_sxx v_syy v_szz v_sxy v_sxz v_syz"
            title2 += " Sxx Syy Szz Sxy Sxz Syz"
        if self.parent.have_dreiding_hbonds:
            thermo_properties += " v_N_hbond v_E_hbond"
            properties += " v_N_hbond v_E_hbond"
            title2 += " N_hbond E_hbond"

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

        # For the heat flux, if requested, we need extra input
        if P["heat flux"] != "never":
            # Unit conversion factor
            if lammps_step.get_lammps_unit_system() == "metal":
                factor = Q_("eV/Å^2/ps")
            else:
                factor = (
                    Q_("kcal/Å^2/fs/mol") / Q_("kcal/mol") * Q_("kcal/mol").to("kJ")
                )
            factor = factor.m_as("W/m^2")
            if ff_form == "class2" or not P["use centroid stress"]:
                # Centroid/stress/atom does not handle class2 ff ... cross-terms?
                lines.append(
                    f"""
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
                )
            else:
                lines.append(
                    f"""
compute             KE all ke/atom
compute             PE all pe/atom

#          centroid doesn't work with kspace, so split into pair and non-pair parts

compute             S_p all stress/atom NULL pair kspace
compute             S_b all centroid/stress/atom NULL bond angle dihedral improper
compute             flux_p all heat/flux KE PE S_p
compute             flux_b all heat/flux KE PE S_b

#          Conversion from kcal/Å^2/fs/mol to W/m^2

variable            factor equal {factor}
variable            Jx equal v_factor*(c_flux_p[1]+c_flux_b[1])/vol
variable            Jy equal v_factor*(c_flux_p[2]+c_flux_b[2])/vol
variable            Jz equal v_factor*(c_flux_p[3]+c_flux_b[3])/vol
"""
                )

        # summary output written 10 times during run so we can see progress
        nevery = 10
        nfreq = int(nsteps / 10)
        nrepeat = int(nfreq / nevery)
        nfreq = nevery * nrepeat
        nfixes += 1
        filename = f"@{self._id[-1]}+nve_summary.trj"
        lines.append(
            f"fix                 {nfixes} all ave/time {nevery} 1 {nfreq} &\n"
            f"                       {properties} &\n"
            "                       off 2 &\n"
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
                "                       off 2 &\n"
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
        lines.append("")

        return {
            "script": lines,
            "postscript": None,
            "use python": False,
        }

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
            title2 += " Sxx Syy Szz Sxy Sxz Syz"
            properties += " v_sxx v_syy v_szz v_sxy v_sxz v_syz"
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
                f"dump_modify         {ndumps} sort id"
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
                f"dump_modify         {ndumps} sort id"
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
        sampling = P["shear stress"]
        if sampling != "never":
            if "interval" in sampling:
                t_s = lammps_step.to_lammps_units(
                    P["shear stress rate"], quantity="time"
                )
                n = max(1, round(t_s / timestep))
            else:
                n = max(1, nsteps // P["shear stress number of samples"])
            filename = f"@{self._id[-1]}+shear_stress.trj"
            nfixes += 1
            dt = (n * P["timestep"]).to_compact()
            text = json.dumps(
                {
                    "code": "LAMMPS",
                    "type": "shear stress",
                    "dt": dt.magnitude,
                    "tunits": str(dt.u),
                    "nsteps": nsteps // n,
                },
                separators=(",", ":"),
            )
            title1 = "!MolSSI trajectory 2.0 " + text
            title2 = "Pxy Pxz Pyz"
            lines.append(
                "\n"
                f"fix                 {nfixes} all ave/time {n} 1 {n} "
                "v_pxy v_pxz v_pyz &\n"
                f"                        title1 '{title1}' &\n"
                f"                        title2 '{title2}' &\n"
                f"                        file {filename}"
            )

        return lines, ncomputes, ndumps, nfixes

    def _to_extxyz(self, data):
        """Create the text of ASE extxyz format for a structure

        Parameters
        ----------
        data : dict(str, any)
            Dictionary containing symbols, coordinates, and other atom properties

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
            line = 'Lattice="'
            line += " ".join([f"{v:.5f}" for v in data["lattice"]])
            line += '"'
        line += " Properties=species:S:1:pos:R:3"
        if "forces" in data:
            line += ":REF_forces:R:3"
            ffact = Q_("kcal/mol/Å").m_as("eV/Å")
        if "velocities" in data:
            line += ":vel:R:3"

        # Add in other properties, like the energy ... currently just energy and stress
        if "Epe" in data:
            # The energy is the potential, not total energy which includes kinetic
            efact = Q_("kcal/mol").m_as("eV")
            line += f" REF_energy={efact*float(data['Epe']):.4f}"
        if "Pxx" in data:
            sfact = Q_("atm").m_as("eV/Å^3")
            Sxx = f"{sfact*float(data['Sxx']):.7f}"
            Syy = f"{sfact*float(data['Syy']):.7f}"
            Szz = f"{sfact*float(data['Szz']):.7f}"
            Syz = f"{sfact*float(data['Syz']):.7f}"
            Sxz = f"{sfact*float(data['Sxz']):.7f}"
            Sxy = f"{sfact*float(data['sxy']):.7f}"
            line += f' REF_stress="{Sxx} {Syy} {Szz} {Syz} {Sxz} {Sxy}" pbc="T T T"'
        else:
            line += ' pbc="F F F"'

        lines.append(line)

        # And the coordinate lines
        for i, symbol in enumerate(data["symbols"]):
            line = f"{symbol:<2}"
            for val in data["xyz"][i]:
                line += f" {val:14.8f}"
            if "forces" in data:
                for val in data["forces"][i]:
                    line += f" {ffact*val:14.8f}"
            if "velocities" in data:
                for val in data["velocities"][i]:
                    line += f" {val:14.8f}"
            lines.append(line)

        # Add empty line so text ends with newline
        lines.append("")

        return "\n".join(lines)
