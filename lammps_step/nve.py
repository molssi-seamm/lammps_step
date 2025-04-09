# -*- coding: utf-8 -*-

"""NVE (microcanonical) dynamics in LAMMPS"""

import json
from pathlib import Path
import traceback

import lammps_step
import logging
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

    def get_input(self, extras=None):
        """Get the input for an NVE dynamics run in LAMMPS"""

        # See what type of forcefield we have and handle it
        ff = self.get_variable("_forcefield")

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
            if ff.ff_form == "class2" or not P["use centroid stress"]:
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
            if ff.ff_form != "class2" and P["use centroid stress"]:
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

        ff = self.get_variable("_forcefield")
        if ff.ff_form == "reaxff":
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
