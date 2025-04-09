# -*- coding: utf-8 -*-

"""NVT (canonical) dynamics in LAMMPS"""

import json

import lammps_step
import logging
import seamm
from seamm_util import units_class, Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __
import random

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("lammps")

thermostat_metadata = {
    "Nose-Hoover": {
        "documentation": "http://lammps.sandia.gov/doc/fix_nvt.html",
        "references": ["Shinoda", "Tuckerman"],
    },
    "Berendsen": {
        "documentation": "http://lammps.sandia.gov/doc/fix_temp_berendsen.html",
        "references": ["Berendsen"],
    },
    "canonical sampling, velocity rescaling (csvr)": {
        "documentation": "http://lammps.sandia.gov/doc/fix_temp_csvr.html",
        "references": ["Bussi1"],
    },
    "canonical sampling, langevin dynamics (csld)": {
        "documentation": "http://lammps.sandia.gov/doc/fix_temp_csvr.html",
        "references": ["Bussi2"],
    },
    "velocity rescaling": {
        "documentation": "http://lammps.sandia.gov/doc/fix_temp_rescale.html",
        "references": [],
    },
    "Langevin": {
        "documentation": "http://lammps.sandia.gov/doc/fix_langevin.html",
        "references": ["Schneider", "Dunweg"],
    },
}


class NVT(lammps_step.NVE):
    def __init__(
        self,
        flowchart=None,
        title="NVT dynamics",
        extension=None,
        logger=logger,
    ):
        """Initialize the node"""

        logger.debug("Creating NVT {}".format(self))

        super().__init__(
            flowchart=flowchart,
            title=title,
            extension=extension,
            logger=logger,
        )

        self.logger.debug("NVT after super init, {}".format(self))

        self.description = "NVT dynamics step in LAMMPS"

        self.logger.debug("NVT.init() creating NVT_Parameters object")

        self._calculation = "nvt"
        self._metadata = lammps_step.metadata
        self.parameters = lammps_step.NVT_Parameters()

        self.logger.debug("NVT.init() completed")

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

        # What will we do?

        if P["T0"] == P["T1"]:
            text = "{time} of canonical (NVT) dynamics at {T0} "
        else:
            text = (
                "{time} of canonical (NVT) dynamics starting "
                " at {T0}, going to {T1}, "
            )
        text += "using a timestep of {timestep} "
        if model is None:
            text += ". The temperature will be controlled "
        else:
            text += f"using the {model} forcefield. The temperature will be controlled "
        if P["thermostat"] == "Nose-Hoover":
            text += "using  a Nose-Hoover thermostat."
            if P["Tchain"] != "3":
                if P["Tloop"] != "1":
                    text += (
                        " The thermostat will use a chain of {Tchain} "
                        "thermostats with {Tloop} subcycles and a "
                    )
                else:
                    text += (
                        " The thermostat will use a chain of {Tchain} "
                        "thermostats and a "
                    )
            elif P["Tloop"] != "1":
                text += " The thermostat will use {Tloop} subcycles and a "
            else:
                text += " The thermostat will use a "
            text += "drag factor of {drag}."
        elif P["thermostat"] == "Berendsen":
            text += "using a Berendsen thermostat with a damping time " "of {Tdamp}"
        elif "csvr" in P["thermostat"]:
            text += (
                "using a canonical sampling thermostat using velocity "
                "rescaling (CSVR) with a damping time of {Tdamp} and "
                "a {seed}."
            )
        elif "csld" in P["thermostat"]:
            text += (
                "using a canonical sampling thermostat using Langevin "
                "dynamics (CSLD) with a damping time of {Tdamp} and "
                "a {seed}."
            )
        elif P["thermostat"] == "velocity rescaling":
            text += (
                "using velocity rescaling every {frequency} with a "
                "temperature window of {window}."
            )
            if P["fraction"] != 1.0:
                text += (
                    " The velocities will only be scaled a fraction "
                    "({fraction}) of the amount needed to fully correct "
                    "the temperature."
                )
        elif P["thermostat"] == "Langevin":
            text += (
                "using a Langevin thermostat with a damping time "
                "of {Tdamp} and a {seed}"
            )
        else:
            text += "using the thermostat given by {thermostat}"

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def describe(self, indent="", json_dict=None):
        """Write out information about what this node will do
        If json_dict is passed in, add information to that dictionary
        so that it can be written out by the controller as appropriate.
        """

        # Can't call super() because it will print too much
        self.visited = True
        job.job("\n" + self.indent + self.header)
        next_node = self.next()

        # Local copies of variables in a dictionary

        P = self.parameters.values_to_dict()
        text = self.description_text(P)
        job.job(__(text, indent=self.indent + "    ", **P))

        return next_node

    def get_input(self, extras=None):
        """Get the input for an NVT dynamics run in LAMMPS"""

        self.description = []

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        _, configuration = self.get_system_configuration()

        # Fix variables with special cases
        timestep, P["timestep"] = self.timestep(P["timestep"])

        if P["seed"] == "random":
            # Apparently the seed must be no larger than 900,000,000
            P["seed"] = int(random.random() * 899999999)

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        self.description.append(
            __(self.description_text(PP), **PP, indent=4 * " ").__str__()
        )

        time = lammps_step.to_lammps_units(P["time"], quantity="time")
        nsteps = round(time / timestep)

        T0 = lammps_step.to_lammps_units(P["T0"], quantity="temperature")
        T1 = lammps_step.to_lammps_units(P["T1"], quantity="temperature")
        Tdamp = lammps_step.to_lammps_units(P["Tdamp"], quantity="time")

        thermo_properties = (
            "time temp press etotal ke pe ebond "
            "eangle edihed eimp evdwl etail ecoul elong"
        )
        properties = "v_time v_temp v_press v_etotal v_ke v_pe v_epair"
        title2 = "tstep t T P Etot Eke Epe Epair"
        if self.parent.have_dreiding_hbonds:
            thermo_properties += " v_N_hbond v_E_hbond"
            properties += " v_N_hbond v_E_hbond"
            title2 += " N_hbond E_hbond"

        lines = []
        lines.append("")
        lines.append(f"# {self.header}")
        lines.append("")
        lines.append("reset_timestep      0")
        lines.append("timestep            {}".format(timestep))
        lines.append("thermo_style        custom {}".format(thermo_properties))
        lines.append("thermo              {}".format(int(nsteps / 100)))

        ncomputes = 0
        ndumps = 0
        nfixes = 0
        if P["thermostat"] == "Nose-Hoover":
            Tchain = P["Tchain"]
            Tloop = P["Tloop"]
            drag = P["drag"]
            nfixes += 1
            lines.append(
                "fix                 {} all nvt ".format(nfixes)
                + "temp {} {} {} ".format(T0, T1, Tdamp)
                + "tchain {} ".format(Tchain)
                + "tloop {} ".format(Tloop)
                + "drag {}".format(drag)
            )
        elif P["thermostat"] == "Berendsen":
            nfixes += 1
            lines.append(
                "fix                 {} ".format(nfixes)
                + "all temp/berendsen "
                + " {} {} {}".format(T0, T1, Tdamp)
            )
            nfixes += 1
            lines.append("fix                 {} ".format(nfixes) + "all nve")
        elif "csvr" in P["thermostat"]:
            seed = P["seed"]
            nfixes += 1
            lines.append(
                "fix                 {} ".format(nfixes)
                + "all temp/csvr "
                + " {} {} {} {}".format(T0, T1, Tdamp, seed)
            )
            nfixes += 1
            lines.append("fix                 {} ".format(nfixes) + "all nve")
        elif "csld" in P["thermostat"]:
            seed = P["seed"]
            nfixes += 1
            lines.append(
                "fix                 {} ".format(nfixes)
                + "all temp/csld "
                + " {} {} {} {}".format(T0, T1, Tdamp, seed)
            )
            nfixes += 1
            lines.append("fix                 {} ".format(nfixes) + "all nve")
        elif P["thermostat"] == "velocity rescaling":
            frequency = lammps_step.to_lammps_units(P["frequency"], quantity="time")

            nevery = max(1, round(nsteps / (frequency / timestep)))
            window = lammps_step.to_lammps_units(P["window"], quantity="temperature")
            fraction = P["fraction"]
            nfixes += 1
            lines.append(
                "fix                 {} ".format(nfixes)
                + "all temp/rescale "
                + "{} {} {} {} {}".format(nevery, T0, T1, window, fraction)
            )
            nfixes += 1
            lines.append("fix                 {} ".format(nfixes) + "all nve")
        elif P["thermostat"] == "Langevin":
            seed = P["seed"]
            nfixes += 1
            lines.append(
                "fix                 {} ".format(nfixes)
                + "all langevin "
                + "{} {} {} {} ".format(T0, T1, Tdamp, seed)
            )
            nfixes += 1
            lines.append("fix                 {} ".format(nfixes) + "all nve")
        else:
            raise RuntimeError(
                "Don't recognize temperature control " + "'{}'".format(P["thermostat"])
            )

        # Add the citation for the thermostat
        for citation in thermostat_metadata[P["thermostat"]]["references"]:
            self.references.cite(
                raw=self._bibliography[citation],
                alias=citation,
                module="lammps_step",
                level=1,
                note="Citation for thermostat.",
            )
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
            lines.append(
                f"""
compute             KE all ke/atom
compute             PE all pe/atom

#          centroid doesn't work with kspace, so split into pair and non-pair parts

compute             S_p all stress/atom NULL pair kspace
compute             S_b all centroid/stress/atom NULL bond angle dihedral improper
compute             flux_p all heat/flux KE PE S_p
compute             flux_b all heat/flux KE PE S_b

#          Conversion from kcal/Å^2/fs/mol to W/m^2")

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
        filename = f"@{self._id[-1]}+nvt_summary.trj"
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
            nrepeat = 1
            nfreq = nevery * nrepeat
            nfixes += 1

            dt = (nevery * P["timestep"]).to_compact()
            if T0 == T1:
                text = json.dumps(
                    {
                        "code": "LAMMPS",
                        "type": "NVT",
                        "dt": dt.magnitude,
                        "tunits": str(dt.u),
                        "nsteps": nsteps // nevery,
                        "T": T0,
                        "Tunits": "K",
                    },
                    separators=(",", ":"),
                )
            else:
                text = json.dumps(
                    {
                        "code": "LAMMPS",
                        "type": "NVT",
                        "dt": dt.magnitude,
                        "units": str(dt.u),
                        "nsteps": nsteps // nevery,
                        "T0": T0,
                        "T1": T1,
                        "Tunits": "K",
                    },
                    separators=(",", ":"),
                )
            title1 = "!MolSSI trajectory 2.0 " + text
            filename = f"@{self._id[-1]}+nvt_state.trj"
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
        filename = f"@{self._id[-1]}+nvt.dump"
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
        lines.append("")

        return {
            "script": lines,
            "postscript": None,
            "use python": False,
        }
