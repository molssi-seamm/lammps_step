# -*- coding: utf-8 -*-

"""NPT (canonical) dynamics in LAMMPS"""

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


class NPT(lammps_step.NVT):
    methods = {
        "Nose-Hoover": {
            "documentation": "https://lammps.sandia.gov/doc/fix_nh.html#fix-npt-command",  # noqa: E501
            "references": ["Shinoda", "Tuckerman"],
        },
        "Berendsen": {
            "documentation": "https://lammps.sandia.gov/doc/fix_press_berendsen.html",  # noqa: E501
            "references": ["Berendsen"],
        },
    }

    def __init__(self, flowchart=None, title="NPT dynamics", extension=None):
        """Initialize the node"""

        logger.debug("Creating NPT {}".format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        logger.debug("NPT after super init, {}".format(self))

        self.description = "NPT dynamics step in LAMMPS"

        logger.debug("NPT.init() creating NPT_Parameters object")

        self._calculation = "npt"
        self._metadata = lammps_step.metadata
        self.parameters = lammps_step.NPT_Parameters()

        logger.debug("NPT.init() completed")

    def description_text(self, P=None):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.
        """

        if not P:
            P = self.parameters.values_to_dict()

        model = self.model

        # What will we do?

        if P["T0"] == P["T1"]:
            text = "{time} of canonical (NPT) dynamics at {T0} "
        else:
            text = (
                "{time} of canonical (NPT) dynamics starting "
                " at {T0}, going to {T1}, "
            )
        text += "using a timestep of {timestep} "
        if model is None:
            text += ". The temperature will be controlled "
        else:
            text += f"using the {model} forcefield. The temperature will be controlled "
        if P["thermostat"] == "Nose-Hoover":
            text += "using a Nose-Hoover thermostat."
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
                "a {random_seed}."
            )
        elif "csld" in P["thermostat"]:
            text += (
                "using a canonical sampling thermostat using Langevin "
                "dynamics (CSLD) with a damping time of {Tdamp} and "
                "a {random_seed}."
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
                "of {Tdamp} and a {random_seed}."
            )
        else:
            text += "using the thermostat given by {thermostat}."

        text += " The pressure will be controlled using a {barostat} barostat."

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
        """Get the input for an NPT dynamics run in LAMMPS"""

        keep_orthorhombic = True

        self.description = []

        _, configuration = self.get_system_configuration()
        if configuration.periodicity == 0:
            raise RuntimeError("Cannot run NPT dynamics on non-periodic system!")

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

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

        self.description.append(__(self.description_text(PP), **PP, indent=4 * " "))

        time = lammps_step.to_lammps_units(P["time"], quantity="time")
        nsteps = round(time / timestep)

        T0 = lammps_step.to_lammps_units(P["T0"], quantity="temperature")
        T1 = lammps_step.to_lammps_units(P["T1"], quantity="temperature")
        Tdamp = lammps_step.to_lammps_units(P["Tdamp"], quantity="time")

        barostat = P["barostat"]
        if barostat == "Berendsen":
            modulus = lammps_step.to_lammps_units(P["modulus"], quantity="pressure")

        # Work out the pressure/stress part of the command
        ptext = self.get_pressure_text(P, keep_orthorhombic)

        thermo_properties = (
            "time temp press etotal ke pe ebond "
            "eangle edihed eimp evdwl etail ecoul elong"
        )

        properties = (
            "v_time v_temp v_press v_vol v_density v_cella v_cellb "
            "v_cellc v_etotal v_ke v_pe v_epair"
        )
        title2 = "tstep t T P V density a b c Etot Eke Epe Epair"
        if self.parent.have_dreiding_hbonds:
            thermo_properties += " v_N_hbond v_E_hbond"
            properties += " v_N_hbond v_E_hbond"
            title2 += " N_hbond E_hbond"

        # and build the LAMMPS script
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
            if barostat == "Nose-Hoover":
                nfixes += 1
                lines.append(
                    "fix                 {} all npt ".format(nfixes)
                    + "temp {} {} {} ".format(T0, T1, Tdamp)
                    + "tchain {} ".format(Tchain)
                    + "tloop {} ".format(Tloop)
                    + "drag {}".format(drag)
                    + ptext
                )

                for citation in NPT.methods["Nose-Hoover"]["references"]:
                    self.references.cite(
                        raw=self._bibliography[citation],
                        alias=citation,
                        module="lammps_step",
                        level=1,
                        note="Citation for NPT barostat.",
                    )
            else:
                nfixes += 1
                lines.append(
                    "fix                 {} all nvt ".format(nfixes)
                    + "temp {} {} {} ".format(T0, T1, Tdamp)
                    + "tchain {} ".format(Tchain)
                    + "tloop {} ".format(Tloop)
                    + "drag {}".format(drag)
                )
                nfixes += 1
                lines.append(
                    "fix                 {} all ".format(nfixes)
                    + "press/berendsen "
                    + ptext
                    + " modulus {}".format(modulus)
                )

                for citation in NPT.methods["Berendsen"]["references"]:
                    self.references.cite(
                        raw=self._bibliography[citation],
                        alias=citation,
                        module="lammps_step",
                        level=1,
                        note="Citation for NPT barostat.",
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
            frequency = P["frequency"]
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
        metadata = lammps_step.thermostat_metadata
        for citation in metadata[P["thermostat"]]["references"]:
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
        filename = f"@{self._id[-1]}+npt_summary.trj"
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
                    indent=4 * " ",
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
            P0tmp = P["Pinitial"].to_compact()
            P1tmp = P["Pfinal"].to(P0tmp.units)
            tmp = {
                "code": "LAMMPS",
                "type": "NPT",
                "dt": dt.magnitude,
                "tunits": str(dt.u),
                "nsteps": nsteps // nevery,
                "Tunits": "K",
                "Punits": str(P0tmp.u),
            }
            if P0tmp == P1tmp:
                tmp["P"] = P0tmp.magnitude
            else:
                tmp["P0"] = P0tmp.magnitude
                tmp["T1"] = P1tmp.magnitude
            if T0 == T1:
                tmp["T"] = T0
            else:
                tmp["T0"] = T0
                tmp["T1"] = T1
            title1 = "!MolSSI trajectory 2.0 " + json.dumps(tmp, separators=(",", ":"))
            filename = f"@{self._id[-1]}+npt_state.trj"
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
        filename = f"@{self._id[-1]}+npt.dump"
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

        for fix in range(1, nfixes + 1):
            lines.append("unfix               {}".format(fix))
        lines.append("")

        return {
            "script": lines,
            "postscript": None,
            "use python": False,
        }

    def get_pressure_text(self, P, keep_orthorhombic):
        """Work out and return the pressure/stress part of the
        'fix npt' or 'fix berendsen' in LAMMPS
        """
        system_type = P["system type"]
        Panneal = P["Panneal"]
        if system_type == "fluid":
            use_stress = "isotropic pressure"
            couple = "x, y and z"
        else:
            use_stress = P["use_stress"]
            couple = P["couple"]

        if use_stress == "isotropic pressure":
            if couple == "x, y and z":
                if keep_orthorhombic:
                    ptext = " iso {P0} {P1} {Pdamp}"
                else:
                    ptext = (
                        " couple xyz "
                        "x {P0} {P1} {Pdamp} "
                        "y {P0} {P1} {Pdamp} "
                        "z {P0} {P1} {Pdamp}"
                    )
                if not keep_orthorhombic:
                    ptext += (
                        " xy 0.0 0.0 {Pdamp} "
                        "xz 0.0 0.0 {Pdamp} "
                        "yz 0.0 0.0 {Pdamp}"
                    )
            elif couple == "x and y":
                ptext = (
                    " couple xy "
                    "x {P0} {P1} {Pdamp} "
                    "y {P0} {P1} {Pdamp} "
                    "z {P0} {P1} {Pdamp}"
                )
                if not keep_orthorhombic:
                    ptext += (
                        " xy 0.0 0.0 {Pdamp} "
                        "xz 0.0 0.0 {Pdamp} "
                        "yz 0.0 0.0 {Pdamp}"
                    )
            elif couple == "x and z":
                ptext = (
                    " couple xz "
                    "x {P0} {P1} {Pdamp} "
                    "y {P0} {P1} {Pdamp} "
                    "z {P0} {P1} {Pdamp}"
                )
                if not keep_orthorhombic:
                    ptext += (
                        " xy 0.0 0.0 {Pdamp} "
                        "xz 0.0 0.0 {Pdamp} "
                        "yz 0.0 0.0 {Pdamp}"
                    )
            elif couple == "y and z":
                ptext = (
                    " couple yz "
                    "x {P0} {P1} {Pdamp} "
                    "y {P0} {P1} {Pdamp} "
                    "z {P0} {P1} {Pdamp}"
                )
                if not keep_orthorhombic:
                    ptext += (
                        " xy 0.0 0.0 {Pdamp} "
                        "xz 0.0 0.0 {Pdamp} "
                        "yz 0.0 0.0 {Pdamp}"
                    )
            else:
                if keep_orthorhombic:
                    ptext = " aniso {P0} {P1} {Pdamp}"
                else:
                    ptext = " tri {P0} {P1} {Pdamp}"

            P0 = lammps_step.to_lammps_units(P["Pinitial"], quantity="pressure")
            if Panneal:
                P1 = lammps_step.to_lammps_units(P["Pfinal"], quantity="pressure")
            else:
                P1 = P0
            Pdamp = lammps_step.to_lammps_units(P["Pdamp"], quantity="time")

            ptext = ptext.format(P0=P0, P1=P1, Pdamp=Pdamp)
        else:
            if couple == "x, y and z":
                ptext = (
                    " couple = xyz "
                    "x {Sxx0} {Sxx1} {Dxx} "
                    "y {Sxx0} {Sxx1} {Dxx} "
                    "z {Sxx0} {Sxx1} {Dxx}"
                )
            elif couple == "x and y":
                ptext = (
                    " couple = xy "
                    "x {Sxx0} {Sxx1} {Dxx} "
                    "y {Sxx0} {Sxx1} {Dxx} "
                    "z {Szz0} {Szz1} {Dzz}"
                )
            elif couple == "x and z":
                ptext = (
                    " couple = xz "
                    "x {Sxx0} {Sxx1} {Dxx} "
                    "y {Syy0} {Syy1} {Dyy} "
                    "z {Sxx0} {Sxx1} {Dxx}"
                )
            elif couple == "y and z":
                ptext = (
                    " couple = yz "
                    "x {Sxx0} {Sxx1} {Dxx} "
                    "y {Syy0} {Syy1} {Dyy} "
                    "z {Syy0} {Syy1} {Dyy}"
                )
            else:
                # elif couple == 'none':
                ptext = (
                    " couple = none "
                    "x {Sxx0} {Sxx1} {Dxx} "
                    "y {Syy0} {Syy1} {Dyy} "
                    "z {Szz0} {Szz1} {Dzz}"
                )

            if not keep_orthorhombic:
                ptext += (
                    " xy {Sxy0} {Sxy1} {Dxy} "
                    "xz {Sxz0} {Sxz1} {Dxz} "
                    "yz {Syz0} {Syz1} {Dyz}"
                )

            Tmp = {}
            Tmp["Sxx0"] = lammps_step.to_lammps_units(
                P["Sxx,initial"], quantity="pressure"
            )
            Tmp["Syy0"] = lammps_step.to_lammps_units(
                P["Syy,initial"], quantity="pressure"
            )
            Tmp["Szz0"] = lammps_step.to_lammps_units(
                P["Szz,initial"], quantity="pressure"
            )
            if Panneal:
                Tmp["Sxx1"] = lammps_step.to_lammps_units(
                    P["Sxx,final"], quantity="pressure"
                )
                Tmp["Syy1"] = lammps_step.to_lammps_units(
                    P["Syy,final"], quantity="pressure"
                )
                Tmp["Szz1"] = lammps_step.to_lammps_units(
                    P["Szz,final"], quantity="pressure"
                )
            else:
                Tmp["Sxx1"] = Tmp["Sxx0"]
                Tmp["Syy1"] = Tmp["Syy0"]
                Tmp["Szz1"] = Tmp["Szz0"]
            Tmp["Dxx"] = lammps_step.to_lammps_units(P["Sxx damp"], quantity="pressure")
            Tmp["Dyy"] = lammps_step.to_lammps_units(P["Syy damp"], quantity="pressure")
            Tmp["Dzz"] = lammps_step.to_lammps_units(P["Szz damp"], quantity="pressure")

            if not keep_orthorhombic:
                Tmp["Sxy0"] = lammps_step.to_lammps_units(
                    P["Sxy,initial"], quantity="pressure"
                )
                Tmp["Sxz0"] = lammps_step.to_lammps_units(
                    P["Sxz,initial"], quantity="pressure"
                )
                Tmp["Syz0"] = lammps_step.to_lammps_units(
                    P["Syz,initial"], quantity="pressure"
                )
                if Panneal:
                    Tmp["Sxy1"] = lammps_step.to_lammps_units(
                        P["Sxy,final"], quantity="pressure"
                    )
                    Tmp["Sxz1"] = lammps_step.to_lammps_units(
                        P["Sxz,final"], quantity="pressure"
                    )
                    Tmp["Syz1"] = lammps_step.to_lammps_units(
                        P["Syz,final"], quantity="pressure"
                    )
                else:
                    Tmp["Sxy1"] = Tmp["Sxy0"]
                    Tmp["Sxz1"] = Tmp["Sxz0"]
                    Tmp["Syz1"] = Tmp["Syz0"]
                Tmp["Dxy"] = lammps_step.to_lammps_units(
                    P["Sxy damp"], quantity="pressure"
                )
                Tmp["Dxz"] = lammps_step.to_lammps_units(
                    P["Sxz damp"], quantity="pressure"
                )
                Tmp["Dyz"] = lammps_step.to_lammps_units(
                    P["Syz damp"], quantity="pressure"
                )

            ptext = ptext.format(**Tmp)

        return ptext
