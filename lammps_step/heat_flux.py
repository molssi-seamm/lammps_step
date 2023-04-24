# -*- coding: utf-8 -*-

"""Non-graphical part of the Heat Flux step in a LAMMPS flowchart
"""

import logging
from pathlib import Path
import pkg_resources

import numpy as np

import lammps_step
from .nve import NVE
import molsystem
import seamm
from seamm_util import Q_, units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("LAMMPS")

# Add this module's properties to the standard properties
path = Path(pkg_resources.resource_filename(__name__, "data/"))
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)

script = """\
#!/usr/bin/env python

from lammps import lammps, LMP_STYLE_ATOM, LMP_TYPE_VECTOR, LMP_TYPE_ARRAY
from mpi4py import MPI
import numpy as np

def J_filter_cb(lmp, ntimestep, nlocal, ids, xyz, fext, args = []):
    # Set forces to 0.0
    fext.fill(0.0)

    # lmp.fix_external_set_energy_global("J_filter", 0.0)
    # lmp.fix_external_set_virial_global("J_filter", [0.0] * 6)

    # update energy
    PE0 = lmp.numpy.extract_fix("PE_ave", LMP_STYLE_ATOM, LMP_TYPE_VECTOR)
    lmp.numpy.fix_external_set_energy_peratom("J_filter", -PE0);

    # and stress
    S_p0 = lmp.numpy.extract_fix("S_p_ave", LMP_STYLE_ATOM, LMP_TYPE_ARRAY)
    S_b0 = lmp.numpy.extract_fix("S_b_ave", LMP_STYLE_ATOM, LMP_TYPE_ARRAY)

    S = np.reshape(S_p0, (nlocal, 6))
    tmp = np.reshape(S_b0, (nlocal, 9))

    S[:, 0] += tmp[:, 0]
    S[:, 1] += tmp[:, 1]
    S[:, 2] += tmp[:, 2]
    S[:, 3] += (tmp[:, 3] + tmp[:, 6]) / 2
    S[:, 4] += (tmp[:, 4] + tmp[:, 7]) / 2
    S[:, 5] += (tmp[:, 5] + tmp[:, 8]) / 2

    lmp.numpy.fix_external_set_virial_peratom("J_filter", -S)

def run_thermal_conductivity():
    lmp = lammps()
    lmp.file("lammps.dat")

    lmp.set_fix_external_callback("J_filter", J_filter_cb, lmp)

    lmp.file("lammps_post.dat")

run_thermal_conductivity()

MPI.Finalize()
"""


class HeatFlux(NVE):
    """
    The non-graphical part of a Heat Flux step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : HeatFluxParameters
        The control parameters for Heat Flux.

    See Also
    --------
    TkHeatFlux,
    HeatFlux, HeatFluxParameters
    """

    def __init__(
        self, flowchart=None, title="Heat Flux", extension=None, logger=logger
    ):
        """A substep for Heat Flux in a subflowchart for LAMMPS.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Heat Flux {self}")

        super().__init__(
            flowchart=flowchart,
            title="Heat Flux",
            extension=extension,
            logger=logger,
        )

        self._calculation = "Heat Flux"
        self._model = None
        self._metadata = lammps_step.metadata
        self.parameters = lammps_step.HeatFluxParameters()

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
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        Returns
        -------
        str
            A description of the current step.
        """
        if not P:
            P = self.parameters.values_to_dict()

        text = (
            "Calculate the heat flux using a microcanonical ensemble (NVE), "
            "running for {time} with a timestep of {timestep}. The "
            "heat flux will be sampled every {sampling}."
        )

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def get_input(self, extras=None):
        """Get the input for a heat flux run in LAMMPS"""

        self.description = []

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        time = lammps_step.to_lammps_units(P["time"], quantity="time")
        timestep, P["timestep"] = self.timestep(P["timestep"])
        timestep_fs = P["timestep"].m_as("fs")
        nsteps = round(time / timestep)

        # Work out the sampling information
        sampling = lammps_step.to_lammps_units(P["sampling"], quantity="time")
        nevery = max(1, round(sampling / timestep))
        nrepeat = 1
        nfreq = nevery * nrepeat

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        self.description.append(
            __(self.description_text(PP), **PP, indent=3 * " ").__str__()
        )

        lines = []
        lines.append("")
        lines.append("#     Heat Flux")
        lines.append("")
        lines.append("reset_timestep      0")
        lines.append("timestep            {}".format(timestep))

        computes = []
        fixes = []
        step_id = "_".join(str(e) for e in self._id)

        # Unit conversion factor
        if lammps_step.get_lammps_unit_system() == "metal":
            factor = Q_("eV/Å^2/ps")
        else:
            factor = Q_("kcal/Å^2/fs/mol") / Q_("kcal/mol") * Q_("kcal/mol").to("kJ")
        factor = factor.m_as("W/m^2")

        lines.append(
            f"""
#          Green-Kubo method

fix                 dynamics all nve

compute             KE all ke/atom
compute             PE all pe/atom

#          centroid doesn't work with kspace, so split into pair and non-pair parts

compute             S_p all stress/atom NULL pair kspace
compute             S_b all centroid/stress/atom NULL bond angle dihedral improper
compute             flux_p all heat/flux KE PE S_p
compute             flux_b all heat/flux KE PE S_b

#          An initial run to settle things down.

run                 1000 post no

#          A run to get the averages for PE and S for scheme 2

fix                 ave all ave/atom 1 1000 1000 c_PE
fix                 p_ave all ave/atom 1 1000 1000 c_S_p[1] c_S_p[2] c_S_p[3] &
                                                   c_S_p[4] c_S_p[5] c_S_p[6]
fix                 b_ave all ave/atom 1 1000 1000 c_S_b[1] c_S_b[2] c_S_b[3] &
                                                   c_S_b[4] c_S_b[5] c_S_b[6] &
                                                   c_S_b[7] c_S_b[8] c_S_b[9]

run                 1000 post no

fix                 PE_ave all store/state 0 f_ave
fix                 S_p_ave all store/state 0 f_p_ave[1] f_p_ave[2] f_p_ave[3] &
                                              f_p_ave[4] f_p_ave[5] f_p_ave[6]
fix                 S_b_ave all store/state 0 f_b_ave[1] f_b_ave[2] f_b_ave[3] &
                                              f_b_ave[4] f_b_ave[5] f_b_ave[6] &
                                              f_b_ave[7] f_b_ave[8] f_b_ave[9]

#          Store the initial values for scheme 1

fix                 PE0 all store/state 0 c_PE
fix                 S_p0 all store/state 0 c_S_p[1] c_S_p[2] c_S_p[3] &
                                           c_S_p[4] c_S_p[5] c_S_p[6]
fix                 S_b0 all store/state 0 c_S_b[1] c_S_b[2] c_S_b[3] &
                                           c_S_b[4] c_S_b[5] c_S_b[6] &
                                           c_S_b[7] c_S_b[8] c_S_b[9]

#          The run is need to capture the store/state fixes

run                 0

unfix               b_ave
unfix               p_ave
unfix               ave

reset_timestep      0

#          Conversion from kcal/Å^2/fs/mol tp W/m^2")

variable            factor equal {factor}
variable            Jx equal v_factor*(c_flux_p[1]+c_flux_b[1])/vol
variable            Jy equal v_factor*(c_flux_p[2]+c_flux_b[2])/vol
variable            Jz equal v_factor*(c_flux_p[3]+c_flux_b[3])/vol

#          The thermo output

thermo              1000
thermo_style        custom time temp press etotal ke pe v_Jx v_Jy v_Jz
thermo_modify       colname v_Jx Jx colname v_Jy Jy colname v_Jz Jz

#          The external fix to adjust the peratom PE and S

fix                 J_filter all external pf/callback 1 1
fix_modify          J_filter energy yes
fix_modify          J_filter virial yes
"""
        )
        # handling of trajectory
        nevery = max(1, round(sampling / timestep))
        nrepeat = 1
        nfreq = nevery * nrepeat

        properties = "v_time v_temp v_press v_etotal v_ke v_pe v_Jx v_Jy v_Jz "
        title1 = (
            "!MolSSI trajectory 1.0 LAMMPS, Heat Flux "
            f" {int(nsteps / nevery)} steps of {timestep_fs * nevery} fs"
        )
        title2 = "tstep t T P Etot Eke Epe Jx Jy Jz "

        lines.append(
            f"""
#          The main trajectory

fix                 trajectory all ave/time {nevery} {nrepeat} {nfreq} &
                       {properties} &
                       off 2 &
                       title1 '{title1}' &
                       title2 '{title2}' &
                       file trajectory_heatflux_{step_id}.seamm_trj
"""
        )

        self.description.append(
            __(
                "The run will be {nsteps:,d} steps of dynamics "
                "sampled every {nevery:n} steps.",
                nsteps=nsteps,
                nevery=nevery,
                indent=7 * " ",
            )
        )

        nevery = 10
        nfreq = int(nsteps / 10)
        nrepeat = int(nfreq / nevery)
        nfreq = nevery * nrepeat

        lines.append(
            f"""
#          summary output written 10 times during run so we can see progress

fix                 summary all ave/time {nevery} {nrepeat} {nfreq} &
                        {properties} &
                        off 2 &
                        title2 '{title2}' &
                        file summary_heatflux_{step_id}.txt
"""
        )

        if extras is not None and "shake" in extras:
            lines.append(extras["shake"].format("constraint"))
            fixes.append("constraint")
        lines.append("")

        post_lines = [
            f"""

run                 {nsteps}

# Reset almost everything
thermo              0
thermo_style        one

variable            factor delete
variable            Jx delete
variable            Jy delete
variable            Jz delete
unfix               summary
unfix               trajectory
unfix               J_filter
unfix               S_b0
unfix               S_p0
unfix               PE0
unfix               PE_ave
uncompute           flux_b
uncompute           flux_p
uncompute           S_b
uncompute           S_p
uncompute           PE
uncompute           KE
unfix               dynamics

"""
        ]
        for compute in computes:
            post_lines.append(f"uncompute           {compute}")
        for fix in fixes:
            post_lines.append(f"unfix               {fix}")

        return {
            "script": lines,
            "postscript": post_lines,
            "use python": True,
            "python script": script,
        }

    def analyze(self, indent="", data={}, properties=None, table=None, **kwargs):
        """Save the heat flux for the thermal conductivity step

        Also print important results to the local step.out file using
        "printer".

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        """
        dir_path = Path(self.directory)
        dir_path.mkdir(parents=True, exist_ok=True)

        # Get the configuration
        _, configuration = self.get_system_configuration(None)

        # We need properties like the temperature and volume.
        dt = Q_(properties["Jx"]["timestep"], "fs")

        Jx = properties["Jx"]["values"]
        Jy = properties["Jy"]["values"]
        Jz = properties["Jz"]["values"]

        # Write out the Heat Flux and related data to a compressed numpy file.
        np.savez_compressed(
            dir_path / "HeatFlux.npz",
            J=np.stack((Jx, Jy, Jz)),
            T=np.double(data["T"]),
            V=np.double(configuration.volume),
            timestep=np.double(dt.magnitude),
        )

        # Pass the property data up to the superclass to handle.
        super().analyze(data=data, table=table)
