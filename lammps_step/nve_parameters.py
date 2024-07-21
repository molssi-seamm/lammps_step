# -*- coding: utf-8 -*-
"""Control parameters for NVE (microcanonical) dynamics"""

import lammps_step
import logging

logger = logging.getLogger(__name__)


class NVE_Parameters(lammps_step.EnergyParameters):
    """The control parameters for NVE dynamics in LAMMPS"""

    parameters = {
        "run_control": {
            "default": "For a fixed length of simulated time.",
            "kind": "enumeration",
            "default_units": None,
            "format_string": "",
            "enumeration": (
                "Until properties converge to the requested accuracy.",
                "For a fixed length of simulated time.",
            ),
            "description": "How long to run? ",
            "help_text": (
                "How to determine when to stop the simulation. "
                "You can give a fixed length of simulation time, "
                "e.g. 15 ps, or you can ask to run long enough to "
                "determine one or more properties to a given "
                "accuracy."
            ),
        },
        "time": {
            "default": 100.0,
            "kind": "float",
            "default_units": "ps",
            "format_string": ".1f",
            "description": "Simulation time:",
            "help_text": ("The time to simulate in the dynamics run."),
        },
        "maximum_time": {
            "default": 1.0,
            "kind": "float",
            "default_units": "ns",
            "format_string": ".1f",
            "description": "Maximum simulation time:",
            "help_text": ("The maximum time to simulate when converging properties."),
        },
        "timestep": {
            "default": "normal",
            "kind": "float",
            "default_units": "fs",
            "enumeration": ("normal", "accurate but slow", "coarse but fast"),
            "format_string": ".1f",
            "description": "Timestep:",
            "help_text": (
                "The time step for the numerical integration in "
                "the dynamics. 1 fs is safe for most systems, "
                "except perhaps at high temperatures. For systems "
                "without hydrogen, helium, lithium or other light "
                "elements, 2-4 fs steps are reasonable. The "
                "timestep needs to be less than 1/10 the highest "
                "frequency. 10^14 Hz is a period if 10 fs, and "
                "corresponds to a frequency of 3,300 wavenumbers or "
                "a wavelength of 3 micrometers.\n"
                "You can enter a value or use the choices, which pick a "
                "reasonable timestep based on the calculation."
            ),
        },
        "control_properties": {
            "default": {},
            "kind": "special",
            "widget": "seamm_widgets.PropertyTable",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": "",
            "description": "Convergence properties",
            "help_text": (
                "The properties to converge when controlling the run " "automatically."
            ),
        },
        "sampling": {
            "default": "50",
            "kind": "float",
            "default_units": "fs",
            "enumeration": ("none",),
            "format_string": ".1f",
            "description": "Sampling frequency:",
            "help_text": (
                "How often to sample the energy, temperature, "
                "etc. during the dynamics run. This controls "
                "writing to the trajectory files. Faster gives "
                "more fidelity to a point, but increases the "
                "file size and slows the calculation down."
                "\nYou can ask for no sampling, give a specific interval, "
                "or allow the system the choose for you."
            ),
        },
    }
    trajectories = {
        "atomic positions": {
            "default": "never",
            "kind": "string",
            "enumeration": ("never", "by number of samples", "by time interval"),
            "format_string": "",
            "description": "Sample the atom positions:",
            "help_text": (
                "How to sample the positions of the atoms, which can be used "
                "for diffusion calculations."
            ),
        },
        "atomic positions rate": {
            "default": 100,
            "kind": "float",
            "default_units": "fs",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Time interval:",
            "help_text": "How often to sample the positions of the atoms",
        },
        "atomic positions number of samples": {
            "default": 1000,
            "kind": "integer",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of samples:",
            "help_text": "How many samples of the positions of the atoms to collect",
        },
        "com positions": {
            "default": "never",
            "kind": "string",
            "enumeration": ("never", "by number of samples", "by time interval"),
            "format_string": "",
            "description": "Sample the positions of the molecules' COM:",
            "help_text": (
                "How to sample the positions of the molecules' COM, which can be used "
                "for diffusion calculations."
            ),
        },
        "com positions rate": {
            "default": 100,
            "kind": "float",
            "default_units": "fs",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Time interval:",
            "help_text": (
                "How often to sample the positions of the centers of mass of the "
                "molecules"
            ),
        },
        "com positions number of samples": {
            "default": 1000,
            "kind": "integer",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of samples:",
            "help_text": (
                "How many samples of the positions of the centers of mass of the "
                "molecules to collect"
            ),
        },
        "atomic velocities": {
            "default": "never",
            "kind": "string",
            "enumeration": ("never", "by number of samples", "by time interval"),
            "format_string": "",
            "description": "Sample the atom velocities:",
            "help_text": (
                "How to sample the velocities of the atoms, which can be used "
                "for diffusion calculations."
            ),
        },
        "atomic velocities rate": {
            "default": 100,
            "kind": "float",
            "default_units": "fs",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Time interval:",
            "help_text": "How often to sample the velocities of the atoms",
        },
        "atomic velocities number of samples": {
            "default": 1000,
            "kind": "integer",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of samples:",
            "help_text": "How many samples of the velocities of the atoms to collect",
        },
        "com velocities": {
            "default": "never",
            "kind": "string",
            "enumeration": ("never", "by number of samples", "by time interval"),
            "format_string": "",
            "description": "Sample the velocities of the molecules' COM:",
            "help_text": (
                "How to sample the velocities of the molecules' COM, which can be used "
                "for diffusion calculations."
            ),
        },
        "com velocities rate": {
            "default": 100,
            "kind": "float",
            "default_units": "fs",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Time interval:",
            "help_text": (
                "How often to sample the velocities of the centers of mass of the "
                "molecules"
            ),
        },
        "com velocities number of samples": {
            "default": 1000,
            "kind": "integer",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of samples:",
            "help_text": (
                "How many samples of the velocities of the centers of mass of the "
                "molecules to collect"
            ),
        },
        "heat flux": {
            "default": "never",
            "kind": "string",
            "enumeration": ("never", "by number of samples", "by time interval"),
            "format_string": "",
            "description": "Sample the heat flux:",
            "help_text": (
                "How to sample the heat flux, usually used for thermal "
                "conductivity calculations. However, it is recommended to use the "
                "heat-flux step instead of using this."
            ),
        },
        "heat flux rate": {
            "default": 100,
            "kind": "float",
            "default_units": "fs",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Time interval:",
            "help_text": "How often to sample the heat flux",
        },
        "heat flux number of samples": {
            "default": 1000,
            "kind": "integer",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of samples:",
            "help_text": "How many samples of the heat flux to collect",
        },
        "use centroid stress": {
            "default": "yes",
            "kind": "boolean",
            "default_units": None,
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Use centroid/stress/atom:",
            "help_text": "Whether to use centroid/stress/atom",
        },
        "shear stress": {
            "default": "never",
            "kind": "string",
            "enumeration": ("never", "by number of samples", "by time interval"),
            "format_string": "",
            "description": "Sample the shear stress:",
            "help_text": (
                "How to sample the shear stress, usually used for viscosity "
                "calculations."
            ),
        },
        "shear stress rate": {
            "default": 100,
            "kind": "float",
            "default_units": "fs",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Time interval:",
            "help_text": "How often to sample the shear stress",
        },
        "shear stress number of samples": {
            "default": 1000,
            "kind": "integer",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of samples:",
            "help_text": "How many samples of the shear stress to collect",
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={
                **NVE_Parameters.parameters,
                **NVE_Parameters.trajectories,
                **defaults,
            },
            data=data,
        )
