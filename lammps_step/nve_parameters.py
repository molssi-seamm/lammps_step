# -*- coding: utf-8 -*-
"""Control parameters for NVE (microcanonical) dynamics"""

import lammps_step
import logging

logger = logging.getLogger(__name__)


class NVE_Parameters(lammps_step.EnergyParameters):
    """The control parameters for NVE dynamics in LAMMPS"""

    parameters = {
        "run_control": {
            "default": "Until properties converge to the requested accuracy.",
            "kind": "enumeration",
            "default_units": None,
            "format_string": "",
            "enumeration": (
                "Until properties converge to the requested accuracy.",
                "For a fixed length of simulated time."
            ),
            "description": "How long to run? ",
            "help_text": (
                "How to determine when to stop the simulation. "
                "You can give a fixed length of simulation time, "
                "e.g. 15 ps, or you can ask to run long enough to "
                "determine one or more properties to a given "
                "accuracy."
            )
        },
        "time": {
            "default": 100.0,
            "kind": "float",
            "default_units": "ps",
            "format_string": ".1f",
            "description": "Simulation time:",
            "help_text": ("The time to simulate in the dynamics run.")
        },
        "maximum_time": {
            "default": 1.0,
            "kind": "float",
            "default_units": "ns",
            "format_string": ".1f",
            "description": "Maximum simulation time:",
            "help_text": (
                "The maximum time to simulate when converging properties."
            )
        },
        "timestep": {
            "default": "normal",
            "kind": "float",
            "default_units": "fs",
            "enumeration": ('normal', 'accurate but slow', 'coarse but fast'),
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
            )
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
                "The properties to converge when controlling the run "
                "automatically."
            )
        },
        "sampling": {
            "default": '50',
            "kind": "float",
            "default_units": "fs",
            "enumeration": ('none',),
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
            )
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**NVE_Parameters.parameters, **defaults},
            data=data
        )
