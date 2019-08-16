# -*- coding: utf-8 -*-
"""Control parameters for NVE (microcanonical) dynamics"""

import lammps_step
import logging

logger = logging.getLogger(__name__)


class NVE_Parameters(lammps_step.EnergyParameters):
    """The control parameters for NVE dynamics in LAMMPS"""

    parameters = {
        "time": {
            "default": 100.0,
            "kind": "float",
            "default_units": "ps",
            "format_string": ".1f",
            "description": "Simulation time:",
            "help_text": ("The time to simulate in the dynamics run.")
        },
        "timestep": {
            "default": "normal",
            "kind": "float",
            "default_units": "fs",
            "enumeration": ('normal', 'accurate but slow', 'coarse but fast'),
            "format_string": ".1f",
            "description": "Timestep:",
            "help_text": ("The time step for the numerical integration in "
                          "the dynamics. 1 fs is safe for most systems, "
                          "except perhaps at high temperatures. For systems "
                          "without hydrogen, helium, lithium or other light "
                          "elements, 2-4 fs steps are reasonable. The "
                          "timestep needs to be less than 1/10 the highest "
                          "frequency. 10^14 Hz is a period if 10 fs, and "
                          "corresponds to a frequency of 3,300 wavenumbers or "
                          "a wavelength of 3 micrometers.")
        },
        "sampling": {
            "default": 20.0,
            "kind": "float",
            "default_units": "fs",
            "enumeration": ('none', ),
            "format_string": ".1f",
            "description": "Sampling frequency:",
            "help_text": ("How often to sample the energy, temperature, "
                          "etc. during the dynamics run. This controls "
                          "writing to the trajectory files. Faster gives "
                          "more fidelity to a point, but increases the "
                          "file size and slows the calculation down.")
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**NVE_Parameters.parameters, **defaults},
            data=data
        )
