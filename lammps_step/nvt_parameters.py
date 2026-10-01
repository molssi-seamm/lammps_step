# -*- coding: utf-8 -*-

"""Control parameters for NVT (canonical) dynamics"""

import lammps_step
import logging

logger = logging.getLogger(__name__)

# The thermostats that use random numbers, and so a seed
random_thermostats = (
    "canonical sampling, velocity rescaling (csvr)",
    "canonical sampling, langevin dynamics (csld)",
    "Langevin",
)


class NVT_Parameters(lammps_step.NVE_Parameters):
    """The control parameters for NVT dynamics in LAMMPS"""

    parameters = {
        "thermostat": {
            "default": "Nose-Hoover",
            "kind": "string",
            "format_string": "s",
            "enumeration": (
                "Nose-Hoover",
                "Berendsen",
                "canonical sampling, velocity rescaling (csvr)",
                "canonical sampling, langevin dynamics (csld)",
                "velocity rescaling",
                "Langevin",
            ),
            "description": "Thermostat:",
            "help_text": ("The thermostat used to control the " "temperature."),
        },
        "T0": {
            "default": 298.15,
            "kind": "float",
            "default_units": "K",
            "format_string": ".2f",
            "description": "Temperature:",
            "help_text": (
                "The temperature, or initial temperature for " "simulated annealing."
            ),
        },
        "T1": {
            "default": 298.15,
            "kind": "float",
            "default_units": "K",
            "format_string": ".2f",
            "description": "Final temperature:",
            "help_text": "The final temperature for simulated annealing.",
        },
        "Tdamp": {
            "applies_when": {"thermostat": {"not": "velocity rescaling"}},
            "default": 100.0,
            "kind": "float",
            "default_units": "fs",
            "format_string": ".1f",
            "description": "Damping time:",
            "help_text": "The damping time constant for thermostat",
        },
        "Tchain": {
            "applies_when": {"thermostat": "Nose-Hoover"},
            "default": 3,
            "kind": "integer",
            "default_units": None,
            "format_string": "d",
            "description": "Thermostat chain:",
            "help_text": "The number of thermostats in the chain.",
        },
        "Tloop": {
            "applies_when": {"thermostat": "Nose-Hoover"},
            "default": 1,
            "kind": "integer",
            "default_units": None,
            "format_string": "d",
            "description": "Thermostat iterations:",
            "help_text": ("The number of sub-iterations for the thermostat."),
        },
        "drag": {
            "applies_when": {"thermostat": "Nose-Hoover"},
            "default": 0.0,
            "kind": "float",
            "default_units": None,
            "format_string": ".1f",
            "description": "Drag:",
            "help_text": (
                "The amount of drag to apply on the thermostat "
                "and barostat if the pressure is controlled. "
                "Typically a value of 0.2 - 2.0 is appropriate."
            ),
        },
        "seed": {
            "applies_when": {"thermostat": random_thermostats},
            "default": "random",
            "kind": "integer",
            "default_units": None,
            "format_string": "",
            "enumeration": ("random",),
            "description": "Random seed:",
            "help_text": (
                "The seed for the random number generator."
                "'random' means to generate a random integer "
                "as the seed."
            ),
        },
        "frequency": {
            "applies_when": {"thermostat": "velocity rescaling"},
            "default": 100.0,
            "kind": "float",
            "default_units": "fs",
            "format_string": ".1f",
            "description": "Rescaling frequency:",
            "help_text": (
                "The frequency (in time units) for the thermostat " "to operate."
            ),
        },
        "window": {
            "applies_when": {"thermostat": "velocity rescaling"},
            "default": 20.0,
            "kind": "float",
            "default_units": "K",
            "format_string": ".1f",
            "description": "Temperature window:",
            "help_text": (
                "The temperature window to use when rescaling."
                "When the instantaneous temperature is outside "
                "the temperature window, the velocities will be "
                "rescaled."
            ),
        },
        "fraction": {
            "applies_when": {"thermostat": "velocity rescaling"},
            "default": 1.0,
            "kind": "float",
            "default_units": None,
            "format_string": ".1f",
            "description": "Fraction to scale:",
            "help_text": (
                "The fraction of the temperature difference "
                "to correct by scaling the velocities."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(defaults={**NVT_Parameters.parameters, **defaults}, data=data)
