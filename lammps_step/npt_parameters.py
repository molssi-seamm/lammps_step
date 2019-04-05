# -*- coding: utf-8 -*-
"""Control parameters for NPT (canonical) dynamics"""

import lammps_step
import logging
import pprint

logger = logging.getLogger(__name__)


class NPT_Parameters(lammps_step.NVT_Parameters):
    """The control parameters for NPT dynamics in LAMMPS"""

    parameters = {
        "system type": {
            "value": 'fluid',
            "default": 'fluid',
            "kind": "string",
            "format_string": "s",
            "enumeration": (
                'fluid',
                'solid'
            ),
            "description": "Type of system:",
            "help_text": ("Whether the system is a fluid, "
                          "which flows under shear stress, "
                          "or a solid, which can resist "
                          "small shear stresses.")
        },
        "barostat": {
            "value": 'Nose-Hoover',
            "default": 'Nose-Hoover',
            "kind": "string",
            "format_string": "s",
            "enumeration": (
                'Nose-Hoover',
                'Berendsen'
            ),
            "description": "Barostat:",
            "help_text": ("The barostat used to control the "
                          "pressure")
        },
        "Panneal": {
            "value": "no",
            "default": "no",
            "kind": "boolean",
            "format_string": "s",
            "enumeration": (
                "no",
                "yes"
            ),
            "description": "Change pressure with time",
            "help_text": ("Change the pressure linearly with time "
                          "during the run, i.e. anneal the system.")
        },
        "use_stress": {
            "value": "isotropic pressure",
            "default": "isotropic pressure",
            "kind": "enumeration",
            "enumeration": (
                "isotropic pressure",
                "general stress"
            ),
            "format_string": "s",
            "description": "Apply",
            "help_text": ("Specify whether to apply (isotropic) "
                          "pressure to the system or explicit stresses "
                          "for each different direction.")
        },
        "couple": {
            "value": "x, y and z",
            "default": "x, y and z",
            "kind": "enumeration",
            "enumeration": (
                "x, y and z",
                "x and y",
                "x and z",
                "y and z",
                "none"
            ),
            "format_string": "s",
            "description": "Directions to couple:",
            "help_text": ("The stress in these directions will be "
                          "averaged and the cell dilated in fixed "
                          "proportions in these directions.")
        },
        "Pinitial": {
            "value": 1.0,
            "default": 1.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "enumeration": tuple(),
            "description": "Initial pressure:",
            "help_text": ("The initial pressure.")
        },
        "Pfinal": {
            "value": 1.0,
            "default": 1.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "enumeration": tuple(),
            "description": "Final pressure:",
            "help_text": ("The final pressure.")
        },
        "Pdamp": {
            "value": 1000.0,
            "default": 1000.0,
            "kind": "float",
            "units": "fs",
            "default_units": "fs",
            "format_string": ".1f",
            "description": "Damping time:",
            "help_text": ("The damping time constant for the barostat. "
                          "Typically a value around 1000 timesteps "
                          "works well.")
        },
        "Sxx,initial": {
            "value": 1.0,
            "default": 1.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Sxx initial:",
            "help_text": "The initial components of the stress tensor."
        },
        "Syy,initial": {
            "value": 1.0,
            "default": 1.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Syy initial:",
            "help_text": "The initial components of the stress tensor."
        },
        "Szz,initial": {
            "value": 1.0,
            "default": 1.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Szz initial:",
            "help_text": "The initial components of the stress tensor."
        },
        "Sxy,initial": {
            "value": 0.0,
            "default": 0.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Sxy initial:",
            "help_text": "The initial components of the stress tensor."
        },
        "Sxz,initial": {
            "value": 0.0,
            "default": 0.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Sxz initial:",
            "help_text": "The initial components of the stress tensor."
        },
        "Syz,initial": {
            "value": 0.0,
            "default": 0.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Syz initial:",
            "help_text": "The initial components of the stress tensor."
        },
        "Sxx,final": {
            "value": 1.0,
            "default": 1.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Sxx final:",
            "help_text": "The final components of the stress tensor."
        },
        "Syy,final": {
            "value": 1.0,
            "default": 1.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Syy final:",
            "help_text": "The final components of the stress tensor."
        },
        "Szz,final": {
            "value": 1.0,
            "default": 1.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Szz final:",
            "help_text": "The final components of the stress tensor."
        },
        "Sxy,final": {
            "value": 0.0,
            "default": 0.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Sxy final:",
            "help_text": "The final components of the stress tensor."
        },
        "Sxz,final": {
            "value": 0.0,
            "default": 0.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Sxz final:",
            "help_text": "The final components of the stress tensor."
        },
        "Syz,final": {
            "value": 0.0,
            "default": 0.0,
            "kind": "float",
            "units": "atm",
            "default_units": "atm",
            "format_string": ".2f",
            "description": "Syz final:",
            "help_text": "The final components of the stress tensor."
        },
        "Sxx damp": {
            "value": 1000.0,
            "default": 1000.0,
            "kind": "float",
            "units": "fs",
            "default_units": "fs",
            "format_string": ".1f",
            "description": "Sxx damping time:",
            "help_text": ("The damping time constant for the barostat. "
                          "Typically a value around 1000 timesteps "
                          "works well.")
        },
        "Syy damp": {
            "value": 1000.0,
            "default": 1000.0,
            "kind": "float",
            "units": "fs",
            "default_units": "fs",
            "format_string": ".1f",
            "description": "Syy damping time:",
            "help_text": ("The damping time constant for the barostat. "
                          "Typically a value around 1000 timesteps "
                          "works well.")
        },
        "Szz damp": {
            "value": 1000.0,
            "default": 1000.0,
            "kind": "float",
            "units": "fs",
            "default_units": "fs",
            "format_string": ".1f",
            "description": "Szz damping time:",
            "help_text": ("The damping time constant for the barostat. "
                          "Typically a value around 1000 timesteps "
                          "works well.")
        },
        "Sxy damp": {
            "value": 1000.0,
            "default": 1000.0,
            "kind": "float",
            "units": "fs",
            "default_units": "fs",
            "format_string": ".1f",
            "description": "Sxy damping time:",
            "help_text": ("The damping time constant for the barostat. "
                          "Typically a value around 1000 timesteps "
                          "works well.")
        },
        "Sxz damp": {
            "value": 1000.0,
            "default": 1000.0,
            "kind": "float",
            "units": "fs",
            "default_units": "fs",
            "format_string": ".1f",
            "description": "Sxz damping time:",
            "help_text": ("The damping time constant for the barostat. "
                          "Typically a value around 1000 timesteps "
                          "works well.")
        },
        "Syz damp": {
            "value": 1000.0,
            "default": 1000.0,
            "kind": "float",
            "units": "fs",
            "default_units": "fs",
            "format_string": ".1f",
            "description": "Syz damping time:",
            "help_text": ("The damping time constant for the barostat. "
                          "Typically a value around 1000 timesteps "
                          "works well.")
        },
        "nreset": {
            "value": "never",
            "default": "never",
            "kind": "integer",
            "units": None,
            "default_units": None,
            "format_string": "d",
            "enumeration": ("never",),
            "description": "Frequency to reset reference cell:",
            "help_text": ("How often, in number of steps, to reset "
                          "the reference cell for the strain energy.")
        },
        "mtk": {
            "value": "yes",
            "default": "yes",
            "kind": "boolean",
            "format_string": "s",
            "enumeration": (
                "no",
                "yes"
            ),
            "description": "Use corrected Hoover barostat",
            "help_text": ("Include the correction terms due to Martyna, "
                          "Tuckerman, and Klein in the equations of "
                          "motion. If not, the original Hoover barostat "
                          "is used, whose volume probability distribution "
                          "function differs from the true NPT and NPH "
                          "ensembles by a factor of 1/V. By default the "
                          "correct equations are used, but in many cases "
                          "the difference is negligible.")
        },
        "modulus": {
            "value": 50.0,
            "default": 50.0,
            "kind": "float",
            "units": "GPa",
            "default_units": "GPa",
            "format_string": ".2f",
            "enumeration": tuple(),
            "description": "Bulk modulus:",
            "help_text": ("The bulk modulus. For liquids a typical "
                          "modulus ranges from a fraction of a GPa "
                          "to a few GPa. For crystalline solids, "
                          "25-500 GPa is a reasonable range.")
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**NPT_Parameters.parameters, **defaults},
            data=data
        )
