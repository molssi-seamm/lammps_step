# -*- coding: utf-8 -*-
"""Control parameters for NVT (canonical) dynamics"""

import logging
import seamm

logger = logging.getLogger(__name__)


class VelocitiesParameters(seamm.Parameters):
    """The control parameters for NVT dynamics in LAMMPS"""

    parameters = {
        "method": {
            "default": "using a random distribution",
            "kind": "enumeration",
            "format_string": "s",
            "enumeration": (
                "using a random distribution",
                "scaling current velocities",
            ),
            "description": "Set the temperature",
            "help_text": (
                "The method to use to adjuust the velocities "
                "and hence the temperature."
            ),
        },
        "T": {
            "default": 298.15,
            "kind": "float",
            "default_units": "K",
            "format_string": ".2f",
            "description": "Temperature:",
            "help_text": "The temperature corresponding to the velocities.",
        },
        "seed": {
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
        "remove_momentum": {
            "default": (
                "remove any translational and, for molecular systems, "
                "rotational momentum (default)."
            ),
            "kind": "enumeration",
            "format_string": "s",
            "enumeration": (
                (
                    "remove any translational and, for molecular systems, "
                    "rotational momentum (default)"
                ),
                "remove translational but not rotational momentum",
                "remove rotational but not translational momentum",
                "remove both translational and rotational momentum",
                "remove neither translational nor rotational momentum",
            ),
            "description": "Momentum:",
            "help_text": (
                "Ensure that there is no overall linear "
                "momentum so the system does not translate."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**VelocitiesParameters.parameters, **defaults}, data=data
        )
