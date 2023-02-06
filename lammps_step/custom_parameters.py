# -*- coding: utf-8 -*-
"""Control parameters for a custom step in LAMMPS"""

import logging
import seamm

logger = logging.getLogger(__name__)


class CustomParameters(seamm.Parameters):
    """The control parameters for Custom dynamics in LAMMPS"""

    parameters = {
        "script": {
            "default": {},
            "kind": "string",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": "",
            "description": "Script:",
            "help_text": "The custom script.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**CustomParameters.parameters, **defaults}, data=data
        )
