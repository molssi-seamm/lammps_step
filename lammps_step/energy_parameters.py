# -*- coding: utf-8 -*-
"""Control parameters for a single-point energy (SPE) in LAMMPS"""

import logging
import seamm

logger = logging.getLogger(__name__)


class EnergyParameters(seamm.Parameters):
    """The control parameters for Energy dynamics in LAMMPS"""

    parameters = {
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": ("The results to save to variables or in " "tables. "),
        },
        "create tables": {
            "default": "yes",
            "kind": "boolean",
            "default_units": None,
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Create tables as needed:",
            "help_text": (
                "Whether to create tables as needed for "
                "results being saved into tables."
            ),
        },
    }

    # Rules shared by the dialog and the flowchart builder (see seamm.Parameters).
    # The simple conditions are the "applies_when" entries in the parameters of each
    # step. The steps that change the cell (Minimization and NPT) apply either a
    # pressure or a general stress: a fluid always takes a pressure, and a solid
    # takes a pressure or a stress, as 'use_stress' says.

    pressure_keys = ()
    """The parameters for the (isotropic) pressure, which apply only when a pressure
    rather than a general stress is applied."""

    not_implemented = ()
    """Parameters that are defined but not yet implemented, so never apply."""

    def _isotropic(self, values):
        """Whether a pressure, rather than a general stress, is applied."""
        system_type = values.get("system type", self["system type"].value)
        if self._is_expr(system_type) or system_type == "fluid":
            return True
        use_stress = values.get("use_stress", self["use_stress"].value)
        return self._is_expr(use_stress) or use_stress == "isotropic pressure"

    def applies(self, key, values=None, _seen=None):
        """As seamm.Parameters.applies, plus: the pressure applies only when a
        pressure, not a general stress, is applied."""
        if values is None:
            values = self.current_values()
        if key in self.not_implemented:
            return False
        if not super().applies(key, values, _seen):
            return False
        if key in self.pressure_keys:
            return self._isotropic(values)
        return True

    def not_applicable_reason(self, key, values=None):
        """Why a parameter does not apply, for the builder's messages."""
        if values is None:
            values = self.current_values()
        if key in self.not_implemented:
            return "it is not implemented yet"
        reason = super().not_applicable_reason(key, values)
        if reason:
            return reason
        if key in self.pressure_keys and not self._isotropic(values):
            return (
                "a general stress, not a pressure, is applied (it applies when "
                "'system type' is 'fluid' or 'use_stress' is 'isotropic pressure')"
            )
        return ""

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={
                **EnergyParameters.parameters,
                **seamm.standard_parameters.structure_handling_parameters,
                **defaults,
            },
            data=data,
        )

        # Do any local editing of defaults
        tmp = self["system name"]
        tmp._data["enumeration"] = ["simulated with {model}", *tmp.enumeration[1:]]
        tmp.default = "keep current name"

        tmp = self["configuration name"]
        tmp._data["enumeration"] = ["simulated with {model}", *tmp.enumeration]
        tmp.default = "simulated with {model}"
