# -*- coding: utf-8 -*-

"""Control parameters for minimization"""
import logging

from .energy_parameters import EnergyParameters


logger = logging.getLogger(__name__)


class MinimizationParameters(EnergyParameters):
    """The control parameters for minimization in LAMMPS"""

    parameters = {
        "convergence": {
            "default": "normal",
            "kind": "string",
            "format_string": "s",
            "enumeration": (
                "crude",
                "loose",
                "normal",
                "tight",
                "custom",
            ),
            "description": "Convergence:",
            "help_text": "The convergence criterion for the minimization.",
        },
        "etol": {
            "default": 0.0,
            "kind": "float",
            "default_units": "kcal/mol",
            "format_string": ".2f",
            "description": "Energy change:",
            "help_text": "Converged when the change in energy is less than this.",
        },
        "ftol": {
            "default": 0.1,
            "kind": "float",
            "default_units": "kcal/mol/Å",
            "format_string": ".2f",
            "description": "Force convergence:",
            "help_text": "Converged when the norm of the forces is less than this.",
        },
        "nsteps": {
            "default": "30*nAtoms",
            "kind": "string",
            "format_string": "",
            "description": "Maximum steps:",
            "help_text": "The maximum number of steps.",
        },
        "nevaluations": {
            "default": "3*nSteps",
            "kind": "string",
            "format_string": "",
            "description": "Maximum evaluations:",
            "help_text": "The maximum number of energy evaluations.",
        },
        "minimizer": {
            "default": "Conjugate Gradient",
            "kind": "string",
            "format_string": "",
            "enumeration": (
                "Conjugate Gradient",
                "Steepest Descent",
                "Hessian-Free truncated Newton",
                "QuickMin",
                "Fire",
            ),
            "description": "Minimizer:",
            "help_text": "Which minimizer to use",
        },
        "timestep": {
            "default": "10",
            "kind": "float",
            "default_units": "fs",
            "enumeration": ("normal", "accurate but slow", "coarse but fast"),
            "format_string": ".1f",
            "description": "Timestep:",
            "help_text": (
                "The time step for the numerical integration in "
                "the dynamics. About 10x that used in MD works well."
            ),
        },
        "optimize cell": {
            "default": "yes",
            "kind": "boolean",
            "format_string": "s",
            "enumeration": ("no", "yes"),
            "description": "Optimize cell",
            "help_text": "Whether to optimize the unit cell.",
        },
        "system type": {
            "default": "fluid",
            "kind": "string",
            "format_string": "s",
            "enumeration": ("fluid", "solid"),
            "description": "Type of system:",
            "help_text": (
                "Whether the system is a fluid, "
                "which flows under shear stress, "
                "or a solid, which can resist "
                "small shear stresses."
            ),
        },
        "allow shear": {
            "default": "no",
            "kind": "boolean",
            "format_string": "s",
            "enumeration": ("no", "yes"),
            "description": "Allow the cell to shear",
            "help_text": "Whether the cell angles can change.",
        },
        "use_stress": {
            "default": "isotropic pressure",
            "kind": "enumeration",
            "enumeration": ("isotropic pressure", "general stress"),
            "format_string": "s",
            "description": "Apply",
            "help_text": (
                "Specify whether to apply (isotropic) "
                "pressure to the system or explicit stresses "
                "for each different direction."
            ),
        },
        "couple": {
            "default": "x, y and z",
            "kind": "enumeration",
            "enumeration": ("x, y and z", "x and y", "x and z", "y and z", "none"),
            "format_string": "s",
            "description": "Directions to couple:",
            "help_text": (
                "The stress in these directions will be "
                "averaged and the cell dilated in fixed "
                "proportions in these directions."
            ),
        },
        "P": {
            "default": 1.0,
            "kind": "float",
            "default_units": "atm",
            "format_string": ".2f",
            "enumeration": tuple(),
            "description": "Pressure:",
            "help_text": "The applied pressure.",
        },
        "Sxx": {
            "default": -1.0,
            "kind": "float",
            "enumeration": ("fixed",),
            "default_units": "atm",
            "format_string": "s",
            "description": "Sxx:",
            "help_text": "The components of the stress tensor.",
        },
        "Syy": {
            "default": -1.0,
            "kind": "float",
            "enumeration": ("fixed",),
            "default_units": "atm",
            "format_string": "s",
            "description": "Syy:",
            "help_text": "The components of the stress tensor.",
        },
        "Szz": {
            "default": -1.0,
            "kind": "float",
            "enumeration": ("fixed",),
            "default_units": "atm",
            "format_string": "s",
            "description": "Szz:",
            "help_text": "The components of the stress tensor.",
        },
        "Syz": {
            "default": "fixed",
            "kind": "float",
            "enumeration": ("fixed",),
            "default_units": "atm",
            "format_string": "s",
            "description": "Syz:",
            "help_text": "The components of the stress tensor.",
        },
        "Sxz": {
            "default": "fixed",
            "kind": "float",
            "enumeration": ("fixed",),
            "default_units": "atm",
            "format_string": "s",
            "description": "Sxz:",
            "help_text": "The components of the stress tensor.",
        },
        "Sxy": {
            "default": "fixed",
            "kind": "float",
            "enumeration": ("fixed",),
            "default_units": "atm",
            "format_string": "s",
            "description": "Sxy:",
            "help_text": "The components of the stress tensor.",
        },
        "nreset": {
            "default": "never",
            "kind": "integer",
            "default_units": None,
            "format_string": "d",
            "enumeration": ("never",),
            "description": "Frequency to reset reference cell:",
            "help_text": (
                "How often, in number of steps, to reset "
                "the reference cell for the strain energy."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={
                **MinimizationParameters.parameters,
                **defaults,
            },
            data=data,
        )
