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
