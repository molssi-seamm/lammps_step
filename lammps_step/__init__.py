# -*- coding: utf-8 -*-

"""
lammps_step
A step for running molecular dynamics using LAMMPS.
"""

# Bring up the classes so that they appear to be directly in
# the lammps_step package.

from lammps_step.lammps_units import set_lammps_unit_system  # noqa: F401
from lammps_step.lammps_units import get_lammps_unit_system  # noqa: F401
from lammps_step.lammps_units import to_lammps_units  # noqa: F401
from lammps_step.lammps_units import from_lammps_units  # noqa: F401

from lammps_step.lammps import bond_style  # noqa: F401
from lammps_step.lammps import angle_style  # noqa: F401
from lammps_step.lammps import dihedral_style  # noqa: F401
from lammps_step.lammps import improper_style  # noqa: F401

from lammps_step.lammps import LAMMPS  # noqa: F401
from lammps_step.lammps_step import LAMMPSStep  # noqa: F401
from lammps_step.tk_lammps import TkLAMMPS  # noqa: F401

from lammps_step.custom import Custom  # noqa: F401
from lammps_step.custom_step import CustomStep  # noqa: F401
from lammps_step.tk_custom import TkCustom  # noqa: F401

from lammps_step.initialization import Initialization  # noqa: F401
from lammps_step.initialization_parameters import InitializationParameters  # noqa: F401, E501
from lammps_step.initialization_parameters import kspace_methods  # noqa: F401
from lammps_step.initialization_step import InitializationStep  # noqa: F401
from lammps_step.tk_initialization import TkInitialization  # noqa: F401

from lammps_step.energy import Energy  # noqa: F401
from lammps_step.energy_parameters import EnergyParameters  # noqa: F401
from lammps_step.energy_step import EnergyStep  # noqa: F401
from lammps_step.tk_energy import TkEnergy  # noqa: F401

from lammps_step.minimization import Minimization  # noqa: F401
from lammps_step.minimization_step import MinimizationStep  # noqa: F401
from lammps_step.tk_minimization import TkMinimization  # noqa: F401

from lammps_step.velocities import Velocities  # noqa: F401
from lammps_step.velocities_parameters import VelocitiesParameters  # noqa: F401, E501
from lammps_step.velocities_step import VelocitiesStep  # noqa: F401
from lammps_step.tk_velocities import TkVelocities  # noqa: F401

from lammps_step.nve import NVE  # noqa: F401
from lammps_step.nve_parameters import NVE_Parameters  # noqa: F401
from lammps_step.nve_step import NVEStep  # noqa: F401
from lammps_step.tk_nve import TkNVE  # noqa: F401

from lammps_step.nvt import NVT  # noqa: F401
from lammps_step.nvt import thermostat_metadata  # noqa: F401
from lammps_step.nvt_parameters import NVT_Parameters  # noqa: F401
from lammps_step.nvt_step import NVTStep  # noqa: F401
from lammps_step.tk_nvt import TkNVT  # noqa: F401

from lammps_step.npt import NPT  # noqa: F401
from lammps_step.npt_parameters import NPT_Parameters  # noqa: F401
from lammps_step.npt_step import NPTStep  # noqa: F401
from lammps_step.tk_npt import TkNPT  # noqa: F401

# Handle versioneer
from ._version import get_versions
__author__ = """Paul Saxe"""
__email__ = 'psaxe@molssi.org'
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

properties = {
    "T":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "temperature",
            "dimensionality": "scalar",
            "type": "float",
            "units": "K"
        },
    "T,stderr":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "stderr of temperature",
            "dimensionality": "scalar",
            "type": "float",
            "units": "K"
        },
    "T,tau":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "autocorrelation time of temperature",
            "dimensionality": "scalar",
            "type": "float",
            "units": "fs"
        },
    "T,inefficiency":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "statistical inefficiency of temperature sampling",
            "dimensionality": "scalar",
            "type": "float",
            "units": ""
        },
    "P":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "pressure",
            "dimensionality": "scalar",
            "type": "float",
            "units": "atm"
        },
    "P,stderr":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "stderr of pressure",
            "dimensionality": "scalar",
            "type": "float",
            "units": "K"
        },
    "P,tau":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "autocorrelation time of pressure",
            "dimensionality": "scalar",
            "type": "float",
            "units": "fs"
        },
    "P,inefficiency":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "statistical inefficiency of pressure sampling",
            "dimensionality": "scalar",
            "type": "float",
            "units": ""
        },
    "density":
        {
            "calculation": ["npt"],
            "description": "pressure",
            "dimensionality": "scalar",
            "type": "float",
            "units": "g/ml"
        },
    "density,stderr":
        {
            "calculation": ["npt"],
            "description": "stderr of density",
            "dimensionality": "scalar",
            "type": "float",
            "units": "K"
        },
    "density,tau":
        {
            "calculation": ["npt"],
            "description": "autocorrelation time of density",
            "dimensionality": "scalar",
            "type": "float",
            "units": "fs"
        },
    "density,inefficiency":
        {
            "calculation": ["npt"],
            "description": "statistical inefficiency of density sampling",
            "dimensionality": "scalar",
            "type": "float",
            "units": ""
        },
    "a":
        {
            "calculation": ["npt"],
            "description": "cell parameter 'a'",
            "dimensionality": "scalar",
            "type": "float",
            "units": "Å"
        },
    "a,stderr":
        {
            "calculation": ["npt"],
            "description": "stderr of cell 'a'",
            "dimensionality": "scalar",
            "type": "float",
            "units": "K"
        },
    "a,tau":
        {
            "calculation": ["npt"],
            "description": "autocorrelation time of cell 'a'",
            "dimensionality": "scalar",
            "type": "float",
            "units": "fs"
        },
    "a,inefficiency":
        {
            "calculation": ["npt"],
            "description": "statistical inefficiency of cell 'a' sampling",
            "dimensionality": "scalar",
            "type": "float",
            "units": ""
        },
    "b":
        {
            "calculation": ["npt"],
            "description": "cell parameter 'b'",
            "dimensionality": "scalar",
            "type": "float",
            "units": "Å"
        },
    "b,stderr":
        {
            "calculation": ["npt"],
            "description": "stderr of cell 'b'",
            "dimensionality": "scalar",
            "type": "float",
            "units": "K"
        },
    "b,tau":
        {
            "calculation": ["npt"],
            "description": "autocorrelation time of cell 'b'",
            "dimensionality": "scalar",
            "type": "float",
            "units": "fs"
        },
    "b,inefficiency":
        {
            "calculation": ["npt"],
            "description": "statistical inefficiency of cell 'b' sampling",
            "dimensionality": "scalar",
            "type": "float",
            "units": ""
        },
    "c":
        {
            "calculation": ["npt"],
            "description": "cell parameter 'c'",
            "dimensionality": "scalar",
            "type": "float",
            "units": "Å"
        },
    "c,stderr":
        {
            "calculation": ["npt"],
            "description": "stderr of cell 'c'",
            "dimensionality": "scalar",
            "type": "float",
            "units": "K"
        },
    "c,tau":
        {
            "calculation": ["npt"],
            "description": "autocorrelation time of cell 'c'",
            "dimensionality": "scalar",
            "type": "float",
            "units": "fs"
        },
    "c,inefficiency":
        {
            "calculation": ["npt"],
            "description": "statistical inefficiency of cell 'c' sampling",
            "dimensionality": "scalar",
            "type": "float",
            "units": ""
        },
    "Etot":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "total energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "kcal/mol"
        },
    "Etot,stderr":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "stderr of total energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "K"
        },
    "Etot,tau":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "autocorrelation time of total energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "fs"
        },
    "Etot,inefficiency":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "statistical inefficiency of total energy sampling",
            "dimensionality": "scalar",
            "type": "float",
            "units": ""
        },
    "Eke":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "kinetic energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "kcal/mol"
        },
    "Eke,stderr":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "stderr of kinetic energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "K"
        },
    "Eke,tau":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "autocorrelation time of kinetic energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "fs"
        },
    "Eke,inefficiency":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description":
                "statistical inefficiency of kinetic energy sampling",
            "dimensionality":
                "scalar",
            "type":
                "float",
            "units":
                ""
        },
    "Epe":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "potential energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "kcal/mol"
        },
    "Epe,stderr":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "stderr of potential energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "K"
        },
    "Epe,tau":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "autocorrelation time of potential energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "fs"
        },
    "Epe,inefficiency":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description":
                "statistical inefficiency of potential energy sampling",
            "dimensionality":
                "scalar",
            "type":
                "float",
            "units":
                ""
        },
    "Epair":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "nonbonded (vdW & electrostatic) energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "kcal/mol"
        },
    "Epair,stderr":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "stderr of nonbond energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "K"
        },
    "Epair,tau":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description": "autocorrelation time of nonbond energy",
            "dimensionality": "scalar",
            "type": "float",
            "units": "fs"
        },
    "Epair,inefficiency":
        {
            "calculation": [
                "nve",
                "nvt",
                "npt",
            ],
            "description":
                "statistical inefficiency of nonbond energy sampling",
            "dimensionality":
                "scalar",
            "type":
                "float",
            "units":
                ""
        },
}
