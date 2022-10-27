# -*- coding: utf-8 -*-

"""
lammps_step
A step for running molecular dynamics using LAMMPS.
"""

# Bring up the classes so that they appear to be directly in
# the lammps_step package.


# The metadata
from .metadata import metadata  # noqa: F401

from .lammps_units import set_lammps_unit_system  # noqa: F401
from .lammps_units import get_lammps_unit_system  # noqa: F401
from .lammps_units import to_lammps_units  # noqa: F401
from .lammps_units import from_lammps_units  # noqa: F401

from .lammps import bond_style  # noqa: F401
from .lammps import angle_style  # noqa: F401
from .lammps import dihedral_style  # noqa: F401
from .lammps import improper_style  # noqa: F401

from .lammps import LAMMPS  # noqa: F401
from .lammps_step import LAMMPSStep  # noqa: F401
from .tk_lammps import TkLAMMPS  # noqa: F401

from .custom import Custom  # noqa: F401
from .custom_step import CustomStep  # noqa: F401
from .tk_custom import TkCustom  # noqa: F401

from .initialization import Initialization  # noqa: F401
from .initialization_parameters import (  # noqa: F401
    InitializationParameters,
)
from .initialization_parameters import kspace_methods  # noqa: F401
from .initialization_step import InitializationStep  # noqa: F401
from .tk_initialization import TkInitialization  # noqa: F401

from .energy import Energy  # noqa: F401
from .energy_parameters import EnergyParameters  # noqa: F401
from .energy_step import EnergyStep  # noqa: F401
from .tk_energy import TkEnergy  # noqa: F401

from .minimization import Minimization  # noqa: F401
from .minimization_step import MinimizationStep  # noqa: F401
from .tk_minimization import TkMinimization  # noqa: F401

from .velocities import Velocities  # noqa: F401
from .velocities_parameters import VelocitiesParameters  # noqa: F401, E501
from .velocities_step import VelocitiesStep  # noqa: F401
from .tk_velocities import TkVelocities  # noqa: F401

from .nve import NVE  # noqa: F401
from .nve_parameters import NVE_Parameters  # noqa: F401
from .nve_step import NVEStep  # noqa: F401
from .tk_nve import TkNVE  # noqa: F401

from .nvt import NVT  # noqa: F401
from .nvt import thermostat_metadata  # noqa: F401
from .nvt_parameters import NVT_Parameters  # noqa: F401
from .nvt_step import NVTStep  # noqa: F401
from .tk_nvt import TkNVT  # noqa: F401

from .npt import NPT  # noqa: F401
from .npt_parameters import NPT_Parameters  # noqa: F401
from .npt_step import NPTStep  # noqa: F401
from .tk_npt import TkNPT  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = """Paul Saxe"""
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
