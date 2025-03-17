# -*- coding: utf-8 -*-

"""Installer for the LAMMPS plug-in.

This handles any further installation needed after installing the Python
package `lammps-step`.
"""

import logging
from pathlib import Path
import pkg_resources
import subprocess

import seamm_installer

logger = logging.getLogger(__name__)


class Installer(seamm_installer.InstallerBase):
    """Handle further installation needed after installing lammps-step.

    The Python package `lammps-step` should already be installed, using `pip`,
    `conda`, or similar. This plug-in-specific installer then checks for the
    LAMMPS executables, installing them if need be, and registers their
    location in seamm.ini.

    .. note:: The LAMMPS step and this assume that the LAMMPS
       executables, if present, are in the same directory. In
       addition, for the parallel version using MPI it is assumed that
       `mpiexec` is either in that same directory, or the correct
       `mpiexec` can be found in the PATH.

    There are a number of ways to determine which are the correct LAMMPS
    executables to use. The aim of this installer is to help the user locate
    the executables. There are a number of possibilities

    #. The correct executables are already available.

        #. If they are already registered in `seamm.ini` there is nothing else
           to do.

        #. They may be in the current path, in which case they need to be added
           to `seamm.ini`.

        #. If a module system is in use, a module may need to be loaded to give
           access to LAMMPS.

        #. They cannot be found automatically, so the user needs to locate the
           executables for the installer.

    #. LAMMPS is not installed on the machine. In this case they can be
       installed in a Conda environment. There is one choice

        #. They can be installed in a separate environment, `seamm-lammps` by
           default.
    """

    def __init__(self, logger=logger):
        # Call the base class initialization, which sets up the commandline
        # parser, amongst other things.
        super().__init__(logger=logger)

        logger.debug("Initializing the LAMMPS installer object.")

        self.environment = "seamm-lammps"
        self.section = "lammps-step"
        self.executables = ["lmp"]
        self.resource_path = Path(pkg_resources.resource_filename(__name__, "data/"))

        # The environment.yaml file for Conda installations.
        logger.debug(f"data directory: {self.resource_path}")
        self.environment_file = self.resource_path / "seamm-lammps.yml"

    def exe_version(self, config):
        """Get the version of the LAMMPS executable.

        Parameters
        ----------
        config : dict
            Dictionary of options for running LAMMPS

        Returns
        -------
        "LAMMPS", str
            The version reported by LAMMPS, or 'unknown'.
        """
        environment = config["conda-environment"]
        conda = config["conda"]
        if environment[0] == "~":
            environment = str(Path(environment).expanduser())
            command = f"'{conda}' run --live-stream -p '{environment}'"
        elif Path(environment).is_absolute():
            command = f"'{conda}' run --live-stream -p '{environment}'"
        else:
            command = f"'{conda}' run --live-stream -n '{environment}'"
        command += " lmp -h"
        try:
            result = subprocess.run(
                command,
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
                shell=True,
            )
        except Exception:
            version = "unknown"
        else:
            version = "unknown"
            lines = result.stdout.splitlines()
            if len(lines) > 1:
                line = lines[1]
                version = " ".join(line.split()[6:])

        return "LAMMPS", version
