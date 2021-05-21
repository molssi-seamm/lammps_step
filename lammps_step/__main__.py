# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""Handle the installation of the LAMMPS step."""

from .installer import Installer


def run():
    """Handle the extra installation needed.

    * Find and/or install the LAMMPS executables.
    * Add or update information in the SEAMM.ini file for LAMMPS
    """

    # print('The is the installer for the LAMMPS step.')
    # Create an installer object
    installer = Installer()
    installer.run()


if __name__ == "__main__":
    run()
