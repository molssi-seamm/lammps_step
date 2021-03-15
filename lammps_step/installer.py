# -*- coding: utf-8 -*-

"""Installer for the LAMMPS plug-in.

This handles any further installation needed after installing the Python
package `lammps-step`.
"""

import logging
import pkg_resources
from pathlib import Path
import shutil
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
    the executables. There are a number of possibilities:

    1. The correct executables are already available.

        1. If they are already registered in `seamm.ini` there is nothing else
           to do.
        2. They may be in the current path, in which case they need to be added
           to `seamm.ini`.
        3. If a module system is in use, a module may need to be loaded to give
           access to LAMMPS.
        3. They cannot be found automatically, so the user needs to locate the
           executables for the installer.

    2. LAMMPS is not installed on the machine. In this case they can be
       installed in a Conda environment. There is one choice:

        1. They can be installed in a separate environment, `seamm-lammps` by
           default.
    """

    def __init__(self, logger=logger):
        # Call the base class initialization, which sets up the commandline
        # parser, amongst other things.
        super().__init__(logger=logger)

        logger.debug('Initializing the LAMMPS installer object.')

        self.executables = ['lmp_serial', 'lmp_mpi']
        # What Conda environment is the default?
        data = self.configuration.get_values('lammps-step')
        if 'conda-environment' in data and data['conda-environment'] != '':
            self.environment = data['conda-environment']
        else:
            self.environment = 'seamm-lammps'

        # The environment.yaml file for Conda installations.
        path = Path(pkg_resources.resource_filename(__name__, 'data/'))
        logger.debug(f"data directory: {path}")
        self.environment_file = path / 'seamm-lammps.yml'

    def check(self):
        """Check the installation and fix errors if requested.

        If the option `yes` is present and True, this method will attempt to
        correct any errors in the configuration file. Use `--yes` on the
        command line to enable this.

        The information in the configuration file is:

            installation
                How LAMMPS is installed. One of `user`, `modules` or `conda`
            conda-environment
                The Conda environment if and only if `installation` = `conda`
            modules
                The environment modules if `installation` = `modules`
            lammps-path
                The path where the LAMMPS executables are. Automatically
                defined if `installation` is `conda` or `modules`, but given
                by the user is it is `user`.

        Returns
        -------
        bool
            True if everything is OK, False otherwise. If `yes` is given as an
            option, the return value is after fixing the configuration.
        """
        self.logger.debug('Entering check method.')
        if not self.configuration.section_exists('lammps-step'):
            if self.options.yes or self.ask_yes_no(
                "There is no section for the LAMMPS step in the configuration "
                f" file ({self.configuration.path}).\nAdd one?",
                default='yes'
            ):
                self.check_configuration_file()
                print("Added the 'lammps-step' section")

        # Get the values from the configuration
        data = self.configuration.get_values('lammps-step')

        # Save the initial values, if any, of the key configuration variables
        if 'lammps-path' in data and data['lammps-path'] != '':
            path = Path(data['lammps-path']).expanduser().resolve()
            initial_lammps_path = path
        else:
            initial_lammps_path = None
        if 'installation' in data and data['installation'] != '':
            initial_installation = data['installation']
        else:
            initial_installation = None
        if 'conda-environment' in data and data['conda-environment'] != '':
            initial_conda_environment = data['conda-environment']
        else:
            initial_conda_environment = None
        if 'modules' in data and data['modules'] != '':
            initial_modules = data['modules']
        else:
            initial_modules = None

        # Is there a valid lammps-path?
        self.logger.debug(
            "Checking for executables in the initial lammps-path "
            f"{initial_lammps_path}."
        )
        if (
            initial_lammps_path is None or
            not self.have_executables(initial_lammps_path)
        ):
            lammps_path = None
        else:
            lammps_path = initial_lammps_path
        self.logger.debug(f"initial-lammps-path = {initial_lammps_path}.")

        # Is there an installation indicated?
        if initial_installation in ('user', 'conda', 'modules'):
            installation = initial_installation
        else:
            installation = None
        self.logger.debug(f"initial-installation = {initial_installation}.")

        if installation == 'conda':
            # Is there a conda environment?
            conda_environment = None
            if (
                initial_conda_environment is None or
                not self.conda.exists(initial_conda_environment)
            ):
                if lammps_path is not None:
                    # see if this path corresponds to a Conda environment
                    for tmp in self.conda.environments:
                        tmp_path = self.conda.path(tmp) / 'bin'
                        if tmp_path == lammps_path:
                            conda_environment = tmp
                            break
                    if conda_environment is not None:
                        if self.options.yes or self.ask_yes_no(
                            "The Conda environment in the config file "
                            "is not correct.\n"
                            f"It should be {conda_environment}. Fix?",
                            default='yes'
                        ):
                            self.configuration.set_value(
                                'lammps-step', 'installation', 'conda'
                            )
                            self.configuration.set_value(
                                'lammps-step', 'conda-environment',
                                conda_environment
                            )
                            self.configuration.set_value(
                                'lammps-step', 'modules', ''
                            )
                            self.configuration.save()
                            print(
                                "Corrected the conda environment to "
                                f"{conda_environment}"
                            )
            else:
                # Have a Conda environment!
                conda_path = self.conda.path(initial_conda_environment) / 'bin'
                self.logger.debug(
                    f"Checking for executables in conda-path: {conda_path}."
                )
                if self.have_executables(conda_path):
                    # All is good!
                    conda_environment = initial_conda_environment
                    if lammps_path is None:
                        if self.options.yes or self.ask_yes_no(
                            "The lammps-path in the config file is not set,"
                            f"but the Conda environment {conda_environment} "
                            "is.\nFix the lammps-path?",
                            default='yes'
                        ):
                            lammps_path = conda_path
                            self.configuration.set_value(
                                'lammps-step', 'lammps-path', lammps_path
                            )
                            self.configuration.set_value(
                                'lammps-step', 'modules', ''
                            )
                            self.configuration.save()
                            print(f"Set the lammps-path to {conda_path}")
                    elif lammps_path != conda_path:
                        if self.options.yes or self.ask_yes_no(
                            f"The lammps-path in the config file {lammps_path}"
                            "is different from that for  the Conda "
                            f"environment {conda_environment} is.\n"
                            "Use the path from the Conda environment?",
                            default='yes'
                        ):
                            lammps_path = conda_path
                            self.configuration.set_value(
                                'lammps-step', 'lammps-path', lammps_path
                            )
                            self.configuration.set_value(
                                'lammps-step', 'modules', ''
                            )
                            self.configuration.save()
                            print(f"Changed the lammps-path to {conda_path}")
                    else:
                        # Everything is fine!
                        pass
        if installation == 'modules':
            print(f"Can't check the actual modules {initial_modules} yet")
            if initial_conda_environment is not None:
                if self.options.yes or self.ask_yes_no(
                    "A Conda environment is given: "
                    f"{initial_conda_environment}.\n"
                    "A Conda environment should not be used when using "
                    "modules. Remove it from the configuration?",
                    default='yes'
                ):
                    self.configuration.set_value(
                        'lammps-step', 'conda-environment', ''
                    )
                    self.configuration.save()
                    print(
                        "Using modules, so removed the conda-environment from "
                        "the configuration"
                    )
        else:
            if lammps_path is None:
                # No path or executables in the path!
                environments = self.conda.environments
                if self.environment in environments:
                    # Make sure it is first!
                    environments.remove(self.environment)
                    environments.insert(0, self.environment)
                for tmp in environments:
                    tmp_path = self.conda.path(tmp) / 'bin'
                    if self.have_executables(tmp_path):
                        if self.options.yes or self.ask_yes_no(
                            "There are no valid executables in the lammps-path"
                            " in the config file, but there are in the Conda "
                            f"environment {tmp}.\n"
                            "Use them?",
                            default='yes'
                        ):
                            conda_environment = tmp
                            lammps_path = tmp_path
                            self.configuration.set_value(
                                'lammps-step', 'lammps-path', lammps_path
                            )
                            self.configuration.set_value(
                                'lammps-step', 'installation', 'conda'
                            )
                            self.configuration.set_value(
                                'lammps-step', 'conda-environment',
                                conda_environment
                            )
                            self.configuration.set_value(
                                'lammps-step', 'modules', ''
                            )
                            self.configuration.save()
                            print(
                                "Will use the conda environment "
                                f"'{conda_environment}'"
                            )
                            break
            if lammps_path is None:
                # Haven't found it. Check in the path.
                lammps_path = self.executables_in_path()
                if lammps_path is not None:
                    if self.options.yes or self.ask_yes_no(
                        "Found LAMMPS executables in the PATH at "
                        f"{lammps_path}\n"
                        "Use them?",
                        default='yes'
                    ):
                        self.configuration.set_value(
                            'lammps-step', 'installation', 'user'
                        )
                        self.configuration.set_value(
                            'lammps-step', 'conda-environment', ''
                        )
                        self.configuration.set_value(
                            'lammps-step', 'modules', ''
                        )
                        self.configuration.save()
                        print("Using the LAMMPS executables at {lammps_path}")

            if lammps_path is None:
                # Can't find LAMMPS
                print(
                    "Cannot find LAMMPS executables. You will need to install "
                    "them."
                )
                if (
                    initial_installation is not None and
                    initial_installation != 'not installed'
                ):
                    if self.options.yes or self.ask_yes_no(
                        "The configuration file indicates that LAMMPS "
                        "is installed, but it can't be found.\n"
                        "Fix the configuration file?",
                        default='yes'
                    ):
                        self.configuration.set_value(
                            'lammps-step', 'installation', 'not installed'
                        )
                        self.configuration.set_value(
                            'lammps-step', 'lammps-path', ''
                        )
                        self.configuration.set_value(
                            'lammps-step', 'conda-environment', ''
                        )
                        self.configuration.set_value(
                            'lammps-step', 'modules', ''
                        )
                        self.configuration.save()
                        print(
                            "Since no LAMMPS executables were found, cleared "
                            "the configuration."
                        )
            else:
                print('The check completed successfully.')

    def check_configuration_file(self):
        """Checks that the lammps-step section is in the configuration file.
        """
        if not self.configuration.section_exists('lammps-step'):
            # Get the text of the data
            path = Path(pkg_resources.resource_filename(__name__, 'data/'))
            path = path / 'configuration.txt'
            text = path.read_text()

            # Add it to the configuration file and write to disk.
            self.configuration.add_section('lammps-step', text)
            self.configuration.save()

    def have_executables(self, path):
        """Check whether the executables are found at the given path.

        Parameters
        ----------
        path : pathlib.Path
            The directory to check.

        Returns
        -------
        bool
            True if at least one of the LAMMPS executables is found.
        """
        for executable in self.executables:
            tmp_path = path / executable
            if tmp_path.exists():
                self.logger.debug(f"Found executables in {path}")
                return True
        self.logger.debug(f"Did not find executables in {path}")
        return False

    def executables_in_path(self):
        """Check whether the executables are found in the PATH.

        Returns
        -------
        pathlib.Path
            The path where the executables are, or None.
        """
        path = None
        for executable in self.executables:
            path = shutil.which(executable)
            if path is not None:
                path = Path(path).expanduser().resolve()
                break
        return path

    def install(self):
        """Install LAMMPS using a Conda environment."""
        print(
            f"Installing Conda environment '{self.environment}'. This "
            "may take a minute or two."
        )
        self.conda.create_environment(
            self.environment_file, name=self.environment
        )
        # Update the configuration file.
        self.check_configuration_file()
        path = self.conda.path(self.environment) / 'bin'
        self.configuration.set_value('lammps-step', 'lammps-path', str(path))
        self.configuration.set_value('lammps-step', 'installation', 'conda')
        self.configuration.set_value(
            'lammps-step', 'conda-environment', self.environment
        )
        self.configuration.set_value('lammps-step', 'modules', '')
        self.configuration.save()
        print('Done!\n')

    def show(self):
        """Show the current installation status."""
        self.logger.debug('Entering show')

        # See if LAMMPS is already registered in the configuration file
        if not self.configuration.section_exists('lammps-step'):
            print(
                "There is no section in the configuration file for the "
                "LAMMPS step (lammps-step)."
            )
        data = self.configuration.get_values('lammps-step')

        # Keep track of where executables are
        serial = None
        mpi = None
        mpiexec = None

        # Is the path in the configuration file?
        if 'lammps-path' in data:
            conf_path = Path(data['lammps-path']).expanduser().resolve()
            if (conf_path / 'lmp_serial').exists():
                serial = conf_path / 'lmp_serial'
                serial_version = self.lammps_version(serial)
            if (conf_path / 'mpiexec').exists():
                mpiexec = conf_path / 'mpiexec'
            else:
                mpiexec = shutil.which('mpiexec')
            if (conf_path / 'lmp_mpi').exists():
                mpi = conf_path / 'lmp_mpi'
                if mpiexec is None:
                    mpi_version = 'unknown'
                else:
                    mpi_version = self.lammps_version(mpi, mpiexec)

            extra = f"from path {conf_path}."
            if 'installation' in data:
                installation = data['installation']
                if installation == 'conda':
                    if (
                        'conda-environment' in data and
                        data['conda-environment'] != ''
                    ):
                        extra = (
                            "from Conda environment "
                            f"{data['conda-environment']}."
                        )
                    else:
                        extra = "from an unknown Conda environment."
                elif installation == 'modules':
                    if 'modules' in data and data['modules'] != '':
                        extra = f"from module(s) {data['modules']}."
                    else:
                        extra = "from unknown modules."
                elif installation == 'user':
                    extra = f"from user-defined path {conf_path}."

            if serial is not None:
                if mpi is not None:
                    if serial_version == mpi_version:
                        print(
                            "LAMMPS serial and mpi executables, version "
                            f"'{serial_version}'"
                        )
                    else:
                        print(
                            "LAMMPS serial executable, version "
                            f"'{serial_version}', and mpi executable, "
                            f"version {mpi_version}"
                        )
                    print(extra)
                else:
                    print(
                        f"LAMMPS serial executable, version {serial_version}"
                    )
                    print(extra)
            elif mpi is not None:
                print(f"LAMMPS mpi executable, version {mpi_version}")
                print(extra)
            else:
                print("LAMMPS is not configured to run.")
        else:
            print("LAMMPS is not configured to run.")

        # Look in the PATH, but only record if not same as in the conf file
        tmp = shutil.which('lmp_serial')
        if tmp is not None:
            tmp = Path(tmp).expanduser().resolve()
            if serial is not None and serial != tmp:
                version = self.lammps_version(tmp)
                print(
                    f"Another serial executable of LAMMPS (version {version}) "
                    "is in the PATH:\n"
                    f"    {tmp}"
                )
        tmp = shutil.which('lmp_mpi')
        if tmp is not None:
            tmp = Path(tmp).expanduser().resolve()
            if mpi is not None and mpi != tmp:
                if mpiexec is None:
                    version = 'unknown'
                else:
                    version = self.lammps_version(tmp, mpiexec)
                print(
                    f"Another mpi executable of LAMMPS (version {version}) "
                    "is in the PATH:\n"
                    f"    {tmp}"
                )

    def uninstall(self):
        """Uninstall the LAMMPS Conda environment."""
        # See if LAMMPS is already registered in the configuration file
        data = self.configuration.get_values('lammps-step')
        if 'installation' in data and data['installation'] == 'conda':
            environment = self.environment
            if 'conda-environment' in data and data['conda-environment'] != '':
                environment = data['conda-environment']
            print(
                f"Uninstalling Conda environment '{environment}'. This "
                "may take a minute or two."
            )
            self.conda.remove_environment(environment)
            # Update the configuration file.
            self.configuration.set_value('lammps-step', 'lammps-path', '')
            self.configuration.set_value('lammps-step', 'modules', '')
            self.configuration.set_value(
                'lammps-step', 'installation', 'not installed'
            )
            self.configuration.set_value(
                'lammps-step', 'conda-environment', ''
            )
            self.configuration.save()
            print('Done!\n')

    def update(self):
        """Update the installation, if possible."""
        # See if LAMMPS is already registered in the configuration file
        data = self.configuration.get_values('lammps-step')
        if 'installation' in data and data['installation'] == 'conda':
            environment = self.environment
            if 'conda-environment' in data and data['conda-environment'] != '':
                environment = data['conda-environment']
            print(
                f"Updating Conda environment '{environment}'. This may "
                "take a minute or two."
            )
            self.conda.update_environment(
                self.environment_file, name=environment
            )
            # Update the configuration file, just in case.
            path = self.conda.path(environment) / 'bin'
            self.configuration.set_value(
                'lammps-step', 'lammps-path', str(path)
            )
            self.configuration.set_value(
                'lammps-step', 'installation', 'conda'
            )
            self.configuration.set_value(
                'lammps-step', 'conda-environment', environment
            )
            self.configuration.set_value('lammps-step', 'modules', '')
            self.configuration.save()
            print('Done!\n')
        else:
            print(
                "Unable to update LAMMPS because it was not installed using "
                "Conda"
            )

    def lammps_version(self, path, mpiexec=None):
        """Get the version of the LAMMPS executable.

        Parameters
        ----------
        path : pathlib.Path
            Path to the executable.

        Returns
        -------
        str
            The version reported by LAMMPS, or 'unknown'.
        """
        if mpiexec is not None:
            try:
                result = subprocess.run(
                    [mpiexec, str(path), '-log', 'none'],
                    stdin=subprocess.DEVNULL,
                    capture_output=True,
                    text=True
                )
            except Exception:
                version = 'unknown'
            else:
                version = 'unknown'
                lines = result.stdout.splitlines()
                if len(lines) > 0:
                    line = lines[0]
                    if line[0:8] == 'LAMMPS (':
                        version = line[8:].rstrip(')')
        else:
            try:
                result = subprocess.run(
                    [str(path), '-log', 'none'],
                    stdin=subprocess.DEVNULL,
                    capture_output=True,
                    text=True
                )
            except Exception:
                version = 'unknown'
            else:
                version = 'unknown'
                lines = result.stdout.splitlines()
                if len(lines) > 0:
                    line = lines[0]
                    if line[0:8] == 'LAMMPS (':
                        version = line[8:].rstrip(')')

        return version
