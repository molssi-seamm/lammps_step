=======
History
=======
2025.5.7 -- Bugfix: error if not saving configuration in minimization
   * The code crashed in the minimization step if the structure was discarded.
     
2025.4.9 -- Adding support for ReaxFF forcefields.
   * Added support for ReaxFF forcefields.
   * Added control for how charges are handled, and options for the various charge
     equilibration methods in LAMMPS.
   * For ReaxFF, include the approximate enthalpy of formation, if available, and added
     it as a property.
   * Improved output for most substeps and included the forcefield used in the output.
     
2025.3.17.1 -- Bugfix: Incorrect handling of Dreiding hydrogen bonds
   * Fixed an issue with the Dreiding forcefield where having no hydrogen bonds in the
     system caused a crash.

2025.3.17 -- Bugfix: LAMMPS installer did not have correct conda environment

2025.3.16 -- Added Dreiding forcefield and support for structure handling
   * Added support for the Dreiding forcefield.
   * Added standard support for structure handling
   * Added the RMSD between initial and final structues in minimization, and
     expanded the possible minimizers and results from minimization.
     
2025.2.7 -- Installation issue
   * Encountered problems with the OpenMPI version of LAMMPS and with the installation
     of OpenKIM, so fixed the CONDA environment file to correct.
   * Changed logging at the INFO level to DEBUG to reduce the output when using geomeTRIC.

2024.8.22 -- Bugfix: error in number of cores.
   * Fixed an error if options for number of cores were strings, not numbers.
     
2024.7.31 -- Bugfix and improvements for parallel runs
   * There was an issue saving the pressure and volume as properties due to a mismatch
     in the names used. This has been corrected.
   * Improved how the code determines the number of processors to use for parallel runs,
     giving salience to the number of atoms in the system, but limiting to the LAMMPS and
     global limits on numbers of cores as well as the hardware available.
     
2024.7.25 -- Bugfix and improvements
   * Bugfix: Fixed issue with the initial seamm.ini file, created if it is missing from
     the installation.
   * Added the ability to set the number of points in the trajectories rather than the
     sampling rate.
   * Added diagnositic information and timings to available results.
     
2024.7.21.1 -- Minor internal change for GUI
   * Switched to new functionality in the SEAMM widgets to simplify the layout of the
     trajectories panel.
     
2024.7.21 -- Improved handling of trajectories and results
   * Improved control over the trajectories of positions, velocities, etc. to allow the
     user to give the number of points in the trajectory rather than the time interval
     of samples
   * Added volume of the cell to properties, and the cell lengths, density, and volume
     for NVE and NVT, where those parameters don't vary but are nonetheless useful in
     subsequent analysis.
     
2024.6.28.1 -- Internal release to fix issue making Docker image.

2024.6.28 -- Added energy and forces to properties
   * Added ability to get the energy and forces from single point calculations to supprt
     e.g. energy scans.
   * Fixed issue with assigning atoms types if they have not been assigned but are None
   * Updated for change in the order of units in the new version of pint
   * Improved analsys based on the output file.
     
2024.3.22 -- Corrected issue with e.g. heat flux calculations
   * Corrected an issue running LAMMPS via Python, introduced in the new scheme for
     executing. It ignored parallelism.
     
2024.3.21 -- Switched to new installation scheme
   * Fully support ~/SEAMM/lammps.ini
   * Updated to new installer
   * Support for Conda and Docker installation.
     
2024.1.18 -- Restructured to support running in containers.

2023.11.7 -- Bugfix: properties that are constant
   * A property, such as the total energy, can be a constant over an MD run due to
     precision of the trajectory. This caused errors because the autocorrelation
     function is not defined. These cases are now detected and the ACF not calculated
     for them.
     
2023.11.6 -- Bugfix in thermal conductivity
   * Due to change in input file name.

2023.9.6 -- Corrected issues with final coordinates; added velocities
   * There was a problem with getting the final coordinates from a dump file. 
   * Added saving and reusing velocities so now a second LAMMPS step will by default use
     the velocities from the previous step, which is what you would expect.

2023.8.31 -- Bugfix: not reading structure correctly after dynamics

2023.8.27 -- Added support for tabulated angle potentials.
   * Support for the CL&P-OPLSAA potential for octahedral PF6-
     
2023.8.21 -- Bugfix: x-axes length on graphics incorrect
   * The size of the x-axes of the trajectory graphs were wrong, often much too large,
     compressing the actual data near the beginning of the graph.
   * Fixed an issue with systems with no non-bonds.

2023.6.17 -- Bugfix: more centroid/stress/atom issues
   * Avoided using centroid/stress/atom for heat flux in standard NVE, NVT, ... dynamics
     with Class 2 forcefield.
   * Added option to not use centroid/stress/atom for any forcefield.
2023.6.16 -- Heat flux with PCFF
   * centroid/stress/atom does not work with Class 2 forcefields, so don't use for PCFF.
2023.5.29 -- Self diffusion and other improvements
   * Added trajectory panel to support diffusion, viscosity and simple thermal
     conductivity.
   * Added support for separate GPU versions of LAMMPS.
   * Added support for command-line arguments to LAMMPS, mainly used for accelerators.
   * Added support for using modules.

2023.4.24 -- Support for thermal conductivity
   * Internal changes to support thermal conductivity with its embedded flowchart.
   * Added the heat flux substep.
   * Now delete output and reference files when rerunning, so the output is clean.
   * Internal changes to support running LAMMPS from a Python driver.
   * Corrected units of properties returned from LAMMPS when e.g. metal units used.
   * Added support for Buckingham potentials
   * Fixed issues with and cleaned up the use of hybrid types for bonds, angles, ....
   * Fixed issues with the alignment of some of the widgets in the GUI.
     
2023.4.9 -- Hid the warning from pymbar
   * Importing pymbar timeseries writes a warning to the terminal about its proper
     usage. SEAMM already handles the warned case, so the message is simply confusing to
     users and hence this release hides it.
     
2023.4.6 -- Better forcefield handling.
   * Added correct molecule numbers for valence forcefields.
   * Correctly handle ReaxFF from OpenKim
   * Updated for some minor changes in OpenKim

2023.2.6 -- Added handling of OPLS-AA forcefield
   * Added handling of the OPLS-AA forcefield
   * Moved documentation to new MolSSI theme and diátaxis layout
   * Cleaned up internale dependencies and workflows for GitHub

2022.10.31 -- Bugfix: properties with commas
  Properties with commas in their name in data/properties.csv need to have quotes to
  protect the property name!

2022.10.27 -- Added properties
  * Added properties to be saved in the database.
  * Updated calls to `pymbar` because the names of methods were changed.
  * Add the missing references for `pymbar`

2021.2.11 (11 February 2021)
  * Updated the README file to give a better description.
  * Updated the short description in setup.py to work with the new installer.
  * Added keywords for better searchability.

2021.2.4.1 (4 February 2021)
  Internal patch to fix CI; no changes for users.

2021.2.4 (4 February 2021)
  Updated for compatibility with the new system classes in MolSystem
  2021.2.2 release.

2020.12.4 (4 December 2020)
  Internal: switching CI from TravisCI to GitHub Actions, and in the
  process moving documentation from ReadTheDocs to GitHub Pages where
  it is consolidated with the main SEAMM documentation.

2020.11.2 (2 November 2020)
  Updated to be compatible with the new command-line argument
  handling.

2020.10.13 (13 October 2020)
  Added capability to run MD until a set of user-selected properties
  are converged to requested accuracy.

2020.9.25 (25 September 2020)
  Updated to be compatible with the new system classes in MolSystem.

2020.8.2.1 (2 August 2020)
  Bugfix: Fixed problem with nonbonds and charges just introduced.

2020.8.2 (2 August 2020)
  Bugfix: Corrected the time units when using `metal` units with
  e.g. EAM potentials.

2020.8.1 (1 August 2020)
  Added support for OpenKIM potentials.

0.9.4 (29 May 2020)
  Cleaned up the output for the statistical analysis.

0.9.3 (29 May 2020)
  Fixed issue with settings for bins in LAMMPS for small nonperiodic
  systems with just a few atoms.

0.9.2 (25 May 2020)
  Switched to using PYMBAR for detecting covergence to equilibrium for
  MD runs. This is a more robust solution than the previous approach.

0.9.1 (24 May 2020)
  Support for rigid water models, such as TIP-3P.

0.9 (15 April 2020)
  Support for plots in the dashboard of properties from MD.
  Added option to produce local HTML for the above plots.

0.8.2 (2020-01-25)
  * No significant changes in functionality.
  * Incorporating changes to the SEAMM infrastructure, which simplify
    the code for plug-ins.
  * Updating the Travis CI to handle incompatible changes in Travis, and
    to use Conda environments in all steps.

0.7.1 (18 December 2019)
  Fixed problem with assigning charges to the system.

0.7.0 (17 December 2019)
  General clean-up of code and output.

0.6 (8 September 2019)
  * Switched to ConfigArgParse for handling command-line arguments.
  * Added the locations of LAMMPS executables to a configuration file
    for easier access.

0.5.2 (31 August 2019)
  Defined the correct requirements for installation.

0.5.1 (30 August 2019)
  Bugfix: corrected the name of the LAMMPS executable.
  
0.5.0 (30 August 2019)
  Added ability to use serial or parallel versions of LAMMPS based on
  an environment variable.

0.3.1 (27 August 2019)
  Added initial, fairly reasonable output.
  
0.2.1 (29 July 2019)
  First release on PyPI.
