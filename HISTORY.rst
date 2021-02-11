=======
History
=======

2021.2.11 (11 February 2021)
----------------------------

* Updated the README file to give a better description.
* Updated the short description in setup.py to work with the new installer.
* Added keywords for better searchability.

2021.2.4.1 (4 February 2021)
----------------------------

* Internal patch to fix CI; no changes for users.

2021.2.4 (4 February 2021)
--------------------------

* Updated for compatibility with the new system classes in MolSystem
  2021.2.2 release.

2020.12.4 (4 December 2020)
---------------------------

* Internal: switching CI from TravisCI to GitHub Actions, and in the
  process moving documentation from ReadTheDocs to GitHub Pages where
  it is consolidated with the main SEAMM documentation.

2020.11.2 (2 November 2020)
---------------------------

* Updated to be compatible with the new command-line argument
  handling.

2020.10.13 (13 October 2020)
----------------------------

* Added capability to run MD until a set of user-selected properties
  are converged to requested accuracy.

2020.9.25 (25 September 2020)
-----------------------------

* Updated to be compatible with the new system classes in MolSystem.

2020.8.2.1 (2 August 2020)
--------------------------

* Bugfix: Fixed problem with nonbonds and charges just introduced.

2020.8.2 (2 August 2020)
------------------------

* Bugfix: Corrected the time units when using `metal` units with
  e.g. EAM potentials.

2020.8.1 (1 August 2020)
------------------------

* Added support for OpenKIM potentials.

0.9.4 (29 May 2020)
-------------------

* Cleaned up the output for the statistical analysis.

0.9.3 (29 May 2020)
-------------------

* Fixed issue with settings for bins in LAMMPS for small nonperiodic
  systems with just a few atoms.

0.9.2 (25 May 2020)
-------------------

* Switched to using PYMBAR for detecting covergence to equilibrium for
  MD runs. This is a more robust solution than the previous approach.

0.9.1 (24 May 2020)
-------------------

* Support for rigid water models, such as TIP-3P.

0.9 (15 April 2020)
-------------------

* Support for plots in the dashboard of properties from MD.
* Added option to produce local HTML for the above plots.

0.8.2 (2020-01-25)
------------------

* No significant changes in functionality.
* Incorporating changes to the SEAMM infrastructure, which simplify
  the code for plug-ins.
* Updating the Travis CI to handle incompatible changes in Travis, and
  to use Conda environments in all steps.

0.7.1 (18 December 2019)
------------------------

* Fixed problem with assigning charges to the system.

0.7.0 (17 December 2019)
------------------------

* General clean-up of code and output.

0.6 (8 September 2019)
----------------------

* Switched to ConfigArgParse for handling command-line arguments.
* Added the locations of LAMMPS executables to a configuration file
  for easier access.

0.5.2 (31 August 2019)
----------------------

* Defined the correct requirements for installation.

0.5.1 (30 August 2019)
----------------------

* Bugfix: corrected the name of the LAMMPS executable.
  
0.5.0 (30 August 2019)
----------------------

* Added ability to use serial or parallel versions of LAMMPS based on
  an environment variable.

0.3.1 (27 August 2019)
----------------------

* Added initial, fairly reasonable output.
  
0.2.1 (29 July 2019)
--------------------

* First release on PyPI.
