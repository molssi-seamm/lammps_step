.. image:: https://img.shields.io/github/issues-pr-raw/molssi-seamm/lammps_step
   :target: https://github.com/molssi-seamm/lammps_step/pulls
   :alt: GitHub pull requests

.. image:: https://github.com/molssi-seamm/lammps_step/workflows/CI/badge.svg
   :target: https://github.com/molssi-seamm/lammps_step/actions
   :alt: Build Status

.. image:: https://codecov.io/gh/molssi-seamm/lammps_step/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/molssi-seamm/lammps_step
   :alt: Code Coverage

.. image:: https://github.com/molssi-seamm/lammps_step/workflows/CodeQL/badge.svg
   :target: https://github.com/molssi-seamm/lammps_step/security/code-scanning
   :alt: Code Quality

.. image:: https://github.com/molssi-seamm/lammps_step/workflows/Release/badge.svg
   :target: https://molssi-seamm.github.io/lammps_step/index.html
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/lammps_step.svg
   :target: https://pypi.python.org/pypi/lammps_step
   :alt: PyPi VERSION

====================
SEAMM LAMMPS Plug-in
====================

A SEAMM plug-in for LAMMPS, a forcefield-based molecular dynamics code.

This plug-in provides a graphical user interface (GUI) for setting up
complex simulations using LAMMPS. It uses a sub-flowchart that
provides steps such as constant pressure and temperature (NPT)
dynamics which give access to the functionality in LAMMPS in a more
consistent and understandable way than the inscrutable fixes that
LAMMPS uses.

These sub-flowcharts mirror the main flowchart in form and function
and can use the same variables such as temperature and pressure that
are accessible anywhere in the flowcharts. This allows "programming" a
LAMMPS workflow in the same familiar way that SEAMM uses to represent
the overall workflow.

* Free software: BSD license
* Documentation: https://molssi-seamm.github.io/lammps_step/index.html
* Code: https://github.com/molssi-seamm/lammps_step


Features
--------

* Use of any forcefield supported by the forcefield plug-in:

  - PCFF
  - OpenKIM: EAM, MEAM, LJ, ReaxFF

* Molecular statics: minimization
* Molecular dynamics: NVE, NVT, and NPT with any of the approaches
  supported in LAMMPS
* Automatic statistical analysis of averages from MD

  - Detection of equilibration
  - Mean and standard error of the mean for the sampling after
    equilibration
  - Autocorrelation function and time
  - Statistical inefficiency
  - Plotting of results in the Dashboard

* Using property values to drive MD. Rather than running MD for a
  length of time, automatically run long enough to determine a set of
  properties within given error bars.

Acknowledgements
----------------

This package was created with Cookiecutter_ and the `molssi-seamm/cookiecutter-seamm-plugin`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`molssi-seamm/cookiecutter-seamm-plugin`: https://github.com/molssi-seamm/cookiecutter-seamm-plugin

Developed by the Molecular Sciences Software Institute (MolSSI_),
which receives funding from the `National Science Foundation`_ under
award ACI-1547580

.. _MolSSI: https://www.molssi.org
.. _`National Science Foundation`: https://www.nsf.gov
