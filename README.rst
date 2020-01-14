===========
LAMMPS step
===========


.. image:: https://img.shields.io/travis/molssi-seamm/lammps_step.svg
           :target: https://travis-ci.org/molssi-seamm/lammps_step
	   :alt: Build Status

.. image:: https://codecov.io/gh/molssi-seamm/lammps_step/branch/master/graph/badge.svg
	   :target: https://codecov.io/gh/molssi-seamm/lammps_step
	   :alt: Code Coverage

.. image:: https://img.shields.io/lgtm/grade/python/g/molssi-seamm/lammps_step.svg?logo=lgtm&logoWidth=18
	   :target: https://lgtm.com/projects/g/molssi-seamm/lammps_step/context:python
	   :alt: Code Quality

.. image:: https://readthedocs.org/projects/lammps-step/badge/?version=latest
           :target: https://lammps-step.readthedocs.io/en/latest/?badge=latest
	   :alt: Documentation Status

.. image:: https://pyup.io/repos/github/molssi-seamm/lammps_step/shield.svg
	   :target: https://pyup.io/repos/github/molssi-seamm/lammps_step/
	   :alt: Updates for Dependencies

.. image:: https://img.shields.io/pypi/v/lammps_step.svg
           :target: https://pypi.python.org/pypi/lammps_step
	   :alt: PyPi VERSION

The LAMMPS step for a SEAMM flowchart


* Free software: BSD license
* Documentation: https://lammps-step.readthedocs.io.


Contributing
------------

This project uses the `git branching model`_ outlined by Vincent
Driessen.  Development is handled by branching feature branches from
the 'develop' branch:

Feature branches 

May branch off from:
  develop

Must merge back into:
  develop

Branch naming convention:
  anything except master, develop, release-*, or hotfix-*

Feature branches (or sometimes called topic branches) are used to
develop new features for the upcoming or a distant future
release. When starting development of a feature, the target release in
which this feature will be incorporated may well be unknown at that
point. The essence of a feature branch is that it exists as long as
the feature is in development, but will eventually be merged back into
develop (to definitely add the new feature to the upcoming release) or
discarded (in case of a disappointing experiment).

Feature branches typically exist in developer repos only, not in origin.

Creating a feature branch 
~~~~~~~~~~~~~~~~~~~~~~~~~~

When starting work on a new feature, branch off from the develop branch::

  $ git checkout -b myfeature develop
  Switched to a new branch "myfeature"

Incorporating a finished feature on develop 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finished features may be merged into the develop branch to definitely
add them to the upcoming release::

  $ git checkout develop
  Switched to branch 'develop'
  $ git merge --no-ff myfeature
  Updating ea1b82a..05e9557
  (Summary of changes)
  $ git branch -d myfeature
  Deleted branch myfeature (was 05e9557).
  $ git push origin develop

The --no-ff flag causes the merge to always create a new commit
object, even if the merge could be performed with a fast-forward. This
avoids losing information about the historical existence of a feature
branch and groups together all commits that together added the
feature.

Features
--------

* TODO

Credits
---------

This package was created with Cookiecutter_ and the `molssi-seamm/cookiecutter-seamm-plugin`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`molssi-seamm/cookiecutter-seamm-plugin`: https://github.com/molssi-seamm/cookiecutter-seamm-plugin
