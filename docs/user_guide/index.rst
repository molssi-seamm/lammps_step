.. _user-guide:

**********
User Guide
**********
The LAMMPS plug-in ...

..
   The following sections cover accessing and controlling this functionality.

   .. toctree::
      :maxdepth: 2
      :titlesonly:

Settings that depend on each other
==================================

Which settings apply depends on others: the k-space settings only with a k-space method,
the charge-equilibration settings only for methods that equilibrate charges (including a
forcefield's default, usually QEq for ReaxFF), and the cell's shear only for a solid.
The step's dialogs show only the settings that apply with the current choices, and the
same rules are used when a flowchart is built or edited without the editor
(``seamm-flowchart`` or SEAMM's MCP server): a setting that would have no effect is
refused, with the reason, and a value that contradicts another is refused too. See
"Flowcharts without the editor" in SEAMM's user guide.

* The damping times (``Pdamp`` and the stress damping times) apply to both barostats:
  LAMMPS's Berendsen barostat needs them too.
* **Allow shear** applies only to a solid with the Nose-Hoover barostat, since the
  Berendsen barostat cannot control a sheared (triclinic) cell. LAMMPS uses one box for
  all the sub-steps in a LAMMPS step, so a Berendsen NPT step is refused before LAMMPS
  starts if the box must be triclinic (a later sub-step allows shear, or the cell is not
  orthorhombic); put such steps in separate LAMMPS steps.
* When LAMMPS itself stops with an error (for example "Lost atoms"), the step stops
  with a message giving the sub-step, LAMMPS's message, the last command, advice for
  common errors, and where ``log.lammps`` is.


Index
=====

* :ref:`genindex`
