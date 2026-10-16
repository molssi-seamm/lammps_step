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

The forcefield
==============

LAMMPS uses the forcefield or model chemistry chosen before it in the flowchart:

* a **Model Chemistry** step with a forcefield -- ``FF@oplsaa``, ``FF@pcff``,
  ``ReaxFF@CHO_cho_2008`` and the others the Forcefield plug-in offers. LAMMPS builds the
  forcefield and types the atoms that do not yet have types for it;
* a **Model Chemistry** step with a machine-learned forcefield, ``xnn:MLFF@<model>``:
  its engine runs alongside LAMMPS, each through ``mdi_bind.sh``, with the Model
  Chemistry step's checkpoint and engine options. It runs on an idle GPU, waiting for
  one if all are busy (``--gpu-wait-time`` in seconds, unlimited by default). When the
  job can never have a GPU -- there is none on the machine, or the scheduler gave it
  none -- it runs on the CPU at once, with the job's cores less one for LAMMPS; on the
  CPU it is much slower except for small systems. Files the engine writes, such as a
  committee's ``committee_spreads.csv``, are kept in the step's directory. How the
  two talk depends on where the engine is installed:

  - in LAMMPS's own environment (``xnn.ini`` and ``lammps.ini`` name the same one, as
    in a combined environment holding both), they run as an MPI pair, and the job
    needs two MPI tasks;
  - in an environment of its own -- a default installation has ``seamm-lammps`` and
    ``seamm-xnn``, whose MPI libraries differ -- they talk over TCP on the same
    machine, starting the engine with ``conda run``.

  ``mlff-transport`` in ``lammps.ini`` (``auto``, ``mpi`` or ``tcp``) overrides the
  choice, and ``mlff-code`` and ``mlff-tcp-code`` change the launch (see the comments
  in ``lammps.ini``). ``gpu-code`` is not used. If either side fails, ``mdi_bind.sh``
  stops the other, so a failed run ends instead of waiting. On Linux, TCP is
  currently slower than MPI: a message can wait up to about 40 ms, which for a fast
  model added a third to the run time. A fix is proposed in the MDI library;
* a **Model Chemistry** step with a quantum method, such as ``MOPAC:SQM@PM6-ORG``:
  LAMMPS runs the dynamics and the program computes the energy and forces over MDI;
* or a **Forcefield** step, as before. A PyTorch model from the Forcefield step still
  needs a GPU on Linux.

If a flowchart has both, the model chemistry is used.

Rigid bonds and angles
======================

Some models are rigid, such as the SPC/E and TIP3P water models in the forcefields
``oplsaa+-spc_e``, ``oplsaa+-tip3p`` and ``nacl_water``. The forcefield says which bonds
and angles are rigid, and LAMMPS holds them at their values with SHAKE (``fix shake``)
during dynamics, so a timestep of 2 fs can be used. LAMMPS cannot minimize with SHAKE,
so a rigid bond or angle is also a stiff harmonic term, 1000 kcal/mol/Å\ :sup:`2` for a
bond and 100 kcal/mol/radian\ :sup:`2` for an angle. It holds the geometry close to
the rigid one in a minimization and is exactly zero when SHAKE holds it in dynamics.

A forcefield that is not rigid, such as ``oplsaa+`` with its flexible SPC water, is
used as it is. The Initialization sub-step's "Fix X-H bond lengths" holds bonds to
hydrogen fixed with SHAKE in any forcefield.

SHAKE is used rather than RATTLE: LAMMPS's RATTLE gave rigid SPC/E water a pressure
about 850 atm too high and a center of mass that drifted.

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


The MD data for analysis: ``md.h5md``
=====================================

At the end of each NVE, NVT, NPT or Heat Flux sub-step, the step writes everything
the run sampled into one file, ``md.h5md``, in the sub-step's directory:

* the time series (temperature, pressure, the energies, the full pressure tensor, the
  volume and density of an NPT run, the heat flux, ...), each on its own time axis
  and in the units LAMMPS produced;
* the atoms' positions, velocities and forces, and the centres of mass of the
  molecules, as far as they were sampled, with the box of each frame;
* the run's parameters and where it came from (the job, step, structure and model).

The file is `H5MD <https://www.nongnu.org/h5md/>`_, an HDF5 format for molecular
dynamics; its layout is described in the documentation of ``seamm_md_analysis``
("The H5MD file"). The analysis steps, such as Diffusivity and Thermal Conductivity,
read it, and it can be read with h5py, ``seamm_md_analysis.open_md`` or MDAnalysis.
LAMMPS still writes its usual files during the run, and they are kept, so a run can
be followed as before and nothing is lost if it stops.

The full pressure tensor is always written with the thermodynamic state, at its
sampling interval. For the viscosity it is needed much more often, every few
femtoseconds: the NVE, NVT and NPT sub-steps can sample it on its own ("Sample the
pressure tensor" on the trajectory tab), and that series is then the one written. It
replaces the earlier choice "Sample the shear stress", which saved only the three
off-diagonal components; flowcharts that use it now sample the full tensor.

Each run is also recorded in the variable ``_md_runs``, which analysis steps later in
the flowchart use to find the data.

Analysis steps check the MD before it runs. When a flowchart is checked (in the
editor before submitting, with ``seamm-flowchart validate``, and when a job starts), an
analysis step that reads the MD before it is told what each NVE, NVT, NPT or Heat Flux
sub-step will sample, and reports anything missing or too coarse. An analysis step
that runs its own MD turns on what is missing itself, and says so in the output. And
just before LAMMPS runs, the analysis steps that will read its MD check the structure,
so that, for example, Diffusivity stops the job before the MD if the system is not
periodic.

Non-periodic systems
====================

The step runs molecules and clusters without a periodic cell, with a shrink-wrapped
box. The pressure, the stress and the heat flux are divided by the volume of the cell,
so for such a system:

* the pressure is not reported, and the stress is not sampled;
* asking for the heat flux or the pressure tensor in an NVE or NVT sub-step gives a
  warning in the output, and they are not sampled;
* the Heat Flux sub-step and NPT dynamics stop with an error before LAMMPS starts.

Index
=====

* :ref:`genindex`
