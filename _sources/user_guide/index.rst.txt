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

The heat flux
=============

The Heat Flux sub-step, and the NVE, NVT and NPT sub-steps when asked to sample it,
give the heat flux for the thermal conductivity,
J = Σ\ :sub:`i` e\ :sub:`i` v\ :sub:`i` − Σ\ :sub:`i` S\ :sub:`i`·v\ :sub:`i`, from
each atom's energy e\ :sub:`i` and stress S\ :sub:`i`. The bonded terms and the
constraint forces of rigid bonds and angles use the centroid form of the stress, which
is the correct one for them (Surblys et al., J. Appl. Phys. 2021), unless "Use
centroid/stress/atom" is turned off; class2 forcefields always use the ordinary form.

The Heat Flux sub-step can also filter the heat flux ("Filter the heat flux"):

``none``
    The default, and the choice for liquids.
``per-atom averages``
    Subtracts each atom's average energy and stress, from a preliminary 6000-step run,
    from the heat flux (Pereverzev and Sewell, 2022). In a crystal, where every atom
    stays on its site, this removes a constant that only adds noise. In a liquid the
    atoms move, and it adds a term that makes the thermal conductivity too large:
    for SPC/E water by about 0.15-0.2 W/m/K, or 20%, which longer runs do not remove.
    Use it only for solids. How much it helps there has not been measured yet.

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
* the run's parameters, including how the electrostatics were handled (the kspace
  style, e.g. ``pppm 1.0e-5``, or ``none``), and where it came from (the job, step,
  structure and model).

The file is `H5MD <https://www.nongnu.org/h5md/>`_, an HDF5 format for molecular
dynamics; its layout is described in the documentation of ``seamm_md_analysis``
("The H5MD file"). The analysis steps, such as Diffusivity and Thermal Conductivity,
read it, and it can be read with h5py, ``seamm_md_analysis.open_md`` or MDAnalysis.
LAMMPS still writes its usual files during the run, so a run can be followed as
before and nothing is lost if it stops. Once ``md.h5md`` is written, the text files
of the sampled positions and velocities of the atoms and molecules
(``atomic_positions.dump_trj``, ``com_velocities.trj``, ...) are removed, since they
are in ``md.h5md`` and are several times larger: 3 ns of 1000 atoms of water, with
the velocities of the molecules every 10 fs, wrote 3.6 GB of them. "Keep the sampled
text files" on the trajectory tab keeps them. The other files, such as the
thermodynamic state and the trajectory, are kept.

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

Extra commands and extra properties
===================================

The NVE, NVT and NPT sub-steps have an **Extras** tab for measuring things that SEAMM
does not sample itself, without leaving SEAMM.

**Extra commands** are LAMMPS ``compute``, ``fix``, ``variable`` and ``dump``
commands for that sub-step only: they are put in after the sub-step's own fixes and
before its sampling, so the sampling can use them, and SEAMM removes everything they
define (``uncompute``, ``unfix``, ``undump``, ``variable ... delete``) when the
sub-step ends. Other commands, such as ``group`` or ``run``, are refused, as are IDs
that are numbers or that SEAMM uses itself: commands that change the system for the
rest of the run belong in a Custom sub-step, whose commands stay in force and whose
output SEAMM does not analyse.

**Extra properties** are quantities to sample, a row each:

* **Name:** letters, digits and ``_``, starting with a letter.
* **Expression:** a LAMMPS equal-style expression, such as ``c_msd[4]``, ``v_myvar``,
  a thermo keyword such as ``density`` or ``pxy``, or arithmetic of these.
* **Units:** a LAMMPS quantity (``energy``, ``pressure``, ``distance``, ``density``,
  ... or ``dimensionless``), recommended, since the units then follow the
  forcefield: a classical forcefield runs in LAMMPS's ``real`` units and a
  machine-learned one in ``metal`` units. Explicit units (``atm``, ``debye``) are
  accepted too, and are then taken as the units of the value LAMMPS reports.
* **Treatment:**

  - *average*, the default: the same analysis as the built-in properties -- the
    mean, its standard error, the equilibration and correlation times, and the
    graphs;
  - *autocorrelation*: also the normalized autocorrelation function with its 95 %
    confidence band, and the correlation time, its integral up to the first zero;
  - *Green–Kubo integral*: prefactor × ∫⟨δA(0) δA(t)⟩ dt, fitted as the Viscosity
    and Thermal Conductivity steps fit theirs. The prefactor is ``V/kT`` (for a
    stress, as for the viscosity), ``V/kT^2`` (for a flux, as for the thermal
    conductivity), ``1/(V kT)`` (for a current, as for an ionic conductivity), or a
    number with units. The error comes from four blocks of the run when they are long
    enough, otherwise from the fit, which is too small, and the output says which.
    The fit smooths the integral over 0.4 ps, so it suits correlations lasting about
    a picosecond or longer, and the property should be sampled often (a few fs for a
    stress).

* **Store:** also save the result (the mean, or the Green–Kubo integral) and its
  standard error as properties of the configuration, ``<name>#LAMMPS#<model>``.

The extra properties are sampled with the built-in ones, at the sub-step's sampling
interval, and appear in its table of properties and graphs; a table *Extra
Properties* gives the result of each treatment, with a graph of the autocorrelation
function (``<name>_acf.graph``) or the running integral (``<name>_gk.graph``). They
are written to ``md.h5md`` as ``observables/user/<name>``, with their units and
expressions, for analysis steps and scripts, and their results are in the step's
result "extra properties". Mistakes are shown in the tab as they are typed, and are
reported when the flowchart is checked, so a job stops before it runs.

Three examples:

* The mean squared displacement: the extra command ``compute msd all msd`` and the
  property ``msd``, expression ``c_msd[4]``, units ``distance``.
* An energy of part of the system: ``compute pe_w water pe/atom`` and ``compute
  sum_w water reduce sum c_pe_w``, with the property ``E_water``, expression
  ``c_sum_w``, units ``energy`` (``water`` being a group defined earlier in a Custom
  sub-step).
* The shear viscosity from one stress component: the property ``Pxy``, expression
  ``pxy``, units ``pressure``, treatment *Green–Kubo integral*, prefactor ``V/kT``,
  sampled every 5 fs. The result is in mPa·s. (The Viscosity step does this
  properly, averaging five components.)

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
