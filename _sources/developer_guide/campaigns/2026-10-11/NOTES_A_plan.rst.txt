The plan
========

Status: **agreed** (2026-10-11): naming and the MLFF launch decided (see
*Decisions*). Phase A's classical and ReaxFF part is implemented (2026-10-12, see
``NOTES_B_phase_A``); OpenKIM, phase B and phase C are not.

This builds on the June sketch, *Phase D -- Forcefield handling folded into
Model Chemistry* (``molssi-seamm.github.io``, campaign ``2026-06-22``,
``NOTES_D``), which it supersedes.

Where things stand
------------------

* **Model chemistries come from the code steps.** ``model_chemistry_step``'s
  ``discover_model_chemistries()`` walks the headless ``org.molssi.seamm``
  namespace and calls ``get_model_chemistry_options()`` on each step that has
  it: today mopac, orca, xtb, vasp and xnn (``xnn:MLFF@<checkpoint>``, one per
  checkpoint in the directories ``xnn.ini`` names). Nothing offers classical
  forcefields.
* **The Forcefield step** sets workspace variables that other steps read:

  - ``_forcefield`` -- a ``seamm_ff_util`` forcefield from a ``.frc`` file
    (``oplsaa``, ``pcff2018``, ``dreiding``, ``lithium_battery``,
    ``nacl_water``, the ``reaxff`` directory, and ``local:``/``personal:``
    files such as LigParGen's), or the strings ``"OpenKIM"`` / ``"PyTorch"``;
  - ``_OpenKIM_Potential`` -- the OpenKIM potential's name;
  - ``_pytorch_model`` -- an MLFF checkpoint's path;

  and its "assign forcefield to structure" task types the atoms.
* **Readers of those variables:** ``lammps_step`` (``lammps.py``,
  ``initialization.py``, ``heat_flux.py``), ``packmol_step`` (types the
  molecules it packs) and ``diffusivity_step`` (its model label).
* **LAMMPS already takes a Model Chemistry**, but only for QM-driven MD
  (``ff_form() == "MDI/QM"``): MOPAC or xTB as an MDI engine over TCP. An MLFF
  goes the Forcefield way: a ``.pt`` checkpoint makes LAMMPS use ``gpu-code``
  from ``lammps.ini``, which launches the xnn engine and LAMMPS as an MPI pair
  (``xnn mdi --ckpt {MODEL} ...`` ``:`` ``lmp -mdi ...``), optionally through
  ``mdi_bind.sh``.

The design
----------

1. **Classical forcefields become model chemistries, offered by
   forcefield_step.** It advertises one model chemistry per forcefield in each
   ``.frc`` file it can find -- the package's and the local and personal ones --
   e.g. ``FF@OPLS-AA``, ``FF@PCFF``, ``ReaxFF@<parameter set>``. It also
   provides the helpers the consumers need: build the ``seamm_ff_util``
   forcefield for a model chemistry, and assign its atom types to a
   configuration (see *Atom types*).

   forcefield_step stays a package **that must be installed** wherever LAMMPS
   or Packmol is: ``lammps_step`` and ``packmol_step`` require it. Like
   ``xnn_step`` it is registered only in the headless ``org.molssi.seamm``
   namespace (under its existing name, ``Forcefield``), which is how the Model
   Chemistry step finds providers and how old flowcharts containing a
   Forcefield step still run; its Tk registration is removed at the end of the
   campaign (see *Retiring the step*).

2. **OpenKIM potentials become model chemistries, offered by lammps_step**,
   since OpenKIM is LAMMPS-specific: ``IP@<potential>`` for the installed
   potentials.

3. **MLFFs come through the Model Chemistry too.** The xnn provider gains an
   MPI mode for ``get_mdi_engine_command()``: it returns the engine's argv
   (checkpoint, device, and xnn's options), and LAMMPS supplies the MPMD launch
   -- the GPU choice, ``mdi_bind.sh`` (2.1 stops the partner when either side
   fails), the MPI slots. This replaces the hard-wired
   ``xnn mdi --ckpt {MODEL} ...`` in ``lammps.ini``'s ``gpu-code``. Two needs
   from the mlff work come with it:

   * a checkpoint given by **path or variable**, not only those in xnn's model
     directories -- checkpoints change hourly during training;
   * **per-run engine options** on the Model Chemistry step (the D4 regime
     ``--dispersion "{regime: large}"``, ``--eeq-reuse``), which must never be
     global in ``gpu-code``.

   **Decided:** the MPI pair stays the launch for MLFFs (today's fast path:
   88 ms/step for 1,554 atoms on an A100); QM-MD keeps MDI over TCP. A
   performance check confirms the new path matches the old.

4. **LAMMPS takes its path from the model chemistry's type:** FF and ReaxFF
   as today (CPU, or Kokkos through ``gpu-cmd-args``); IP through
   ``pair kim``; MLFF the MPI pair; QM MDI over TCP. The timing record's class
   (``ff_form``, ``model``) comes from the model chemistry.

5. **Packmol and diffusivity read the Model Chemistry:** Packmol types its
   molecules through forcefield_step's helper when the model chemistry is a
   classical forcefield; diffusivity takes its model label from it.

6. **Backward compatibility.** With no ``_model_chemistry`` the consumers keep
   reading the Forcefield step's variables, unchanged, for some time.

7. **GUI** (prevent invalid combinations): a Model Chemistry feeding LAMMPS
   lists only what LAMMPS can run -- FF, ReaxFF, IP, MLFF and MDI-capable QM --
   and shows the MLFF options (checkpoint path, engine options) only for an
   MLFF. The run-time checks stay as the backstop.

Type names
----------

Decided 2026-10-11, and recorded on the site's model-chemistry naming page
(``docs/background/model_chemistry_naming.rst``):

* ``SQM`` -- semiempirical (``MOPAC:SQM@PM7``, ``xTB:SQM@GFN2-xTB``);
* ``DFT`` -- Kohn-Sham DFT, including double hybrids;
* ``QC`` -- wavefunction methods: ``QC@HF``, ``QC@MP2``, CCSD(T), DLPNO;
* ``FF`` -- classical (valence) forcefields: ``FF@OPLS-AA``, ``FF@PCFF``;
* ``ReaxFF`` -- reactive forcefields, the method being the parameter set;
* ``IP`` -- interatomic potentials (OpenKIM): ``IP@<OpenKIM ID>``;
* ``MLFF`` -- machine-learned forcefields: ``xnn:MLFF@<model>``.

There is no catch-all type. The earlier names -- ``VFF``, ``rxnFF`` /
``ReactiveFF``, ``ML`` -- are not used. One change reaches existing code:
**orca_step advertises its wavefunction methods as ``HF@...`` and ``MP2@...``
today and must change to ``QC@HF`` / ``QC@MP2``.** Model-chemistry strings in
existing flowcharts (``ORCA:HF@...``, ``ORCA:MP2@...``) are rewritten by the
converter below, and orca_step accepts the old form, with a warning, until the
conversions are done.

Atom types
----------

The June sketch left open when atom types are assigned, since typing changes the
configuration and has no QM analogue. The answer: **the consumer types the
atoms when it needs them**, through forcefield_step's helper -- LAMMPS while
writing its input, Packmol for each molecule before packing -- with the same
rule the consumers already use: type the structure if it has no types for this
forcefield (``atom_types_<forcefield>``), keep existing ones otherwise, so that
re-running is harmless and hand-assigned types survive. Choosing a model
chemistry types nothing by itself, and the Forcefield step's separate "assign
forcefield to structure" task is no longer needed.

Retiring the step
-----------------

Removing the Tk registration removes the step from the editor's palette, but
the editor also needs a Tk plug-in for every step of a flowchart it opens: an
old flowchart with a Forcefield step could still be *run* headless, but not
*opened*. So, in order:

1. Release the provider and the LAMMPS/Packmol/diffusivity changes with the Tk
   interface still in place, the Forcefield step's description marking it
   deprecated and the step warning when it runs.
2. Release a flowchart converter (as for the ``1-end`` index syntax) that
   rewrites a Forcefield step as a Model Chemistry step, and ``HF@`` /
   ``MP2@`` model chemistries as ``QC@HF`` / ``QC@MP2``, and run it over the
   flowcharts we know of: the installations' (``~/SEAMM``, ``~/SEAMM_DEV``,
   ChemAI, ARC, MolSSI10), Dropbox and the production directories on ARC and
   ChemAI.
3. Then remove the Tk registration. A flowchart missed by the converter still
   runs; it has to be converted before it can be edited.

The compatibility code in the consumers (6) can go a release or two later.

Phases
------

A. **Classical forcefields and OpenKIM through the Model Chemistry** -- design
   1, 2, 4 (non-MLFF), 5. Lowest risk: the classical code paths are unchanged,
   only where they get the forcefield from. Releases: forcefield_step,
   lammps_step, packmol_step, diffusivity_step, model_chemistry_step if its GUI
   needs the new types.
B. **MLFFs through the Model Chemistry** -- design 3 and 4 (MLFF): xnn_step
   (MPI mode, path/variable checkpoints, engine options) and lammps_step, with
   the performance check against today's ``gpu-code`` path.
C. **Retirement** -- design 7 (GUI narrowing), the deprecation, the converter
   and the conversions, then the Tk registration's removal.

Tests
-----

* Run-path tests (``seamm_exec.testing``) for each type: a classical forcefield
  (OPLS-AA water), ReaxFF, an OpenKIM potential, an MLFF (on a GPU node), and
  QM-MD as today -- each also run through the old Forcefield step to show the
  two give the same LAMMPS input.
* Headless GUI tests (``reset_dialog`` over every choice) for the Model
  Chemistry step's narrowed lists and the MLFF-only options.
* The converter: every Forcefield variant in, a valid Model Chemistry
  flowchart out; ``seamm-flowchart validate`` on all converted flowcharts.

Decisions
---------

Made by Paul on 2026-10-11:

1. **Naming:** ``FF@OPLS-AA`` style (amended 2026-10-12: the forcefields keep the
   names their files give them, ``FF@oplsaa+``, ``FF@pcff``), and the type vocabulary above (``ReaxFF``,
   ``MLFF``, ``QC``).
2. **The MLFF launch:** the MPI pair through ``mdi_bind.sh``.
3. **The provider:** forcefield_step, installed and headless-only like
   ``xnn_step``, its Tk registration removed after the conversions.

Still open:

4. **When the Tk registration goes**, i.e. how long after the converter.
5. **When the compatibility code goes.**

Related
-------

* ``2026-09-15``: how an MLFF run reaches its engine today (``gpu-code``,
  engine selection, GPUs and tasks).
* ``mdi_bind.sh`` 2.1 (lammps-step 2026.10.11.1): each side of an MDI pair
  stops the other when it fails; tested in TinkerCliffs jobs 7944051-54.
