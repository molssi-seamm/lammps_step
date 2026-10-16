=======
History
=======
2026.10.16.4 -- The electrostatics method in the MD data
    * md.h5md records the run's kspace style (e.g. 'pppm 1.0e-5', or 'none'), for
      analyses whose formulas depend on it: the dielectric constant assumes Ewald or
      PPPM with conducting boundaries.

2026.10.16.3 -- Bugfix: rigid water was not held rigid in MD
    * In 2026.10.16.1 and 2026.10.16.2, NVE, NVT and NPT dynamics dropped the SHAKE
      constraints, so rigid models such as SPC/E water ran with flexible bonds and
      angles (the stiff harmonic terms that hold them in minimizations). Results of
      MD with rigid models from those two versions should be rerun; the shear
      viscosity of SPC/E water came out about twice too large. 2026.10.16 and earlier
      are not affected.

2026.10.16.2 -- Bugfix: LAMMPS said it used more cores than it did; summaries are averages
    * In a batch job, LAMMPS can only use as many MPI processes as the job has tasks.
      A job with one task and several CPUs per task ran LAMMPS on one process, but
      the output said it was using all the cores the system's size allowed. It now
      gives the number of processes LAMMPS actually runs on, and when the job's
      allocation is what limits it, explains that MPI processes need tasks (e.g.
      --ntasks) rather than CPUs per task.
    * The summary of an NVE, NVT, NPT or Heat Flux run (<ensemble>_summary.trj), ten
      lines over the run, now gives on each line the average over that tenth of the
      run, sampled every 10 steps, rather than the values at that one step.

2026.10.16.1 -- Extra commands and extra properties in the MD sub-steps
    * The NVE, NVT and NPT sub-steps have an Extras tab. Extra commands are LAMMPS
      compute, fix, variable and dump commands for that sub-step only, removed when
      it ends. Extra properties are LAMMPS expressions, each with a name and units,
      sampled and analysed with the built-in properties, as an average, an
      autocorrelation function with its correlation time, or a Green-Kubo integral
      with a chosen prefactor. They are printed, graphed, written to md.h5md as
      observables/user/<name>, and can be stored as properties of the configuration.
    * Mistakes in the extra commands or properties are shown in the dialog as they
      are typed, and stop a job before it runs.

2026.10.16 -- Rigid models from the forcefield, and the heat flux of liquids corrected
    * Which bonds and angles are rigid now comes from the forcefield, e.g. the SPC/E
      and TIP3P water in the new forcefields oplsaa+-spc_e and oplsaa+-tip3p, and
      SHAKE holds exactly those. Before, any water was made rigid, whatever the
      forcefield; a flexible water model such as oplsaa+'s SPC is now flexible.
    * A rigid bond or angle is also a stiff harmonic term, so that minimizations,
      which cannot use SHAKE, keep the geometry.
    * Bugfix: the Heat Flux sub-step turned SHAKE on only after the runs that settle
      the system and average each atom's energy and stress, so rigid water ran
      unconstrained there and the run failed after a few femtoseconds.
    * Bugfix: when LAMMPS stopped with an error in the Heat Flux sub-step, the job
      hung instead of stopping.
    * Bugfix: the Heat Flux sub-step subtracted each atom's average energy and stress
      from the heat flux. That is right for a crystal, but in a liquid, where the
      atoms move, it made the thermal conductivity too large, by about 20% for SPC/E
      water. By default the heat flux is no longer filtered; the new "Filter the heat
      flux" can still subtract the per-atom averages, for solids. Without the filter,
      Heat Flux is a plain LAMMPS run with no preliminary runs.
    * Bugfix: the heat flux of the NVE, NVT and NPT sub-steps left out the forces
      holding rigid bonds and angles, and in NVT and NPT counted the convective term
      twice. Both are corrected, and the computes they define are removed afterwards,
      so a later sub-step can sample the heat flux again.

2026.10.15.4 -- Bugfix: rigid water gave the wrong pressure
    * Rigid water molecules, and fixed X-H bond lengths, were held with LAMMPS's
      RATTLE, which for rigid SPC/E water gave a pressure about 850 atm too high, an
      anisotropic pressure tensor and a center of mass that drifted ever faster, so
      the viscosity and diffusion coefficients were wrong. SHAKE is now used, with a
      tolerance of 1e-6, which gives the right pressure and conserves momentum.

2026.10.15.3 -- Sampling the full pressure tensor, for the viscosity
    * The NVE, NVT and NPT sub-steps can sample the full pressure tensor often, as
      the viscosity needs, with "Sample the pressure tensor". It replaces "Sample the
      shear stress", which saved only the three off-diagonal components; flowcharts
      that used it now sample the full tensor.
    * The sampled pressure tensor is written to ``md.h5md`` as ``pressure_tensor``,
      in place of the coarser one from the thermodynamic state, and analysis steps
      such as Viscosity can ask for it to be turned on.
    * Bugfix: "Keep the internal geometry of water molecules fixed", on by default,
      did nothing since 2020, so water was always flexible. Water molecules (an O
      bonded to two H's) are now held rigid with RATTLE.
    * Bugfix: "Fix X-H bond lengths" in the Initialization sub-step looked up the
      elements of the wrong atoms, so it stopped with an error or constrained the
      wrong bonds.
    * Bugfix: the nacl_water forcefield, with SPC/E water and ions, could not be used:
      LAMMPS stopped with "Special bonds not implemented".
    * Requires seamm-md-analysis 2026.10.15.2.

2026.10.15.2 -- Analysis steps can check and adjust the MD before it runs
    * The NVE, NVT, NPT and Heat Flux sub-steps tell analysis steps what each run will
      sample, before it runs, and can turn on sampling an analysis needs.
    * Just before LAMMPS runs, the analysis steps that will read its MD check the
      structure, so that, for example, Diffusivity stops the job before the MD if the
      system is not periodic.

2026.10.15.1 -- The MD data in H5MD files; non-periodic systems
    * At the end of each NVE, NVT, NPT and Heat Flux sub-step, everything the run sampled
      is written to ``md.h5md`` in the sub-step's directory: the time series, the full
      pressure tensor, the trajectory and centres of mass with the box of each frame,
      and the run's parameters and provenance. The analysis steps read it. LAMMPS's own
      files are still written and kept.
    * Each MD run is recorded in the variable ``_md_runs``.
    * Bugfix: sampling the trajectory of a non-periodic system (a molecule or cluster)
      failed with "Variable name sxx ... does not exist". For a non-periodic system the
      pressure is no longer reported, the heat flux and shear stress are not sampled
      (with a warning), and the Heat Flux sub-step stops with a clear message.

2026.10.15 -- Bugfix: the heat flux was far too large
    * The Heat Flux sub-step removes each atom's average energy and stress from the
      heat flux, to reduce its noise. The stress was handed to LAMMPS in the wrong
      units and with the wrong sign, which changed it by 68,568 times its average
      instead of removing it. For liquid argon the heat flux was about 10,000 times
      too large and the thermal conductivity about 4x10^7 W/m/K; it is now about 0.1
      W/m/K, as measured. Every thermal conductivity computed with the Heat Flux
      sub-step for a system without bonds (e.g. with class II forcefields such as
      PCFF) was affected.
    * For molecular systems the average stress was not removed at all, and the
      convective part of the heat flux was counted twice. Both are corrected, also
      for the heat flux sampled during NVE dynamics.
    * The Heat Flux sub-step sampled the heat flux 1000 times per run, e.g. every 20
      fs for a 20 ps run. That is too coarse for the fastest vibrations, which then
      alias into the Green-Kubo integral: for liquid water the thermal conductivity
      came out 3-4 times too large. The heat flux is now sampled every two timesteps,
      and the meaningless "50" for the heat flux sampling is no longer shown in the
      dialog. For water the thermal conductivity is now about 1 W/m/K.
    * After the Heat Flux sub-step the final coordinates and velocities were not read
      back into the configuration, because the step looked for a file of the wrong
      name. They are now.

2026.10.14.1 -- Bugfix: MLFF dynamics on the CPU, NVE run length and short runs
    * A machine-learned forcefield from the Model Chemistry step now runs on the CPU
      when the job can have no GPU, using the job's cores less one for LAMMPS. Before,
      the step either waited forever for a GPU or stopped with "PyTorch forcefield
      requested, but no GPUs!". The Forcefield step's PyTorch models still need a GPU.
    * With busy GPUs the step waits for one for ``--gpu-wait-time`` seconds; a finite
      wait never ended before.
    * NVE dynamics with a machine-learned forcefield, or any model in LAMMPS's metal
      units, ran 1000 times the steps asked for (e.g. 40,000 steps for 0.02 ps at
      0.5 fs). NVT and NPT were right.
    * Dynamics of fewer than 100 steps (NVE, NVT, NPT, Heat Flux) stopped with "Illegal
      fix ave/time nfreq value: 0".
    * Files the MLFF engine writes, such as a committee's spreads, are now kept rather
      than removed when the run finishes.

2026.10.14 -- Bugfix: machine-learned forcefields with the engine in its own environment
    * A machine-learned forcefield from the Model Chemistry step failed on a default
      installation ("xnn: not found"): its engine and LAMMPS were always run as an MPI
      pair in LAMMPS's environment, which works only when one environment holds both,
      as on a cluster set up that way. When ``xnn.ini`` names another environment, the
      engine now starts there with ``conda run`` and the two talk over TCP on the same
      machine. ``mlff-transport`` in ``lammps.ini`` (``auto``, ``mpi`` or ``tcp``)
      overrides the choice and ``mlff-tcp-code`` changes the LAMMPS command. On Linux
      TCP is currently slower than MPI, by about a third for a fast model, until the
      MDI library sets TCP_NODELAY (proposed upstream).
    * mdi_bind.sh 2.3 runs a pair started without mpirun: ``SEAMM_MDI_ROLE`` names each
      side and ``SEAMM_MDI_PAIR`` the pair. Stopping a partner also stops its child
      processes, such as an engine started by ``conda run``.
    * The helper scripts in ``<root>/bin`` (mdi_bind.sh, cpu_bind.sh, gpu_bind.sh,
      mdi_monitor.sh) are now upgraded when the plug-in ships a newer version, keeping
      the old copy as, e.g., ``mdi_bind.sh.2.2``; before, the step only logged a
      warning that was easy to miss, and old copies failed (mdi_bind.sh 1.3 does not
      run with macOS's bash). A job running the old script is not affected. A script
      without a ``#MolSSI`` header line is left alone.

2026.10.12.1 -- Machine-learned forcefields from the Model Chemistry step
    * LAMMPS can take a machine-learned forcefield (``xnn:MLFF@...``) from a Model
      Chemistry step, with its checkpoint and engine options, and runs its engine and
      LAMMPS as an MPI pair on a GPU, each through ``mdi_bind.sh`` -- as ``gpu-code`` in
      ``lammps.ini`` did for the Forcefield step's PyTorch models, which still works.
      ``mlff-code`` in ``lammps.ini`` can change the launch. Requires xnn-step
      2026.10.12 and model-chemistry-step 2026.10.12.1.

2026.10.12 -- Forcefields from the Model Chemistry step
    * LAMMPS can take its forcefield from a Model Chemistry step: ``FF@oplsaa``,
      ``ReaxFF@CHO_cho_2008`` and the others the Forcefield plug-in offers. It builds the
      forcefield and types the atoms as it does with the Forcefield step, which still
      works as before. A model chemistry with a quantum method still runs the dynamics
      with the program over MDI.
    * Requires forcefield-step 2026.10.12.1 and model-chemistry-step 2026.10.12.

2026.10.11.2 -- Bugfix: mdi_bind.sh 2.2 really stops a failed MDI pair
    * 2.1 did not end a failed MLFF run. LAMMPS does not exit after an error in an MDI
      pair -- it waits in MPI_Finalize for the engine -- so neither side ever exited;
      and stopping a partner stopped only the shell around it, leaving the program
      running. 2.2 starts each program with exec, so that stopping it stops the
      program, and watches a copy of LAMMPS's output for ERROR lines: if LAMMPS has
      not exited SEAMM_ERROR_GRACE seconds (default 30) after one, it is stopped, and
      with it the engine. Tested on TinkerCliffs with LAMMPS and an xnn model: a
      failed run now ends within a minute instead of holding its GPU until cancelled.
    * The driver's output passes through unchanged, and line-buffered as it is on
      mpirun's terminal, so a hung LAMMPS no longer loses its output.

2026.10.11.1 -- Bugfix: a failed LAMMPS run no longer leaves its MDI engine holding the GPU
    * In an MDI pair (e.g. a machine-learned potential driving LAMMPS) a LAMMPS error
      left the engine waiting for its driver forever: the job held its GPU and stayed
      RUNNING until cancelled. mdi_bind.sh (version 2.1) now stops the engine when the
      driver fails, and the driver when the engine fails; after a normal finish the
      engine has SEAMM_ENGINE_GRACE seconds (default 120) to exit.
    * Installing the step does not replace an existing ~/SEAMM/bin/mdi_bind.sh; the
      step warns that a newer version is available. Replace it by writing the new file
      beside it and renaming it over the old one, so that running jobs keep their copy.

2026.10.11 -- The step declares what its cost model is made of
    * ``lammps.TIMING_SPEC`` -- atoms as the size variable, the forcefield as the class,
      and MD steps as the unit -- is passed with each timing record (seamm-exec
      2026.10.6.1 and later write it beside the records), so the cost model is fitted
      from the step's own description.
    * Each timing record names the forcefield (its form -- PyTorch, OpenKIM, MDI/QM or
      a classical form -- and the model or checkpoint), so that machine-learned and
      classical forcefields, orders of magnitude apart in cost, are fitted
      separately; a run on a GPU also records the GPU's model.
    * Removed the Docker image (and its recipe and release job) and the Docker settings
      in the default lammps.ini, which were no longer used or maintained.
    * Requires Python 3.12 or later.

2026.10.6.1 -- Bugfix: forcefield runs failed after LAMMPS finished
    * 2026.10.6's timing record used a variable bound only for QM-MD over MDI, so
      every ordinary forcefield run raised ``UnboundLocalError`` once LAMMPS had
      finished and its results were never analyzed (ChemAI jobs 5177-5188). The
      record now takes only what every run has.

2026.10.6 -- Timing records that a cost model can be fitted to
    * Each LAMMPS run appends a record to ``~/.seamm.d/timing/lammps.csv`` through
      ``seamm_exec.timing``: the machine class, processes, wall time and outcome, and
      the descriptors of the calculation: the sub-steps run, the structure, and from the log the ``Loop time`` lines: MD steps, atoms, processes and LAMMPS's own time (steps times atoms is the unit of cost). This replaces the step's own CSV
      (SMILES, formula and the whole parameter dictionary as JSON), which grew
      without bound. See seamm_exec's campaign of 2026-10-05.
    * Requires seamm-exec 2026.10.6.
2026.10.5 -- Job-level files for parallel loops
    * Files named ``/name`` are written to the job's own directory (``Node.job_path``)
      and read from it or, in an iteration of a parallel loop (loop_step 2026.10.5),
      from the enclosing iterations and the job. Nothing changes for other jobs.
    * Requires seamm 2026.10.5.
2026.10.3 -- Bugfix: QM/MD over MDI used the program's method name and default basis
    * The MDI engine for a QM-driven LAMMPS run was launched with the model
      chemistry's method name alone, so ORCA ran def2-SVP whatever basis was
      chosen, and a functional whose name the Model Chemistry step had to alter
      was not recognized. It now gets the program's own keyword and the chosen
      basis (with model_chemistry_step 2026.10.3).

2026.10.1 -- Settings that depend on each other; LAMMPS errors reported
    * The dialogs show only the settings that apply with the current choices, and SEAMM's
      flowchart tools use the same rules. Settings that were hidden but used now apply
      where they matter: the charged-atom threshold whenever there is a k-space method;
      the QEq settings for the forcefield's default charges; the damping times for the
      Berendsen barostat too. Allow shear applies only to a solid with the Nose-Hoover
      barostat.
    * A Berendsen NPT step whose box must be triclinic failed inside LAMMPS; it is now
      refused before LAMMPS starts, saying why.
    * Bugfix: when LAMMPS itself stopped with an error, the step failed later with an
      unrelated IndexError. It now stops with LAMMPS's message, the sub-step, the last
      command, advice for common errors, and where the log is.
    * Bugfix: an explicitly chosen k-space method raised KeyError; Minimization read
      Sxy from Sxz; NPT's stress damping times were converted as pressures, and with y
      and z coupled z's final stress was wrong; NPT's mtk and nreset were not used;
      Velocities' default momentum choice was not one of its own choices.
    * Documented in the user guide. Needs seamm 2026.10.1.
    * Internal: CI now installs the package's declared dependencies with uv rather than
      a conda test environment; model-chemistry-step, which it uses, is now declared.

2026.9.29 -- Bugfix: the CPU/GPU binding scripts work on any machine
    * ``mdi_bind.sh``, ``gpu_bind.sh`` and ``cpu_bind.sh`` bound the codes to fixed CPU
      cores laid out for one machine (an AMD EPYC 7763 with two GPUs). Elsewhere they
      failed: under SLURM or another scheduler, when the job's cores were not those
      ("taskset: failed to set pid's affinity: Invalid argument"); on a third or
      later GPU ("no CPU binding defined for GPU 2"); and on machines with fewer
      cores. On other layouts they silently pinned to unrelated cores.
    * Version 2.0 of the scripts works the binding out from the machine. Under a
      scheduler (SLURM, PBS, LSF) it leaves the binding to the scheduler. Otherwise
      it uses the CPUs next to the GPU, from ``nvidia-smi topo -m``, limited to the
      CPUs the job may use; GPUs sharing those CPUs split them, and the engine and
      LAMMPS take halves. If none of that can be determined the code runs unbound
      rather than failing. ``SEAMM_BIND=none`` turns binding off, and
      ``SEAMM_ENGINE_CPUS``/``SEAMM_DRIVER_CPUS`` (``SEAMM_CPUS`` for
      ``cpu_bind.sh``) choose the CPUs by hand.
    * The scripts also find the MPI rank under MPICH and srun, not only OpenMPI, and
      ``gpu_bind.sh`` no longer overrides a GPU allocation the scheduler made.
    * The warning that a newer version of a helper script is available never
      appeared, because its version header was misread; it does now. The copies in
      ``~/SEAMM/bin`` are not replaced automatically, since they may have been
      edited: delete an older ``mdi_bind.sh``, ``gpu_bind.sh`` or ``cpu_bind.sh`` to
      get version 2.0, unless you tuned it for your machine.
    * The documentation and the ``lammps.ini`` examples now give the plain command
      for an MDI engine (``mpirun --bind-to none -np 1 <engine> : -np 1 lmp``), which
      works on workstations and clusters alike, with ``mdi_bind.sh`` as an option.

2026.9.25 -- Bugfix: the LAMMPS step failed to load on Python 3.12 outside conda
    * The step used the GPUtil package to find idle GPUs. GPUtil is unmaintained and
      imports ``distutils``, which was removed from Python 3.12, so in any
      environment without setuptools (a plain venv, for example) the LAMMPS step
      silently failed to load and was missing from the flowchart editor. It worked
      in conda environments only because conda installs setuptools everywhere. GPU
      detection now queries ``nvidia-smi`` directly and GPUtil is no longer a
      dependency.

2026.9.20 -- Bugfix: a large system could be given a net charge
    * The charges were written to the structure file for LAMMPS with three
      decimals, while a forcefield assigns four or more. For a small molecule the
      rounding is harmless, but over a large system it accumulates: 500 molecules
      of ethylene carbonate, 5000 atoms whose charges add up to zero, were handed
      to LAMMPS as a system with a charge of -1. They are now written with six
      decimals, which leaves the total as the forcefield set it.
    * Warnings from assigning the forcefield -- in particular that the charges did
      not add up and have been adjusted -- now appear in the step's output and the
      job's output. They went only to the log, which for a job in a queue is the
      queueing system's output file if it is kept at all, so a change to the
      charges could pass unnoticed.
    * Internal: pinned bibtexparser to 1.x and added Pmw in the test environment,
      both of which the tests need and neither of which was there.

2026.9.16 -- Bugfix: trajectory output failed with machine-learned potentials
    * Writing a trajectory (dump or extxyz) from a run using a PyTorch/MDI potential
      such as MACE died in LAMMPS with "Unknown dump_modify keyword: H". The element
      list given to LAMMPS had one entry per atom rather than one per atom type, so
      for anything beyond a handful of atoms the surplus entries were read as
      keywords. The list is now per type, as it already was for OpenKIM models.
    * The Python engine for ML potentials has been renamed from xnns to xnn. The
      examples in the default lammps.ini and comments now use the new name. If your
      ~/SEAMM/lammps.ini has a gpu-code line using `xnns mdi`, change it to `xnn mdi`
      and reinstall the engine, since checkpoints saved with the new package can only
      be loaded by it.

2026.9.15.2 -- Bugfix: concurrent GPU jobs bound to the same cores
    * mdi_bind.sh chooses which cores to run the engine and the driver on from a map
      of the GPUs on the machine, so that each is near its card. It was looking that
      map up by the GPU's index within the job's own allocation, which is 0 for every
      job given a single GPU whichever card that is, so two jobs running at once both
      took the cores belonging to GPU 0 and shared them. It now looks up by the
      physical index, and jobs on different cards get different cores.
    * Note that ~/SEAMM/bin/mdi_bind.sh is not replaced once it exists, since the core
      layout in it is specific to the machine and meant to be edited. Delete it to
      pick up this version.

2026.9.15.1 -- Machine-learned forcefields are served by an engine, whichever one
    * Which engine evaluates a machine-learned forcefield is now set entirely by
      gpu-code in lammps.ini. The model file is offered to it two ways -- as the
      SEAMM_FF environment variable, and as {MODEL} for substitution into the command
      -- so an engine taking the model as an argument, such as 'xnns mdi --ckpt
      {MODEL}', works alongside one reading the environment, such as mace-mdi, with
      nothing in the plug-in needing to know which is configured.
    * Every PyTorch forcefield except an ML-IAP one now goes to that engine. The
      choice used to depend on the model file being named '.mace.pt', with anything
      else falling back to 'pair_style mace'. A file name says nothing about how a
      model is evaluated, and reaching that pair style needs a LAMMPS built with it --
      which is the build the MDI route exists to avoid needing. A checkpoint named for
      what it was rather than for that convention failed with "Unrecognized pair style
      'mace'". The pair style is no longer written.
    * The comm_modify line an MDI run needs was applied on the same file-name test,
      so it now applies to every MDI run.

2026.9.15 -- Bugfix: a classical forcefield was sent to the machine-learning engine
    * Asking for a GPU made any LAMMPS job take the GPU command from lammps.ini,
      whatever the forcefield. Where that command launches an MDI engine -- the usual
      setup for machine-learned potentials -- a classical run such as OPLS-AA was
      handed to that engine, which has no model to serve and died with
      "FileNotFoundError: 'Unknown'" before LAMMPS ran at all. The GPU count comes
      from the scheduler, so it says nothing about what the calculation needs; the
      choice now depends on whether the configured command suits the forcefield.
      A Kokkos GPU command, which is how a classical forcefield is meant to use a
      GPU, is unaffected.
    * This only became reachable when queues began offering a GPU count at
      submission, so the first classical job submitted with one hit it.

2026.9.14.3 -- Say plainly when a GPU job has too few tasks
    * A GPU command that drives an MDI engine runs the engine and LAMMPS as separate
      MPI ranks, and those counts are written into the command rather than derived from
      the allocation. A job given fewer tasks than the command asks for died inside
      mpirun with "All nodes which are allocated for this job are already filled",
      which mentions neither MDI, nor tasks, nor what to change -- and arrived wrapped
      in a conda error, obscuring it further. The step now checks before launching and
      says how many tasks the command needs, how many the job was given, and that an
      MDI run needs one task for the engine and one for the driver.
    * This is easy to hit now that queues can offer a GPU count at submission:
      requesting one GPU and one task looks reasonable but cannot work.

2026.9.14.2 -- Bugfix: could run on a GPU the scheduler had not allocated
    * When choosing GPUs the step asked GPUtil, which reads nvidia-smi and therefore
      reports physical GPU indices, knowing nothing about CUDA_VISIBLE_DEVICES. A
      scheduler that allocates a subset of a node's GPUs sets that variable, which also
      renumbers the devices: given physical GPU 1, a job's own device 0 is that GPU. The
      step would report GPU 1, so the calculation could run on a card it had not been
      given, alongside whichever job actually held it. GPUs are now selected in the
      numbering CUDA itself uses, and a GPU outside the allocation is never chosen.
    * mdi_bind.sh no longer overwrites CUDA_VISIBLE_DEVICES when a scheduler has
      already set it. That variable is not composable -- setting it again is
      interpreted against the machine's full set of devices, not against the
      allocation -- so overwriting it was another way to land on the wrong GPU. The GPU
      monitor still logs the physical card, mapped back for the purpose.
    * Both only affected machines where a scheduler hands out individual GPUs. Where
      every job sees every GPU the behaviour is unchanged.

2026.9.14.1 -- Internal: the MACE MDI engine is no longer duplicated here
    * This plug-in shipped its own copy of the MACE MDI engine as mace_mdi.py and
      installed it into ~/SEAMM/bin. The same engine is maintained in the lammps-mdi
      package, and the two had drifted well apart: the bundled copy was missing, among
      other things, the fix for loading a GPU-saved model on a machine without a usable
      GPU. It has been removed, and the example configurations now use 'mace-mdi', the
      console script from lammps-mdi.
    * An existing ~/SEAMM/bin/mace_mdi.py is left in place, so a configuration pointing
      at it keeps working. The step now warns when gpu-code still refers to it, since
      that copy no longer receives fixes.
    * Note that lammps-mdi is deliberately NOT a dependency of this plug-in: it pulls in
      PyTorch and MACE, which belong in the LAMMPS environment rather than the SEAMM
      one. The Getting Started guide shows how to install it there.

2026.9.14 -- Bugfix: machine-learned forcefields could not run at all
    * The conda environment asked for conda-forge's LAMMPS, which is built without
      the MDI package. Since machine-learned forcefields are run by driving a separate
      engine over MDI, they could not work: LAMMPS rejected the '-mdi' option and then
      hung until the job hit its time limit rather than reporting the error. The
      environment now uses an MDI-enabled LAMMPS, and OpenMPI, so that LAMMPS and the
      engine share one MPI and one MDI library.
    * PyTorch and the MACE stack are deliberately not part of that environment, as the
      correct PyTorch build depends on the machine's NVIDIA driver. Install them with
      'lammps-mdi install-ml', which works this out for you. The Getting Started guide
      now covers setting up a machine-learned forcefield end to end.
    * Fixed mdi_bind.sh reporting failure after a run that had in fact succeeded, which
      made a completed MDI calculation look like a failed one. Note that an existing
      ~/SEAMM/bin/mdi_bind.sh is not replaced automatically, since it is meant to be
      edited; delete it to pick up the fixed version.

2026.3.26 -- Bugfix: MDI stress not calculated and error catching
    * MDI requires an option in fix mdi/qm for it to request the virial; it is not
      automatic for e.g. NPT calculations. The code now adds the keywords for periodic
      calculations but not molecular ones.
    * Added code to catch errors in the analysis of the GPU performance in cases where
      data is missing or cannot be interpreted correctly, so that calculation continues
      without error.
    * Fixed an issue returning the stdout and stderr from the LAMMPS/MDI
      calculation. LAMMPS logs to stdout and MDI to stderr to keep the streams separate.

2026.3.24 -- Bugfix: Issues calculating the stress and its sign.
    * Added an option to the fix mdi/qm line to ensure that MDI calculates the stress.
    * Fixed the sign issue in the stress between LAMMPS and MDI.

2026.3.23 -- Improved the performance output and results for MD

2026.3.22.1 -- Fixed issue running LAMMPS using a conda environment
    * Fixed a problem running LAMMPS directly in a conda environment, due to the
      construction of the command line having an error.

2026.3.22 -- Added support for MDI and GPUs, mainly for MLFFs
    * Added monitoring of GPU utilization and memory utilization
    * Added more output about the performance of the run to the job output.
    * Added the ability to run MACE MLFFs using MDI. This avoids needing a GPU version
      of LAMMPS that supports ML-IAP, making it much easier to run these MLFFs.
    
2026.2.25 -- Added support for ML-IAP forcefields.
    * Handle forcefields with the extension .pt and containing 'mliap' in their name using
      ML-IAP.

2026.1.28 -- Bugfix: Fixed error in analysis of non-periodic molecules
    * Recent additions to the analysis introduced a bug for non-periodic systems, where
      in some circumstances the code crashed due to interpreting the cell info in LAMMPS
      dump files.

2026.1.8 -- Improved handling of atom position and velocity trajectories.
    * Added the time and units for the time to the atomic poisition and velocity
      trajectories used in e.g. diffusivity to ensure the correct times are used.

2025.12.9 -- Enhancements: Running on GPUs and output from minimization
    * Ensured that at least 1 GPU is used when running using PyTorch potentials. This
      allows running from the commmand-line or locally in addition to using SLURM, which
      already worked.
    * Added the cell parameters, volume, and density, plus the change in each, whne
      running minimizations. These are now results that can be stores in tables and
      saved in variables.

2025.11.22 -- Bugfixes: volume units and converting pressure units.
    * Fixed the incorrect units for volume in the metadata for LAMMPS.
    * Fixed a crash caused by LAMMPS storing large pressures as integers, not floating
      point numbers.
      
2025.11.19.1 -- Bugfix: Properties not handled in MD runs
    * Due to an inadvertant error in the previous release, properties were not placed in
      local variables, tables or as properties on the configuration if requested.
      
2025.11.19 -- Enhancements to properties and storing trajectories in the database

    This release enhances the LAMMPS step with property calculations and trajectory storage
    capabilities. The main focus is adding formation energy calculations (DfE0) using
    tabulated element energies, storing stress as a 6-component Voigt vector, and
    enabling trajectory import into the database as system configurations.

    * Added formation energy (DfE0) and atomization energy calculations using tabulated
      reference data for atomic energies
    * Implemented trajectory save functionality to store MD trajectories as
      configurations in the database 
    * Reorganized stress tensor output to follow Voigt ordering (Sxx, Syy, Szz, Syz,
      Sxz, Sxy) and added stress vector property 
      
2025.11.9 -- Bugfix: Error with units in trajectories when using OpenKIM

     This release adds support for volume units across all LAMMPS unit systems to fix
     a bug in trajectory analysis when using OpenKIM potentials that use "metal" units.
     The changes enable proper unit conversion for volume quantities during trajectory
     analysis. Key changes are

     * Added "volume" definitions to all seven LAMMPS unit systems
     * Added volume dimensionality mappings in the unit conversion dictionaries
     * Implemented unit conversion logic in trajectory analysis to convert from LAMMPS
       units to SEAMM units
      
2025.11.6 -- Bugfix: Error in writing extxyz trajectory files.
    * Fixed the condition check from "Pxx" to "Sxx" to properly detect stress data
    * Corrected the case of the stress component variable from 'sxy' to 'Sxy' for consistency
     
2025.11.5 -- Bugfix: Fixed errors with trajectories using OpenKIM potentials
   * Fixed an error writing the dumpfile for trajectories for e.g. MEAM potentials where
     there are no charges, and the Open KIM atom type does not support charges on the
     atoms.
   * Fixed an error with the element labels for OpenKIM potentials.

2025.9.24 -- Bugfix: Error with cell angles for sheared cells
   * Fixed an error with sheared cells that actualy had 90 degree angles, which were
     incorrectly made to be 0 degress :-)
     
2025.9.21.1 -- Added cell angle parameters to properties
   * The cell lengths, but not angles, were available as properties. This update adds
     the angles.
     
2025.9.21 -- Bugfix: error in NPT dynamics
   * Fixed an error that slipped through and caused the code to crash.
     
2025.9.20 -- Completed handling of stress in NPT dynamics
   * Added the ability to define the offdiagonal stresses XY, XZ, and YZ, as well as to
     fix any of the directions. This is now handled in the same way as the Minimization
     step handles the stress.

2025.9.2 -- Bugfix: error in the stress parameters in the GUI
   * There was an error in the setup of the stress parameters, which caused a
     crash. This is now fixed.
     
2025.8.27 -- Added cell optimization to minimization step
   * Added the ability to optimize the cell during minimization, giving control over
     applied pressure or stress as well as fixing lattice directions.
   * Allowing "=" as the first character of a parameter to signal that it is a variable
     or expression that needs to be evaluated before use.

2025.8.21 -- Added stress and improved graph handling
   * Added stress for periodic systems to the state trajectory, and added to the
     properties and results so that it can be saved and used elsewhere.
   * Improved the formatting of graphs and added options for writing the graphs to a
     variety of formats such as HTML, PDF, and JPG.

2025.8.2 -- Initial pass at trajectories.
   * Added option to save the geometry & energies for a trajectory, and converting to
     the extended xy format that e.g. ASE uses, and the MACE ML code uses as input.
   * Bugfix: an error slipped into the tabulated angle term used for PF6- when
     formatting the forcefield parameters and tracking the provenance. This is fixed.
     
2025.7.23 -- Added a complete description of number of cores used.
   * Added a detailed description of how the codes decided on the number of cores to
     use, so that users can understand and change the limits in seamm.ini to their
     desired values.
   * Changed the default number of atoms per core (natoms-per-core) to 100.

2025.6.25 -- Bugfix: ReaxFF memory allocation
   * Apparently the memory handling for reafxx in LAMMPS is kludgy, and there is an
     empirical fudge-factor that helps solve errors with e.g. numbers of bonds. This
     version increases the factor to 1.6 until we get more experience.

2025.6.23 -- Bugfix: handling of GPUs and KOKKOS arguments
   * Fixed issue detecting GPUs caused by differences in SLURM enviroment variables
   * Fixed inconsistencies of dash vs underscore in configuration arguments cmd-args and
     gpu-cmd-args, which are used for e.g. Kokkos.

2025.6.21 -- Added support for MACE PyTorch potentials.

2025.5.27 -- Bugfix: error in formatting forcefield parameters
   * The cleanup of the format of 'structure.dat' accidentally converted some floating
     point numbers into integers, with loss of precision. This is now fixed.

2025.5.26 -- Added more metadata to the structure.dat file.
   * Added the original parameters, when available, plus the version of each parameter
     to the 'structure.dat' file input to LAMMPS to aid in understanding which
     forcefield parameters were used and where they came from.
   * Cleaned up the format of the 'structure.dat' file to make it more readable.

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
