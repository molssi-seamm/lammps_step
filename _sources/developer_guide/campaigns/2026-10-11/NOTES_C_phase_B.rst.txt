Phase B: machine-learned forcefields (2026-10-12)
=================================================

What was built
--------------

* **xnn_step**: ``get_mdi_engine_command(..., mdi_method="MPI", checkpoint=,
  device=)`` returns the engine's argv for an MPI pair: ``xnn mdi --ckpt <file>
  --device <device> <engine options> -mdi "-role ENGINE -name XNN -method MPI"``,
  with no ``conda run`` and no thread-count prefix. The engine runs in the
  driver's (LAMMPS's) environment, since the pair must share MPI and the MDI
  library. ``model_chemistry_option()`` describes one checkpoint;
  ``model_chemistry_for_checkpoint(path)`` does so for a checkpoint given by path,
  named after its stem, with the same refusal of dispersion removed twice.

* **model_chemistry_step**: two parameters that apply only to an MLFF (a rule
  shared by the dialog and the builder):

  - ``checkpoint``: a path, relative to the job's directory unless absolute, or a
    ``$variable``. At run time it is described by the step that offers the
    owner's MLFFs, and becomes the published model chemistry
    (``xnn:MLFF@<stem>``). With a checkpoint the model chosen need not be one of
    those listed; the builder leaves the check to the run.
  - ``engine options``: split with ``shlex`` into ``options["engine_args"]``.

* **lammps_step**:

  - ``mlff_model_chemistry()`` and ``pytorch_model()``. ``ff_form()`` is
    ``PyTorch`` for an MLFF model chemistry, so the existing PyTorch machinery
    applies unchanged: GPU polling, one pair per GPU, the ``fix mdi/qm`` input
    and the timing record.
  - ``_mlff_launch()`` writes ``mdi_launch.sh``, which runs as
    ``bash mdi_launch.sh``, as the QM path does. The script holds the pair
    (``MLFF_PAIR``, or ``mlff-code`` from lammps.ini): ``mpirun --bind-to none
    -np 1 {bind} {engine} : -np 1 {bind} lmp -mdi "..." -in input.dat``, where
    ``{bind}`` is ``<root>/bin/mdi_bind.sh``.
  - ``{bind}`` and ``{engine}`` are filled in by Python, not by the executor's
    ``str.format`` loop, because engine options can contain braces
    (``--dispersion "{regime: large}"``).
  - The device is ``cuda:0``, the GPU that mdi_bind.sh makes visible from
    ``SEAMM_GPUS``. Without a GPU it is xnn.ini's ``device``.
  - A job allocated fewer than two MPI tasks is refused before launching.
  - The Forcefield step's ``.pt`` route and ``gpu-code`` are unchanged.

Validation
----------

mlff's 1,554-atom NaCl(aq) box (512 H2O + 9 NaCl), with the ``ew10_s2`` MACE-LES
model, its subtracted D4 added back and ``--eeq-reuse``, was run as 1,000 steps of
NVT at 298.15 K with a 1 fs step. Falcon job 634309 (T4, fal144) ran both routes
one after the other on the same GPU:

* old: a Forcefield step with the ``.pt`` and lammps.ini's ``gpu-code``;
* new: a Model Chemistry step with ``xnn:MLFF@ew10_s2``, the checkpoint by path,
  and ``--eeq-reuse`` as an engine option.

Results:

* Both exited 0, and ``input.dat`` was identical.
* pe(0) was -1722.884 eV in both, matching mlff's reference of -1722.88403 eV
  (xnn 0.11.1).
* Energies and temperatures were identical to the printed precision through
  0.4 ps, then diverged slowly, as float32 MD does: at most 0.02 eV at the
  printed 0.1 ps points by 1 ps.
* Loop time for 1,000 steps: 412.0 s (old) and 412.7 s (new), i.e. 2.427 and
  2.423 steps/s. The two launches are the same.

Left for later
--------------

* GUI narrowing (design 7, phase C): a Model Chemistry feeding LAMMPS lists only
  what LAMMPS can run.
* The converter and the retirement of the Forcefield step's dialog (phase C).
* OpenKIM (``IP@``), deferred by Paul.
