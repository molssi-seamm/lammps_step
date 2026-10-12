Phase A: classical forcefields and ReaxFF (2026-10-12)
======================================================

What was built
--------------

* **forcefield_step** (``model_chemistry.py``; classmethods on ``ForcefieldStep``):

  - ``get_model_chemistry_options()`` offers one model chemistry per ``#define`` in
    every ``.frc`` file found, with the file's own names: ``FF@<name>``
    (``FF@oplsaa+``, ``FF@pcff``; ``dreiding/X6`` is ``FF@dreiding-X6``, the
    grammar's reserved characters replaced by ``-``) and ``ReaxFF@<set>`` for
    ``reaxff/<set>``. (Paul, 2026-10-12: the names are the files', not a table of
    display names such as ``FF@OPLS-AA``.) Files another ``.frc`` file includes
    (``#include``) are parts of that forcefield and are not offered: e.g.
    ``local:OPLS-AA/ligpargen.frc`` extends ``oplsaa+``. There is no owner, since LAMMPS evaluates them itself. They
    are not MDI engines (``mdi_only`` gives nothing).
  - The files are searched in the same order as data files: the user's
    (``~/.seamm.d/data/Forcefields``), then the installation's
    (``<root>/data/Forcefields``, ``~/SEAMM/data/Forcefields``), then the package's.
    The first forcefield of a name wins, so a user's file shadows a shipped one.
  - ``forcefield(mc, node)`` builds the ``seamm_ff_util`` forcefield as the
    Forcefield step does, with the node's references for the citations.
  - ``model_chemistry_data_files(mc, node)`` returns the ``local:``/``personal:``
    file and its includes, for a submitted job.

* **model_chemistry_step**:

  - ``discover_model_chemistries(..., data_path=)`` passes the flowchart's data
    path (``Node.data_path``: in a JobServer the job's own ``data/`` first) to
    providers whose ``get_model_chemistry_options`` accepts it. That is how a job
    finds the user's forcefield sent with it, which is not in the run host's own
    directories. Providers without the argument are called as before.
  - ``ModelChemistry.list_data_files()`` asks the offering step for its files.
  - The dialog has no Program for a level offered without an owner, and composes
    ``FF@oplsaa``. ``_programs()`` drops ``None``, because ``sorted()`` of a mix
    raises.

* **lammps_step**: ``forcefield_model_chemistry()`` and ``forcefield()``. The
  forcefield is built through the provider (``plugin_manager.get(mc["step"])``)
  and cached per level and file. ``ff_form()`` returns the forcefield's form for
  FF/ReaxFF and ``MDI/QM`` for other model chemistries only. ``initialization``
  and ``heat_flux`` call ``parent.forcefield()``.

* **packmol_step** types with the forcefield model chemistry when there is one,
  else with ``_forcefield``. **diffusivity_step** names the model chemistry's level.

Verification
------------

Packmol (12 ethanol, periodic) and then LAMMPS Initialization + Energy were run
both ways, through a Forcefield step and through a Model Chemistry step, for
``FF@oplsaa+`` (the Forcefield step's default for ``oplsaa.frc``) and
``ReaxFF@CHO_cho_2008``:

* ``input.dat`` (and ReaxFF's ``forcefield.dat``) were identical;
* ``structure.dat`` was identical apart from the coordinates (Packmol packs
  randomly): the same atom types, charges, bonded terms and coefficients;
* the citations were the same, including the forcefield's (Jorgensen).

Notes for later phases
----------------------

* The Forcefield step's default for ``oplsaa.frc`` is its first forcefield,
  ``oplsaa+``: the converter (phase C) must map it to ``FF@oplsaa+``.
* **OpenKIM (IP@) is not done.** The plan offers *installed* potentials, but
  none are installed in practice: this Mac has none, and ChemAI's LAMMPS envs
  have only kim-api's examples (``LennardJones_Ar``, ``ex_model_*``). Listing
  them means running ``kim-api-collections-management list`` in LAMMPS's
  environment (about 0.3 s).
* Releases: forcefield_step and model_chemistry_step first; then lammps_step and
  packmol_step require them (the released dialog would compose ``None:FF@...``,
  and the old forcefield_step has no ``forcefield()``).
