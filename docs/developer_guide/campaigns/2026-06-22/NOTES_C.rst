.. _notes-mdi-lammps-phase-c:

==========================================================================
Phase C -- LAMMPS driver: read ``_model_chemistry``, launch the MDI engine
==========================================================================

:Author: Paul Saxe (with Claude)
:Date: 2026-06-24
:Status: Implemented (2026-06-24) -- ``lammps.py`` + ``initialization.py`` +
         unit tests done (14 pass). Manual end-to-end run still pending.
:Campaign: LAMMPS + MOPAC/xTB QM-MD via MDI

.. contents:: Contents
   :depth: 2
   :local:


Phase numbering
===============

This campaign grew a step. ``NOTES_A`` (written 2026-06-23, in
``molssi-seamm.github.io``) calls the LAMMPS driver work "Phase B", because at
that point the only two halves were *MOPAC engine* and *LAMMPS driver*. We then
inserted the **Model Chemistry step** as a first-class producer (``NOTES_B``, in
``model_chemistry_step``), which took the name "Phase B". The LAMMPS driver is
therefore **Phase C**. Where ``NOTES_A`` says "Phase B (lammps_step)", read
"Phase C". The contract did not change, only the label.


Scope
=====

Phase C delivers the **LAMMPS side** of QM-MD over MDI: ``lammps_step`` reads the
``_model_chemistry`` workspace variable published by the Model Chemistry step
(Phase B), validates it can be driven via MDI, asks the owning program step for
its engine launch command (Phase A's ``get_mdi_engine_command``), allocates a
TCP rendezvous, composes a launch script that runs engine + LAMMPS together, and
runs it -- closing the thin line "select MOPAC PM6-ORG, then drive LAMMPS MD with
it via MDI".

In scope:

* Detect the QM/MDI path: ``_model_chemistry`` present (mirrors how the existing
  code branches on ``_forcefield`` being ``"OpenKIM"`` / ``"PyTorch"``).
* Validate ``options["mdi_capable"]`` (and ``options["periodic_mdi"]`` for a
  periodic configuration) with a helpful error if the selection cannot be driven.
* Resolve the owning step via the Stevedore handle
  (``self.flowchart.plugin_manager.manager[mc["step"]].obj``) and call
  ``get_mdi_engine_command(...)`` on it.
* Allocate a free TCP port (``NOTES_A`` ``_free_tcp_port`` helper) and choose the
  hostname; pass the *same* values to engine and driver.
* Build the LAMMPS (driver) argv for TCP transport from ``lammps.ini`` -- the
  same ini-reading machinery already in ``LAMMPS.run()``.
* Compose the bash launch script: engine backgrounded first (it is the
  listener), driver second, ``wait`` last (``NOTES_A`` "Launch script shape").
* Generate the LAMMPS input deck: ``fix ... all mdi/qm ... elements <types>``,
  reusing the existing MACE/PyTorch path in ``initialization.py``.
* Charge / multiplicity from the ``configuration`` object
  (``configuration.charge`` / ``configuration.spin_multiplicity``), per SEAMM
  convention -- never from the model-chemistry step.

Deferred:

* xTB (revisit once MOPAC works end to end -- it just needs its own Phase A).
* Non-conda engine launches (local / modules / docker); Phase A's
  ``get_mdi_engine_command`` already raises ``NotImplementedError`` for these.
* Concurrent-job port-collision retry (``NOTES_A`` "Port allocation").
* GPU binding for the QM engine. MOPAC is CPU-only, so the MACE ``mdi_bind.sh`` /
  ``mdi_monitor.sh`` GPU-pinning wrappers are *not* on the MOPAC path.


What we consume (the Phase B contract)
======================================

``self.get_variable("_model_chemistry")`` returns::

    {
        "model_chemistry": "MOPAC:SQM@PM6-ORG",   # canonical string
        "program": "MOPAC", "type": "SQM", "method": "PM6-ORG",
        "basis": None, "cutoff": None,
        "step": "<stevedore plugin name>",         # resolution handle
        "options": {                               # full get_model_chemistry_options() entry
            "mdi_capable": True,
            "periodic_mdi": False,
            "mdi_method_arg": "PM6-ORG",
            "elements": "...", "description": "...", ...
        },
    }

Phase C reads ``mc["step"]`` to resolve the owner, ``mc["method"]`` (or
``options["mdi_method_arg"]``) for the engine, and the two ``options`` booleans
to gate MDI launch. It needs nothing else from the program step except the argv
that ``get_mdi_engine_command`` returns.


The existing MACE path, and why MOPAC is a *new* path
=====================================================

``lammps_step`` already drives an MDI engine -- MACE -- but the launch model is
materially different from MOPAC's, so Phase C adds a parallel path rather than
extending the MACE one.

.. list-table::
   :header-rows: 1
   :widths: 18 41 41

   * - Aspect
     - MACE (existing)
     - MOPAC (Phase C)
   * - Transport
     - MDI ``-method MPI``
     - MDI ``-method TCP -port -hostname``
   * - Launch
     - Single ``mpirun ... : ...`` multiple-program line, *hardwired in*
       ``lammps.ini``
     - bash script: engine ``&`` then driver, then ``wait``
   * - Engine env
     - same conda env as LAMMPS, plain ``python``
     - **separate** conda env (``seamm-mopac``) via ``conda run``
   * - Engine argv
     - fixed string in ``lammps.ini``
     - returned by the program step's ``get_mdi_engine_command()``
   * - Engine binary
     - ``data/mace_mdi.py``
     - ``mopac_step/data/mopac_mdi.py``
   * - GPU pinning
     - ``mdi_bind.sh`` / ``mdi_monitor.sh``
     - none (MOPAC is CPU-only)

Shared, and reused as-is: the input-deck generation. Both emit
``fix mdi_fix all mdi/qm [virial yes] elements <atom types>`` and *no*
``pair_style`` (``initialization.py``, ``PyTorch_input``). The QM path differs
only in *who is launched and how*, not in what LAMMPS is told to do.

Why TCP, not MPI, for MOPAC: the engine lives in a different conda environment.
You cannot put two different conda environments into one ``mpirun ... : ...``
MPMD launch and have MDI's MPI transport rendezvous across them. TCP decouples
the two processes -- each is launched however its own step wants (``conda run``
for MOPAC, plain for LAMMPS) and they meet on a socket. This is exactly the split
of responsibility ``NOTES_A`` fixed: driver owns the rendezvous (port/host),
engine owns its own launch line.


Proposed design
===============

Detection -- ``ff_form()`` / ``run()`` (``lammps.py``)
------------------------------------------------------

Today ``ff_form()`` branches on ``_forcefield`` being ``"OpenKIM"`` /
``"PyTorch"`` / a forcefield object. Add a QM branch *before* touching
``_forcefield``, because a QM-MD flowchart has a Model Chemistry step and **no**
Forcefield step::

    if self.variable_exists("_model_chemistry"):
        return "MDI/QM"           # new form
    ff = self.get_variable("_forcefield")
    ...

``run()`` then, for ``ff_form == "MDI/QM"``, sets ``self.model = "MDI/QM/" +
mc["model_chemistry"]`` (so summaries/labels read sensibly) and skips the
forcefield-assignment machinery.

Engine launch -- new helper in ``lammps.py``
---------------------------------------------

A focused method, called from ``run()`` on the QM path, that returns the engine
argv and the rendezvous it picked::

    def _mdi_engine_launch(self, configuration):
        mc = self.get_variable("_model_chemistry")
        options = mc["options"]

        periodic = configuration.periodicity != 0
        if not options.get("mdi_capable"):
            raise ValueError(
                f"The model chemistry '{mc['model_chemistry']}' cannot be "
                "driven via MDI."
            )
        if periodic and not options.get("periodic_mdi"):
            raise ValueError(
                f"The model chemistry '{mc['model_chemistry']}' is not "
                "validated for periodic systems via MDI."
            )

        port = _free_tcp_port()                 # NOTES_A helper, host "localhost"
        step = self.flowchart.plugin_manager.manager[mc["step"]].obj
        engine_argv = step.get_mdi_engine_command(
            self.flowchart.executor,
            self.global_options,
            method=mc["method"],
            port=port,
            hostname="localhost",
            charge=configuration.charge,
            multiplicity=configuration.spin_multiplicity,
        )
        return engine_argv, port

Driver argv (LAMMPS, TCP DRIVER side)
-------------------------------------

Built from the existing ``lammps.ini`` ``code`` / ``cmd-args`` keys (the CPU
path; GPU/Kokkos is irrelevant for a CPU QM engine), plus the MDI driver flag
and the input file::

    driver_argv = [config["code"], *config["cmd-args"].split(),
                   "-mdi", f"-role DRIVER -name LAMMPS -method TCP "
                           f"-port {port} -hostname localhost",
                   "-in", "input.dat"]

Launch script (``NOTES_A`` shape)
---------------------------------

::

    #!/bin/bash
    set -e
    <shlex.join(engine_argv)> &
    ENGINE_PID=$!
    <shlex.join(driver_argv)>
    wait $ENGINE_PID

The engine is backgrounded first because it is the listener
(``MDI_Accept_Communicator``); the driver connects and retries. We then run this
script through ``executor.run(..., shell=True)``, the way the existing path runs
its command. ``conda run --live-stream`` (already in the engine argv from
Phase A) keeps the engine's stdout/stderr interleaved into the job output.

Input deck -- reuse ``initialization.py``
-----------------------------------------

Route ``ff_form == "MDI/QM"`` to the same code as ``"PyTorch"`` for the ``.pt``
branch: emit no ``pair_style`` and::

    fix    mdi_fix all mdi/qm [virial yes] elements <atom types>

(``virial yes`` only when periodic). The element list comes from the system, as
it does on the MACE path. This is the one place the two MDI paths converge.


Decisions (resolved 2026-06-24)
===============================

#. **D1 -- New ``ff_form`` value ``"MDI/QM"``, discriminated by**
   ``_model_chemistry`` **being present.** *Resolved: yes.* The branch is checked
   *before* ``_forcefield`` (a QM-MD flowchart has no Forcefield step, so
   ``_forcefield`` does not exist). No GUI toggle -- variable-presence mirrors the
   existing ``_forcefield`` / ``_OpenKIM_Potential`` precedent.

#. **D2 -- TCP + bash launch script.** *Resolved: TCP + script.* The engine runs
   in a separate conda env (``seamm-mopac``) via ``conda run``; two conda envs
   cannot share one ``mpirun ... : ...`` MPMD launch for MDI's MPI transport, so
   the MACE MPMD line is **not** reused. A second launch path is the accepted
   cost; it preserves Phase A's "engine owns its own env" design.

#. **D3 -- ``_mdi_engine_launch`` helper; rest inline.** *Resolved, with a small
   amendment during implementation.* The ``_mdi_engine_launch`` method picks the
   port and builds the engine argv. The driver-argv + script assembly were
   pulled out of ``run()`` into a pure module-level ``_mdi_launch_script(engine_argv,
   port, config, ce)`` (rather than left inline) **purely so the script shape and
   the driver-line template resolution are unit-testable** -- ``run()`` now just
   calls the two helpers. No new ``get_executor_config``-style method on
   ``lammps_step``; the inline ini-reading already there suffices.

#. **D4 -- Reuse ``code`` / ``cmd-args``.** *Resolved.* The driver argv is built
   from the existing CPU ``lammps.ini`` keys (GPU/Kokkos keys ignored -- MOPAC is
   CPU-only), with the MDI DRIVER flag and ``-in input.dat`` appended. No
   dedicated ``mdi-code`` key.

#. **D5 -- Simple ``_free_tcp_port``, retry deferred.** *Resolved.* Accept the
   TOCTOU window for the one-job-per-node case; defer the bind-failure retry
   (``NOTES_A`` open question 3) to a later pass.

#. **D6 -- Charge / multiplicity from the configuration.** *Resolved.* Read
   ``configuration.charge`` and ``configuration.spin_multiplicity`` and pass them
   into ``get_mdi_engine_command`` (what ``mopac.py`` itself does). Note: LAMMPS's
   ``fix mdi/qm`` does **not** send ``>TOTCHARGE`` / ``>ELEC_MULT``, so these are
   fixed for the whole run at launch -- acceptable for a single-species MD box.


Tests
=====

* Unit: ``_free_tcp_port`` returns a bindable integer.
* Unit: the QM detection branch in ``ff_form`` (mock ``_model_chemistry``
  present / absent).
* Unit: launch-script composition -- given a stub ``_model_chemistry`` and a
  stub program step whose ``get_mdi_engine_command`` returns a known argv, assert
  the script backgrounds the engine, runs the driver with the matching
  ``-port`` / ``-hostname``, and waits on the engine. No real engine, no real
  socket.
* Unit: MDI-capability gating raises the right ``ValueError`` for a
  non-``mdi_capable`` selection and for periodic-without-``periodic_mdi``.
* (Manual / integration, out of CI) the end-to-end thin line: MOPAC PM6-ORG
  driving a few MD steps on a small molecule.


References
==========

* Phase A engine contract: ``molssi-seamm.github.io`` ..
  ``campaigns/2026-06-22/NOTES_A.rst`` (``get_mdi_engine_command``,
  ``get_executor_config``, ``_free_tcp_port``, launch-script shape, port race).
* Phase B producer + contract: ``model_chemistry_step`` ..
  ``campaigns/2026-06-22/NOTES_B.rst`` and ``model_chemistry_step/grammar.py``.
* Existing MDI path (reused input deck): ``lammps_step/lammps_step/lammps.py``
  (``ff_form`` ~629, model selection ~649, ini/command building ~1001-1149) and
  ``initialization.py`` (``PyTorch_input``, ``fix mdi/qm`` ~807).
* LAMMPS ``fix mdi/qm``: https://docs.lammps.org/fix_mdi_qm.html
* MDI Library (TCP transport, roles): https://molssi-mdi.github.io/MDI_Library/
