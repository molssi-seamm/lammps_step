.. _notes-mdi-lammps-phase-c-runbook:

==================================================================
Phase C runbook -- manual end-to-end verification of QM-MD via MDI
==================================================================

:Author: Paul Saxe (with Claude)
:Date: 2026-06-24
:Status: Procedure -- not yet executed
:Campaign: LAMMPS + MOPAC/xTB QM-MD via MDI

.. contents:: Contents
   :depth: 2
   :local:


Purpose
=======

The Phase C unit tests (``tests/test_mdi_qm.py``) cover detection, MDI-capability
gating, engine resolution, and launch-script composition with stubs -- but
nothing in CI actually starts a MOPAC engine or runs LAMMPS. This runbook is the
manual "thin line" check: **select MOPAC PM6-ORG, then drive a few LAMMPS MD
steps with it via MDI**, end to end, and confirm it works. Run it once after the
Phase C code lands, and again whenever the engine contract or launch path
changes.


Prerequisites
=============

Three conda environments, the standard SEAMM split (the launch script relies on
exactly this -- the driver runs in seamm-lammps, the engine in seamm-mopac):

* **seamm** -- the environment SEAMM itself runs in. Must have, importable:
  ``lammps_step`` (with the Phase C code), ``model_chemistry_step``,
  ``mopac_step`` (with the Phase A engine contract). Install editable if testing
  local changes::

      conda activate seamm
      pip install -e /Users/psaxe/Work/SEAMM/lammps_step
      pip install -e /Users/psaxe/Work/SEAMM/model_chemistry_step
      pip install -e /Users/psaxe/Work/SEAMM/mopac_step

* **seamm-lammps** -- the LAMMPS binary, **built with the MDI package** so
  ``fix mdi/qm`` exists. Verify::

      conda run -n seamm-lammps lmp -h | grep -i mdi      # expect mdi/qm listed

* **seamm-mopac** -- MOPAC plus the engine's imports (``mopactools``,
  ``pymdi``/``mdi``, ``numpy``, ``seamm_util``). The engine script
  ``mopac_mdi.py`` is shipped in ``mopac_step/data/`` and run from there; it must
  import *only* packages present in this environment (NOTES_A, Option C). Verify::

      conda run -n seamm-mopac python -c "import mopactools, mdi, numpy, seamm_util"

The ini files must select the conda installs:

* ``~/SEAMM/mopac.ini`` -- the active executor section has
  ``installation = conda`` and ``conda-environment = seamm-mopac``. (Non-conda
  raises ``NotImplementedError`` by design -- NOTES_A.)
* ``~/SEAMM/lammps.ini`` -- the active executor section has a ``code`` key (e.g.
  ``mpirun -np {NTASKS} lmp``) and ``installation = conda`` /
  ``conda-environment = seamm-lammps``.


Pre-flight smoke checks (cheap, isolate failures early)
=======================================================

Do these before a full job; each isolates one layer.

#. **Discovery** -- the program-step sweep offers MOPAC PM6-ORG. Check the
   producing classmethod directly (no flowchart needed)::

       conda run -n seamm python -c "
       from mopac_step.mopac_step import MOPACStep
       opts = MOPACStep.get_model_chemistry_options(mdi_only=True)
       print(sorted(opts))
       assert any(o.get('model_chemistry') == 'MOPAC:SQM@PM6-ORG'
                  and o['mdi_capable'] for o in opts.values())
       print('PM6-ORG is MDI-capable: OK')
       "

   Or, in the SEAMM GUI, open a Model Chemistry step and confirm
   ``MOPAC:SQM@PM6-ORG`` appears in the combobox.

#. **Engine stands alone** -- the MOPAC MDI engine accepts a TCP connection
   without LAMMPS. In one shell start the engine::

       conda run -n seamm-mopac python \
           /Users/psaxe/Work/SEAMM/mopac_step/mopac_step/data/mopac_mdi.py \
           -mdi "-role ENGINE -name MOPAC -method TCP -port 8021 -hostname localhost" \
           --method PM6-ORG --charge 0 --multiplicity 1

   It should bind the port and block in ``MDI_Accept_Communicator`` (log line
   "MDI connection established" appears only once a driver connects). Ctrl-C to
   stop. A bind error here means the port is taken -- pick another.

#. **MDI TCP semantics** -- confirm MDI's TCP transport needs the *same explicit*
   ``-port`` on both sides and has no port-0/auto option (NOTES_A "Unverified"
   warning). If a port-file handshake exists, prefer it and revisit ``NOTES_C``
   (it would remove the D5 race).


The thin-line run
=================

#. **Build the flowchart** (SEAMM GUI or a ``.flow`` file):

   #. **Model Chemistry** step -- select ``MOPAC:SQM@PM6-ORG``. (Leave
      "periodic" off for the first run.)
   #. **LAMMPS** step -- inside it: an **Initialization** sub-step, then a short
      **molecular dynamics** sub-step (NVE or NVT, ~50-100 steps, ~0.5 fs).

#. **Pick a tiny non-periodic system** -- e.g. a single water molecule or a
   small organic molecule, charge 0, multiplicity 1. Set charge/multiplicity on
   the configuration (these flow to the engine via
   ``configuration.charge`` / ``configuration.spin_multiplicity`` -- D6).

#. **Run the job.**


What to check
=============

In the LAMMPS step's working directory:

* ``mdi_launch.sh`` exists and looks like::

      #!/bin/bash
      set -e
      '<conda>' run --live-stream -n seamm-mopac python .../mopac_mdi.py -mdi '-role ENGINE ... -port NNNNN -hostname localhost' --method PM6-ORG ... &
      ENGINE_PID=$!
      mpirun -np 1 lmp -mdi "-role DRIVER -name LAMMPS -method TCP -port NNNNN -hostname localhost" -in input.dat
      wait $ENGINE_PID

  -- the **same port NNNNN** on both lines (the rendezvous), engine backgrounded
  first.

* ``input.dat`` contains ``fix mdi_fix all mdi/qm elements ...`` and **no**
  ``pair_style`` (the QM engine provides the forces).

* ``stdout.txt`` / ``stderr.txt`` -- the engine's "MDI connection established"
  appears, then LAMMPS runs the requested number of steps with finite energies
  and forces. No "address already in use".

* ``log.lammps`` -- the MD loop completed; temperature/energy are physically
  sane (not NaN, not exploding on step 1).

Success = the MD trajectory advances the requested steps with forces that
clearly came from PM6-ORG (e.g. reasonable bond lengths / energies), and the
engine exits cleanly when LAMMPS sends ``EXIT``.


Troubleshooting
===============

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Symptom
     - Likely cause
     - Fix / where to look
   * - ``ValueError: ... cannot be driven via MDI``
     - Method not in ``_MDI_CAPABLE_METHODS`` (gating in ``_mdi_engine_launch``)
     - Pick an MDI-capable method (PM7, PM6-D3H4, PM6-ORG, PM6, AM1, RM1).
   * - ``ValueError: ... not validated for periodic``
     - Periodic system + method not in ``_MDI_PERIODIC_VALIDATED``
     - Use a non-periodic system, or a method in ``_MDI_PERIODIC_VALIDATED``
       (mopac_step).
   * - ``NotImplementedError`` from the engine command
     - ``mopac.ini`` installation != conda
     - Set the executor section to ``installation = conda`` (non-conda is a
       Phase B+ item, NOTES_A).
   * - Engine dies "address already in use"
     - Port race (TOCTOU, D5)
     - Re-run. Persisting under concurrency -> add the bind-failure retry
       (NOTES_A open Q3).
   * - Hangs forever at connection
     - ``-port`` / ``-hostname`` mismatch, or MDI TCP port semantics (NOTES_A)
     - Confirm both lines in ``mdi_launch.sh`` share the port; verify MDI TCP
       needs matching ports.
   * - ``fix mdi/qm`` unknown in LAMMPS
     - LAMMPS built without the MDI package
     - Rebuild seamm-lammps with ``-DPKG_MDI=on``.
   * - Engine ImportError under conda run
     - Engine imported a seamm / mopac_step module not in seamm-mopac
     - Keep ``mopac_mdi.py`` importing only seamm-mopac packages (NOTES_A,
       Option C).
   * - Ghost-atom / comm cutoff warnings from LAMMPS
     - ``comm_modify`` intentionally dropped on the QM path (NOTES_C)
     - Revisit ``MDI_QM_input`` -- the MACE path set ``comm_modify cutoff 2.0``;
       add if needed.
   * - Wrong charge / spin in the QM energies
     - configuration charge / multiplicity unset
     - Set them on the system; they pass through D6.


After a successful run
======================

* Flip ``NOTES_C`` status to "Verified end-to-end (<date>)".
* Then the deferred items become live work: periodic systems, xTB (its own Phase
  A), non-conda launches, and the concurrent-job port retry.


References
==========

* Phase C design + decisions: ``NOTES_C.rst`` (this directory).
* Phase A engine contract / port race / MDI TCP warning:
  ``molssi-seamm.github.io`` .. ``campaigns/2026-06-22/NOTES_A.rst``.
* LAMMPS ``fix mdi/qm``: https://docs.lammps.org/fix_mdi_qm.html
* MDI Library (TCP transport): https://molssi-mdi.github.io/MDI_Library/
