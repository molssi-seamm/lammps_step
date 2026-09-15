NOTES A -- How a model reaches its engine
=========================================

The question
------------

A machine-learned forcefield is not evaluated by LAMMPS. LAMMPS drives a
separate engine process over MDI, ships it the atoms, and takes back energy and
forces. Three things therefore have to line up, and they belong to different
owners:

=========================  ==============  ===================================
question                   owned by        where it lives
=========================  ==============  ===================================
*which* model              the flowchart   ``_pytorch_model`` (a path)
*which engine* serves it   the model       -- nothing, before this campaign
*how* to launch it here    the machine     ``gpu-code`` in ``lammps.ini``
=========================  ==============  ===================================

The middle row was missing, and its absence was filled by guessing from the
model's file name.

What the code used to do
------------------------

``initialization.py`` chose between an MDI run and a native pair style on the
name of the model file:

.. code-block:: python

   if "mliap" in model:
       pair_style mliap ...
   elif model.endswith(".mace.pt"):
       fix mdi/qm ...                 # the MDI engine
   else:
       pair_style mace ...            # LAMMPS evaluates it itself

and ``lammps.py`` hard-coded one engine per machine, because ``gpu-code`` was a
single template and ``SEAMM_FF`` -- the variable ``mace-mdi`` reads -- was the
only way a model reached an engine.

That works for exactly as long as there is one engine. A file name says nothing
about how a model is evaluated: it is chosen by whoever trained it, and it
changes when someone renames a file. An ``xnns`` checkpoint named for what it
is, ``2026-08-08.mace_les_xnns.pt``, does not end in ``.mace.pt``, so it took
the ``else`` branch and died with::

    ERROR: Unrecognized pair style 'mace' (src/force.cpp:275)

Reaching that pair style needs a LAMMPS built with a MACE pair style -- which
is the build the MDI route exists so as not to need. Nobody was using it.

The same test also gated the ``comm_modify`` line that an MDI run needs, so it
was missing from exactly the runs that had just been routed to MDI. That would
have been the next failure.

What it does now
----------------

**The plug-in decides whether a model goes to an engine, not which engine.**
Everything except an ML-IAP model -- which LAMMPS really does evaluate itself
-- is served over MDI:

.. code-block:: python

   if "mliap" in model:
       pair_style mliap ...
   else:
       fix mdi/qm ...

**The machine decides which engine.** ``gpu-code`` is the only place that
names one. To make that possible the model is offered both ways engines take
it:

``SEAMM_FF``
    An environment variable. ``mace-mdi`` reads it.

``{MODEL}``
    Substituted into the command like ``{NTASKS}``, by the executor's existing
    mechanism (``command.format(**config, **ce)``). For engines taking the
    model as an argument: ``xnns mdi --ckpt {MODEL}``.

Both are always set, so nothing in ``lammps_step`` knows or needs to know which
engine is configured. Adding a third engine is a ``lammps.ini`` edit.

Why not detect the model type
-----------------------------

Two alternatives were considered and rejected.

*Keep using the file name, but look for ``xnns`` in it.* It would have worked
that day -- the files are named that way -- and it is one line. It was rejected
because the failure mode is silent: rename a file and the run starts the wrong
engine, which gives ``KeyError: 'hidden_irreps'`` if you are lucky and wrong
numbers if you are not.

*Identify the model by reading it.* This is feasible without torch, which the
SEAMM environment does not have and must not have: a ``.pt`` file is a zip, and
the pickled module references distinguish the formats. A prototype identified
an ``xnns`` checkpoint, a ``mace-torch`` model and a deployed TorchScript
archive correctly using only the standard library.

It was rejected because it answers a question nobody needs answered. Detection
matters only while more than one engine is in play, and the decision below
means there is one.

The decision: ``xnns`` is the mechanism
---------------------------------------

``xnns`` serves every model family it knows -- MACE, NequIP, Allegro, CACE,
SchNet, ANI, PhysNet, HDNNP, BAMBOO -- through one interface, and its MDI
engine speaks the same protocol ``lammps_step`` already drives. It is therefore
the single mechanism for machine-learned forcefields, and where it does not yet
cover something, that is a gap to close in ``xnns`` rather than a reason to add
a second path here.

One consequence is worth recording, because it would otherwise look like an
oversight: ``xnns`` cannot load a ``mace-torch`` model. ``MDIEngine.from_checkpoint``
expects the trainer's ``{"model": ..., "cfg": ...}`` and raises
``TypeError: 'ScaleShiftMACE' object is not subscriptable`` on a bare model
object. This did not need solving, because the older MACE models were being
retrained -- they carried a density error from the DFT functional used and were
not scientifically meaningful. Had they been kept, the choice would have been
between a transitional fallback to ``mace-mdi`` and asking for the loader to be
generalised.

``mace-mdi`` remains usable: it is a ``gpu-code`` line, and the shipped
``lammps.ini`` documents it as the MACE-specific alternative.

What this leaves
----------------

* File names no longer decide behaviour anywhere in this path. Checkpoints can
  be named for what they are.
* A new engine is a configuration change, not a code change.
* ``ff_form()`` still returns ``"PyTorch"``. The engine is a sub-distinction of
  it, and widening ``ff_form`` would have rippled into the classical/Kokkos/MDI
  logic for no gain.
