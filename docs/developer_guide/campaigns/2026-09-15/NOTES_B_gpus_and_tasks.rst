NOTES B -- GPUs, tasks, and what the scheduler is telling us
============================================================

Making GPUs requestable per job turned out to expose three assumptions in this
plug-in that had held only because nothing had ever asked for a GPU explicitly.
They are recorded together because they share a root: **an allocation describes
what a job was given, not what the calculation needs.**

A GPU in the allocation is not a reason to use ``gpu-code``
-----------------------------------------------------------

``NGPUS`` comes from the scheduler. ``seamm_exec`` reads ``SLURM_JOB_GPUS``:

.. code-block:: python

   if "JOB_GPUS" in ce:
       ce["NGPUS"] = len(ce["JOB_GPUS"].split(","))

and this plug-in used its presence alone to select ``gpu-code``:

.. code-block:: python

   if "NGPUS" in ce:
       cmd = ["{gpu-code}"]

So the first classical job submitted with ``--gres=gpu:1`` -- an OPLS-AA water
equilibration -- was handed to the MACE engine, which has no model to serve.
``SEAMM_FF`` is ``"Unknown"`` for a non-PyTorch forcefield, and the engine died
with ``FileNotFoundError: 'Unknown'`` before LAMMPS ran at all.

``gpu-code`` serves two unrelated purposes: a Kokkos command is how a
*classical* forcefield uses a GPU, and an MDI command serves a machine-learned
one. The test is therefore not "are there GPUs" but "does this command suit
this forcefield", which is decided by whether it drives MDI. A classical
forcefield with an MDI ``gpu-code`` now runs on the CPU and says why.

GPU indices are not the scheduler's indices
-------------------------------------------

GPU selection went through ``GPUtil``, which reads ``nvidia-smi`` and therefore
reports *physical* indices, knowing nothing about ``CUDA_VISIBLE_DEVICES``. A
scheduler that allocates a subset of a node's GPUs sets that variable, and doing
so **renumbers** the devices: given physical GPU 1, the job's own device 0 *is*
that GPU.

The step reported 1. The calculation would have run on a card it had not been
allocated, beside whichever job actually held it, while the reserved one sat
idle. Selection now happens in the numbering CUDA itself uses.

``mdi_bind.sh`` had the same bug from the other end: it exported
``CUDA_VISIBLE_DEVICES`` unconditionally. That variable is **not composable** --
setting it again is interpreted against the machine's full set of devices, not
against the allocation -- so writing an index there was itself a way to land on
the wrong GPU. It now leaves a scheduler's value alone, and maps back to a
physical index only for the ``nvidia-smi`` monitor, which ignores the variable.

Neither bug appears where every job sees every GPU, which is why they went
unnoticed for so long: they need a scheduler handing out individual GPUs.

An MDI run needs two tasks
--------------------------

An MDI ``gpu-code`` launches the engine and the LAMMPS driver as separate MPI
ranks, and those counts are written into the command rather than derived from
the allocation. A job given fewer tasks than the command asks for dies inside
``mpirun`` with::

    All nodes which are allocated for this job are already filled.

which mentions neither MDI, nor tasks, nor anything to change -- and arrives
wrapped in ``ERROR conda.cli.main_run:execute(127)``.

This is easy to hit precisely because the obvious reading is wrong: asking for
**one GPU and one task** looks correct and cannot work. The step now counts the
ranks the command asks for, compares them with ``SLURM_NTASKS``, and says so
before launching.

Configuration consequences
--------------------------

On a machine with GPUs, ``gres`` should be *overridable* rather than defaulted.
A default makes every job request a GPU, including CPU-only work, which
serialises the queue to the number of cards regardless of any concurrency
limit. Jobs that need one ask for it.

Nothing requests a GPU automatically. A machine-learned forcefield run gets one
only if the submission includes it. Whether that should be inferred is left
open: the forcefield is often only known at run time, from a variable, so the
information is not reliably available at submission -- which is the same reason
a step cannot declare its own resource needs.
