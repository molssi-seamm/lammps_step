NOTES C -- The LAMMPS environment, and what belongs in it
=========================================================

Machine-learned forcefields put a large, machine-specific stack next to LAMMPS:
PyTorch, an engine, neighbour-list and equivariance libraries. Where that stack
lives, and who installs it, caused more trouble during this campaign than any
single bug.

LAMMPS must be built with MDI
-----------------------------

conda-forge's LAMMPS recipe does not enable ``PKG_MDI``. The switch is still
listed in ``lmp -h``, so it looks supported, but it is rejected at startup::

    ERROR: Invalid command-line argument: -mdi (src/lammps.cpp:486)

and the rejection happens inside ``Error::universe_all``, whose error path calls
``MPI_Barrier`` over ``MPI_COMM_WORLD``. The engine is a separate MPMD rank that
never joins that barrier, so the job does not fail -- it **hangs**, until the
wall clock runs out. A four-hour hang, from a one-line configuration difference.

``seamm-lammps.yml`` therefore takes LAMMPS from a channel that builds it with
MDI, pinned to an OpenMPI variant, for two reasons beyond MDI itself:

* ``pymdi`` must share LAMMPS's MPI *and* its MDI library. A conda build links
  ``$PREFIX/lib/libmdi.so``, which the ``pymdi`` package provides, so driver and
  engine use one library. The PyPI ``pymdi`` wheel bundles its own, which gives
  two and ``Invalid communicator`` at run time. ``pip show pymdi`` reports it
  *absent* when conda provides it, so pip will install the wheel over the top
  given the chance -- do not use an extra that depends on ``pymdi`` in this
  environment.
* ``mdi_bind.sh`` reads ``OMPI_COMM_WORLD_LOCAL_RANK``, which MPICH never sets.
  Under MPICH every rank took the rank-0 branch and bound to the same cores.

PyTorch is not in the environment file
--------------------------------------

The correct PyTorch build depends on the machine's NVIDIA driver, so it cannot
be pinned portably. A CUDA 13 wheel on a CUDA 12.2 driver imports without
complaint and then reports no GPU, so the first thing to touch the device fails
-- for a machine-learned forcefield that is a ``torch.load`` of a GPU-serialised
model, which reads as a model problem and is not.

``lammps-mdi install-ml`` exists for this: it picks the wheel from the driver,
then installs the rest in an order that stops anything re-resolving torch from
PyPI. The environment file says so and points at it.

The ordering matters and is easy to get wrong by hand:

#. torch, from the CUDA-specific index
#. the neighbour-list and equivariance libraries
#. MACE and the engine

``mace-torch`` depends on torch, so installing it first pulls whatever the
default index calls newest. Some packages cap the torch version they support,
so they must be resolved *together* with torch and from the same index, or they
drag it back down off PyPI.

What does **not** belong here
-----------------------------

``lammps_step`` must not depend on the engine package. Its dependencies are
PyTorch and the ML stack, which belong beside LAMMPS, not in the SEAMM
environment this plug-in installs into. The command template only names the
engine, so no Python dependency is needed.

This is not hypothetical: a ``make install`` that pulled a GPU extra put torch
into a SEAMM development environment as a side effect of cutting a release, and
broke that environment's test run on an unrelated torch/numpy clash. The
plug-in's own CI had it right -- ``pip install --no-deps -e .`` -- while its
Makefile did not.

Scripts are deployed, not imported
----------------------------------

``mdi_bind.sh`` and friends are copied into ``~/SEAMM/bin`` the first time the
step runs, and deliberately **not** replaced afterwards, since the CPU topology
in them is machine-specific and meant to be edited. A newer version only
produces a message. That is the right default, but it means a fix in a release
does not reach an existing installation: the file has to be deleted for the new
one to be installed.
