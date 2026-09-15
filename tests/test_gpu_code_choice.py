"""A GPU in the allocation is not by itself a reason to use ``gpu-code``.

NGPUS comes from the scheduler, so it is set for any job that asked for a GPU,
whatever the forcefield. An MDI ``gpu-code`` launches a machine-learned engine,
which only a PyTorch forcefield has a model for; a Kokkos one is how a
classical forcefield uses a GPU.
"""

from lammps_step.lammps import _gpu_code_is_usable

MDI = (
    'mpirun -np 1 ~/SEAMM/bin/mdi_bind.sh mace-mdi -mdi "-role ENGINE" '
    ': -np 1 ~/SEAMM/bin/mdi_bind.sh lmp -mdi "-role DRIVER"'
)
KOKKOS = "mpirun -np {NTASKS} --bind-to none ~/SEAMM/bin/gpu_bind.sh lmp"


def test_classical_forcefield_cannot_use_an_mdi_command():
    """The regression: OPLS-AA was launched into the MACE engine, which got
    SEAMM_FF='Unknown' and died with FileNotFoundError."""
    assert not _gpu_code_is_usable("OPLS-AA", MDI)


def test_pytorch_forcefield_can():
    assert _gpu_code_is_usable("PyTorch", MDI)


def test_classical_forcefield_can_use_kokkos():
    """A Kokkos gpu-code is exactly how a classical forcefield uses a GPU, so
    this must not be caught by the same guard."""
    assert _gpu_code_is_usable("OPLS-AA", KOKKOS)


def test_pytorch_with_kokkos_is_left_alone():
    assert _gpu_code_is_usable("PyTorch", KOKKOS)


def test_empty_gpu_code_is_not_mdi():
    assert _gpu_code_is_usable("OPLS-AA", "")


def test_other_classical_forcefields():
    for ff in ("OPLS-AA", "AMBER", "CHARMM", "Unknown", "MDI/QM"):
        assert not _gpu_code_is_usable(ff, MDI), ff
