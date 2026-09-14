"""An MDI gpu-code declares its own rank count; the allocation must cover it.

The engine and the LAMMPS driver are separate MPI ranks written literally into
the command, not derived from the allocation. Given fewer tasks, mpirun dies
with "All nodes which are allocated for this job are already filled", which
mentions neither MDI nor ntasks.
"""

import pytest

from lammps_step.lammps import _allocated_mpi_slots, _required_mpi_slots

MDI = (
    'mpirun -np 1 ~/SEAMM/bin/mdi_bind.sh mace-mdi -mdi "-role ENGINE" '
    ': -np 1 ~/SEAMM/bin/mdi_bind.sh lmp -mdi "-role DRIVER"'
)


def test_mdi_pair_needs_the_sum_not_the_max():
    """Both segments run at once, so 1 + 1 is 2 -- the bug that bit job 4778."""
    assert _required_mpi_slots(MDI) == 2


def test_colons_need_no_surrounding_spaces():
    assert _required_mpi_slots("mpirun -np 1 a: -np 2 b") == 3


def test_single_segment():
    assert _required_mpi_slots("mpirun -np 4 lmp") == 4


def test_short_and_long_flags():
    assert _required_mpi_slots("mpirun -n 2 a : --np 3 b") == 5


def test_placeholder_is_not_guessed():
    """{NTASKS} is resolved later from the allocation, so it cannot conflict
    with it; returning a number here would be inventing one."""
    assert _required_mpi_slots("mpirun -np {NTASKS} lmp") is None


def test_no_rank_count_at_all():
    assert _required_mpi_slots("lmp -in input.dat") is None


@pytest.mark.parametrize("variable", ["SLURM_NTASKS", "SLURM_NPROCS"])
def test_allocation_read_from_either_variable(monkeypatch, variable):
    monkeypatch.delenv("SLURM_NTASKS", raising=False)
    monkeypatch.delenv("SLURM_NPROCS", raising=False)
    monkeypatch.setenv(variable, "2")
    assert _allocated_mpi_slots() == 2


def test_no_batch_system(monkeypatch):
    monkeypatch.delenv("SLURM_NTASKS", raising=False)
    monkeypatch.delenv("SLURM_NPROCS", raising=False)
    assert _allocated_mpi_slots() is None


def test_junk_value_is_ignored(monkeypatch):
    monkeypatch.delenv("SLURM_NPROCS", raising=False)
    monkeypatch.setenv("SLURM_NTASKS", "not-a-number")
    assert _allocated_mpi_slots() is None


def test_the_failing_and_working_cases():
    """Job 4778 had ntasks=1 and failed; 4779 had 2 and ran."""
    assert _required_mpi_slots(MDI) > 1
    assert _required_mpi_slots(MDI) <= 2
