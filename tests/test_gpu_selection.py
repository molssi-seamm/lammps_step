"""GPU selection must respect CUDA_VISIBLE_DEVICES.

GPUtil reads nvidia-smi and so reports *physical* GPU indices, knowing nothing
about CUDA_VISIBLE_DEVICES. A scheduler that allocates a subset of a node's
GPUs sets that variable, which also renumbers the devices -- so the physical
index GPUtil reports is not the index the job should use.
"""

import pytest

from lammps_step.lammps import _cuda_visible_devices, _usable_gpus


@pytest.fixture
def cvd(monkeypatch):
    def _set(value):
        if value is None:
            monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
        else:
            monkeypatch.setenv("CUDA_VISIBLE_DEVICES", value)

    return _set


def test_unset_means_every_gpu_visible(cvd):
    cvd(None)
    assert _cuda_visible_devices() is None


def test_parses_indices_in_order(cvd):
    cvd("2,0")
    assert _cuda_visible_devices() == [2, 0]


def test_empty_means_no_gpus(cvd):
    """SLURM sets this to empty for a job allocated no GPUs."""
    cvd("")
    assert _cuda_visible_devices() == []
    assert _usable_gpus([0, 1], []) == []


def test_uuid_entries_are_unmappable(cvd):
    cvd("GPU-3221d296,GPU-5882380f")
    assert _cuda_visible_devices() == [None, None]


def test_without_the_variable_physical_indices_pass_through():
    assert _usable_gpus([0, 1], None) == [0, 1]


def test_allocated_gpu_is_renumbered():
    """The regression: allocated physical GPU 1, GPUtil reports it as 1, but
    the job must use device 0 -- returning 1 would run on a GPU it does not
    have."""
    assert _usable_gpus([0, 1], [1]) == [0]


def test_busy_allocated_gpu_is_excluded():
    """Our GPU is not in GPUtil's available list, so nothing is usable -- we
    must not fall back to a physical index we were not given."""
    assert _usable_gpus([0], [1]) == []


def test_order_follows_the_variable_not_the_hardware():
    assert _usable_gpus([0, 1], [1, 0]) == [0, 1]
    assert _usable_gpus([1], [1, 0]) == [0]
    assert _usable_gpus([0], [1, 0]) == [1]


def test_uuid_allocation_is_used_whole():
    """Unmappable, but discarding the allocation would be worse."""
    assert _usable_gpus([0, 1], [None, None]) == [0, 1]
