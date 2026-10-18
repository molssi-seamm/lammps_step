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


# --------------------------------------------------------------------------- #
# Waiting for an idle GPU
# --------------------------------------------------------------------------- #


class _Clock:
    def __init__(self):
        self.now = 1.0e9  # an epoch time, as time.time() gives

    def time(self):
        return self.now

    def sleep(self, seconds):
        self.now += seconds


def _wait(monkeypatch, idle, every, visible=None, wait=None):
    """_wait_for_gpus with nvidia-smi reporting ``idle`` GPUs at the given load
    and ``every`` GPU at any load."""
    from lammps_step import lammps

    if visible is None:
        monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    else:
        monkeypatch.setenv("CUDA_VISIBLE_DEVICES", visible)
    monkeypatch.setattr(
        lammps,
        "_available_gpus",
        lambda maxload, timeout=10: list(every if maxload >= 1.0 else idle),
    )
    clock = _Clock()
    gpus = lammps._wait_for_gpus(0.1, wait, sleep=clock.sleep, clock=clock.time)
    return gpus, clock.now - 1.0e9


def test_an_idle_gpu_is_used_at_once(monkeypatch):
    assert _wait(monkeypatch, idle=[1], every=[0, 1]) == ([1], 0)


def test_no_gpu_on_the_machine_runs_on_the_cpu_at_once(monkeypatch):
    assert _wait(monkeypatch, idle=[], every=[]) == ([], 0)


def test_an_allocation_without_gpus_runs_on_the_cpu_at_once(monkeypatch):
    assert _wait(monkeypatch, idle=[], every=[0, 1], visible="") == ([], 0)


def test_busy_gpus_are_waited_for_as_long_as_asked(monkeypatch):
    """A finite wait ends (it once compared the time elapsed with the epoch)."""
    assert _wait(monkeypatch, idle=[], every=[0, 1], wait=5) == ([], 5)


# ---------------------------------------------------------------------------
# Claims: runs started together spread over the idle GPUs
# ---------------------------------------------------------------------------
def test_runs_started_together_take_different_gpus(tmp_path, monkeypatch):
    from lammps_step.lammps import _claim_gpus, _release_gpus

    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    runs = [_claim_gpus([0, 1, 2, 3], 1, directory=tmp_path) for _ in range(4)]
    assert [gpus for gpus, _ in runs] == [[0], [1], [2], [3]]

    # A fifth run shares a GPU, as before there were claims
    fifth, claims = _claim_gpus([0, 1, 2, 3], 1, directory=tmp_path)
    assert fifth == [0] and claims == []

    # A finished run frees its GPU for the next
    _release_gpus(runs[2][1])
    assert _claim_gpus([0, 1, 2, 3], 1, directory=tmp_path)[0] == [2]
    for _, claims in runs:
        _release_gpus(claims)


def test_several_gpus_for_one_run(tmp_path, monkeypatch):
    from lammps_step.lammps import _claim_gpus, _release_gpus

    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    first, a = _claim_gpus([0, 1, 2, 3], 2, directory=tmp_path)
    second, b = _claim_gpus([0, 1, 2, 3], 2, directory=tmp_path)
    assert first == [0, 1] and second == [2, 3]
    _release_gpus(a + b)


def test_claims_are_on_the_physical_gpu(tmp_path, monkeypatch):
    """Two allocations numbering their GPUs from 0 do not block each other."""
    from lammps_step.lammps import _claim_gpus, _gpu_claim_name, _release_gpus

    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "2,3")
    assert _gpu_claim_name(0).endswith("-gpu2.lock")
    held, claims = _claim_gpus([0], 1, directory=tmp_path)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "3")
    other, more = _claim_gpus([0], 1, directory=tmp_path)
    assert other == [0] and len(more) == 1  # physical 3 is free
    assert (tmp_path / _gpu_claim_name(0)).exists()
    _release_gpus(claims + more)
