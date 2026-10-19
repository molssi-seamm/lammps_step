"""The environment for runs on GPUs must not bind the processes to CPUs."""

from lammps_step.lammps import _gpu_environment


def test_gpu_environment_has_the_gpus_and_one_thread():
    env = _gpu_environment([0, 2])
    assert env["SEAMM_GPUS"] == "0,2"
    assert env["OMP_NUM_THREADS"] == "1"


def test_gpu_environment_does_not_bind_cpus():
    # OMP_PROC_BIND/OMP_PLACES pinned every engine and driver in a job to one core
    env = _gpu_environment([0, 1, 2, 3])
    assert "OMP_PROC_BIND" not in env
    assert "OMP_PLACES" not in env


def test_every_launch_is_single_threaded():
    """An inherited OMP_NUM_THREADS, e.g. from a parallel Loop, is overridden: each
    MPI rank must have one thread."""
    from lammps_step.lammps import _launch_environment

    for gpus in ([], [1]):
        env = _launch_environment(gpus)
        assert env["OMP_NUM_THREADS"] == "1"
        assert env["MKL_NUM_THREADS"] == "1"
    assert "SEAMM_GPUS" not in _launch_environment([])
    assert _launch_environment([1])["SEAMM_GPUS"] == "1"
