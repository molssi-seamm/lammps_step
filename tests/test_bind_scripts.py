#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The CPU/GPU binding helper scripts must work on any machine, not just one."""

import importlib.resources
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from lammps_step.lammps import _get_script_version

DATA = importlib.resources.files("lammps_step") / "data"
SCRIPTS = ("mdi_bind.sh", "gpu_bind.sh", "cpu_bind.sh", "mdi_monitor.sh")
SCHEDULER = ("SLURM_JOB_ID", "PBS_JOBID", "LSB_JOBID")
LINUX = sys.platform.startswith("linux") and Path("/proc/self/status").exists()


def run(script, args, tmp_path, **env):
    """Run a helper script with a clean scheduler environment plus ``env``."""
    environ = {k: v for k, v in os.environ.items() if k not in SCHEDULER}
    environ.pop("CUDA_VISIBLE_DEVICES", None)
    environ.update(env)
    return subprocess.run(
        ["bash", str(DATA / script), *args],
        cwd=tmp_path,
        env=environ,
        capture_output=True,
        text=True,
        timeout=60,
    )


@pytest.mark.parametrize("script", SCRIPTS)
def test_headers_are_read(script):
    file_type, version = _get_script_version((DATA / script).read_text())
    assert file_type == f"lammps_step:{script[:-3]}"
    assert version


@pytest.mark.parametrize("script", SCRIPTS)
def test_syntax(script):
    assert subprocess.run(["bash", "-n", str(DATA / script)]).returncode == 0


@pytest.mark.parametrize("rank", ["0", "1"])
def test_mdi_bind_leaves_binding_to_the_scheduler(tmp_path, rank):
    result = run(
        "mdi_bind.sh",
        ["sh", "-c", "echo ran; exit 3"],
        tmp_path,
        SLURM_JOB_ID="123",
        OMPI_COMM_WORLD_LOCAL_RANK=rank,
    )
    assert result.stdout.strip() == "ran"
    assert result.returncode == 3  # the command's own status
    assert "CPUs not bound" in result.stderr


def test_mdi_bind_rank_from_other_launchers(tmp_path):
    """MPICH/Hydra sets MPI_LOCALRANKID, not OpenMPI's variable."""
    result = run(
        "mdi_bind.sh", ["true"], tmp_path, SEAMM_BIND="none", MPI_LOCALRANKID="1"
    )
    assert result.returncode == 0
    assert result.stderr.startswith("Driver (rank 1)")


def test_mdi_bind_can_be_turned_off(tmp_path):
    result = run("mdi_bind.sh", ["true"], tmp_path, SEAMM_BIND="none")
    assert result.returncode == 0
    assert "CPUs not bound" in result.stderr


@pytest.mark.skipif(not LINUX or not shutil.which("taskset"), reason="needs taskset")
def test_mdi_bind_explicit_cpus(tmp_path):
    result = run(
        "mdi_bind.sh",
        ["sh", "-c", "grep Cpus_allowed_list /proc/self/status"],
        tmp_path,
        SEAMM_ENGINE_CPUS="0",
        OMPI_COMM_WORLD_LOCAL_RANK="0",
    )
    assert result.returncode == 0
    assert result.stdout.split()[-1] == "0"


def _expand(text):
    cpus = set()
    for part in text.split(","):
        a, _, b = part.partition("-")
        cpus.update(range(int(a), int(b or a) + 1))
    return cpus


@pytest.mark.skipif(not LINUX or not shutil.which("taskset"), reason="needs taskset")
@pytest.mark.parametrize("rank", ["0", "1"])
def test_mdi_bind_stays_within_the_allowed_cpus(tmp_path, rank):
    """Without a scheduler it binds, but only to CPUs this process may use."""
    allowed = _expand(
        next(
            line.split()[-1]
            for line in Path("/proc/self/status").read_text().splitlines()
            if line.startswith("Cpus_allowed_list")
        )
    )
    result = run(
        "mdi_bind.sh",
        ["sh", "-c", "grep Cpus_allowed_list /proc/self/status"],
        tmp_path,
        OMPI_COMM_WORLD_LOCAL_RANK=rank,
    )
    assert result.returncode == 0, result.stderr
    bound = _expand(result.stdout.split()[-1])
    assert bound and bound <= allowed


def test_gpu_bind_leaves_binding_to_the_scheduler(tmp_path):
    result = run(
        "gpu_bind.sh",
        ["sh", "-c", "echo $CUDA_VISIBLE_DEVICES"],
        tmp_path,
        SLURM_JOB_ID="123",
        SEAMM_GPUS="0",
        CUDA_VISIBLE_DEVICES="3",
    )
    assert result.returncode == 0
    assert result.stdout.strip() == "3"  # the scheduler's allocation is kept
    assert "CPUs not bound" in result.stderr


def test_gpu_bind_needs_seamm_gpus(tmp_path):
    assert run("gpu_bind.sh", ["true"], tmp_path).returncode == 1


def test_cpu_bind_leaves_binding_to_the_scheduler(tmp_path):
    result = run(
        "cpu_bind.sh",
        ["sh", "-c", "echo ran"],
        tmp_path,
        SLURM_JOB_ID="123",
        SEAMM_NP="4",
        SEAMM_DEBUG="1",
    )
    assert result.returncode == 0
    assert result.stdout.strip() == "ran"
    assert "not bound" in result.stderr
