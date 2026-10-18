#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The CPU/GPU binding helper scripts must work on any machine, not just one."""

import importlib.resources
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

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


def start_pair(tmp_path, engine, driver, namespace="test-pair", path=None, **env):
    """Start the engine (rank 0) and driver (rank 1) of one MDI pair, as mpirun would,
    sharing a working directory and the launch's namespace (each mpirun has its
    own)."""
    procs = []
    for rank, cmd in (("0", engine), ("1", driver)):
        environ = {k: v for k, v in os.environ.items() if k not in SCHEDULER}
        environ.pop("CUDA_VISIBLE_DEVICES", None)
        environ.update(
            SEAMM_BIND="none",
            OMPI_COMM_WORLD_LOCAL_RANK=rank,
            PMIX_NAMESPACE=namespace,
            **env,
        )
        if path:
            environ["PATH"] = f"{path}:{environ['PATH']}"
        procs.append(
            subprocess.Popen(
                ["bash", str(DATA / "mdi_bind.sh"), *cmd],
                cwd=tmp_path,
                env=environ,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        )
    return procs


def finish(procs, timeout=40):
    import time

    start = time.monotonic()
    codes = [p.wait(timeout=timeout) for p in procs]
    return codes, time.monotonic() - start


def test_a_failed_driver_stops_the_engine(tmp_path):
    """A LAMMPS error must not leave the engine holding the GPU (mlff, ChemAI 11977)."""
    engine, driver = start_pair(
        tmp_path, ["sleep", "300"], ["bash", "-c", "sleep 1; exit 1"]
    )
    (e_rc, d_rc), seconds = finish([engine, driver])
    assert d_rc == 1 and e_rc != 0
    assert seconds < 20  # not the engine's 300 s
    assert not list(tmp_path.glob(".mdi_bind.*.pid"))


def test_a_lingering_engine_is_stopped_after_the_grace(tmp_path):
    """After a normal finish the engine gets SEAMM_ENGINE_GRACE seconds to exit."""
    engine, driver = start_pair(
        tmp_path,
        ["sleep", "300"],
        ["bash", "-c", "sleep 1; exit 0"],
        SEAMM_ENGINE_GRACE="2",
    )
    (e_rc, d_rc), seconds = finish([engine, driver])
    assert d_rc == 0 and e_rc != 0
    assert 2 <= seconds < 20


def test_a_failed_engine_stops_the_driver(tmp_path):
    engine, driver = start_pair(
        tmp_path, ["bash", "-c", "sleep 1; exit 3"], ["sleep", "300"]
    )
    (e_rc, d_rc), seconds = finish([engine, driver])
    assert e_rc == 3 and d_rc != 0
    assert seconds < 20


def test_a_normal_pair_is_untouched(tmp_path):
    engine, driver = start_pair(
        tmp_path, ["bash", "-c", "sleep 1; exit 0"], ["bash", "-c", "sleep 2; exit 0"]
    )
    (e_rc, d_rc), seconds = finish([engine, driver])
    assert (e_rc, d_rc) == (0, 0)


def alive(pattern):
    """Whether a process whose full command line is `pattern` is running."""
    return (
        subprocess.run(["pgrep", "-f", "-x", pattern], capture_output=True).returncode
        == 0
    )


def test_the_partners_program_is_stopped_not_just_its_wrapper(tmp_path):
    """2.1 recorded the PID of a shell around the program, so the engine itself
    survived the stop (an orphan holding the GPU)."""
    engine, driver = start_pair(
        tmp_path, ["sleep", "303"], ["bash", "-c", "sleep 1; exit 1"]
    )
    (e_rc, d_rc), seconds = finish([engine, driver])
    assert d_rc == 1 and seconds < 20
    assert not alive("sleep 303")


def test_lammps_hanging_after_an_error_is_stopped(tmp_path):
    """LAMMPS in an MDI pair prints its ERROR and then waits in MPI_Finalize for
    the engine, so neither side ever exits (TinkerCliffs 7944053)."""
    engine, driver = start_pair(
        tmp_path,
        ["sleep", "304"],
        ["bash", "-c", 'echo "ERROR: Unknown command: bogus (src/input.cpp:314)"; '
                       "exec sleep 305"],  # fmt: skip
        SEAMM_ERROR_GRACE="2",
    )
    (e_rc, d_rc), seconds = finish([engine, driver])
    assert d_rc != 0 and e_rc != 0
    assert seconds < 30
    assert not alive("sleep 304") and not alive("sleep 305")
    out = driver.stdout.read()
    assert "ERROR: Unknown command" in out  # the output still gets through


def test_output_without_errors_is_untouched(tmp_path):
    engine, driver = start_pair(
        tmp_path,
        ["bash", "-c", "sleep 1; exit 0"],
        ["bash", "-c", "echo line one; echo WARNING: fine; echo line three; exit 0"],
    )
    (e_rc, d_rc), seconds = finish([engine, driver])
    assert (e_rc, d_rc) == (0, 0)
    assert driver.stdout.read().splitlines() == [
        "line one",
        "WARNING: fine",
        "line three",
    ]


def test_an_error_on_stderr_is_caught(tmp_path):
    engine, driver = start_pair(
        tmp_path,
        ["sleep", "307"],
        ["bash", "-c", 'echo "ERROR: on stderr" >&2; exec sleep 308'],
        SEAMM_ERROR_GRACE="2",
    )
    (e_rc, d_rc), seconds = finish([engine, driver])
    assert seconds < 30
    assert not alive("sleep 307") and not alive("sleep 308")
    assert "ERROR: on stderr" in driver.stderr.read()


def test_the_driver_output_is_byte_for_byte(tmp_path):
    """mlff parses thermo blocks from the job's stdout: nothing may change."""
    text = "Step Temp\n  1  300.5\ttab\\slash\r\nno final newline"
    script = f"printf '%s' '{text}'; printf 'warn\\nend' >&2"
    engine, driver = start_pair(
        tmp_path, ["bash", "-c", "sleep 1"], ["bash", "-c", script]
    )
    finish([engine, driver])
    out = driver.stdout.read()
    err = driver.stderr.read()
    assert (
        out
        == subprocess.run(
            ["bash", "-c", f"printf '%s' '{text}'"], capture_output=True, text=True
        ).stdout
    )
    # The wrapper adds its own status lines to stderr; the program's bytes are intact
    assert "warn\nend" in err


def test_the_gpu_monitor_is_cleaned_up(tmp_path):
    """The engine side's nvidia-smi monitor must not outlive the engine."""
    fake = tmp_path / "bin"
    fake.mkdir()
    (fake / "nvidia-smi").write_text(
        "#!/bin/bash\n"
        'case "$*" in *--query-gpu=timestamp*) exec sleep 309;; esac\n'
        "exit 0\n"
    )
    (fake / "nvidia-smi").chmod(0o755)
    for driver in (["bash", "-c", "sleep 1; exit 1"], ["bash", "-c", "sleep 1"]):
        engine, drv = start_pair(
            tmp_path, ["sleep", "310"], driver, path=str(fake), SEAMM_ENGINE_GRACE="2"
        )
        finish([engine, drv])
        assert not alive("sleep 309") and not alive("sleep 310")


def test_two_pairs_in_one_directory_keep_apart(tmp_path):
    """md2.slurm runs two mpiruns on GPU 0 in one step: one pair's failure must not
    stop the other's engine."""
    failing = start_pair(
        tmp_path,
        ["sleep", "311"],
        ["bash", "-c", "sleep 1; exit 1"],
        namespace="prterun-node-1001_1",
    )
    healthy = start_pair(
        tmp_path,
        ["bash", "-c", "sleep 4; exit 0"],
        ["bash", "-c", "sleep 5; exit 0"],
        namespace="prterun-node-1002_1",
    )
    codes, _ = finish(failing)
    assert codes[1] == 1 and not alive("sleep 311")
    codes, _ = finish(healthy)
    assert codes == [0, 0]  # not stopped by the other pair


def test_mdi_bind_role_without_mpirun(tmp_path):
    """A pair over TCP has no MPI rank, so it names its roles."""
    result = run(
        "mdi_bind.sh", ["true"], tmp_path, SEAMM_BIND="none", SEAMM_MDI_ROLE="driver"
    )
    assert result.returncode == 0
    assert result.stderr.startswith("Driver (rank 1)")
    result = run(
        "mdi_bind.sh",
        ["true"],
        tmp_path,
        SEAMM_BIND="none",
        SEAMM_MDI_ROLE="engine",
        OMPI_COMM_WORLD_LOCAL_RANK="1",
    )
    assert result.stderr.startswith("Engine (rank 0)")
    result = run("mdi_bind.sh", ["true"], tmp_path, SEAMM_MDI_ROLE="both")
    assert result.returncode == 2
    assert "SEAMM_MDI_ROLE" in result.stderr


def _tcp_pair(tmp_path, engine, driver, **env):
    """Run an engine and a driver through mdi_bind.sh as a pair not started by
    mpirun, the way lammps_step's TCP launch script does."""
    from lammps_step.lammps import _mlff_tcp_script

    bind = DATA / "mdi_bind.sh"
    script = _mlff_tcp_script(engine, 1, {"mlff-tcp-code": driver}, bind)
    (tmp_path / "launch.sh").write_text(script)
    environ = {k: v for k, v in os.environ.items() if k not in SCHEDULER}
    environ.update(SEAMM_BIND="none", SEAMM_ENGINE_GRACE="2", **env)
    return subprocess.run(
        ["bash", "launch.sh"],
        cwd=tmp_path,
        env=environ,
        capture_output=True,
        text=True,
        timeout=60,
    )


def test_tcp_pair_engine_failure_stops_the_driver(tmp_path):
    """An engine that dies starting must not leave LAMMPS waiting on its port."""
    t0 = time.monotonic()
    result = _tcp_pair(tmp_path, ["sh", "-c", "sleep 1; exit 4"], "sh -c 'sleep 50' x")
    assert time.monotonic() - t0 < 30
    assert result.returncode != 0
    assert "Stopping the partner process" in result.stderr


def test_tcp_pair_driver_failure_stops_the_engine_and_its_children(tmp_path):
    """The engine runs under "conda run", which does not pass on signals, so its
    child must be stopped too."""
    marker = "317.5"  # an unusual sleep, to find any leftover
    t0 = time.monotonic()
    result = _tcp_pair(
        tmp_path, ["sh", "-c", f"sleep {marker} & wait"], "sh -c 'sleep 1; exit 3' x"
    )
    assert time.monotonic() - t0 < 30
    assert result.returncode == 3
    time.sleep(1)
    leftover = subprocess.run(["pgrep", "-f", f"sleep {marker}"], capture_output=True)
    assert leftover.stdout == b""


def test_tcp_pair_driver_that_never_starts_stops_the_engine(tmp_path):
    """If the driver fails before mdi_bind.sh runs (a broken wrapper, say), nothing
    stops an engine waiting for its connection but the launch script itself."""
    from lammps_step.lammps import _mlff_tcp_script

    marker = "318.5"  # an unusual sleep, to find any leftover
    bind = DATA / "mdi_bind.sh"
    script = _mlff_tcp_script(
        ["sh", "-c", f"sleep {marker} & wait"], 1, {"mlff-tcp-code": "x"}, bind
    )
    lines = script.splitlines()
    i = next(n for n, line in enumerate(lines) if "SEAMM_MDI_ROLE=driver" in line)
    lines[i] = "no-such-wrapper-xyz"  # the driver fails to start: rc 127
    (tmp_path / "launch.sh").write_text("\n".join(lines) + "\n")
    environ = {k: v for k, v in os.environ.items() if k not in SCHEDULER}
    environ.update(SEAMM_BIND="none", SEAMM_ENGINE_GRACE="2")
    t0 = time.monotonic()
    result = subprocess.run(
        ["bash", "launch.sh"],
        cwd=tmp_path,
        env=environ,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert time.monotonic() - t0 < 30
    assert result.returncode == 127
    time.sleep(1)
    leftover = subprocess.run(["pgrep", "-f", f"sleep {marker}"], capture_output=True)
    assert leftover.stdout == b""


def test_tcp_pair_success(tmp_path):
    result = _tcp_pair(tmp_path, ["sh", "-c", "sleep 1"], "sh -c 'sleep 2' x")
    assert result.returncode == 0, result.stderr


def test_install_helper_script(tmp_path):
    from lammps_step.lammps import _install_helper_script

    dest = tmp_path / "mdi_bind.sh"
    v1 = "#!/bin/bash\n#MolSSI lammps_step:mdi_bind 1.3\necho old\n"
    v2 = "#!/bin/bash\n#MolSSI lammps_step:mdi_bind 2.10\necho new\n"
    assert "Installed" in _install_helper_script(dest, v1)
    assert os.access(dest, os.X_OK)
    message = _install_helper_script(dest, v2)
    assert "from version 1.3 to 2.10" in message
    assert dest.read_text() == v2
    assert (tmp_path / "mdi_bind.sh.1.3").read_text() == v1
    assert os.access(dest, os.X_OK)
    # The same or an older version, no header, or another script: left alone
    assert _install_helper_script(dest, v2) is None
    assert _install_helper_script(dest, v1) is None
    assert dest.read_text() == v2
    for text in (
        "#!/bin/bash\necho mine\n",
        "#!/bin/bash\n#MolSSI lammps_step:gpu_bind 1.0\n",
    ):
        dest.write_text(text)
        assert _install_helper_script(dest, v2) is None
        assert dest.read_text() == text
    assert not list(tmp_path.glob(".*.tmp"))


# ----------------------------------------------------------------------
# TCP_NODELAY on MDI's TCP connections
# ----------------------------------------------------------------------
CC = shutil.which(os.environ.get("CC", "cc").split()[0])
NEEDS_CC = pytest.mark.skipif(
    not (LINUX and CC), reason="the shim needs Linux and a C compiler"
)


def test_tcp_scripts_carry_the_nodelay_preamble_unless_turned_off():
    from lammps_step.lammps import (
        NODELAY_SOURCE,
        _mdi_launch_script,
        _mlff_pair_script,
        _mlff_tcp_script,
    )

    qm = {"code": "lmp", "installation": "conda"}
    for config in ({}, {"tcp-nodelay": "yes"}):
        script = _mlff_tcp_script(["xnn"], 5555, config, "b")
        assert NODELAY_SOURCE in script
        # before either side starts
        assert script.index(NODELAY_SOURCE) < script.index("SEAMM_MDI_ROLE=engine")
        script = _mdi_launch_script(["mopac"], 5555, {**qm, **config}, {})
        assert script.index(NODELAY_SOURCE) < script.index("mopac &")
    for value in ("no", "No", "false", "off", "0"):
        config = {"tcp-nodelay": value}
        assert NODELAY_SOURCE not in _mlff_tcp_script(["xnn"], 5555, config, "b")
        assert NODELAY_SOURCE not in _mdi_launch_script(
            ["mopac"], 5555, {**qm, **config}, {}
        )
    # An MPI pair has no TCP connection
    assert NODELAY_SOURCE not in _mlff_pair_script(["xnn"], {}, "b")


def test_the_nodelay_source_is_packaged():
    from lammps_step.lammps import _nodelay_source

    text = _nodelay_source()
    for name in ("TCP_NODELAY", "int connect(", "int accept(", "int accept4("):
        assert name in text


SOCKET_PAIR = """
import socket
server = socket.socket()
server.bind(("127.0.0.1", 0))
server.listen(1)
client = socket.create_connection(server.getsockname())
accepted, _ = server.accept()
option = (socket.IPPROTO_TCP, socket.TCP_NODELAY)
print(client.getsockopt(*option) != 0, accepted.getsockopt(*option) != 0)
"""


@NEEDS_CC
def test_the_shim_sets_tcp_nodelay_on_both_ends(tmp_path):
    from lammps_step.lammps import _nodelay_source

    (tmp_path / "mdi_nodelay.c").write_text(_nodelay_source())
    library = tmp_path / "mdi_nodelay.so"
    subprocess.run(
        [CC, "-shared", "-fPIC", "-O2", "-o", str(library), "mdi_nodelay.c", "-ldl"],
        cwd=tmp_path,
        check=True,
    )

    def sockets(**env):
        environ = {k: v for k, v in os.environ.items() if k != "LD_PRELOAD"}
        environ.update(env)
        result = subprocess.run(
            [sys.executable, "-c", SOCKET_PAIR],
            env=environ,
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.split()

    assert sockets() == ["False", "False"]
    assert sockets(LD_PRELOAD=str(library)) == ["True", "True"]


def _preload_seen(tmp_path, **env):
    """Run the TCP launch script with an engine that records its LD_PRELOAD."""
    from lammps_step.lammps import NODELAY_SOURCE, _nodelay_source

    (tmp_path / NODELAY_SOURCE).write_text(_nodelay_source())
    engine = ["sh", "-c", 'echo "$LD_PRELOAD" > preload.txt']
    result = _tcp_pair(tmp_path, engine, "true", **env)
    return result, (tmp_path / "preload.txt").read_text().strip()


@NEEDS_CC
def test_the_tcp_script_preloads_the_shim(tmp_path):
    result, preload = _preload_seen(tmp_path, LD_PRELOAD="")
    assert result.returncode == 0, result.stderr
    assert preload.endswith("/mdi_nodelay.so")
    assert not (tmp_path / "mdi_nodelay.err").exists()
    # the compiled library is removed when the script ends
    assert not Path(preload).exists()


@pytest.mark.skipif(not LINUX, reason="the shim is only used on Linux")
def test_the_run_goes_on_without_a_compiler(tmp_path):
    result, preload = _preload_seen(tmp_path, LD_PRELOAD="", CC="no-such-cc-xyz")
    assert result.returncode == 0
    assert preload == ""
    assert "Warning: could not compile mdi_nodelay.c" in result.stderr
    # with the compiler's own message
    assert "no-such-cc-xyz" in result.stderr


@pytest.mark.skipif(LINUX, reason="checks the behaviour off Linux")
def test_nothing_is_preloaded_off_linux(tmp_path):
    result, preload = _preload_seen(tmp_path, LD_PRELOAD="")
    assert result.returncode == 0
    assert preload == ""
    assert "Warning" not in result.stderr
