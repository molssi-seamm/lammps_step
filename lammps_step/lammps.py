# -*- coding: utf-8 -*-

"""A node or step for LAMMPS in a flowchart"""

import configparser
from contextlib import contextmanager
import copy
import importlib
import inspect
import json
import logging
from math import sqrt, exp, degrees, radians, cos, acos
from pathlib import Path
import os
import os.path
import platform
import pprint
import re
import shlex
import shutil
import socket
import subprocess
import string
import sys
import tempfile
import textwrap
import time
import traceback
import warnings

import bibtexparser
import numpy as np
import pandas
import statsmodels.tsa.stattools as stattools
from tabulate import tabulate

import lammps_step
import molsystem
import seamm
import seamm_exec
from seamm_ff_util import tabulate_angle
import seamm_util
import seamm_util.printing as printing
from seamm_util import CompactJSONEncoder, Configuration, units_class
from seamm_util.printing import FormattedText as __

# from pymbar import timeseries

logger = logging.getLogger("lammps")
job = printing.getPrinter()
printer = printing.getPrinter("lammps")

# statsmodels 0.15 warns that acf will return an object instead of a tuple; keep the
# tuple where the option exists.
_acf_options = (
    {"result_object": False}
    if "result_object" in inspect.signature(stattools.acf).parameters
    else {}
)

# The header of the helper scripts, e.g. "#MolSSI lammps_step:mdi_bind 2.0".
VERSION_PATTERN = re.compile(r"^#\s?MolSSI\s+(lammps_step:\S+)\s+([\d.]+)")


# Advice for common LAMMPS errors: (pattern in the message, advice)
_lammps_error_advice = (
    (
        "Lost atoms",
        "Atoms moved too far in one timestep. This usually means a poor starting "
        "structure, or too long a timestep: try minimizing the structure first, "
        "or a shorter timestep.",
    ),
    (
        "Out of range atoms",
        "Atoms moved too far for the k-space grid. This usually means a poor "
        "starting structure, or too long a timestep: try minimizing the structure "
        "first, or a shorter timestep.",
    ),
    (
        "Bond atoms missing",
        "Bonded atoms moved too far apart. This usually means a poor starting "
        "structure, or too long a timestep.",
    ),
    (
        "Non-numeric",
        "The simulation blew up (the energy or pressure is not a number). Check "
        "the starting structure and the timestep.",
    ),
    (
        "triclinic",
        "The cell is triclinic (not orthorhombic), which this command cannot "
        "handle.",
    ),
)

_lammps_error_re = re.compile(r"^ERROR(?: on proc \d+)?: (.*?)(?: \(src/\S+\))?\s*$")


def lammps_error(log, screen=""):
    """Find the error, if any, in LAMMPS's log (log.lammps) or screen output.

    Parameters
    ----------
    log : str
        The text of log.lammps.
    screen : str
        The screen output (stdout), which has errors raised by processes other
        than the first, which are not in the log.

    Returns
    -------
    dict or None
        The error "message", the "step" running (from the "# Step ..." comments
        in the input), the "last command", and any "advice"; or None if there
        was no error.
    """
    for text in (log, screen):
        lines = text.splitlines()
        for i, line in enumerate(lines):
            match = _lammps_error_re.match(line)
            if match is None:
                continue
            message = match.group(1).strip()
            last_command = ""
            for line2 in lines[i + 1 : i + 3]:  # noqa: E203
                if line2.startswith("Last command:"):
                    last_command = " ".join(line2[len("Last command:") :].split())
                    break
            step = ""
            for line2 in reversed(lines[:i]):
                if line2.startswith("# Step "):
                    step = line2[2:].strip()
                    break
            advice = ""
            for pattern, text2 in _lammps_error_advice:
                if pattern in message:
                    advice = text2
                    break
            return {
                "message": message,
                "step": step,
                "last command": last_command,
                "advice": advice,
            }
    return None


def lammps_error_message(error, directory=None):
    """A readable message for an error found by lammps_error."""
    where = f" in {error['step']}" if error["step"] != "" else ""
    text = f"LAMMPS failed{where}: {error['message']}"
    if not text.endswith("."):
        text += "."
    if error["last command"] != "":
        text += f" The last command was '{error['last command']}'."
    if error["advice"] != "":
        text += " " + error["advice"]
    if directory is not None:
        text += f" See log.lammps in {directory} for details."
    return text


def _get_script_version(text):
    """Extract file type and version from a MolSSI header line.

    Returns:
        (file_type, version) e.g. ("lammps_step:mdi_bind", "1.0")
        or (None, None) if no header found.
    """
    for line in text.splitlines()[:3]:
        m = VERSION_PATTERN.match(line)
        if m:
            return m.group(1), m.group(2)
    return None, None


def _version_tuple(version):
    """ "2.10" -> (2, 10), or None if it is not dotted integers."""
    try:
        return tuple(int(part) for part in version.split("."))
    except (AttributeError, ValueError):
        return None


def _install_helper_script(dest, text):
    """Install a helper script, or upgrade an older copy of it, in <root>/bin.

    A copy is upgraded only when its MolSSI header names the same script at a
    lower version; a copy without a header, of another script, or of the same or
    a newer version is left alone. The old copy is kept as ``<name>.<version>``.
    The new file is written beside it and renamed over it, so a job running the
    old script -- bash reads a script as it goes -- keeps reading the old one.

    Returns
    -------
    str or None
        What was done, or None if nothing was.
    """
    dest = Path(dest)
    backup = None
    if dest.exists():
        file_type, version = _get_script_version(dest.read_text())
        new_type, new_version = _get_script_version(text)
        old, new = _version_tuple(version), _version_tuple(new_version)
        if file_type != new_type or old is None or new is None or new <= old:
            return None
        backup = dest.with_name(f"{dest.name}.{version}")
        shutil.copy2(dest, backup)
    tmp = dest.with_name(f".{dest.name}.{os.getpid()}.tmp")
    tmp.write_text(text)
    os.chmod(tmp, 0o755)
    os.replace(tmp, dest)
    if backup is None:
        return f"Installed {dest.name} in {dest.parent}"
    return (
        f"Upgraded {dest.name} in {dest.parent} from version {version} to "
        f"{new_version}; the old one is kept as {backup.name}."
    )


def _gpu_code_is_usable(ff_form, gpu_code):
    """Whether a forcefield can actually use the configured ``gpu-code``.

    A GPU in the allocation is not on its own a reason to use ``gpu-code``.
    ``NGPUS`` comes from the scheduler -- seamm_exec reads ``SLURM_JOB_GPUS``
    -- so it is set for any job that asked for a GPU, whatever the forcefield.

    ``gpu-code`` serves two unrelated purposes. A Kokkos command is how a
    *classical* forcefield uses a GPU. An MDI command launches a
    machine-learned engine alongside LAMMPS, and only a PyTorch forcefield has
    a model for it to serve: run a classical one that way and the engine is
    handed ``SEAMM_FF="Unknown"`` and dies with FileNotFoundError before
    LAMMPS gets anything.

    So the question is not "are there GPUs" but "does this command suit this
    forcefield", which is decided by whether it drives MDI.
    """
    return ff_form == "PyTorch" or "-mdi" not in gpu_code


def _required_mpi_slots(command):
    """Total MPI ranks an mpirun command line asks for, or None if unknown.

    An MDI run is MPMD -- ``mpirun -np 1 <engine> : -np 1 <driver>`` -- so what
    it needs is the *sum* of the -np values across the colon-separated
    segments, not the largest of them. Returns None when any count is a
    template placeholder such as {NTASKS}, which is only resolved later and is
    derived from the allocation anyway, or when there is no -np at all.
    """
    total = 0
    found = False
    tokens = command.replace(":", " : ").split()
    for i, token in enumerate(tokens[:-1]):
        if token in ("-np", "-n", "--np", "--n"):
            try:
                total += int(tokens[i + 1])
            except ValueError:
                return None
            found = True
    return total if found else None


def _allocated_mpi_slots():
    """MPI tasks the batch system allocated, or None if not under one."""
    for variable in ("SLURM_NTASKS", "SLURM_NPROCS"):
        value = os.environ.get(variable)
        if value:
            try:
                return int(value)
            except ValueError:
                pass
    return None


def _available_gpus(maxload, timeout=10):
    """Physical indices of the NVIDIA GPUs whose load and memory use are at
    or below ``maxload`` (a fraction, 0-1), as ``nvidia-smi`` reports them.

    Replaces GPUtil, which is unmaintained and imports ``distutils``, so it no
    longer imports on Python 3.12 outside a conda environment (conda ships
    setuptools, which provides a shim; a plain venv does not). Returns an
    empty list when ``nvidia-smi`` is missing, fails, or times out, i.e. when
    there are no usable GPUs. A GPU whose numbers cannot be parsed is skipped,
    as GPUtil did with ``includeNan=False``.
    """
    smi = shutil.which("nvidia-smi")
    if smi is None:
        return []
    try:
        result = subprocess.run(
            [
                smi,
                "--query-gpu=index,utilization.gpu,memory.used,memory.total",
                "--format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except (OSError, subprocess.SubprocessError):
        return []
    if result.returncode != 0:
        return []
    return _parse_nvidia_smi(result.stdout, maxload)


def _parse_nvidia_smi(text, maxload):
    """Parse ``nvidia-smi --query-gpu=index,utilization.gpu,memory.used,
    memory.total --format=csv,noheader,nounits`` output into the indices of
    the GPUs at or below ``maxload`` in both utilization and memory use."""
    gpus = []
    for line in text.splitlines():
        fields = [f.strip() for f in line.split(",")]
        if len(fields) != 4:
            continue
        try:
            index = int(fields[0])
            load = float(fields[1]) / 100
            used = float(fields[2])
            total = float(fields[3])
        except ValueError:
            continue
        if total <= 0:
            continue
        memory = used / total
        if load <= maxload and memory <= maxload:
            gpus.append(index)
    return gpus


def _cuda_visible_devices():
    """Parse CUDA_VISIBLE_DEVICES into the physical GPU indices it exposes.

    Returns None when the variable is unset, meaning every GPU on the machine
    is visible and a CUDA program's device *i* is physical GPU *i*. Otherwise
    returns one entry per visible GPU, in order: the physical index, or None
    for an entry given as a UUID rather than an index.

    The ordering is the point. CUDA_VISIBLE_DEVICES renumbers devices, so with
    CUDA_VISIBLE_DEVICES=1 the process sees exactly one GPU and it is device
    *0*, not device 1. A scheduler sets this when it allocates a subset of a
    node's GPUs.
    """
    raw = os.environ.get("CUDA_VISIBLE_DEVICES")
    if raw is None:
        return None

    devices = []
    for entry in (e.strip() for e in raw.split(",")):
        if entry == "":
            continue
        try:
            devices.append(int(entry))
        except ValueError:
            devices.append(None)
    return devices


def _usable_gpus(available, visible):
    """GPUs this process may use, in the numbering CUDA itself uses.

    ``available`` is what ``nvidia-smi`` reports (via ``_available_gpus``),
    which are *physical* indices -- nvidia-smi does not know about
    CUDA_VISIBLE_DEVICES. ``visible`` is the parsed variable.

    Without this, a job allocated GPU 1 by the scheduler would be told to use
    GPU 0, because that is the physical index nvidia-smi reports as idle: it would
    run on a GPU it was never given, alongside whoever actually holds it.
    """
    if visible is None:
        return list(available)

    if any(device is None for device in visible):
        # UUID form: these cannot be lined up with nvidia-smi's indices, so treat
        # the whole allocation as usable rather than silently discarding it.
        return list(range(len(visible)))

    return [
        logical for logical, physical in enumerate(visible) if physical in available
    ]


#: A leading NAME=value item of an engine's argv
_ASSIGNMENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=")


def _engine_cpu_threads(ncores="available"):
    """The CPU threads for an MLFF engine running without a GPU.

    The CPUs this process may use (the job's, under a scheduler) less one for the
    LAMMPS driver, and no more than ``ncores`` unless that is ``"available"``.
    """
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError:  # macOS
        n = os.cpu_count() or 1
    n = max(1, n - 1)
    if str(ncores) != "available":
        n = max(1, min(n, int(ncores)))
    return n


#: Why an MLFF model chemistry runs without a GPU
_MLFF_ON_CPU = (
    " No GPU is available, so the engine of the machine-learned forcefield runs "
    "on the CPU, which is much slower except for small systems."
)


def _wait_for_gpus(maxload, wait=None, sleep=time.sleep, clock=time.time):
    """The usable GPUs that are idle, waiting up to ``wait`` seconds for one.

    Parameters
    ----------
    maxload : float
        The load and memory use (a fraction) at or below which a GPU is idle.
    wait : float or None
        The seconds to wait for an idle GPU; ``None`` waits as long as it takes.
    sleep, clock : callable
        ``time.sleep`` and ``time.time``, replaceable for testing.

    Returns
    -------
    [int]
        The idle GPUs, in CUDA's numbering; empty to run on the CPU. That is at
        once when this process can never have a GPU -- there is none on the
        machine, or CUDA_VISIBLE_DEVICES gives it none -- and otherwise when no
        GPU became idle within ``wait``.
    """
    visible = _cuda_visible_devices()
    t_end = None if wait is None else clock() + wait
    while True:
        gpus = _usable_gpus(_available_gpus(maxload), visible)
        if len(gpus) > 0:
            return gpus
        if len(_usable_gpus(_available_gpus(1.0), visible)) == 0:
            return []
        if t_end is not None and clock() >= t_end:
            return []
        sleep(1)


#: Where the runs on a machine claim the GPUs they use, so that runs started
#: together -- in one allocation, or by the task layer -- spread over the idle GPUs
#: rather than all taking the first one
GPU_CLAIM_DIR = Path(tempfile.gettempdir()) / "seamm-gpu-claims"


def _gpu_claim_name(logical):
    """The name of the claim on a GPU given in CUDA's numbering: the physical GPU,
    which is the same for every process on the machine."""
    raw = [e.strip() for e in os.environ.get("CUDA_VISIBLE_DEVICES", "").split(",")]
    raw = [e for e in raw if e != ""]
    if "CUDA_VISIBLE_DEVICES" in os.environ and logical < len(raw):
        physical = raw[logical]
    else:
        physical = str(logical)
    physical = re.sub(r"[^A-Za-z0-9._-]", "_", physical)
    return f"{socket.gethostname()}-gpu{physical}.lock"


def _claim_gpus(gpus, n, directory=None):
    """Choose ``n`` of the idle GPUs, claiming each so other runs choose others.

    A claim is a lock (``flock``) on a file per physical GPU, held until
    :func:`_release_gpus`, i.e. for the run; the system drops it if the process
    dies. Unclaimed GPUs are taken first, in order. If there are not ``n`` of
    them, GPUs claimed by other runs make up the number, as before claims existed.

    Parameters
    ----------
    gpus : [int]
        The idle GPUs, in CUDA's numbering, from :func:`_wait_for_gpus`.
    n : int
        How many to use.
    directory : pathlib.Path, optional
        Where the claims are kept, by default :data:`GPU_CLAIM_DIR`.

    Returns
    -------
    ([int], [int])
        The GPUs chosen, and the open files holding the claims, for
        :func:`_release_gpus`.
    """
    try:
        import fcntl
    except ImportError:  # Windows: no claims
        return list(gpus[:n]), []
    directory = GPU_CLAIM_DIR if directory is None else Path(directory)
    try:
        directory.mkdir(parents=True, exist_ok=True)
        directory.chmod(0o1777)  # shared by everyone on the machine
    except OSError:
        pass
    chosen, claims, taken = [], [], []
    for gpu in gpus:
        if len(chosen) >= n:
            break
        try:
            fd = os.open(
                directory / _gpu_claim_name(gpu), os.O_RDWR | os.O_CREAT, 0o666
            )
        except OSError:
            taken.append(gpu)
            continue
        try:
            os.fchmod(fd, 0o666)
        except OSError:
            pass
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            os.close(fd)
            taken.append(gpu)
            continue
        chosen.append(gpu)
        claims.append(fd)
    chosen.extend(taken[: n - len(chosen)])
    return chosen, claims


def _release_gpus(claims):
    """Release the claims :func:`_claim_gpus` made."""
    for fd in claims:
        try:
            os.close(fd)
        except OSError:
            pass


def _launch_environment(gpus):
    """The environment variables for launching LAMMPS.

    LAMMPS runs in parallel with MPI ranks, one per core, so each rank must be
    single-threaded. The environment may say otherwise: an iteration of a parallel
    Loop sets OMP_NUM_THREADS to all of its cores for the codes it runs, and every
    rank then started that many threads (56 ranks x 14 threads on 56 cores, with no
    progress). It applies to the whole launch script, so a QM engine over MDI (MOPAC,
    a small serial calculation) is single-threaded too; an MLFF engine on the CPU is
    given its own number of threads on its command line.

    Parameters
    ----------
    gpus : [int]
        The GPUs of the run, if any.

    Returns
    -------
    {str: str}
    """
    env = {"OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if len(gpus) != 0:
        env.update(_gpu_environment(gpus))
    return env


def _gpu_environment(gpus):
    """The environment variables for a run on the given GPUs.

    One OpenMP thread is usually best when the GPU does the work. OMP_PROC_BIND and
    OMP_PLACES are deliberately not set: with them the OpenMP runtime binds each
    process's first thread to the first CPU it may use as torch or LAMMPS starts, and
    the programs it starts inherit that, so every engine and driver in a job shared
    one core. mdi_bind.sh and gpu_bind.sh bind the CPUs themselves where they can.
    """
    return {
        "SEAMM_GPUS": ",".join(str(i) for i in gpus),
        "OMP_NUM_THREADS": "1",
    }


def _free_tcp_port(hostname="localhost"):
    """Ask the OS for a free TCP port for the MDI rendezvous.

    Bind a throwaway socket to port 0, read back the port the OS assigned,
    release it, and return that integer for both the engine and the LAMMPS
    driver to use. There is an unavoidable race between releasing the socket
    and the engine re-binding it (classic TOCTOU); acceptable for the typical
    one-job-per-node case. See campaigns/2026-06-22/NOTES_C.rst (D5).
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        s.bind((hostname, 0))
        return s.getsockname()[1]
    finally:
        s.close()


#: The launch of an MLFF engine and the LAMMPS driver as an MPI pair. {bind} is
#: mdi_bind.sh and {engine} the engine's command; "-in input.dat" is appended.
#: lammps.ini's "mlff-code" replaces it, e.g. to give mpirun other options.
MLFF_PAIR = (
    "mpirun --bind-to none -np 1 {bind} {engine} "
    ': -np 1 {bind} lmp -mdi "-role DRIVER -name LAMMPS -method MPI"'
)


def _mlff_pair_script(engine_argv, config, bind):
    """The bash script launching an MLFF engine and LAMMPS as an MPI pair.

    The template is ``mlff-code`` from lammps.ini, else :data:`MLFF_PAIR`; its
    ``{bind}`` and ``{engine}`` are filled in here, so the engine's options may
    contain braces (``--dispersion "{regime: large}"``), which the executor's
    substitution in a command line would take for placeholders. An ini value on
    several lines is joined into one.
    """
    template = config.get("mlff-code", "").strip() or MLFF_PAIR
    template = re.sub(r"[ \t]*\\?[ \t]*\n[ \t]*", " ", template)
    line = template.replace("{bind}", shlex.quote(str(bind)))
    line = line.replace("{engine}", shlex.join(engine_argv))
    return "#!/bin/bash\n" + line + " -in input.dat\n"


#: The launch of an MLFF engine in another environment than LAMMPS: the two meet
#: over TCP on localhost. The engine's command carries its own ``conda run``. Each
#: side runs through mdi_bind.sh ({bind}), which stops the other when one fails --
#: otherwise LAMMPS would wait forever for an engine that died starting. The
#: engine dials until LAMMPS listens, so it may start first. lammps.ini's
#: "mlff-tcp-code" replaces the LAMMPS line, e.g. to name another executable.
MLFF_TCP_DRIVER = 'lmp -mdi "-role DRIVER -name LAMMPS -method TCP -port {port}"'


#: The shim that sets TCP_NODELAY on MDI's TCP connections (data/mdi_nodelay.c). It
#: is written into the job's directory, and the launch script compiles and preloads it.
NODELAY_SOURCE = "mdi_nodelay.c"


def _tcp_nodelay(config):
    """Whether to set TCP_NODELAY on MDI's TCP connections: lammps.ini's
    "tcp-nodelay", yes by default."""
    value = str(config.get("tcp-nodelay", "yes")).strip().lower()
    return value not in ("no", "false", "off", "0")


def _nodelay_source():
    """The text of the TCP_NODELAY shim."""
    path = importlib.resources.files("lammps_step") / "data" / NODELAY_SOURCE
    return path.read_text()


def _nodelay_preamble():
    """Bash that compiles mdi_nodelay.c and preloads it into what follows.

    MDI sets no TCP_NODELAY, so each message can wait up to ~40 ms for the other
    side's delayed acknowledgement. The shim is compiled where the job runs, so it
    matches that machine; it needs a C compiler ($CC, else cc) and Linux, where
    LD_PRELOAD works. Without them, the run goes on without it, with a warning.
    """
    return (
        "# MDI sets no TCP_NODELAY, so a message can wait ~40 ms for the other side\n"
        f"# to acknowledge the last. On Linux, compile and preload {NODELAY_SOURCE},\n"
        "# which sets it.\n"
        f'if [ "$(uname -s)" = Linux ] && [ -f {NODELAY_SOURCE} ]; then\n'
        "    NODELAY_DIR=$(mktemp -d 2>/dev/null || true)\n"
        '    if [ -n "$NODELAY_DIR" ] && ${CC:-cc} -shared -fPIC -O2 \\\n'
        f'            -o "$NODELAY_DIR/mdi_nodelay.so" {NODELAY_SOURCE} -ldl \\\n'
        "            2> mdi_nodelay.err; then\n"
        "        trap 'rm -rf \"$NODELAY_DIR\"' EXIT\n"
        '        export LD_PRELOAD="${LD_PRELOAD:+$LD_PRELOAD:}'
        '$NODELAY_DIR/mdi_nodelay.so"\n'
        "        rm -f mdi_nodelay.err\n"
        "    else\n"
        '        [ -n "$NODELAY_DIR" ] && rm -rf "$NODELAY_DIR"\n'
        f'        echo "Warning: could not compile {NODELAY_SOURCE}"\\\n'
        '            "(see mdi_nodelay.err), so MDI runs over TCP without"\\\n'
        '            "TCP_NODELAY and may be slower." >&2\n'
        "    if [ -s mdi_nodelay.err ]; then sed 's/^/    /' mdi_nodelay.err >&2; fi\n"
        "    fi\n"
        "fi\n"
    )


def _mlff_tcp_script(engine_argv, port, config, bind):
    """The bash script launching an MLFF engine and LAMMPS as a pair over TCP.

    Leading ``NAME=value`` items of the engine's argv (thread counts) are passed
    through ``env``, since mdi_bind.sh runs the command with ``exec``. The script
    keeps the environment, such as LD_PRELOAD, adding only the TCP_NODELAY shim
    (unless lammps.ini's "tcp-nodelay" is "no"; see :func:`_nodelay_preamble`), and
    exits with LAMMPS's status, else the engine's.

    When LAMMPS fails, mdi_bind.sh normally stops the engine. If the driver never
    got as far as mdi_bind.sh -- a wrapper in "mlff-tcp-code" that cannot start,
    for instance -- the engine would wait for a connection for ever, so the script
    stops the engine and everything under it once its grace period has passed.
    """
    assignments = []
    argv = list(engine_argv)
    while argv and re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", argv[0]):
        assignments.append(argv.pop(0))
    engine = shlex.join(["env", *assignments, *argv] if assignments else argv)

    driver = config.get("mlff-tcp-code", "").strip() or MLFF_TCP_DRIVER
    driver = re.sub(r"[ \t]*\\?[ \t]*\n[ \t]*", " ", driver)
    driver = driver.replace("{port}", str(port))
    bind = shlex.quote(str(bind))
    return (
        "#!/bin/bash\n"
        "# The engine is in another environment than LAMMPS, so they meet over TCP.\n"
        f"export SEAMM_MDI_PAIR=tcp{port}\n"
        + (_nodelay_preamble() if _tcp_nodelay(config) else "")
        + f"SEAMM_MDI_ROLE=engine {bind} {engine} &\n"
        "ENGINE_PID=$!\n"
        f"SEAMM_MDI_ROLE=driver {bind} {driver} -in input.dat\n"
        "rc=$?\n"
        "stop_tree() {\n"
        "    local child\n"
        '    for child in $(pgrep -P "$1" 2>/dev/null); do stop_tree "$child"; done\n'
        '    kill -TERM "$1" 2>/dev/null\n'
        "}\n"
        'if [ "$rc" -ne 0 ]; then\n'
        "    # An engine still running after the grace is waiting for a driver that\n"
        "    # never came: stop it rather than waiting for ever.\n"
        "    for (( i=0; i<${SEAMM_ENGINE_GRACE:-120}; i++ )); do\n"
        '        kill -0 "$ENGINE_PID" 2>/dev/null || break\n'
        "        sleep 1\n"
        "    done\n"
        '    kill -0 "$ENGINE_PID" 2>/dev/null && stop_tree "$ENGINE_PID"\n'
        "fi\n"
        "wait $ENGINE_PID\n"
        "engine_rc=$?\n"
        '[ "$rc" -eq 0 ] && rc=$engine_rc\n'
        "exit $rc\n"
    )


def _summarize_run_result(result, max_text=2000):
    """The executor's result for the log, without the contents of the files.

    The result holds every returned file -- trajectories of gigabytes among them --
    and formatting it whole took about nine times the size of the files, enough to
    exhaust a machine running several MD runs at once. Files are given by their
    size, and other long text is shortened.
    """
    summary = {}
    for key, value in result.items():
        if isinstance(value, dict) and ("data" in value or "path" in value):
            if "path" in value:
                path = Path(value["path"])
                entry = {
                    "path": str(path),
                    "size": path.stat().st_size if path.exists() else None,
                }
            else:  # seamm-exec before 2026.10.18 read the files
                data = value["data"]
                entry = {"size": None if data is None else len(data)}
            if value.get("exception") is not None:
                entry["exception"] = value["exception"]
            summary[key] = entry
        elif isinstance(value, (str, bytes)) and len(value) > max_text:
            summary[key] = value[:max_text] + f"... ({len(value)} characters)"
        else:
            summary[key] = value
    return summary


def _conda_environment(config):
    """The conda environment an ini section names, as a path if it gives one."""
    environment = config.get("conda-environment", "").strip()
    if environment.startswith("~") or "/" in environment:
        return str(Path(environment).expanduser())
    return environment


def _same_environment(config, engine_config):
    """Whether LAMMPS and an MDI engine run in the same environment.

    Each is described by its ini section (lammps.ini, and e.g. xnn.ini). Two
    "conda" installations are the same when they name the same environment; an
    environment given by name matches one given by a path ending in that name.
    Two "local" installations share the PATH, so they count as the same.
    """
    installation = config.get("installation", "conda")
    if installation != engine_config.get("installation", "conda"):
        return False
    if installation != "conda":
        return True
    a = _conda_environment(config)
    b = _conda_environment(engine_config)
    if a == b:
        return True
    if "/" in a and "/" not in b:
        return Path(a).name == b
    if "/" in b and "/" not in a:
        return Path(b).name == a
    return False


def _mdi_launch_script(engine_argv, port, config, ce, hostname="localhost"):
    """Compose the bash launch script for QM-MD over MDI.

    The QM engine is the listener (it binds the port and blocks in
    MDI_Accept_Communicator), so it is backgrounded first; the LAMMPS driver
    runs second and connects; then we wait on the engine. ``executor.run``
    wraps this whole script in the seamm-lammps conda env, while the engine
    line carries its own ``conda run -n seamm-mopac`` so each code stays in its
    own environment. See campaigns/2026-06-22/NOTES_C.rst.

    Parameters
    ----------
    engine_argv : list of str
        The fully resolved engine launch argv (from the program step's
        get_mdi_engine_command); rendered with shlex.join.
    port : int
        The TCP port both sides use.
    config : dict
        The lammps.ini section. The driver command comes from its "code" /
        "cmd-args" keys (D4).
    ce : dict
        The computational environment, for resolving ini templates such as
        {NTASKS}.
    hostname : str
        The rendezvous host (default "localhost").

    Returns
    -------
    str
        The launch-script text.
    """
    # Driver (LAMMPS) command: resolve the ini templates the same way
    # executor.run does (nested .format), then add the MDI DRIVER flag.
    driver = config["code"]
    while True:
        resolved = driver.format(**config, **ce)
        if resolved == driver:
            break
        driver = resolved
    if config.get("cmd-args", "") != "":
        driver += " " + config["cmd-args"]
    driver += (
        f' -mdi "-role DRIVER -name LAMMPS -method TCP '
        f'-port {port} -hostname {hostname}" -in input.dat'
    )

    return (
        "#!/bin/bash\n"
        "set -e\n"
        + (_nodelay_preamble() if _tcp_nodelay(config) else "")
        + f"{shlex.join(engine_argv)} &\n"
        "ENGINE_PID=$!\n"
        f"{driver}\n"
        "wait $ENGINE_PID\n"
    )


# Temporarily used here to stop pymbar's annoying warning.
@contextmanager
def logging_disabled(highest_level=logging.CRITICAL):
    """
    A context manager that will prevent any logging messages
    triggered during the body from being processed.

    :param highest_level: the maximum logging level in use.
      This would only need to be changed if a custom level greater than CRITICAL
      is defined.

    From Simon Weber https://gist.github.com/simon-weber/7853144
    """
    # two kind-of hacks here:
    #    * can't get the highest logging level in effect => delegate to the user
    #    * can't get the current module-level override => use an undocumented
    #       (but non-private!) interface

    previous_level = logging.root.manager.disable

    logging.disable(highest_level)

    try:
        yield
    finally:
        logging.disable(previous_level)


with logging_disabled(highest_level=logging.WARNING):
    from pymbar import timeseries


# Add LAMMPS's properties to the standard properties
path = importlib.resources.files("lammps_step") / "data"
csv_file = path / "properties.csv"
molsystem.add_properties_from_file(csv_file)

# Rigid bonds and angles are held by SHAKE during dynamics. LAMMPS cannot minimize
# with SHAKE (its "kbond" restraints make the line search fail), so they are also
# stiff harmonic terms, which hold the geometry while minimizing and are exactly
# zero when SHAKE holds the bond or angle at its value. In kcal/mol/Å^2 and
# kcal/mol/radian^2.
rigid_bond_k = 1000.0
rigid_angle_k = 100.0

bond_style = {
    "rigid_bond": "harmonic",
    "quadratic_bond": "harmonic",
    "quartic_bond": "class2",
    "fene": "fene",
    "morse": "morse",
}

angle_style = {
    "rigid_angle": "harmonic",
    "quadratic_angle": "harmonic",
    "quartic_angle": "class2",
    "cosine": "cosine",
    "cosine/squared": "cosine/squared",
    "simple_fourier_angle": "fourier/simple",
    "tabulated_angle": "table",
}

dihedral_style = {
    "torsion_1": "harmonic",
    "torsion_3": "class2",
    "torsion_opls": "opls",
    "torsion_charmm": "charmm",
}

improper_style = {
    "wilson_out_of_plane": "class2",
    "improper_opls": "cvff",
    "dreiding_out_of_plane": "umbrella",
}


def pvalue(values, key, fmt="9.4f", default="n/a"):
    """Return the formatted value and original value of a parameter.

    Parameters
    ----------
    values : [str: str]
        The dictionary of parameter values.

    key : str
        The parameter name

    fmt : str = "9.4f"
        The format for the parameter

    default : str = "n/a"
        The default value for the original parameter, if not find.

    Returns
    -------
    (str, str)
        The formatted parameter value and the original form it came from.
    """
    val = float(values[key])
    original = values["original " + key] if "original " + key in values else default

    return f"{val:{fmt}}", original


def pline(
    i,
    function,
    values,
    keys,
    types,
    rtypes,
    fmt="10.4f",
    default="n/a",
):
    """Create a line for the parameters in structure.dat"""
    line = f"{i:6d} {function}"
    for key in keys:
        value = values[key]
        if isinstance(value, float):
            line += f" {value:{fmt}}"
        elif isinstance(value, int):
            fmt2 = fmt.split(".")[0] + "d"
            line += f" {value:{fmt2}}"
        elif isinstance(value, str):
            if "." in value:
                line += f" {float(value):{fmt}}"
            else:
                fmt2 = fmt.split(".")[0] + "s"
                line += f" {value:{fmt2}}"
        else:
            raise ValueError(
                f"Cannot handle forcefield parameter of type {type(value).__name__}:"
                f" {value}"
            )

    line += " #"

    for typ in types:
        line += f" {typ:6s}"
    line += " -->"
    for typ in rtypes:
        line += f" {typ:6s}"

    if "version" in values:
        line += f" v{values['version']}"

    for key in keys:
        oval = values["original " + key] if "original " + key in values else default
        line += f" {oval:7s}"

    return line


def _output_text(directory, names):
    """The text of the first of ``names`` found in ``directory``, or None."""
    for name in names:
        path = Path(directory) / name
        if path.exists():
            try:
                return path.read_text(errors="replace")
            except OSError:
                return None
    return None


#: What the cost model is made of (seamm_exec.timing_model.Spec as plain data): atoms
#: as the size variable, MD steps as the unit, and the forcefield as the class -- its
#: form (PyTorch, OpenKIM, MDI/QM, or a classical form) and the model itself, since an
#: MLFF and a classical forcefield differ by orders of magnitude per atom-step. Written
#: beside the records when a run is recorded.
TIMING_SPEC = {
    "size": ["n_atoms"],
    "klass": ["ff_form", "model"],
    "task": None,
    "units": "md_steps",
    "multiplier": None,
    "default_alpha": 0.8,
}


def _record_kwargs():
    """``spec=`` for seamm-exec releases that take it (2026.10.6.1 on)."""
    return {"spec": TIMING_SPEC} if hasattr(seamm_exec, "TimingSpec") else {}


def _gpu_model():
    """The model of this node's first GPU (e.g. "NVIDIA A100-SXM4-80GB"), so that the
    timing records tell an A30 from an A100; None if it cannot be found."""
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except Exception:
        return None
    lines = [x.strip() for x in out.stdout.splitlines() if x.strip()]
    return lines[0] if out.returncode == 0 and lines else None


def timing_descriptors(control, log_text=None, configuration=None):
    """The descriptors of a LAMMPS run for its timing record (seamm_exec's
    campaign of 2026-10-05): the structure, the kinds of calculation in the
    run (the sub-steps' names), and from the log the ``Loop time`` lines --
    MD steps, atoms, processes and LAMMPS's own time summed over the runs --
    whose product of steps and atoms is the unit of cost.
    """
    d = {}
    names = []
    for item in control or ():
        try:
            names.append(str(item[0]) if isinstance(item, (list, tuple)) else str(item))
        except Exception:
            pass
    d["calculations"] = " + ".join(names)[:200]
    d["n_calculations"] = len(names)
    if configuration is not None:
        d.update(seamm_exec.structure_descriptors(configuration))
    if log_text:
        loops = re.findall(
            r"Loop time of ([\d.eE+-]+) on (\d+) procs for (\d+) steps"
            r" with (\d+) atoms",
            log_text,
        )
        d["n_runs"] = len(loops)
        d["md_steps"] = sum(int(x[2]) for x in loops) if loops else None
        d["code_seconds"] = sum(float(x[0]) for x in loops) if loops else None
        d["procs"] = int(loops[0][1]) if loops else None
        m = re.search(r"Total wall time:\s+(\d+):(\d+):(\d+)", log_text)
        if m:
            d["total_wall"] = (
                int(m.group(1)) * 3600 + int(m.group(2)) * 60 + int(m.group(3))
            )
        d["terminated_normally"] = "Total wall time" in log_text
    return d


def allocation_note(wanted, ntasks):
    """Why LAMMPS runs on fewer processes than it could: the batch job allocated
    fewer MPI tasks. CPUs per task do not count, since each MPI process is a
    task."""
    return (
        f" However, this job was allocated only {ntasks} MPI "
        f"task{'s' if ntasks != 1 else ''}, so LAMMPS runs on {ntasks} "
        f"process{'es' if ntasks != 1 else ''}. Each MPI process needs a task: ask "
        f"the batch system for {wanted} tasks (e.g. --ntasks={wanted}) rather than "
        "CPUs per task."
    )


class LAMMPS(seamm.Node):
    display_units = {
        "T": "K",
        "P": "atm",
        "t": "fs",
        "V": "Å^3",
        "density": "g/mL",
        "a": "Å",
        "b": "Å",
        "c": "Å",
        "alpha": "º",
        "beta": "º",
        "gamma": "º",
        "Etot": "kcal/mol",
        "DfH0_reax": "kcal/mol",
        "Eke": "kcal/mol",
        "Epe": "kcal/mol",
        "Emol": "kcal/mol",
        "Epair": "kcal/mol",
        "Jx": "W/m^2",
        "Jy": "W/m^2",
        "Jz": "W/m^2",
        "Kappa_x": "W/K/m",
        "Kappa_y": "W/K/m",
        "Kappa_z": "W/K/m",
        "Kappa": "W/K/m",
        "u": "kcal/mol",
        "N_hbond": "",
        "E_hbond": "kcal/mol",
        "Sxx": "atm",
        "Syy": "atm",
        "Szz": "atm",
        "Syz": "atm",
        "Sxz": "atm",
        "Sxy": "atm",
    }
    display_title = {
        "T": "Temperature",
        "P": "Pressure",
        "t": "Time",
        "V": "Volume",
        "density": "Density",
        "a": "a lattice parameter",
        "b": "b lattice parameter",
        "c": "c lattice parameter",
        "Etot": "Total Energy",
        "DfH0_reax": "Reax \N{GREEK CAPITAL LETTER DELTA}fH\N{SUPERSCRIPT ZERO}",
        "Eke": "Kinetic Energy",
        "Epe": "Potential Energy",
        "Emol": "Molecular Energy, Valence Terms",
        "Epair": "Pair (Nonbond) Energy",
        "Jx": "heat flux in x",
        "Jy": "heat flux in y",
        "Jz": "heat flux in z",
        "Kappa_x": "thermal conductivity in x",
        "Kappa_y": "thermal conductivity in y",
        "Kappa_z": "thermal conductivity in z",
        "Kappa": "thermal conductivity",
        "u": "atom PE",
        "N_hbond": "Number of H-bonds",
        "E_hbond": "Energy of H-bonds",
        "Sxx": "Sxx",
        "Syy": "Syy",
        "Szz": "Szz",
        "Syz": "Syz",
        "Sxz": "Sxz",
        "Sxy": "Sxy",
    }

    def __init__(
        self, flowchart=None, namespace="org.molssi.seamm.lammps", extension=None
    ):
        """Setup the main LAMMPS step

        Keyword arguments:
        """
        logger.debug("Creating LAMMPS {}".format(self))

        # The subflowchart
        self.subflowchart = seamm.Flowchart(
            parent=self, name="LAMMPS", namespace=namespace
        )
        self._initialization_node = None
        self._trajectory = []
        self._data = {}
        self._have_dreiding_hbonds = False
        self._atomic_energy_sum = 0.0
        self._eex = None
        self._force_triclinic = False  # Whether to insist on a triclinic cell
        self._triclinic_box = False  # Whether the box written is triclinic

        self._results = {}  # Storage for computational and timing results

        super().__init__(
            flowchart=flowchart, title="LAMMPS", extension=extension, logger=logger
        )

    @property
    def version(self):
        """The semantic version of this module."""
        return lammps_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return lammps_step.__git_revision__

    @property
    def eex(self):
        """The energy expression."""
        return self._eex

    @eex.setter
    def eex(self, value):
        self._eex = value

    @property
    def have_dreiding_hbonds(self):
        """Whether the system has Dreiding hydrogen bonds."""
        return self._have_dreiding_hbonds

    @have_dreiding_hbonds.setter
    def have_dreiding_hbonds(self, value):
        self._have_dreiding_hbonds = value

    @property
    def results(self):
        """The storage for results."""
        return self._results

    @property
    def force_triclinic(self):
        """Whether to force a triclinic cell."""
        return self._force_triclinic

    @force_triclinic.setter
    def force_triclinic(self, value):
        self._force_triclinic = value

    @staticmethod
    def box_to_cell(lx, ly, lz, xy, xz, yz):
        """Convert the LAMMPS box definition to cell parameters."""
        if xy == 0 and xz == 0 and yz == 0:
            a = lx
            b = ly
            c = lz
            alpha = 90.0
            beta = 90.0
            gamma = 90.0
        else:
            a = lx
            b = sqrt(ly**2 + xy**2)
            c = sqrt(lz**2 + xz**2 + yz**2)
            alpha = degrees(acos((xy * xz + lx * yz) / (b * c)))
            beta = degrees(acos(xz / c))
            gamma = degrees(acos(xy / b))

        return (a, b, c, alpha, beta, gamma)

    @staticmethod
    def cell_to_box(a, b, c, alpha, beta, gamma):
        """Convert cell parameters to the LAMMPS box."""
        if alpha == 90 and beta == 90 and gamma == 90:
            lx = a
            ly = b
            lz = c
            xy = xz = yz = 0.0
        else:
            lx = a
            xy = b * cos(radians(gamma))
            xz = c * cos(radians(beta))
            ly = sqrt(b**2 - xy**2)
            yz = (b * c * cos(radians(alpha)) - xy * xz) / ly
            lz = sqrt(c**2 - xz**2 - yz**2)

        return (lx, ly, lz, xy, xz, yz)

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = self.step_type
        parser = seamm_util.getParser()

        # Remember if the parser exists ... this type of step may have been
        # found before
        parser_exists = parser.exists(parser_name)

        # Create the standard options, e.g. log-level
        result = super().create_parser(name=parser_name)

        if parser_exists:
            return result

        # LAMMPS specific options
        parser.add_argument(
            parser_name,
            "--ncores",
            default="available",
            help=(
                "The maximum number of cores to use for LAMMPS. "
                "Default: all available cores."
            ),
        )
        parser.add_argument(
            parser_name,
            "--ngpus",
            default="available",
            help=(
                "The maximum number of gpus to use for LAMMPS. "
                "Default: all available gpus."
            ),
        )
        parser.add_argument(
            parser_name,
            "--gpu-load",
            default="10",
            help=(
                "The percent maximum cpu or memory usage on GPU considered idle. "
                "Default: 10"
            ),
        )
        parser.add_argument(
            parser_name,
            "--gpu-wait-time",
            default="unlimited",
            help=(
                "How long to wait for a GPU to become available, in seconds. "
                "Default: unlimited"
            ),
        )
        parser.add_argument(
            parser_name,
            "--atoms-per-core",
            type=int,
            default="100",
            help="the optimal number of atoms per core for LAMMPS",
        )
        parser.add_argument(
            parser_name,
            "--graph-formats",
            default=tuple(),
            choices=("html", "png", "jpeg", "webp", "svg", "pdf"),
            nargs="+",
            help="extra formats to write for graphs",
        )
        parser.add_argument(
            parser_name,
            "--graph-fontsize",
            default=15,
            help="Font size in graphs, defaults to 15 pixels",
        )
        parser.add_argument(
            parser_name,
            "--graph-width",
            default=1024,
            help="Width of graphs in formats that support it, defaults to 1024",
        )
        parser.add_argument(
            parser_name,
            "--graph-height",
            default=1024,
            help="Height of graphs in formats that support it, defaults to 1024",
        )
        if False:
            parser.add_argument(
                parser_name,
                "--modules",
                nargs="*",
                default=None,
                help="the environment modules to load for LAMMPS",
            )
            parser.add_argument(
                parser_name,
                "--gpu-modules",
                nargs="*",
                default=None,
                help="the environment modules to load for the GPU version of LAMMPS",
            )
            parser.add_argument(
                parser_name,
                "--lammps-path",
                default=None,
                help="the path to the LAMMPS executables",
            )
            parser.add_argument(
                parser_name,
                "--lammps-serial",
                default="lmp_serial",
                help="the serial version of LAMMPS",
            )
            parser.add_argument(
                parser_name,
                "--lammps-mpi",
                default="lmp_mpi",
                help="the mpi version of LAMMPS",
            )
            parser.add_argument(
                parser_name,
                "--cmd-args",
                default="",
                help="the command-line arguments for LAMMPS, e.g. '-k on'",
            )
            parser.add_argument(
                parser_name,
                "--gpu-cmd-args",
                default="",
                help=(
                    "the command-line arguments for GPU version of LAMMPS, e.g. '-k on'"
                ),
            )
            parser.add_argument(
                parser_name, "--mpiexec", default="mpiexec", help="the mpi executable"
            )

        return result

    def set_id(self, node_id):
        """Set the id for node to a given tuple"""
        self._id = node_id

        # and set our subnodes
        self.subflowchart.set_ids(self._id)

        return self.next()

    def description_text(self, P=None):
        """Return a short description of this step.

        Return a nicely formatted string describing what this step will
        do.

        Keyword arguments:
            P: a dictionary of parameter values, which may be variables
                or final values. If None, then the parameters values will
                be used as is.
        """

        self.subflowchart.root_directory = self.flowchart.root_directory

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        text = self.header + "\n\n"

        while node is not None:
            try:
                text += __(node.description_text(), indent=4 * " ").__str__()
            except Exception as e:
                print(
                    "Error describing LAMMPS flowchart: {} in {}".format(
                        str(e), str(node)
                    )
                )
                self.logger.critical(
                    "Error describing LAMMPS flowchart: {} in {}".format(
                        str(e), str(node)
                    )
                )
                raise
            except:  # noqa: E722
                print(
                    "Unexpected error describing LAMMPS flowchart: {} in {}".format(
                        sys.exc_info()[0], str(node)
                    )
                )
                self.logger.critical(
                    "Unexpected error describing LAMMPS flowchart: {} in {}".format(
                        sys.exc_info()[0], str(node)
                    )
                )
                raise
            text += "\n"
            node = node.next()

        return text

    #: Model-chemistry types that are forcefields LAMMPS evaluates itself
    FORCEFIELD_TYPES = ("FF", "ReaxFF")

    def forcefield_model_chemistry(self):
        """The model chemistry, if it is a forcefield (``FF`` or ``ReaxFF``),
        else None."""
        if not self.variable_exists("_model_chemistry"):
            return None
        mc = self.get_variable("_model_chemistry")
        if mc.get("type") in self.FORCEFIELD_TYPES:
            return mc
        return None

    def forcefield(self):
        """The forcefield for this run.

        From a forcefield model chemistry (``FF@oplsaa``, ``ReaxFF@...``), built
        by the step that offers it (forcefield_step) and kept for the next run with
        the same model chemistry; otherwise the Forcefield step's ``_forcefield``:
        a seamm_ff_util forcefield, or "OpenKIM" or "PyTorch".
        """
        mc = self.forcefield_model_chemistry()
        if mc is None:
            return self.get_variable("_forcefield")
        key = (mc["level"], (mc.get("options") or {}).get("path"))
        cached = getattr(self, "_mc_forcefield", None)
        if cached is None or cached[0] != key:
            provider = self.flowchart.plugin_manager.get(mc["step"])
            cached = (key, provider.forcefield(mc, node=self))
            self._mc_forcefield = cached
        return cached[1]

    def mlff_model_chemistry(self):
        """The model chemistry, if it is a machine-learned forcefield (MLFF),
        else None."""
        if not self.variable_exists("_model_chemistry"):
            return None
        mc = self.get_variable("_model_chemistry")
        return mc if mc.get("type") == "MLFF" else None

    def pytorch_model(self):
        """The machine-learned forcefield's checkpoint: an MLFF model
        chemistry's, else the Forcefield step's ``_pytorch_model``."""
        mc = self.mlff_model_chemistry()
        if mc is not None:
            return (mc.get("options") or {}).get("path")
        return self.get_variable("_pytorch_model")

    def ff_form(self):
        """Return the form of the forcefield."""
        # A forcefield model chemistry is evaluated by LAMMPS itself, and an MLFF
        # by its engine run with LAMMPS as an MPI pair, as the Forcefield step's
        # PyTorch models are; any other model chemistry is QM driven over MDI. A
        # QM-MD flowchart has no Forcefield step, so check this *before*
        # "_forcefield", which would not exist. See campaigns/2026-06-22/NOTES_C.rst
        # (D1) and campaigns/2026-10-11.
        if self.variable_exists("_model_chemistry"):
            if self.mlff_model_chemistry() is not None:
                return "PyTorch"
            if self.forcefield_model_chemistry() is None:
                return "MDI/QM"
        ff = self.forcefield()
        if ff == "OpenKIM":
            ff_form = "OpenKIM"
        elif ff == "PyTorch":
            ff_form = "PyTorch"
        else:
            try:
                ff_form = ff.ff_form
            except Exception:
                ff_form = "unknown"
        return ff_form

    def _mdi_engine_launch(self, configuration):
        """Resolve the "_model_chemistry" selection into an MDI engine launch.

        Reads the "_model_chemistry" workspace variable published by the Model
        Chemistry step (Phase B), validates that the selection can be driven via
        MDI, resolves the owning program step, and asks it for the engine launch
        argv (Phase A's get_mdi_engine_command). LAMMPS is the driver and owns
        the rendezvous, so the TCP port is chosen here and passed to the engine;
        the driver side must use the same port and hostname.

        Parameters
        ----------
        configuration : molsystem._Configuration
            The current configuration, for periodicity and the charge /
            multiplicity (the SEAMM convention -- never from the model
            chemistry step).

        Returns
        -------
        (list of str, int)
            The engine launch argv and the TCP port both sides must use.
        """
        mc = self.get_variable("_model_chemistry")
        options = mc["options"]

        periodic = configuration.periodicity != 0
        if not options.get("mdi_capable", False):
            raise ValueError(
                f"The model chemistry '{mc['level']}' cannot be driven "
                "via MDI; choose an MDI-capable model chemistry."
            )
        if periodic and not options.get("periodic_mdi", False):
            raise ValueError(
                f"The model chemistry '{mc['level']}' is not validated "
                "for periodic systems via MDI."
            )

        port = _free_tcp_port()
        step = self.flowchart.plugin_manager.get(mc["step"])
        # The engine's real keyword (an ORCA functional's alias undone) and the
        # user's basis, for programs that take one; MOPAC, xTB and MLFFs take
        # a method alone.
        kwargs = {}
        if options.get("mdi_basis_arg") is not None:
            kwargs["basis"] = options["mdi_basis_arg"]
        engine_argv = step.get_mdi_engine_command(
            self.flowchart.executor,
            self.global_options,
            method=options.get("mdi_method_arg") or mc["method"],
            port=port,
            hostname="localhost",
            charge=configuration.charge,
            multiplicity=configuration.spin_multiplicity,
            n_atoms=configuration.n_atoms,
            **kwargs,
        )
        return engine_argv, port

    def _mlff_transport(self, config, step):
        """How an MLFF engine and LAMMPS communicate: "MPI" or "TCP".

        As an MPI pair the engine runs in LAMMPS's environment, so the two must
        share it -- and with it MPI and the MDI library -- as in a combined
        environment holding both. When the engine's step is configured for
        another environment (a default installation: seamm-lammps and, e.g.,
        seamm-xnn) the two meet over TCP instead. lammps.ini's ``mlff-transport``
        ("auto", "mpi" or "tcp") overrides the choice.
        """
        choice = config.get("mlff-transport", "auto").strip().lower() or "auto"
        if choice in ("mpi", "tcp"):
            return choice.upper()
        if choice != "auto":
            raise ValueError(
                f"mlff-transport in lammps.ini is '{choice}', but must be 'auto', "
                "'mpi' or 'tcp'."
            )
        if not hasattr(step, "get_executor_config"):
            return "MPI"
        engine_config = step.get_executor_config(
            self.flowchart.executor, self.global_options
        )
        return "MPI" if _same_environment(config, engine_config) else "TCP"

    def _mlff_launch(self, config, bin_dir, gpus):
        """The launch script for an MLFF model chemistry: its MDI engine and the
        LAMMPS driver, each through ``mdi_bind.sh``.

        They run as an MPI pair when the engine shares LAMMPS's environment, and
        over TCP when it does not (see :meth:`_mlff_transport`). The engine's argv
        comes from the step offering the model chemistry
        (``get_mdi_engine_command``): the checkpoint, the device -- ``cuda:0``,
        the GPU that mdi_bind.sh makes visible from ``SEAMM_GPUS``, if one was
        chosen -- and the run's engine options.
        """
        mc = self.mlff_model_chemistry()
        options = mc.get("options") or {}
        _, configuration = self.get_system_configuration()
        step = self.flowchart.plugin_manager.get(mc["step"])
        transport = self._mlff_transport(config, step)
        if transport == "MPI":
            allocated = _allocated_mpi_slots()
            if allocated is not None and allocated < 2:
                raise RuntimeError(
                    "A machine-learned forcefield runs its engine and LAMMPS as two "
                    f"MPI tasks, but this job was allocated {allocated}. Resubmit "
                    "asking for at least 2 tasks (ntasks)."
                )
            port = None
        else:
            port = _free_tcp_port()
        engine = step.get_mdi_engine_command(
            self.flowchart.executor,
            self.global_options,
            method=options.get("mdi_method_arg") or mc["method"],
            port=port,
            hostname="localhost",
            charge=configuration.charge,
            multiplicity=configuration.spin_multiplicity,
            n_atoms=configuration.n_atoms,
            mdi_method=transport,
            checkpoint=options.get("path"),
            device="cuda:0" if gpus else None,
            extra_args=options.get("engine_args") or None,
        )
        if not gpus and not (engine and _ASSIGNMENT.match(engine[0])):
            # mdi_bind.sh keeps an engine beside a GPU on one thread; on the CPU
            # it needs the job's cores, less the one LAMMPS uses
            n = _engine_cpu_threads(self.global_options.get("ncores", "available"))
            engine = ["env", f"OMP_NUM_THREADS={n}", f"MKL_NUM_THREADS={n}", *engine]
        where = f" on GPU {gpus[0]}." if gpus else " without a GPU."
        if transport == "MPI":
            text = f"The engine for '{mc['level']}' runs with LAMMPS as an MPI pair"
        else:
            text = (
                f"The engine for '{mc['level']}' is in another environment than "
                "LAMMPS, so the two communicate over TCP"
            )
        printer.important(__(text + where, indent=4 * " "))
        if transport != "MPI" and _tcp_nodelay(config):
            printer.important(
                __(
                    "MDI does not yet set TCP_NODELAY on the connection, without "
                    "which each step can wait up to about 40 ms, so on Linux SEAMM "
                    "compiles and preloads a small library that sets it.",
                    indent=4 * " ",
                )
            )
        if hasattr(step, "engine_notes"):
            # e.g. whether the engine has the fast paths it needs on a GPU
            try:
                notes = step.engine_notes(
                    self.flowchart.executor,
                    self.global_options,
                    "cuda:0" if gpus else None,
                )
            except Exception as e:
                self.logger.debug(f"Could not get the engine's notes: {e}")
                notes = []
            for note, warning in notes:
                prefix = "Warning: " if warning else ""
                printer.important(__(prefix + note, indent=4 * " "))
        if transport == "MPI":
            return _mlff_pair_script(engine, config, bin_dir / "mdi_bind.sh")
        return _mlff_tcp_script(engine, port, config, bin_dir / "mdi_bind.sh")

    def model_chemistry(self, task):
        """Compose the full ``driver:task|level`` provenance label for a task.

        ``task`` is the LAMMPS operation -- ``"SP"``, ``"OPT"``, ``"MD"``. For a
        QM-MD-over-MDI run the level (owner/type/method/...) comes from the
        ``_model_chemistry`` variable, giving e.g.
        ``"LAMMPS:MD|MOPAC:SQM@PM6-ORG"`` (or, classical, the owner collapses to
        ``"LAMMPS:MD|FF@oplsaa"``). For classical / MLFF / OpenKIM runs the
        level is not yet classified into ``type@method``, so the legacy bare
        label (``self.model``) is returned unchanged. See
        ``model_chemistry_naming.rst``.
        """
        if self.variable_exists("_model_chemistry"):
            from model_chemistry_step.grammar import compose_model_chemistry

            mc = self.get_variable("_model_chemistry")
            return compose_model_chemistry(
                driver="LAMMPS",
                task=task,
                owner=mc["owner"],
                type=mc["type"],
                method=mc["method"],
                basis=mc["basis"],
                cutoff=mc["cutoff"],
            )
        return self.model

    def run(self):
        """Run a LAMMPS simulation"""

        next_node = super().run(printer)

        # Set the model
        try:
            if self.variable_exists("_model_chemistry"):
                # Task-agnostic provenance base; the per-operation substeps
                # compose the full "LAMMPS:<task>|<level>" label via
                # model_chemistry(task). See model_chemistry_naming.rst.
                mc = self.get_variable("_model_chemistry")
                self.model = mc["level"]
            else:
                ff = self.get_variable("_forcefield")
                if ff == "OpenKIM":
                    self.model = "OpenKIM/" + self.get_variable("_OpenKIM_Potential")
                elif ff == "PyTorch":
                    path = Path(self.get_variable("_pytorch_model"))
                    self.model = "PyTorch/" + path.stem
                else:
                    self.model = ff.current_forcefield
        except Exception:
            self.model = None

        # Print headers
        printer.important(self.header)

        system_db = self.get_variable("_system_db")
        configuration = system_db.system.configuration

        n_atoms = configuration.n_atoms
        if n_atoms == 0:
            printer.important("    There is no structure! Must exit")
            self.logger.error("LAMMPS run(): there is no structure!")
            raise RuntimeError("LAMMPS run(): there is no structure!")

        # Let the analysis steps that will read this MD check the structure first
        self.check_following_analyses(configuration)

        # Initialize storage
        self._results = {}

        # Get the options
        o = self.options
        global_options = self.global_options

        # Whether to run parallel and if so, how many mpi processes
        ff_form = self.ff_form()

        # Only worry about GPU's when using PyTorch, for the moment
        if ff_form == "PyTorch" and platform.system() != "Darwin":
            maxload = float(o["gpu_load"]) / 100
            wait = o["gpu_wait_time"]
            gpus = _wait_for_gpus(maxload, None if wait == "unlimited" else float(wait))
        else:
            gpus = []
            ngpus = 0

        if global_options["parallelism"] in ("any", "mpi"):
            np = n_atoms // o["atoms_per_core"] + 1
            np_text = (
                f"The calculation can use {np} cores given {n_atoms} atoms and "
                f"the requested {o['atoms_per_core']} atoms per core."
            )
            if o["ngpus"] == "available":
                ngpus = len(gpus)
            else:
                ngpus = min(len(gpus), int(o["ngpus"]))
            if o["ncores"] != "available":
                if np > int(o["ncores"]):
                    np_text += (
                        " The maximum number of cores for LAMMPS limits the "
                        f"number of cores used to  {o['ncores']}."
                    )
                np = min(np, int(o["ncores"]))
            if global_options["ncores"] != "available":
                if np > int(global_options["ncores"]):
                    np_text += (
                        " The actual number of cores is limited by the global "
                        f"limit of {int(global_options['ncores'])} cores."
                    )
                np = min(np, int(global_options["ncores"]))
            if ff_form == "PyTorch" and platform.system() != "Darwin":
                if ngpus < 1 and self.mlff_model_chemistry() is not None:
                    # The engine of an MLFF model chemistry runs on the CPU
                    np = 1
                    np_text += _MLFF_ON_CPU
                elif ngpus < 1:
                    printer.important(
                        "    Using PyTorch forcefield, but no GPUs! Must exit."
                    )
                    raise RuntimeError(
                        "LAMMPS: PyTorch forcefield requested, but no GPUs!"
                        " Must stop."
                    )
                else:
                    np = min(np, ngpus)
                    ngpus = np
                    np_text += (
                        " Because LAMMPS is using a PyTorch forcefield, the number "
                        f"of GPUs must be the same as the number of cores = {ngpus}."
                    )
        else:
            np_text = "No parallelism requested: {global_options['parallelism']}"
            np = 1
            if ff_form == "PyTorch":
                if ngpus < 1 and self.mlff_model_chemistry() is not None:
                    np_text += _MLFF_ON_CPU
                elif ngpus < 1:
                    printer.important(
                        "    Using PyTorch forcefield, but no GPUs! Must exit."
                    )
                    raise RuntimeError(
                        "LAMMPS: PyTorch forcefield requested, but no GPUs!"
                        " Must stop."
                    )
                else:
                    ngpus = 1
                    np_text += (
                        " Because LAMMPS is using a PyTorch forcefield, the number "
                        f"of GPUs must be the same as the number of cores = {ngpus}."
                    )
        all_gpus = gpus
        if ngpus > 0:
            gpus = gpus[0:ngpus]

        # Get to the actual work
        self.subflowchart.root_directory = self.flowchart.root_directory

        files = {}

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        extras = {}

        history_nodes = []

        # Create overall directory for the lammps step

        os.makedirs(self.directory, exist_ok=True)

        files = {}

        control = []
        if self.variable_exists("_model_chemistry"):
            mc = self.get_variable("_model_chemistry")
            control.append(["model_chemistry", mc["level"]])
        else:
            ff = self.get_variable("_forcefield")
            if ff == "OpenKIM":
                potential = self.get_variable("_OpenKIM_Potential")
                control.append(["forcefield", "OpenKIM " + potential])
            elif ff == "PyTorch":
                model = self.get_variable("_pytorch_model")
                control.append(["forcefield", "PyTorch " + model])
            else:
                control.append(["forcefield", ff.current_forcefield])
        while node is not None:
            P = node.parameters.current_values_to_dict(
                context=seamm.flowchart_variables._data
            )

            # For the timing data, get the parameters used
            tmp = [node.title]
            for key, value in P.items():
                if isinstance(value, units_class):
                    tmp.append((key, f"{value:~P}"))
                else:
                    tmp.append((key, value))
            control.append(tmp)

            if isinstance(node, lammps_step.Initialization):
                # Reset whether the cell needs to be triclinic
                self.force_triclinic = False

                initialization_header, eex = self._get_node_input(
                    node=node, extras={"read_data": True}
                )
                self.eex = eex

                # Need the masses for e.g. timestep
                self._data["masses"] = [float(v[0]) for v in eex["masses"]]

                self._initialization_node = node

                files["input.dat"] = copy.deepcopy(initialization_header)

                # Find the bond & angle types as needed for shake/rattle
                shake = self.shake_fix(P, eex)
                if shake != "":
                    extras["shake"] = shake
            else:
                history_nodes.append(node)
            node = node.next()

        files = self._prepare_input(
            files,
            nodes=history_nodes,
            read_restart=False,
            write_restart=True,
            extras=extras,
        )

        # Procrastinate making the strucuture file as long as possible,
        # because it may need to be forced to be triclinic by substeps

        (
            structure_data,
            pair_table,
            bond_table,
            angle_table,
            dihedral_table,
        ) = self.structure_data()

        files["structure.dat"] = structure_data
        self.check_barostats(history_nodes)
        if "forcefield" in self.eex:
            files["forcefield.dat"] = self.eex["forcefield"]

        if bond_table != "":
            files["tabulated_bonds.dat"] = bond_table
        if angle_table != "":
            files["tabulated_angles.dat"] = angle_table
        if dihedral_table != "":
            files["tabulated_dihedrals.dat"] = dihedral_table

        self.logger.debug("structure.dat:\n" + files["structure.dat"])

        # And run
        self._timing_control = list(control)
        control = []

        # Claim the GPUs only now, for the run: other runs then choose others
        gpu_claims = []
        if ngpus > 0:
            gpus, gpu_claims = _claim_gpus(all_gpus, ngpus)
        try:
            files = self._execute_single_sim(
                files, np=np, np_text=np_text, ngpus=ngpus, gpus=gpus
            )
        finally:
            _release_gpus(gpu_claims)

        # Analyze any GPU logs
        values = {}
        gpu_logs = list(self.wd.glob("gpu_*.log"))
        if len(gpu_logs) > 0:
            with open(self.wd / "gpu_statistics.log", "w") as fd:
                stats = lammps_step.print_multi_gpu_summary(gpu_logs, fd=fd)

            if stats is not None:
                # Add the values to the data to store
                values["average gpu utilization"] = stats["Average GPU utilization"]
                values["gpu memory used"] = stats["Maximum GPU memory used"]
                values["gpu memory free"] = stats["Minimum free GPU memory"]
                # Print a table with the gpu statistics to the output.
                table = {
                    "Statistic": [*stats.keys()],
                    "Value": [*stats.values()],
                    "Units": [],
                }
                tmp = table["Units"]
                for key in stats.keys():
                    if "%" in key:
                        tmp.append("%")
                    elif "memory" in key:
                        tmp.append("GiB")
                    elif "utilization" in key:
                        tmp.append("%")
                    elif "time" in key:
                        tmp.append("s")
                    else:
                        tmp.append("")

                max_atoms = int(n_atoms * (1 + 100 / stats["% GPU memory used"]))
                values["gpu maximum atoms"] = max_atoms
                table["Statistic"].append("# of atoms")
                table["Value"].append(n_atoms)
                table["Units"].append("")
                table["Statistic"].append("Maximum # of atoms possible")
                table["Value"].append(max_atoms)
                table["Units"].append("")

                tmp = tabulate(
                    table,
                    headers="keys",
                    tablefmt="rounded_outline",
                    colalign=("center", "decimal", "left"),
                    disable_numparse=True,
                )
                length = len(tmp.splitlines()[0])
                text_lines = []
                text_lines.append("GPU Usage Statistics".center(length))
                text_lines.append(tmp)
                printer.important(
                    textwrap.indent("\n".join(text_lines), self.indent + 7 * " ")
                )
                printer.important("")

        # And analyze the results
        self.analyze(nodes=history_nodes, values=values)

        self._trajectory = []

        printer.normal("")

        return next_node

    def _execute_single_sim(
        self, files, np=1, np_text="", ngpus=0, gpus=[], return_files=None
    ):
        """
        Step #1: Execute input file
        """

        return_files = [
            "stdout.txt",
            "stderr.txt",
            "summary_*.txt",
            "trajectory_*.seamm_trj",
            "*.json",
            "*.trj",
            "*.dump_trj",
            "*.dump_trj.gz",
            "*.restart.*",
            "*.dump",
            "*.dump.*",
            "*.log",
            "*.dat",
            "log.lammps",
            "log.cite",
            "run_lammps",
        ]
        # and what an MLFF model chemistry's engine writes, such as a committee's
        # spreads, which the cleanup of the run would otherwise remove
        mc = self.mlff_model_chemistry()
        if mc is not None:
            return_files += list((mc.get("options") or {}).get("engine_files") or [])

        # Check for already having run
        path = Path(self.directory) / "success.dat"
        if path.exists():
            result = {}
            path = Path(self.directory) / "stdout.txt"
            if path.exists():
                result["stdout"] = path.read_text()
            path = Path(self.directory) / "stderr.txt"
            if path.exists():
                result["stderr"] = path.read_text()
            else:
                result["stderr"] = ""
        else:
            # For QM-MD over MDI the QM engine does essentially all the work and
            # is the bottleneck; LAMMPS is just the driver, so run it on a single
            # core (the OMP/MKL thread caps below keep the engine single-threaded
            # too). See campaigns/2026-06-22/NOTES_C.rst.
            if self.ff_form() == "MDI/QM":
                np = 1

            # Set up the computational limits and get the computational enviroment
            cl = {"NTASKS": np}
            ce = seamm_exec.computational_environment(cl)

            # For PyTorch need at least 1 GPU and only 1 MPI process per GPU
            ff_form = self.ff_form()
            if ff_form == "PyTorch" and ngpus > 0:
                if "NGPUS" not in ce:
                    ce["NGPUS"] = ngpus
                if "NTASKS" not in ce:
                    ce["NTASKS"] = ngpus
                else:
                    if ce["NTASKS"] > ngpus:
                        ce["NTASKS"] = ngpus

            executor = self.flowchart.executor

            # Read configuration file for LAMMPS if it exists
            executor_type = executor.name
            full_config = configparser.ConfigParser()
            ini_dir = Path(self.global_options["root"]).expanduser()
            path = ini_dir / "lammps.ini"

            # If the config file doesn't exists, get the default
            resources = importlib.resources.files("lammps_step") / "data"
            if not path.exists():
                ini_text = (resources / "lammps.ini").read_text()
                txt_config = Configuration(path)
                txt_config.from_string(ini_text)

                # Work out the conda info needed
                conda_exe = os.environ.get(
                    "CONDA_EXE", "<please fill in the path to conda>"
                )
                txt_config.set_value("local", "conda", conda_exe)
                txt_config.set_value("local", "conda-environment", "seamm-lammps")
                txt_config.save()
                if conda_exe == "<please fill in the path to conda>":
                    msg = (
                        f"Cannot locate the conda executable. Please edit {path}\n"
                        "look for the <conda> key in the local section, and put in\n"
                        "the path, which you can get using 'echo $CONDA_EXE'."
                    )
                    printer.important(msg)
                    raise RuntimeError(msg)
                else:
                    printer.normal(f"Wrote the LAMMPS configuration file to {path}")
                    printer.normal("")

            # Ensure the helper scripts for GPU/MDI execution are installed
            # These are copied from the package data directory
            # to ~/SEAMM/bin/ so they can be referenced in lammps.ini commands.
            bin_dir = ini_dir / "bin"
            bin_dir.mkdir(exist_ok=True)

            helper_scripts = [
                "cpu_bind.sh",
                "gpu_bind.sh",
                "mdi_bind.sh",
                "mdi_monitor.sh",
            ]

            for script_name in helper_scripts:
                src_text = (resources / script_name).read_text()
                message = _install_helper_script(bin_dir / script_name, src_text)
                if message is not None:
                    self.logger.info(message)
                    if "Upgraded" in message:
                        printer.normal(__(message, indent=4 * " "))

            full_config.read(ini_dir / "lammps.ini")

            if executor_type not in full_config:
                path = shutil.which("lmp")
                mpi_path = shutil.which("mpirun")
                if path is None or mpi_path is None:
                    raise RuntimeError(
                        f"No section for '{executor_type}' in LAMMPS ini file "
                        f"({ini_dir / 'lammps.ini'}), nor in the defaults, nor "
                        "in the path!"
                    )
                else:
                    txt_config = Configuration(path)
                    txt_config.add_section(executor_type)
                    txt_config.set_value(executor_type, "installation", "local")
                    txt_config.set_value(
                        executor_type,
                        "code",
                        f"{mpi_path} -np {{NTASKS}} {path}",
                    )
                    txt_config.set_value(
                        executor_type,
                        "python",
                        f"mpirun -np {{NTASKS}} {shutil.which('python')}",
                    )
                    txt_config.save()
                    printer.normal(f"Wrote the LAMMPS configuration file to {path}")
                    printer.normal("")
                    full_config.read(ini_dir / "lammps.ini")
                    full_config.set(executor_type, "code", str(path))

            config = dict(full_config.items(executor_type))

            # Environment variables
            env = _launch_environment(gpus if "NGPUS" in ce else [])

            # The forcefield file for PyTorch/MDI runs (keyed off the forcefield
            # form, not the model label, which is now a grammar string).
            #
            # Two ways in, because engines differ in how they take the model.
            # SEAMM_FF is the environment variable mace-mdi reads. {MODEL} is
            # for engines that take it as an argument -- `xnn mdi --ckpt
            # {MODEL}` -- and is substituted into the command by the executor
            # along with {NTASKS} and the rest. Setting both leaves the choice
            # of engine entirely to gpu-code in lammps.ini, with nothing here
            # needing to know which one is configured.
            if ff_form == "PyTorch":
                model_path = self.pytorch_model()
                env["SEAMM_FF"] = model_path
                ce["MODEL"] = model_path
            else:
                env["SEAMM_FF"] = "Unknown"

            executor_type = executor.name
            if executor_type not in full_config:
                raise RuntimeError(
                    f"No section for '{executor_type}' in LAMMPS ini file "
                    f"({ini_dir / 'lammps.ini'})"
                )
            config = dict(full_config.items(executor_type))

            # Setup the command lines
            cmd = []

            # A GPU in the allocation is not on its own a reason to use
            # gpu-code. NGPUS comes from the scheduler (seamm_exec reads
            # SLURM_JOB_GPUS), so it is set for any job that asked for a
            # GPU, whatever the forcefield. But an MDI gpu-code launches a
            # machine-learned engine, which only makes sense for a PyTorch
            # forcefield: run a classical one that way and the engine is
            # handed SEAMM_FF="Unknown" and dies with FileNotFoundError.
            # A Kokkos-style gpu-code has no such restriction -- that is
            # how a classical forcefield is meant to use a GPU -- so key
            # off whether the command drives MDI, not on GPUs alone.
            gpu_code = config.get("gpu-code", "")
            use_gpu_code = "NGPUS" in ce
            # MDI/QM builds its own launch script below and never consults
            # gpu-code, so it needs the decision but not the advice.
            if (
                use_gpu_code
                and ff_form != "MDI/QM"
                and not _gpu_code_is_usable(ff_form, gpu_code)
            ):
                printer.important(
                    __(
                        f"The '{ff_form}' forcefield cannot use the MDI "
                        "engine that gpu-code in lammps.ini launches, which "
                        "serves machine-learned potentials. Running on the "
                        "CPU instead. To use a GPU for this forcefield, set "
                        "gpu-code to a Kokkos command and gpu-cmd-args to "
                        "the matching accelerator flags.",
                        indent=4 * " ",
                    )
                )
                use_gpu_code = False
            elif use_gpu_code and not _gpu_code_is_usable(ff_form, gpu_code):
                use_gpu_code = False

            if ff_form == "MDI/QM":
                # QM-MD over MDI: launch the QM engine (the listener, in its own
                # conda environment) and the LAMMPS driver together, rendezvousing
                # over TCP. See campaigns/2026-06-22/NOTES_C.rst.
                _, configuration = self.get_system_configuration()
                engine_argv, port = self._mdi_engine_launch(configuration)
                files["mdi_launch.sh"] = _mdi_launch_script(
                    engine_argv, port, config, ce
                )
                cmd = ["bash", "mdi_launch.sh"]
                if NODELAY_SOURCE in files["mdi_launch.sh"]:
                    files[NODELAY_SOURCE] = _nodelay_source()
            elif self.mlff_model_chemistry() is not None:
                # An MLFF from the Model Chemistry step: its engine and LAMMPS run
                # as an MPI pair, launched as lammps.ini's gpu-code used to.
                files["mdi_launch.sh"] = self._mlff_launch(config, bin_dir, gpus)
                cmd = ["bash", "mdi_launch.sh"]
                if NODELAY_SOURCE in files["mdi_launch.sh"]:
                    files[NODELAY_SOURCE] = _nodelay_source()
            elif "run_lammps" in files:
                cmd = ["{python}", "run_lammps"]
                if (
                    not use_gpu_code
                    and "cmd-args" in config
                    and config["cmd-args"] != ""
                ):
                    cmd.extend(["--cmd-args", config["cmd-args"]])
                if (
                    use_gpu_code
                    and "gpu-cmd-args" in config
                    and config["gpu-cmd-args"] != ""
                ):
                    cmd.extend(["--cmd-args", config["gpu-cmd-args"]])
            else:
                if use_gpu_code:
                    # This plug-in used to ship its own copy of the MACE MDI
                    # engine as mace_mdi.py. It is now maintained once, in the
                    # lammps-mdi package, and reached through its `mace-mdi`
                    # console script. An existing ~/SEAMM/bin/mace_mdi.py is
                    # left alone, so a configuration pointing at it keeps
                    # working -- on a copy that predates several fixes.
                    # An MDI gpu-code launches the engine and the driver as
                    # separate MPI ranks. Those counts are written into the
                    # command, not derived from the allocation, so a job given
                    # fewer tasks than the command asks for dies inside mpirun
                    # with "All nodes which are allocated for this job are
                    # already filled" -- which says nothing about MDI, ranks or
                    # how to fix it. Say it plainly instead, before launching.
                    needed = _required_mpi_slots(config.get("gpu-code", ""))
                    allocated = _allocated_mpi_slots()
                    if (
                        needed is not None
                        and allocated is not None
                        and needed > allocated
                    ):
                        raise RuntimeError(
                            f"The GPU command in lammps.ini needs {needed} MPI "
                            f"tasks, but this job was allocated {allocated}. An "
                            "MDI run needs one task for the engine and one for "
                            "the LAMMPS driver, so asking for a single GPU is "
                            "not enough on its own. Resubmit requesting at "
                            f"least {needed} tasks (ntasks)."
                        )

                    if "mace_mdi.py" in config.get("gpu-code", ""):
                        printer.important(
                            __(
                                "Warning: gpu-code in lammps.ini runs the bundled "
                                "mace_mdi.py. That copy is no longer shipped or "
                                "updated with this plug-in and is missing later "
                                "fixes. Replace 'python ~/SEAMM/bin/mace_mdi.py' "
                                "with 'mace-mdi', from the lammps-mdi package "
                                "installed in the LAMMPS environment.",
                                indent=4 * " ",
                            )
                        )
                    cmd = ["{gpu-code}"]
                    if "gpu-cmd-args" in config and config["gpu-cmd-args"] != "":
                        cmd.extend(config["gpu-cmd-args"].split())
                else:
                    cmd = ["{code}"]
                    if "cmd-args" in config and config["cmd-args"] != "":
                        cmd.extend(config["cmd-args"].split())
                cmd.extend(["-in", "input.dat"])

            # Say what is producing the forces -- a model chemistry (QM via
            # MDI) or a forcefield -- before reporting the parallelism.
            if ff_form == "MDI/QM":
                printer.important(
                    __(
                        f"Using the model chemistry '{self.model}', with the "
                        "energy and forces evaluated by a QM engine via MDI.",
                        indent=4 * " ",
                    )
                )
            elif self.model:
                printer.important(f"    Using the forcefield '{self.model}'.")

            if use_gpu_code:
                printer.important(
                    f"    LAMMPS running with {np} processes and {ce['NGPUS']} gpus."
                )
            elif ff_form == "MDI/QM":
                printer.important(
                    __(
                        "LAMMPS is running on a single core: it only drives the "
                        "dynamics and does very little work, while the QM engine "
                        "(reached over MDI) computes the energy and forces.",
                        indent=4 * " ",
                    )
                )
            else:
                # The batch system's allocation may allow fewer MPI tasks than
                # the step asked for: say how many LAMMPS actually gets
                ntasks = int(ce.get("NTASKS", np))
                printer.important(f"    LAMMPS using MPI with {ntasks} processes.")
                if ntasks < np:
                    np_text += allocation_note(np, ntasks)
                if np_text != "":
                    printer.important("")
                    printer.important(__(np_text, indent=4 * " "))
            printer.important("")

            cmd.extend([">", "stdout.txt", "2>", "stderr.txt"])

            t0 = time.time_ns()

            result = executor.run(
                cmd=cmd,
                config=config,
                directory=self.directory,
                files=files,
                return_files=return_files,
                in_situ=True,
                shell=True,
                ce=ce,
                env=env,
            )

            t = (time.time_ns() - t0) / 1.0e9
            # The record must not depend on this function's locals: the
            # configuration is only bound in the MDI branch above (2026.10.6
            # raised UnboundLocalError here for every forcefield run).
            self.record_timing(ce, t, result)

            if not result:
                self.logger.error("There was an error running LAMMPS")
                return None

            if self.logger.isEnabledFor(logging.DEBUG):
                self.logger.debug("\n" + pprint.pformat(_summarize_run_result(result)))

            # Did LAMMPS fail? Say why here, rather than failing later on missing
            # output. (No success.dat is written, so a rerun runs LAMMPS again.)
            self.check_for_lammps_error(result)

            # f = os.path.join(self.directory, "stdout.txt")
            # with open(f, mode="w") as fd:
            #     fd.write(result["stdout"])
            # f = os.path.join(self.directory, "stderr.txt")
            # with open(f, mode="w") as fd:
            #     fd.write(result["stderr"])

        # Add the citations, getting the version from stdout and any citations.
        # The returned files are in the directory; the executor does not read them.
        path = Path(self.directory) / "log.cite"
        cite = path.read_text() if path.exists() else None
        self._add_lammps_citations(result["stdout"], cite=cite)

        initialization_header, eex = self._get_node_input(
            node=self._initialization_node, extras={"read_data": False}
        )
        self.eex = eex
        files["input.dat"] = copy.deepcopy(initialization_header)

        # Write a small file to say that LAMMPS ran successfully, so cancel
        # skip if rerunning.
        path = Path(self.directory) / "success.dat"
        path.write_text("success")

        return files

    def _prepare_input(
        self,
        files,
        nodes=None,
        iteration=0,
        read_restart=False,
        write_restart=False,
        extras=None,
    ):
        _, configuration = self.get_system_configuration()
        python_script = None
        postscript = None
        if isinstance(nodes, list) is False:
            node_ids = [nodes._id[1]]
            todo = self._get_node_input(node=nodes, extras=extras)
            new_input_data = todo["script"]
            if todo["postscript"] is not None:
                postscript = todo["postscript"]
            if todo["use python"] and "python script" in todo:
                python_script = todo["python script"]
        else:
            node_ids = []
            new_input_data = []
            for n in nodes:
                node_ids.append(n._id[1])
                todo = self._get_node_input(node=n, extras=extras)
                new_input_data += todo["script"]
                if todo["postscript"] is not None:
                    postscript = todo["postscript"]
                if todo["use python"] and "python script" in todo:
                    python_script = todo["python script"]

        if postscript is None:
            if configuration.periodicity == 0:
                new_input_data.append(
                    "write_dump         all custom  final.dump id xu yu zu vx vy vz"
                    " modify flush yes sort id"
                )
            else:
                new_input_data.append(
                    "write_dump         all custom  final.dump id xsu ysu zsu vx vy vz"
                    " modify flush yes sort id"
                )
            new_input_data.append("")
            new_input_data.append("")
            new_input_data.append("info               computes fixes dumps out log")

        files["input.dat"] += new_input_data

        files["input.dat"] = "\n".join(files["input.dat"])
        self.logger.debug("input.dat:\n" + files["input.dat"])

        if postscript is not None:
            if configuration.periodicity == 0:
                postscript.append(
                    "write_dump         all custom  final.dump id xu yu zu vx vy vz"
                    " modify flush yes sort id"
                )
            else:
                postscript.append(
                    "write_dump         all custom  final.dump id xsu ysu zsu vx vy vz"
                    " modify flush yes sort id"
                )
            postscript.append("")
            postscript.append("info               computes fixes dumps out log")
            files["input_post.dat"] = "\n".join(postscript)
        if python_script is not None:
            files["run_lammps"] = python_script

        return files

    def _get_node_input(self, node=None, extras=None):
        try:
            ret = node.get_input(extras=extras)
        except Exception as e:
            print("Error running LAMMPS flowchart: {} in {}".format(str(e), str(node)))
            self.logger.critical(
                "Error running LAMMPS flowchart: {} in {}".format(str(e), str(node))
            )
            raise
        except:  # noqa: E722
            print(
                "Unexpected error running LAMMPS flowchart: {} in {}".format(
                    sys.exc_info()[0], str(node)
                )
            )
            self.logger.critical(
                "Unexpected error running LAMMPS flowchart: {} in {}".format(
                    sys.exc_info()[0], str(node)
                )
            )
            raise
        return ret

    def following_nodes(self):
        """The steps after this one in execution order, up to the next LAMMPS step.

        Loop bodies after this step are not entered; a loop around this step is
        left through its exit. Only the flowchart's edges are followed, so no
        step's state changes.
        """
        flowchart = self.flowchart
        result = []
        seen = {id(self)}
        node = self
        while True:
            following = None
            for edge in flowchart.edges(node, direction="out"):
                if edge.edge_type != "execution":
                    continue
                if edge.edge_subtype in ("next", "exit"):
                    following = edge.node2
            if following is None or id(following) in seen:
                break
            if isinstance(following, LAMMPS):
                break
            seen.add(id(following))
            result.append(following)
            node = following
        return result

    def check_following_analyses(self, configuration):
        """Ask the analysis steps that will read this MD whether the structure suits.

        Each step after this one, up to the next LAMMPS step, that has a
        ``check_configuration`` method is asked. One that cannot use the structure,
        e.g. an analysis needing a periodic system, raises an exception, so the job
        stops before the MD rather than after it.

        A check may also return warnings, e.g. that the system is small for the
        analysis, which are printed before the MD runs.

        Parameters
        ----------
        configuration : molsystem._Configuration
            The structure about to be simulated.
        """
        for node in self.following_nodes():
            check = getattr(node, "check_configuration", None)
            if check is not None:
                warnings = check(configuration)
                for warning in warnings or []:
                    printer.important(
                        __(f"Warning: {warning}", indent=self.indent + 4 * " ")
                    )

    def check_for_lammps_error(self, result=None):
        """Raise an error, with a readable message, if LAMMPS failed."""
        directory = Path(self.directory)
        path = directory / "log.lammps"
        log = path.read_text(errors="replace") if path.exists() else ""
        screen = ""
        if result is not None and isinstance(result.get("stdout"), str):
            screen = result["stdout"]
        else:
            path = directory / "stdout.txt"
            if path.exists():
                screen = path.read_text(errors="replace")
        error = lammps_error(log, screen)
        if error is None:
            return
        text = lammps_error_message(error, directory)
        # FormattedText formats the text, so protect any braces in it
        safe = text.replace("{", "{{").replace("}", "}}")
        printer.important(__(safe, indent=self.indent + 4 * " ", wrap=True))
        printer.important("")
        raise RuntimeError(text)

    def check_barostats(self, nodes):
        """Check that the box suits the barostats, before running LAMMPS.

        LAMMPS writes one box for all the steps, and press/berendsen cannot work
        with a triclinic box, which a later step allowing shear or a cell that is
        not orthorhombic requires.
        """
        if not self._triclinic_box:
            return
        for node in nodes:
            if getattr(node, "uses_berendsen", False):
                if self.force_triclinic:
                    why = "another step in this LAMMPS step allows the cell to shear"
                else:
                    why = "the cell is not orthorhombic"
                raise RuntimeError(
                    f"{node.title}: the Berendsen barostat cannot be used with a "
                    f"triclinic cell, which is needed because {why}. Use the "
                    "Nose-Hoover barostat, or put the steps in separate LAMMPS steps."
                )

    def structure_data(self):
        """Create the LAMMPS structure file from the energy expression"""
        lines = []
        pair_table = []
        bond_table = []
        angle_table = []
        dihedral_table = []

        eex = self.eex

        lines.append("Structure file for LAMMPS generated by a MolSSI flowchart")
        lines.append("{:10d} atoms".format(eex["n_atoms"]))
        lines.append("{:10d} atom types".format(eex["n_atom_types"]))
        if "n_bonds" in eex and eex["n_bonds"] > 0:
            lines.append("{:10d} bonds".format(eex["n_bonds"]))
            lines.append("{:10d} bond types".format(eex["n_bond_types"]))
        if "n_angles" in eex and eex["n_angles"] > 0:
            lines.append("{:10d} angles".format(eex["n_angles"]))
            lines.append("{:10d} angle types".format(eex["n_angle_types"]))
        if "n_torsions" in eex and eex["n_torsions"] > 0:
            lines.append("{:10d} dihedrals".format(eex["n_torsions"]))
            lines.append("{:10d} dihedral types".format(eex["n_torsion_types"]))
        if "n_oops" in eex and eex["n_oops"] > 0:
            lines.append("{:10d} impropers".format(eex["n_oops"]))
            lines.append("{:10d} improper types".format(eex["n_oop_types"]))

        # Find the box limits
        self._triclinic_box = False
        periodicity = eex["periodicity"]
        if periodicity == 3:
            a, b, c, alpha, beta, gamma = eex["cell"]
            lx, ly, lz, xy, xz, yz = LAMMPS.cell_to_box(a, b, c, alpha, beta, gamma)

            lines.append("{} {} xlo xhi".format(0.0, lx))
            lines.append("{} {} ylo yhi".format(0.0, ly))
            lines.append("{} {} zlo zhi".format(0.0, lz))

            xy = xy if abs(xy) > 1.0e-06 else 0.0
            xz = xz if abs(xz) > 1.0e-06 else 0.0
            yz = yz if abs(yz) > 1.0e-06 else 0.0

            self._triclinic_box = (
                self.force_triclinic or xy != 0.0 or xz != 0.0 or yz != 0.0
            )
            if self._triclinic_box:
                lines.append("{} {} {} xy xz yz".format(xy, xz, yz))
        else:
            x, y, z, index = eex["atoms"][0]
            xlo = xhi = x
            ylo = yhi = y
            zlo = zhi = z
            for x, y, z, index in eex["atoms"]:
                xlo = x if x < xlo else xlo
                xhi = x if x > xhi else xhi
                ylo = y if y < ylo else ylo
                yhi = y if y > yhi else yhi
                zlo = z if z < zlo else zlo
                zhi = z if z > zhi else zhi

            # Some extra space....
            xlo -= 10.0
            xhi += 10.0
            ylo -= 10.0
            yhi += 10.0
            zlo -= 10.0
            zhi += 10.0

            lines.append(f"{xlo:12.6f} {xhi:12.6f} xlo xhi")
            lines.append(f"{ylo:12.6f} {yhi:12.6f} ylo yhi")
            lines.append(f"{zlo:12.6f} {zhi:12.6f} zlo zhi")

        # the atoms and their masses, etc.
        lines.append("")
        lines.append("Atoms")
        lines.append("")

        if "charges" in eex:
            if "molecule" in eex:
                for i, mol, xyz_index, q in zip(
                    range(1, eex["n_atoms"] + 1),
                    eex["molecule"],
                    eex["atoms"],
                    eex["charges"],
                ):
                    x, y, z, index = xyz_index
                    lines.append(
                        f"{i:6d} {mol + 1:6d} {index:6d} {q:12.6f} {x:12.7f} "
                        f"{y:12.7f} {z:12.7f}"
                    )
            else:
                for i, xyz_index, q in zip(
                    range(1, eex["n_atoms"] + 1), eex["atoms"], eex["charges"]
                ):
                    x, y, z, index = xyz_index
                    lines.append(
                        f"{i:6d} {index:6d} {q:12.6f} {x:12.7f} {y:12.7f} {z:12.7f}"
                    )
        else:
            for i, xyz_index in enumerate(eex["atoms"]):
                x, y, z, index = xyz_index
                lines.append(f"{i+1:6d} {index:6d} {x:12.7f} {y:12.7f} {z:12.7f}")

        _, configuration = self.get_system_configuration()
        if configuration.atoms.have_velocities:
            lines.append("")
            lines.append("Velocities")
            lines.append("")
            for i, vxyz in enumerate(
                configuration.atoms.get_velocities(fractionals=False), start=1
            ):
                vx, vy, vz = vxyz
                lines.append(f"{i:6d} {vx:12.7f} {vy:12.7f} {vz:12.7f}")

        lines.append("")
        lines.append("Masses")
        lines.append("")
        self._data["masses"] = []
        for i, parameters in enumerate(eex["masses"], start=1):
            mass, itype = parameters
            lines.append("{:6d} {} # {}".format(i, mass, itype))
            self._data["masses"].append(float(mass))

        if "extra molecule data" in eex:
            lines.append("")
            lines.append("Molecules")
            lines.append("")
            for i, molecule in enumerate(eex["extra molecule data"], start=1):
                lines.append(f"{i:6d} {molecule + 1:6d}")

        # nonbonds
        if "nonbond parameters" in eex:
            # If using a hybrid/overlay form for e.g. Dreiding nonbonds, cannot use the
            # section in the structure file because it must have N or N*(N+1) lines,
            # which is not what we want
            forms = set([v[0] for v in eex["nonbond parameters"]])
            use_hybrid = len(forms) > 1

            if not use_hybrid:
                lines.append("")
                if len(eex["nonbond parameters"]) > eex["n_atom_types"]:
                    lines.append("PairIJ Coeffs")
                else:
                    lines.append("Pair Coeffs")
                lines.append("")
                i = 1
                j = 1
                for parameters in eex["nonbond parameters"]:
                    form, values, types, parameters_type, real_types = parameters
                    if form == "nonbond(9-6)":
                        lines.append(
                            pline(i, "", values, ["eps", "rmin"], types, real_types)
                        )
                        i += 1
                    elif form == "nonbond(12-6)":
                        lines.append(
                            pline(i, "", values, ["eps", "sigma"], types, real_types)
                        )
                        i += 1
                    elif form == "buckingham":
                        lines.append(
                            pline(
                                j,
                                f"{i:6d}",
                                values,
                                ["A", "rho", "C"],
                                types,
                                real_types,
                            )
                        )
                        if j == i:
                            i += 1
                            j = 1
                        else:
                            j += 1
        # bonds
        if "n_bonds" in eex and eex["n_bonds"] > 0:
            lines.append("")
            lines.append("Bonds")
            lines.append("")
            for counter, tmp in zip(range(1, eex["n_bonds"] + 1), eex["bonds"]):
                i, j, index = tmp
                lines.append("{:6d} {:6d} {:6d} {:6d}".format(counter, index, i, j))

            lines.append("")
            lines.append("Bond Coeffs")
            lines.append("")

            styles = set([bond_style[v[0]] for v in eex["bond parameters"]])
            use_hybrid = len(styles) > 1

            for counter, parameters in zip(
                range(1, eex["n_bond_types"] + 1), eex["bond parameters"]
            ):
                form, values, types, parameters_type, real_types = parameters
                if form == "rigid_bond":
                    function = "harmonic" if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        {**values, "K2": rigid_bond_k},
                        ["K2", "R0"],
                        types,
                        real_types,
                    )
                elif form == "quadratic_bond":
                    function = "harmonic" if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        values,
                        ["K2", "R0"],
                        types,
                        real_types,
                    )
                elif form == "quartic_bond":
                    function = "class2" if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        values,
                        ["R0", "K2", "K3", "K4"],
                        types,
                        real_types,
                    )
                lines.append(line)

        # angles
        if "n_angles" in eex and eex["n_angles"] > 0:
            lines.append("")
            lines.append("Angles")
            lines.append("")
            for counter, tmp in enumerate(eex["angles"], start=1):
                i, j, k, index = tmp
                lines.append(
                    "{:6d} {:6d} {:6d} {:6d} {:6d}".format(counter, index, i, j, k)
                )

            lines.append("")
            lines.append("Angle Coeffs")
            lines.append("")

            quartic_function = "class2" if "n_bond-bond_types" in eex else "quartic"
            styles = set([angle_style[v[0]] for v in eex["angle parameters"]])
            use_hybrid = len(styles) > 1

            for counter, parameters in enumerate(eex["angle parameters"], start=1):
                form, values, types, parameters_type, real_types = parameters
                if form == "rigid_angle":
                    function = "harmonic" if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        {**values, "K2": rigid_angle_k},
                        ["K2", "Theta0"],
                        types,
                        real_types,
                    )
                elif form == "quadratic_angle":
                    function = "harmonic" if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        values,
                        ["K2", "Theta0"],
                        types,
                        real_types,
                    )
                elif form == "quartic_angle":
                    function = quartic_function if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        values,
                        ["Theta0", "K2", "K3", "K4"],
                        types,
                        real_types,
                    )
                elif form == "cosine":
                    function = "cosine" if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        values,
                        ["K2"],
                        types,
                        real_types,
                    )
                elif form == "cosine/squared":
                    function = "cosine/squared" if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        values,
                        ["K2", "Theta0"],
                        types,
                        real_types,
                    )
                elif form == "simple_fourier_angle":
                    function = "fourier/simple" if use_hybrid else ""
                    values["phase"] = "-1"
                    values["original phase"] = "-1"
                    line = pline(
                        counter,
                        function,
                        values,
                        ["K", "phase", "n"],
                        types,
                        real_types,
                    )
                elif form == "tabulated_angle":
                    function = "table" if use_hybrid else ""
                    key = f"{types[0]}-{types[1]}-{types[2]}"
                    line = f"{counter:6d} {function} tabulated_angles.dat {key}"
                    angle_table.extend(self.angle_table(key, values))
                    line += (
                        f" # {types[0]:6} {types[1]:6} {types[2]:6}"
                        f" --> {real_types[0]:6} {real_types[1]:6} {real_types[2]:6} "
                    )
                lines.append(line)

            # bond-bond coefficients, which must match angles in order & number
            if "n_bond-bond_types" in eex:
                lines.append("")
                lines.append("BondBond Coeffs")
                lines.append("")
                for i, parameters, angles in zip(
                    range(1, eex["n_bond-bond_types"] + 1),
                    eex["bond-bond parameters"],
                    eex["angle parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    angle_form = angles[0]
                    if angle_form == "quartic_angle":
                        fn = "class2" if use_hybrid else ""
                        lines.append(
                            pline(i, fn, values, ["K", "R10", "R20"], types, real_types)
                        )
                    else:
                        lines.append(pline(i, "skip", [], types, real_types))

                # bond-angles coefficients, which must match angles in order &
                # number
                lines.append("")
                lines.append("BondAngle Coeffs")
                lines.append("")
                for i, parameters, angles in zip(
                    range(1, eex["n_bond-angle_types"] + 1),
                    eex["bond-angle parameters"],
                    eex["angle parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    angle_form = angles[0]
                    if angle_form == "quartic_angle":
                        fn = "class2" if use_hybrid else ""
                        lines.append(
                            pline(
                                i,
                                fn,
                                values,
                                ["K12", "K23", "R10", "R20"],
                                types,
                                real_types,
                            )
                        )
                    else:
                        lines.append(pline(i, "skip", [], types, real_types))

        # torsions
        if "n_torsions" in eex and eex["n_torsions"] > 0:
            lines.append("")
            lines.append("Dihedrals")
            lines.append("")
            for counter, tmp in zip(range(1, eex["n_torsions"] + 1), eex["torsions"]):
                i, j, k, l, index = tmp
                lines.append(
                    "{:6d} {:6d} {:6d} {:6d} {:6d} {:6d}".format(
                        counter, index, i, j, k, l
                    )
                )

            lines.append("")
            lines.append("Dihedral Coeffs")
            lines.append("")

            forms = set([v[0] for v in eex["torsion parameters"]])
            use_hybrid = len(forms) > 1

            for counter, parameters in zip(
                range(1, eex["n_torsion_types"] + 1), eex["torsion parameters"]
            ):
                form, values, types, parameters_type, real_types = parameters
                if form == "torsion_1":
                    Phi0 = values["Phi0"]

                    # Discover form is
                    #  KPhi * [1 + cos(n*Phi - Phi0)]
                    #  with trans = 180
                    #
                    #  For ethane, Phi0 = 0 so at Phi=180 E is min. Correct

                    # Lammps for is
                    #  KPhi * [1 + d*cos(n*Phi)]
                    #  with trans = 180
                    #
                    # Again for ethane, d=+1 and at Phi=180, E is min.
                    #
                    # Phi0 = 0   ==> d=+1
                    # Phi0 = 180 ==> d=-1

                    if float(Phi0) == 0.0:
                        d = "-1"
                    elif float(Phi0) == 180.0:
                        d = "+1"
                    else:
                        raise RuntimeError(
                            "LAMMPS cannot handle Phi0 = {}".format(Phi0)
                        )
                    values["d"] = d
                    values["d"] = values["original Phi0"]

                    function = "harmonic" if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        values,
                        ["Kphi", "d", "n"],
                        types,
                        real_types,
                    )
                elif form == "torsion_3":
                    function = "class2" if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        values,
                        ["V1", "Phi0_1", "V2", "Phi0_2", "V3", "Phi0_3"],
                        types,
                        real_types,
                    )
                elif form == "torsion_opls":
                    function = "opls" if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        values,
                        ["V1", "V2", "V3", "V4"],
                        types,
                        real_types,
                    )
                elif form == "torsion_charmm":
                    function = "charmm" if use_hybrid else ""
                    line = pline(
                        counter,
                        function,
                        values,
                        ["K", "n", "Phi0", "weight"],
                        types,
                        real_types,
                    )
                lines.append(line)

            # middle bond-torsion_3 coefficients, which must match torsions
            # in order & number
            if "n_middle_bond-torsion_3_types" in eex:
                lines.append("")
                lines.append("MiddleBondTorsion Coeffs")
                lines.append("")
                for counter, parameters, torsions in zip(
                    range(1, eex["n_middle_bond-torsion_3_types"] + 1),
                    eex["middle_bond-torsion_3 parameters"],
                    eex["torsion parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    torsion_form = torsions[0]
                    if torsion_form == "torsion_3":
                        function = "class2" if use_hybrid else ""
                        line = pline(
                            counter,
                            function,
                            values,
                            ["V1", "V2", "V3", "R0"],
                            types,
                            real_types,
                        )
                    else:
                        line = pline(
                            counter,
                            "skip",
                            values,
                            [],
                            types,
                            real_types,
                        )
                    lines.append(line)

                # end bond-torsion_3 coefficients, which must match torsions
                # in order & number
                lines.append("")
                lines.append("EndBondTorsion Coeffs")
                lines.append("")
                for counter, parameters, torsions in zip(
                    range(1, eex["n_end_bond-torsion_3_types"] + 1),
                    eex["end_bond-torsion_3 parameters"],
                    eex["torsion parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    torsion_form = torsions[0]
                    if torsion_form == "torsion_3":
                        function = "class2" if use_hybrid else ""
                        line = pline(
                            counter,
                            function,
                            values,
                            [
                                "V1_L",
                                "V2_L",
                                "V3_L",
                                "V1_R",
                                "V2_R",
                                "V3_R",
                                "R0_L",
                                "R0_R",
                            ],
                            types,
                            real_types,
                        )
                        lines.append(line)
                    else:
                        lines.append(pline(i, "skip", [], types, real_types))

                # angle-torsion_3 coefficients, which must match torsions
                # in order & number
                lines.append("")
                lines.append("AngleTorsion Coeffs")
                lines.append("")
                for counter, parameters, torsions in zip(
                    range(1, eex["n_angle-torsion_3_types"] + 1),
                    eex["angle-torsion_3 parameters"],
                    eex["torsion parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    torsion_form = torsions[0]
                    if torsion_form == "torsion_3":
                        function = "class2" if use_hybrid else ""
                        line = pline(
                            counter,
                            function,
                            values,
                            [
                                "V1_L",
                                "V2_L",
                                "V3_L",
                                "V1_R",
                                "V2_R",
                                "V3_R",
                                "Theta0_L",
                                "Theta0_R",
                            ],
                            types,
                            real_types,
                        )
                        lines.append(line)
                    else:
                        lines.append(pline(i, "skip", [], types, real_types))

                # angle-angle-torsion_1 coefficients, which must match torsions
                # in order & number
                lines.append("")
                lines.append("AngleAngleTorsion Coeffs")
                lines.append("")
                for counter, parameters, torsions in zip(
                    range(1, eex["n_angle-angle-torsion_1_types"] + 1),
                    eex["angle-angle-torsion_1 parameters"],
                    eex["torsion parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    torsion_form = torsions[0]
                    if torsion_form == "torsion_3":
                        function = "class2" if use_hybrid else ""
                        line = pline(
                            counter,
                            function,
                            values,
                            ["K", "Theta0_L", "Theta0_R"],
                            types,
                            real_types,
                        )
                        lines.append(line)
                    else:
                        lines.append(pline(i, "skip", [], types, real_types))

                # bond-bond_1_3 coefficients, which must match torsions
                # in order & number
                lines.append("")
                lines.append("BondBond13 Coeffs")
                lines.append("")
                for counter, parameters, torsions in zip(
                    range(1, eex["n_bond-bond_1_3_types"] + 1),
                    eex["bond-bond_1_3 parameters"],
                    eex["torsion parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    torsion_form = torsions[0]
                    if torsion_form == "torsion_3":
                        function = "class2" if use_hybrid else ""
                        line = pline(
                            counter,
                            function,
                            values,
                            ["K", "R10", "R30"],
                            types,
                            real_types,
                        )
                        lines.append(line)
                    else:
                        lines.append(pline(i, "skip", [], types, real_types))

        # out-of-planes
        if "n_oops" in eex and eex["n_oops"] > 0:
            lines.append("")
            lines.append("Impropers")
            lines.append("")

            # Need forms to reorder impropers....
            form = {i: f[0] for i, f in enumerate(eex["oop parameters"], start=1)}

            for counter, tmp in enumerate(eex["oops"], start=1):
                i, j, k, l, index = tmp
                if form[index] in ("improper_opls", "dreiding_out_of_plane"):
                    lines.append(
                        "{:6d} {:6d} {:6d} {:6d} {:6d} {:6d}".format(
                            counter, index, j, i, k, l
                        )
                    )
                else:
                    lines.append(
                        "{:6d} {:6d} {:6d} {:6d} {:6d} {:6d}".format(
                            counter, index, i, j, k, l
                        )
                    )

            lines.append("")
            lines.append("Improper Coeffs")
            lines.append("")
            for counter, parameters in zip(
                range(1, eex["n_oop_types"] + 1), eex["oop parameters"]
            ):
                form, values, types, parameters_type, real_types = parameters
                if form == "wilson_out_of_plane":
                    lines.append(
                        pline(
                            counter,
                            "",
                            values,
                            ["K", "Chi0"],
                            types,
                            real_types,
                        )
                    )
                elif form == "dreiding_out_of_plane":
                    # divide by number of oops from atom (3 for all at the moment!)
                    lines.append(
                        pline(
                            counter,
                            "",
                            values,
                            ["K2", "Psi0"],
                            types,
                            real_types,
                        )
                    )
                elif form == "improper_opls":
                    # divide by two because OPLS uses V2/2 and CVS V.
                    values["V2"] = f"{float(values['V2']) / 2:9.4f}"
                    values["sign"] = "-1"
                    values["original sign"] = "-1"
                    values["n"] = "2"
                    values["original n"] = "2"
                    lines.append(
                        pline(
                            counter,
                            "",
                            values,
                            ["V2", "sign", "n"],
                            types,
                            real_types,
                        )
                    )
                else:
                    raise RuntimeError(f"Can't handle oop form '{form}'")

            # angle-angle
            if "n_angle-angle_types" in eex:
                lines.append("")
                lines.append("AngleAngle Coeffs")
                lines.append("")
                for counter, parameters in zip(
                    range(1, eex["n_angle-angle_types"] + 1),
                    eex["angle-angle parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    line = pline(
                        counter,
                        function,
                        values,
                        [
                            "K1",
                            "K2",
                            "K3",
                            "Theta10",
                            "Theta20",
                            "Theta30",
                        ],
                        types,
                        real_types,
                    )
                    lines.append(line)
        lines.append("")
        lines.append("")
        return (
            "\n".join(lines),
            "\n".join(pair_table),
            "\n".join(bond_table),
            "\n".join(angle_table),
            "\n".join(dihedral_table),
        )

    def record_timing(self, ce, wall, result):
        """Append this run's timing record (``~/.seamm.d/timing/lammps.csv``)
        with :func:`timing_descriptors`; never raises."""
        try:
            _, configuration = self.get_system_configuration()
            text = _output_text(
                self.directory, ("lammps.out", "log.lammps", "stdout.txt", "output.txt")
            )
            descriptors = timing_descriptors(
                getattr(self, "_timing_control", None), text, configuration
            )
            # What the run computed with: the class of its cost
            try:
                descriptors["ff_form"] = str(self.ff_form() or "")
            except Exception:
                descriptors["ff_form"] = ""
            descriptors["model"] = str(getattr(self, "model", None) or "")
            ngpus = ce.get("NGPUS", 0) or 0
            seamm_exec.record_timing(
                "lammps",
                wall,
                descriptors,
                ntasks=ce.get("NTASKS"),
                ngpus=ngpus,
                gpu_model=_gpu_model() if ngpus else None,
                state="finished" if result else "failed",
                in_situ=True,
                **_record_kwargs(),
            )
        except Exception as e:  # pragma: no cover - must never stop the step
            self.logger.warning(f"Could not record the timing of the LAMMPS run: {e}")

    def analyze(self, indent="", nodes=None, values={}, **kwargs):
        """Analyze the output of the calculation"""
        if isinstance(nodes, list) is False:
            nodes = [nodes]

        ret = {node._id[1]: None for node in nodes}

        run_dir = Path(self.directory)

        # Divide the log file into sections for the steps
        log_file = run_dir / "log.lammps"
        if not log_file.is_file():
            raise RuntimeError(f"Log file {log_file} is missing!")
        lines = log_file.read_text().split("\n")
        sections = {}
        section = None
        for line in lines:
            if line.startswith("# Step "):
                section = line[2:]
                sections[section] = []
            elif section is not None:
                sections[section].append(line)
            if line.startswith("Total wall time:"):
                try:
                    h, m, s = line.split()[3].split(":")
                    _time = 3600 * float(h) + 60 * float(m) + float(s)
                    self.results["t_lammps_wall"] = _time
                except Exception as _e:
                    print(f"Wall time exception {_e}")

        for node in nodes:
            for value in node.description:
                printer.important(value)
                printer.important(" ")

            subdir = run_dir / str(node._id[-1])

            # Move any old style files
            id_str = "_".join(str(e) for e in node._id)
            paths = sorted(run_dir.glob(f"*trajectory*{id_str}*.seamm_trj"))
            for path in paths:
                ensemble = str(path).split("_")[1]
                path.rename(subdir / f"{ensemble}_state.trj")

            paths = sorted(run_dir.glob(f"*summary*{id_str}*.seamm_trj"))
            for path in paths:
                ensemble = str(path).split("_")[1]
                path.rename(subdir / f"{ensemble}_summary.trj")

            # Find the state trajectory, if any
            paths = sorted(subdir.glob("*state.trj"))
            if len(paths) > 1:
                raise RuntimeError(f"More than one state.trj file for {id_str}.")
            elif len(paths) == 1:
                control_properties = lambda x: x not in ["tstep"]  # noqa: E731
                node_data, table = self.analyze_trajectory(
                    paths[0], control_properties=control_properties, node=node
                )
                # Save the state data as JSON
                data_file = subdir / "state.json"
                with data_file.open("w") as fd:
                    json.dump(
                        node_data, fd, indent=4, cls=CompactJSONEncoder, sort_keys=True
                    )

                # Get the values from the node data
                for k, v in node_data.items():
                    values[k] = v["mean"]
                    for key in ("stderr", "tau", "inefficiency", "n_samples"):
                        if key in v:
                            values[f"{k},{key}"] = v[key]

                # Get the stress vector is the stresses exist
                stress = []
                tmp = {
                    "stderr": [],
                    "tau": [],
                    "inefficiency": [],
                    "n_samples": [],
                }
                for key in ("Sxx", "Syy", "Szz", "Syz", "Sxz", "Sxy"):
                    if key in values:
                        stress.append(values[key])
                        for item in tmp.keys():
                            if f"{key},{item}" in values:
                                tmp[item].append(values[f"{key},{item}"])
                if len(stress) == 6:
                    values["stress"] = stress
                for item, vector in tmp.items():
                    if len(vector) == 6:
                        values[f"stress,{item}"] = vector
            else:
                node_data = None
                table = None

            values["model"] = self.model
            node.analyze(
                data=values,
                properties=node_data,
                table=table,
                output=sections[node.header],
            )

            ret[node._id[1]] = node_data

        return ret

    def analyze_trajectory(
        self,
        path,
        sampling_rate=20,
        control_properties=None,
        node=None,
    ):
        """Read a trajectory file and do the statistical analysis"""
        results = {}

        table = {
            "Property": [],
            "Value": [],
            " ": [],
            "StdErr": [],
            "Units": [],
            "convergence": [],
            "tau": [],
            "inefficiency": [],
        }

        # Process the trajectory data
        # temporary until we sort out multiple runs
        self._trajectory = []
        with path.open() as fd:
            file_data = pandas.read_csv(
                fd,
                sep=" ",
                header=0,
                comment="!",
                usecols=control_properties,
                index_col="t",
            )
            self._trajectory.append(file_data.iloc[:-1])

        dt = lammps_step.from_lammps_units(
            file_data.index[1] - file_data.index[0], "fs"
        )
        dt_fs = dt.m_as("fs")
        data = pandas.concat(self._trajectory)
        data = data.reset_index(drop=True)
        data.index *= dt_fs

        self.logger.debug("Columns: {}".format(data.columns))
        self.logger.debug("  Types:\n{}".format(data.dtypes))

        printer.normal(f"       Analysis of {path.name}\n")

        # Work out the time step, rather than give the whole vector
        t = data.index
        t_units = "fs"
        len_trj = (len(t) - 1) * dt_fs
        if len_trj >= 4000000000:
            t_units = "ms"
        elif len_trj >= 4000000:
            t_units = "ns"
        elif len_trj >= 4000:
            t_units = "ps"
        t_max = float((len(t) - 1) * dt.m_as(t_units))

        # The user's extra properties: their names, units, and whether the values
        # are in LAMMPS's units (a quantity) or the user's own (explicit units)
        extra_meta = {}
        if node is not None and hasattr(node, "extra_definitions"):
            extra_meta = lammps_step.extras.column_metadata(node.extra_definitions())

        for column in data.columns:
            if "Unnamed:" in column:
                continue
            meta_column = column.rstrip("0123456789")
            convert = True

            if column in extra_meta:
                meta_title, meta_units, convert = extra_meta[column]
            elif meta_column in LAMMPS.display_title:
                meta_title = LAMMPS.display_title[meta_column]
            elif column in LAMMPS.display_title:
                meta_title = LAMMPS.display_title[column]
            else:
                meta_title = f"unknown: {column}"
            if column in extra_meta:
                pass
            elif meta_column in LAMMPS.display_units:
                meta_units = LAMMPS.display_units[meta_column]
            elif column in LAMMPS.display_units:
                meta_units = LAMMPS.display_units[column]
            else:
                meta_units = "???"

            have_warning = False
            have_acf_warning = False
            # Ignore first point, t=0, cause it might not be right.
            yy = data[column].to_numpy()[1:]

            self.logger.info("Analyzing {}, nsamples = {}".format(column, len(yy)))

            # Change from the units LAMMPS is using to the ones that we want
            if convert and meta_units not in ("???", "º", ""):
                try:
                    factor = lammps_step.from_lammps_units(1.0, meta_units)
                    # Do not use inplace operator *= because yy might be int
                    yy = factor.magnitude * yy
                except Exception as e:
                    raise RuntimeError(
                        f"Cannot convert LAMMPS units for column '{column}' to "
                        f"{meta_units}.\n"
                        f"{type(yy)=} {type(factor.magnitude)=}\n"
                        "Please notify support of this error!"
                        f"\n{str(e)}"
                    )

            # compute indices of uncorrelated timeseries using pymbar
            # Their algorithm is quadratic in length of Y unless you
            # use 'nskip'. I set it so there are about 100 time origins, so
            # the convergence time is accurate to about 1%.
            nskip = yy.size // 100 + 1
            conv, inefficiency, Neff_max = timeseries.detect_equilibration(
                yy, nskip=nskip
            )

            self.logger.info(
                "  converged in {} steps, inefficiency = {}, Neff_max = {}".format(
                    conv, inefficiency, Neff_max
                )
            )

            if np.isnan(inefficiency) or np.isnan(Neff_max):
                # Apparently didn't converge!
                table["Property"].append(
                    extra_meta[column][0] if column in extra_meta else column
                )
                table["Value"].append("")
                table[" "].append("")
                table["StdErr"].append("")
                table["Units"].append("")
                table["convergence"].append("unconverged")
                table["tau"].append("")
                table["inefficiency"].append("")

                have_acf = False
                is_converged = False
            else:
                is_converged = True
                tau = dt_fs * (inefficiency - 1) / 2
                if tau < dt_fs / 2:
                    tau = dt_fs / 2
                t0 = conv * dt_fs
                y_t_equil = yy[conv:]
                indices = timeseries.subsample_correlated_data(
                    y_t_equil, g=inefficiency
                )

                if len(indices) == 0:
                    print("Problem with column " + column)
                    print("yy")
                    print(yy)
                    print("y_t_equil")
                    print(y_t_equil)
                    print("indices")
                    print(indices)
                    continue

                y_n = y_t_equil[indices]
                n_samples = len(y_n)
                mean = float(y_n.mean())
                std = float(y_n.std())
                sem = std / sqrt(n_samples)

                # Get the autocorrelation function
                if len(y_t_equil) < 8:
                    have_acf = False
                    have_acf_warning = True
                    acf_warning = "^"
                elif all(y_t_equil == y_t_equil[0]):
                    # A constant value, so no ACF
                    have_acf = False
                    acf_warning = ""
                else:
                    have_acf = True
                    acf_warning = ""
                    nlags = 4 * int(round(inefficiency + 0.5))

                    if nlags > int(len(y_t_equil) / 2):
                        nlags = int(len(y_t_equil) / 2)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        acf, confidence = stattools.acf(
                            y_t_equil,
                            nlags=nlags,
                            alpha=0.05,
                            fft=True,
                            adjusted=False,
                            **_acf_options,
                        )

                results[column] = {}
                results[column]["mean"] = mean
                results[column]["stderr"] = sem
                results[column]["n_sample"] = n_samples
                results[column]["short_production"] = have_acf_warning
                results[column]["tau"] = float(tau)
                results[column]["inefficiency"] = float(inefficiency)
                results[column]["timestep"] = float(dt_fs)
                results[column]["rootname"] = path.stem
                if have_acf:
                    results[column]["acf"] = acf.tolist()
                    results[column]["acf_confidence"] = confidence.tolist()

                # Don't print or graph some properties, like heat flux
                if column in ("Jx", "Jy", "Jz"):
                    results[column]["values"] = data[column].to_numpy()

                # The production part of an extra property, in its display
                # units, for its treatment (lammps_step.extras)
                if column in extra_meta:
                    results[column]["units"] = meta_units
                    if not hasattr(node, "_extra_series"):
                        node._extra_series = {}
                    node._extra_series[column] = y_t_equil

                # Work out units on convergence time
                conv_units = "fs"
                t_conv = t0
                if t0 >= 1000000000:
                    conv_units = "ms"
                    t_conv = t0 / 1000000000
                elif t0 >= 1000000:
                    conv_units = "ns"
                    t_conv = t0 / 1000000
                elif t0 >= 1000:
                    conv_units = "ps"
                    t_conv = t0 / 1000

                # Work out units on autocorrelation time
                tau_units = "fs"
                t_tau = tau
                if tau >= 1000000000:
                    tau_units = "ms"
                    t_tau = tau / 1000000000
                elif tau >= 1000000:
                    tau_units = "ns"
                    t_tau = tau / 1000000
                elif tau >= 1000:
                    tau_units = "ps"
                    t_tau = tau / 1000

                if n_samples < 100:
                    have_warning = True
                    warn = "*"
                else:
                    warn = " "

                results[column]["few_neff"] = have_warning

                label = extra_meta[column][0] if column in extra_meta else column
                table["Property"].append(f"{label}{warn}")
                table["Value"].append(f"{mean:.3f}")
                table[" "].append("±")
                table["StdErr"].append(f"{sem:.3f}")
                table["Units"].append(meta_units)
                table["convergence"].append(f"{t_conv:.2f} {conv_units}")
                table["tau"].append(f"{t_tau:.1f} {tau_units}{acf_warning}")
                table["inefficiency"].append(f"{inefficiency:.1f}")

                if column == "Etot" and self.ff_form() == "reaxff":
                    Eat = self._atomic_energy_sum
                    if Eat != 0.0:
                        results["DfH0_reax"] = {}
                        results["DfH0_reax"]["mean"] = mean + Eat
                        results["DfH0_reax"]["stderr"] = sem
                        results["DfH0_reax"]["n_sample"] = n_samples
                        results["DfH0_reax"]["short_production"] = have_acf_warning
                        results["DfH0_reax"]["tau"] = float(tau)
                        results["DfH0_reax"]["inefficiency"] = float(inefficiency)
                        results["DfH0_reax"]["timestep"] = float(dt_fs)
                        results["DfH0_reax"]["rootname"] = path.stem
                        if have_acf:
                            results["DfH0_reax"]["acf"] = acf.tolist()
                            results["DfH0_reax"]["acf_confidence"] = confidence.tolist()
                        results["DfH0_reax"]["few_neff"] = have_warning

                        table["Property"].append(f"DfH0_reax{warn}")
                        table["Value"].append(f"{mean + Eat:.3f}")
                        table[" "].append("±")
                        table["StdErr"].append(f"{sem:.3f}")
                        table["Units"].append(meta_units)
                        table["convergence"].append(f"{t_conv:.2f} {conv_units}")
                        table["tau"].append(f"{t_tau:.1f} {tau_units}{acf_warning}")
                        table["inefficiency"].append(f"{inefficiency:.1f}")

            # Create graphs of the property
            figure = self.create_figure(
                module_path=(self.__module__.split(".")[0], "seamm"),
                template="line.graph_template",
                fontsize=self.options["graph_fontsize"],
                title=meta_title,
            )

            # The autocorrelation function
            if have_acf:
                plot_acf = figure.add_plot("acf")

                t_acf_units = "fs"
                len_acf = (len(acf) - 1) * dt_fs
                if len_acf >= 2000000000:
                    t_acf_units = "ms"
                elif len_acf >= 2000000:
                    t_acf_units = "ns"
                elif len_acf >= 2000:
                    t_acf_units = "ps"

                x_acf_axis = plot_acf.add_axis(
                    "x", label="Time ({})".format(t_acf_units)
                )
                y_acf_axis = plot_acf.add_axis("y", label="acf", anchor=x_acf_axis)
                x_acf_axis.anchor = y_acf_axis

                # Put the fit to the autocorrelation time in first so the
                # subsequent trajectory trace sits in top
                ts = 0.0
                fit = [1.0]
                for step in range(len(acf) - 1):
                    ts += dt_fs
                    fit.append(exp(-ts / tau))

                plot_acf.add_trace(
                    x_axis=x_acf_axis,
                    y_axis=y_acf_axis,
                    name="fit",
                    x0=0,
                    dx=dt.m_as(t_acf_units),
                    xlabel="t",
                    xunits=t_acf_units,
                    y=fit,
                    ylabel="fit",
                    yunits="",
                    color="gray",
                )

                # the partly transparent error band
                yplus = []
                yminus = []
                t_acf = []
                tmp = 0.0
                for lower, upper in confidence:
                    t_acf.append(tmp)
                    yplus.append(upper)
                    yminus.append(lower)
                    tmp += dt.m_as(t_acf_units)

                plot_acf.add_trace(
                    x_axis=x_acf_axis,
                    y_axis=y_acf_axis,
                    name="stderr",
                    x=t_acf + t_acf[::-1],
                    xlabel="t",
                    xunits=t_acf_units,
                    y=yplus + yminus[::-1],
                    ylabel="stderr",
                    yunits=meta_units,
                    color="rgba(211,211,211,0.5)",
                    fill="toself",
                )

                # And the acf plot last
                plot_acf.add_trace(
                    x_axis=x_acf_axis,
                    y_axis=y_acf_axis,
                    name="acf",
                    x0=0,
                    dx=dt.m_as(t_acf_units),
                    xlabel="t",
                    xunits=t_acf_units,
                    y=list(acf),
                    ylabel="acf",
                    yunits="",
                    color="red",
                )

            # The property data over the trajectory
            y = list(data[column])

            plot = figure.add_plot("trj")

            ylabel = meta_title
            if meta_units != "":
                ylabel += f" ({meta_units})"

            x_axis = plot.add_axis("x", label="Time ({})".format(t_units))
            y_axis = plot.add_axis("y", label=ylabel, anchor=x_axis)
            x_axis.anchor = y_axis

            # Add the trajectory, error band and median value in that order so
            # stack in a nice order.

            # Add the trajectory
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=column,
                x0=0,
                dx=dt.m_as(t_units),
                xlabel="t",
                xunits=t_units,
                y=list(y),
                ylabel=column,
                yunits=meta_units,
                color="#4dbd74",
            )

            if is_converged:
                # the partly transparent error band
                t_min = t0 / dt_fs * dt.m_as(t_units)
                plot.add_trace(
                    x_axis=x_axis,
                    y_axis=y_axis,
                    name="sem",
                    x=[t_min, t_max, t_max, t_min],
                    xlabel="t",
                    xunits=t_units,
                    y=[mean + sem, mean + sem, mean - sem, mean - sem],
                    ylabel="sem",
                    yunits=meta_units,
                    color="rgba(211,211,211,0.5)",
                    fill="toself",
                )

                # and finally the median value so it is on top
                plot.add_trace(
                    x_axis=x_axis,
                    y_axis=y_axis,
                    name="average",
                    x=[t_min, t_max],
                    xlabel="t",
                    xunits=t_units,
                    y=[mean, mean],
                    ylabel="average",
                    yunits=meta_units,
                    color="black",
                )

            if have_acf:
                figure.grid_plots("trj - acf")
            else:
                figure.grid_plots("trj")
            if node is None:
                path = Path(f"{path.stem}_{column}.graph")
            else:
                node_path = Path(node.directory)
                node_path.mkdir(parents=True, exist_ok=True)
                path = node_path / f"{column}.graph"

            figure.write_file(path)

            # Other requested formats
            if "graph_formats" in self.options:
                formats = self.options["graph_formats"]
                # If from seamm.ini, is a single string so parse.
                if isinstance(formats, str):
                    formats = shlex.split(formats)
                for _format in formats:
                    figure.write_file(
                        path.with_suffix("." + _format),
                        width=int(self.options["graph_width"]),
                        height=int(self.options["graph_height"]),
                    )

        # Add citations for pymbar
        self.references.cite(
            raw=self._bibliography["mbar"],
            alias="pymbar-1",
            module="lammps_step",
            level=1,
            note="The main reference for pymbar.",
        )
        self.references.cite(
            raw=self._bibliography["histogram"],
            alias="pymbar-2",
            module="lammps_step",
            level=2,
            note="The second citation for pymbar",
        )
        self.references.cite(
            raw=self._bibliography["convergence"],
            alias="pymbar-3",
            module="lammps_step",
            level=2,
            note="The third citation for pymbar",
        )

        return results, table

    def shake_fix(self, P, eex):
        """Create the 'fix shake' line for rigid bonds and angles and X-H bonds.

        Parameters
        ----------
        P : dict
            The parameters for the initialization step as a dict.
        eex : dict
            The energy expression for this calculation

        Returns
        -------
        line : str
            The correct fix line for LAMMPS
        """

        bond_types = {}
        angle_types = {}

        # The element of each atom, 1-based. eex["elements"] is by atom type.
        type_elements = eex["elements"]
        elements = [""] + [type_elements[atom[3] - 1] for atom in eex["atoms"]]
        bonds = eex["bonds"] if eex.get("n_bonds", 0) > 0 else []

        # Bonds and angles the forcefield makes rigid, e.g. in water models
        for index, parameters in enumerate(eex.get("bond parameters", []), start=1):
            if parameters[0] == "rigid_bond":
                bond_types[index] = 1
        for index, parameters in enumerate(eex.get("angle parameters", []), start=1):
            if parameters[0] == "rigid_angle":
                angle_types[index] = 1

        # Fixing bond lengths of X-H bonds...
        fix_bonds = P["fix_XH_bond_lengths"]
        if fix_bonds == "CH":
            for i, j, index in bonds:
                if {elements[i], elements[j]} == {"C", "H"}:
                    bond_types[index] = 1
        elif fix_bonds == "all":
            for i, j, index in bonds:
                if elements[i] == "H" or elements[j] == "H":
                    bond_types[index] = 1

        # And the result is ....
        if len(bond_types) > 0:
            # SHAKE, not RATTLE: with LAMMPS's fix rattle, rigid SPC/E water had a
            # pressure about 850 atm too high, an anisotropic pressure tensor and a
            # drifting center of mass, whatever its tolerance. SHAKE gives the
            # right pressure and conserves momentum.
            result = "fix                 {} all shake 1.0e-6 50 1000 b "
            for bond_type in bond_types.keys():
                result += " " + str(bond_type)
            if len(angle_types) > 0:
                result += " a "
                for angle_type in angle_types.keys():
                    result += " " + str(angle_type)
        else:
            result = ""

        return result

    def read_dump(self, dumpfile):
        """Read the LAMMPS dumpfile and update the system.

        Parameters
        ----------
        dumpfile : str
            The filename (or path) to the dumpfile.
        """
        self.logger.debug("Reading dump file '{}'".format(dumpfile))

        _, configuration = self.get_system_configuration()
        n_atoms = configuration.n_atoms

        cell = None
        xyz = None
        vxyz = None

        section = ""
        section_lines = []
        with open(dumpfile, "r") as fd:
            lineno = 0
            for line in fd:
                line = line.strip()
                lineno += 1
                if lineno == 1:
                    if line[0:5] != "ITEM:":
                        raise RuntimeError(
                            "Error reading dump file '" + dumpfile + "': The "
                            "first line is incorrect! (" + line + ")"
                        )
                    section = line[6:].strip()
                    section_lines = []
                    self.logger.debug("   section = " + section)
                    continue

                if line[0:5] == "ITEM:":
                    # end a section
                    self.logger.debug("  processing section '{}'".format(section))
                    if "BOX BOUNDS" in section:
                        if len(section.split()) == 8:
                            xlo_bound, xhi_bound, xy = section_lines[0].split()
                            ylo_bound, yhi_bound, xz = section_lines[1].split()
                            zlo, zhi, yz = section_lines[2].split()

                            xlo_bound = float(xlo_bound)
                            xhi_bound = float(xhi_bound)
                            ylo_bound = float(ylo_bound)
                            yhi_bound = float(yhi_bound)
                            zlo = float(zlo)
                            zhi = float(zhi)
                            xy = float(xy)
                            xz = float(xz)
                            yz = float(yz)

                            xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
                            xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
                            ylo = ylo_bound - min(0.0, yz)
                            yhi = yhi_bound - max(0.0, yz)
                            cell = LAMMPS.box_to_cell(
                                xhi - xlo, yhi - ylo, zhi - zlo, xy, xz, yz
                            )
                        else:
                            xlo, xhi = section_lines[0].split()
                            ylo, yhi = section_lines[1].split()
                            zlo, zhi = section_lines[2].split()

                            xlo = float(xlo)
                            xhi = float(xhi)
                            ylo = float(ylo)
                            yhi = float(yhi)
                            zlo = float(zlo)
                            zhi = float(zhi)

                            cell = (xhi - xlo, yhi - ylo, zhi - zlo, 90, 90, 90)
                    elif section == "NUMBER OF ATOMS":
                        if int(section_lines[0]) != n_atoms:
                            raise RuntimeError(
                                "Number of atoms has changed! {} to {}".format(
                                    n_atoms, section_lines[0]
                                )
                            )
                    elif "ATOMS" in section:
                        xyz = []
                        vxyz = []
                        keys = section.split()[1:]
                        if keys[1:4] == ["x", "y", "z"] or keys[1:4] == [
                            "xu",
                            "yu",
                            "zu",
                        ]:
                            fractional = False
                        elif keys[1:4] == ["xs", "ys", "zs"] or keys[1:4] == [
                            "xsu",
                            "ysu",
                            "zsu",
                        ]:
                            fractional = True
                        else:
                            logger.error(f"Can't handle dump file, {keys=}")
                        if len(keys) >= 7 and keys[4:7] == ["vx", "vy", "vz"]:
                            have_velocities = True
                            factor = lammps_step.from_lammps_units(1.0, "fs").magnitude
                            factor = 1 / factor
                            for tmp in section_lines:
                                x, y, z, vx, vy, vz = tmp.split()[1:7]
                                xyz.append([float(x), float(y), float(z)])
                                vxyz.append(
                                    [
                                        factor * float(vx),
                                        factor * float(vy),
                                        factor * float(vz),
                                    ]
                                )
                        else:
                            have_velocities = False
                            for tmp in section_lines:
                                x, y, z = tmp.split()[1:4]
                                xyz.append([float(x), float(y), float(z)])
                    section = line[6:].strip()
                    section_lines = []
                else:
                    section_lines.append(line)

        # Clean up the last section
        xyz = []
        vxyz = []
        if "ATOMS" in section:
            self.logger.debug("  processing section '{}'".format(section))
            self.logger.debug("  handling the atoms")
            keys = section.split()[1:]
            if keys[1:4] == ["x", "y", "z"] or keys[1:4] == ["xu", "yu", "zu"]:
                fractional = False
            elif keys[1:4] == ["xs", "ys", "zs"] or keys[1:4] == ["xsu", "ysu", "zsu"]:
                fractional = True
            else:
                logger.error(f"Can't handle dump file, {keys=}")
            if len(keys) >= 7 and keys[4:7] == ["vx", "vy", "vz"]:
                have_velocities = True
                factor = 1.0 / lammps_step.from_lammps_units(1.0, "fs").magnitude
                for tmp in section_lines:
                    x, y, z, vx, vy, vz = tmp.split()[1:7]
                    xyz.append([float(x), float(y), float(z)])
                    vxyz.append(
                        [factor * float(vx), factor * float(vy), factor * float(vz)]
                    )
            else:
                have_velocities = False
                for tmp in section_lines:
                    x, y, z = tmp.split()[1:4]
                    xyz.append([float(x), float(y), float(z)])

        if not have_velocities:
            vxyz = None

        if configuration.periodicity == 0:
            cell = None

        return xyz, fractional, cell, vxyz

    def _add_lammps_citations(self, text, cite=None):
        """Add the two main citations for LAMMPS, getting the version from stdout
        text.

        Parameters
        ----------
        text : str
            The standard output from LAMMPS

        Returns
        -------
        None
        """
        # Add the JCP paper
        self.references.cite(
            raw=self._bibliography["PLIMPTON19951"],
            alias="lammps-jcp",
            module="lammps_step",
            level=1,
            note="The principle LAMMPS citation.",
        )

        # And the citation to the LAMMPS code itself
        lines = text.splitlines()

        if len(lines) == 0:
            return

        line = lines[0]
        tmp = line.split(" ")
        if len(tmp) != 4:
            self.logger.info(f"Cannot get LAMMPS version: '{line}'")
            return

        month = tmp[2]
        year = tmp[3].rstrip(")")
        version = " ".join(tmp[1:])

        try:
            template = string.Template(self._bibliography["lammps"])

            citation = template.substitute(month=month, version=version, year=year)

            self.references.cite(
                raw=citation,
                alias="lammps-exe",
                module="lammps_step",
                level=1,
                note="The principle citation for the LAMMPS executable.",
            )

        except Exception as e:
            printer.important(f"Exception in citation {type(e)}: {e}")
            printer.important(traceback.format_exc())

        # If there is a log.cite file, process it
        if cite is not None:
            self.logger.debug("log.cite\n" + cite + "\n")
            bibliography = {}
            tmp = bibtexparser.loads(cite).entries_dict
            writer = bibtexparser.bwriter.BibTexWriter()
            for key, data in tmp.items():
                self.logger.info(f"      {key}")
                bibliography[key] = writer._entry_to_bibtex(data)
            self.logger.debug("Bibliography\n" + pprint.pformat(bibliography))

            for entry in bibliography:
                if entry.lower() in ("commment",):
                    continue
                self.references.cite(
                    raw=bibliography[entry],
                    alias=entry,
                    module="lammps_step",
                    level=1,
                    note="LAMMPS citations from log.cite.",
                )

    def angle_table(self, key, values):
        """Create a section of the tabulated angle file.

        Parameters
        ----------
        key : str
            The name for this section in the file.

        values : {str: int, float, or str}
            The dictionary of the constants for this angle term. It looks like
            this::

                {
                    'reference': '5',
                    'Eqn': 'K/8*(1-cos(n*Theta)) + A/(2*Rb*sin(Theta/2))**12',
                    'K': 278.4416826003824,
                    'n': 4,
                    'Rb': 1.606,
                    'A': 11.950286806883364,
                    'zero-shift': 0.001
                }

        Returns
        -------
        [str]
            A list of lines of the tabulated angle file.
        """
        lines = []

        data = {**values}
        if "reference" in data:
            del data["reference"]
        eqn = data.pop("Eqn")
        thetas, Es, dEs = tabulate_angle(eqn, data)

        lines.append(key)
        lines.append(f"N {len(thetas)}")
        lines.append("")
        for i, data in enumerate(zip(thetas, Es, dEs), start=1):
            theta, E, dE = data
            lines.append(f"{i:5d} {theta:9.4f} {E:40}{-dE:40}")
        lines.append("")

        return lines

    def get_dump(self, dumpfile):
        """Read the LAMMPS dumpfile and return the data.

        Parameters
        ----------
        dumpfile : str
            The filename (or path) to the dumpfile.
        """
        self.logger.debug("Reading dump file '{}'".format(dumpfile))

        # system_db = self.get_variable("_system_db")
        # configuration = system_db.system.configuration
        # periodicity = configuration.periodicity
        # n_atoms = configuration.n_atoms

        result = {
            "timestep": [],
            "n_atoms": [],
            "cell": [],
            "data": [],
        }
        section = ""
        section_lines = []
        with open(dumpfile, "r") as fd:
            lineno = 0
            for line in fd:
                line = line.strip()
                lineno += 1
                if lineno == 1:
                    if line[0:5] != "ITEM:":
                        raise RuntimeError(
                            "Error reading dump file '" + dumpfile + "': The "
                            "first line is incorrect! (" + line + ")"
                        )
                    section = line[6:].strip()
                    section_lines = []
                    self.logger.debug("   section = " + section)
                    continue

                if line[0:5] == "ITEM:":
                    # end a section
                    self.logger.debug("  processing section '{}'".format(section))
                    if section == "TIMESTEP":
                        result["timestep"].append(int(section_lines[0]))
                    elif section == "NUMBER OF ATOMS":
                        n_atoms = int(section_lines[0])
                        result["n_atoms"].append(n_atoms)
                    elif "BOX BOUNDS" in section:
                        if len(section.split()) == 8:
                            xlo_bound, xhi_bound, xy = section_lines[0].split()
                            ylo_bound, yhi_bound, xz = section_lines[1].split()
                            zlo, zhi, yz = section_lines[2].split()

                            xlo_bound = float(xlo_bound)
                            xhi_bound = float(xhi_bound)
                            ylo_bound = float(ylo_bound)
                            yhi_bound = float(yhi_bound)
                            zlo = float(zlo)
                            zhi = float(zhi)
                            xy = float(xy)
                            xz = float(xz)
                            yz = float(yz)

                            xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
                            xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
                            ylo = ylo_bound - min(0.0, yz)
                            yhi = yhi_bound - max(0.0, yz)
                            cell = LAMMPS.box_to_cell(
                                xhi - xlo, yhi - ylo, zhi - zlo, xy, xz, yz
                            )
                        else:
                            xlo, xhi = section_lines[0].split()
                            ylo, yhi = section_lines[1].split()
                            zlo, zhi = section_lines[2].split()

                            xlo = float(xlo)
                            xhi = float(xhi)
                            ylo = float(ylo)
                            yhi = float(yhi)
                            zlo = float(zlo)
                            zhi = float(zhi)

                            cell = (xhi - xlo, yhi - ylo, zhi - zlo, 90, 90, 90)
                        result["cell"].append(cell)
                    elif "ATOMS" in section:
                        keys = section.split()[2:]
                        if "fields" not in result:
                            result["fields"] = keys
                        elif keys != result["fields"]:
                            raise ValueError(
                                f"Error reading dump file '{dumpfile}': The fields in "
                                "the ATOMS section do not match the previous fields!\n"
                                f"  timestep: {result['timestep'][-1]}\n"
                                f"  previous: {result['fields']=}\n"
                                f"   current: {keys=}"
                            )
                        tmp = [[None] * n_atoms for key in keys]
                        for line in section_lines:
                            values = line.split()
                            _id = int(values[0]) - 1
                            for index, value in enumerate(values[1:]):
                                key = keys[index]
                                if key in ("mol", "proc", "procp1", "type"):
                                    value = int(value)
                                elif key in ("element",):
                                    value = value.strip()
                                else:
                                    value = float(value)
                                tmp[index][_id] = value
                        result["data"].append(tmp)
                    section = line[6:].strip()
                    section_lines = []
                else:
                    section_lines.append(line)

        # Clean up the last section
        if "ATOMS" in section:
            keys = section.split()[2:]
            if "fields" not in result:
                result["fields"] = keys
            elif keys != result["fields"]:
                raise ValueError(
                    f"Error reading dump file '{dumpfile}': The fields in "
                    "the ATOMS section do not match the previous fields!\n"
                    f"  timestep: {result['timestep'][-1]}\n"
                    f"  previous: {result['fields']=}\n"
                    f"   current: {keys=}"
                )
            tmp = [[None] * n_atoms for key in keys]
            for line in section_lines:
                values = line.split()
                _id = int(values[0]) - 1
                for index, value in enumerate(values[1:]):
                    key = keys[index]
                    if key in ("mol", "proc", "procp1", "type"):
                        value = int(value)
                    elif key in ("element",):
                        value = value.strip()
                    else:
                        value = float(value)
                    tmp[index][_id] = value
            result["data"].append(tmp)

        return result
