#!/bin/bash
#MolSSI lammps_step:mdi_bind 2.2
# mdi_bind.sh — Resource binding for an MDI engine (e.g. an ML potential) + LAMMPS
#
# Runs one rank of an MPMD pair: the engine (rank 0) on the selected GPU, and the
# LAMMPS driver (rank 1) with no GPU. Also starts an nvidia-smi monitor for the
# engine's GPU.
#
# CPU binding is worked out from the machine itself, never assumed:
#   * Under a batch scheduler (SLURM, PBS, LSF) the job already has its own cores
#     and GPU, so nothing is bound.
#   * Otherwise the engine and driver share the CPUs next to the GPU (from
#     `nvidia-smi topo -m`), limited to the CPUs this process may use. GPUs that
#     share those CPUs split them evenly, and the engine and driver take halves.
#   * If that cannot be determined, or taskset is missing, nothing is bound.
#
# Environment:
#   SEAMM_GPUS          Comma-separated list of GPU IDs (default: "0")
#   SEAMM_MEMORY_LOG    Path for the GPU monitoring log (optional)
#   SEAMM_BIND          "none" to never bind CPUs
#   SEAMM_ENGINE_CPUS   CPUs for the engine, e.g. "0-7" (overrides the choice)
#   SEAMM_DRIVER_CPUS   CPUs for the LAMMPS driver, e.g. "8-15"
#   SEAMM_ENGINE_GRACE  Seconds the engine may take to exit after the driver
#                       finishes normally before it is stopped (default 120)
#   SEAMM_ERROR_GRACE   Seconds LAMMPS may take to exit after printing an ERROR
#                       before it is stopped (default 30)
#
# Each side stops the other when it fails: a LAMMPS error otherwise leaves the
# engine waiting for its driver forever, holding the GPU while the job sits in
# the queue as RUNNING. LAMMPS itself does not exit after an error in an MDI pair
# -- it waits in MPI_Finalize for the engine -- so a copy of the driver's stdout
# and stderr is watched for "ERROR" lines (the output itself passes through tee
# unchanged), and LAMMPS is stopped if it has not exited a little later. The two
# record their programs' PIDs in files named for the pair (the MPI launch's
# namespace, which each mpirun has its own of, and the GPU), so pairs sharing a
# directory keep apart.
#
# Usage:
#   mpirun --bind-to none \
#     -np 1 mdi_bind.sh xnn mdi --ckpt model.pt -mdi "-role ENGINE ... -method MPI" \
#     : -np 1 mdi_bind.sh lmp -mdi "-role DRIVER ... -method MPI" -in input.dat
#
# Or for multi-GPU (simultaneous independent runs):
#   SEAMM_GPUS=0 mpirun ... -np 1 mdi_bind.sh <engine> ... : -np 1 mdi_bind.sh lmp ...  &
#   SEAMM_GPUS=1 mpirun ... -np 1 mdi_bind.sh <engine> ... : -np 1 mdi_bind.sh lmp ...  &

set -uo pipefail

SEAMM_GPUS="${SEAMM_GPUS:-0}"
# The rank within this pair: 0 is the engine, 1 the driver. Each MPI launcher names
# it differently (OpenMPI, MPICH/Hydra, PMIx, SLURM's srun).
LOCAL_RANK="${OMPI_COMM_WORLD_LOCAL_RANK:-${MPI_LOCALRANKID:-${PMIX_RANK:-${PMI_RANK:-${SLURM_LOCALID:-0}}}}}"

IFS=',' read -ra GPU_ARRAY <<< "$SEAMM_GPUS"
GPU_ID="${GPU_ARRAY[0]}"  # Use first GPU for this simulation

# The physical index of that GPU, which is not the same number.  SEAMM_GPUS
# counts within the allocation, so a job given one GPU always sees 0 whichever
# card that is.  nvidia-smi uses the machine's own numbering, so it needs
# translating -- otherwise every concurrent job picks the CPUs belonging to GPU 0
# and they land on top of each other.
if [ -n "${CUDA_VISIBLE_DEVICES+set}" ] && [ -n "$CUDA_VISIBLE_DEVICES" ]; then
    PHYSICAL_GPU=$(echo "$CUDA_VISIBLE_DEVICES" | cut -d, -f$((GPU_ID + 1)))
    [ -z "$PHYSICAL_GPU" ] && PHYSICAL_GPU="$GPU_ID"
else
    PHYSICAL_GPU="$GPU_ID"
fi

# ---------------------------------------------------------------------------
# Choosing CPUs
# ---------------------------------------------------------------------------

# "0-3,8,10-11" -> "0 1 2 3 8 10 11"
expand_cpus() {
    local part a b i out=()
    IFS=',' read -ra parts <<< "$1"
    for part in "${parts[@]}"; do
        if [[ $part =~ ^([0-9]+)-([0-9]+)$ ]]; then
            a=${BASH_REMATCH[1]}; b=${BASH_REMATCH[2]}
            for (( i=a; i<=b; i++ )); do out+=("$i"); done
        elif [[ $part =~ ^[0-9]+$ ]]; then
            out+=("$part")
        fi
    done
    echo "${out[*]:-}"
}

# The CPUs this process may run on (Linux), or nothing.
allowed_cpus() {
    awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status 2>/dev/null
}

# The CPU affinity nvidia-smi reports for a physical GPU, e.g. "0-63", or nothing.
gpu_affinity() {
    command -v nvidia-smi > /dev/null 2>&1 || return 0
    nvidia-smi topo -m 2>/dev/null | sed 's/\x1b\[[0-9;]*m//g' | awk -v g="GPU$1" '
        $1 == g {
            for (i = 2; i <= NF; i++)
                if ($i ~ /^[0-9]+([-,][0-9]+)+$/ || ($i ~ /^[0-9]+$/ && i > 2)) {
                    print $i; exit
                }
        }'
}

# The physical GPUs whose CPU affinity is the same as this one's, in order.
gpus_sharing() {
    local mine="$1" g
    command -v nvidia-smi > /dev/null 2>&1 || return 0
    for g in $(nvidia-smi --query-gpu=index --format=csv,noheader 2>/dev/null); do
        [ "$(gpu_affinity "$g")" = "$mine" ] && echo "$g"
    done
}

under_scheduler() {
    [ -n "${SLURM_JOB_ID:-}" ] || [ -n "${PBS_JOBID:-}" ] || [ -n "${LSB_JOBID:-}" ]
}

# The CPUs for "engine" or "driver", as a taskset list, or nothing for no binding.
choose_cpus() {
    local role="$1"
    if [ "$role" = engine ] && [ -n "${SEAMM_ENGINE_CPUS:-}" ]; then
        echo "$SEAMM_ENGINE_CPUS"; return
    fi
    if [ "$role" = driver ] && [ -n "${SEAMM_DRIVER_CPUS:-}" ]; then
        echo "$SEAMM_DRIVER_CPUS"; return
    fi
    [ "${SEAMM_BIND:-}" = none ] && return
    under_scheduler && return
    command -v taskset > /dev/null 2>&1 || return

    local allowed affinity cpus=() c
    allowed=" $(expand_cpus "$(allowed_cpus)") "
    [ "$allowed" = "  " ] && return
    affinity="$(gpu_affinity "$PHYSICAL_GPU")"
    for c in $(expand_cpus "$affinity"); do
        [[ $allowed == *" $c "* ]] && cpus+=("$c")
    done
    [ ${#cpus[@]} -eq 0 ] && read -ra cpus <<< "$allowed"

    # GPUs that share these CPUs split them; this GPU takes its own share.
    local peers=() k=0 m=1 i
    if [ -n "$affinity" ]; then
        read -ra peers <<< "$(gpus_sharing "$affinity" | tr '\n' ' ')"
        m=${#peers[@]}
        [ "$m" -lt 1 ] && m=1
        for i in "${!peers[@]}"; do
            [ "${peers[$i]}" = "$PHYSICAL_GPU" ] && k=$i
        done
    fi
    local n=${#cpus[@]} share start
    share=$(( n / m )); [ "$share" -lt 1 ] && share=1
    start=$(( (k * share) % n ))
    local mine=("${cpus[@]:$start:$share}")

    # The engine takes the first half of the share, the driver the rest.
    local half=$(( (${#mine[@]} + 1) / 2 )) chosen
    if [ "$role" = engine ] || [ ${#mine[@]} -lt 2 ]; then
        chosen=("${mine[@]:0:$half}")
    else
        chosen=("${mine[@]:$half}")
    fi
    local IFS=,
    echo "${chosen[*]}"
}

# Become the command, bound to the CPUs if any were chosen. Run in a background
# subshell, so that the PID recorded for it is the program's own: stopping it then
# stops the program, not just a shell around it.
exec_bound() {
    local cpus="$1"; shift
    if [ -n "$cpus" ]; then
        exec taskset -c "$cpus" "$@"
    else
        exec "$@"
    fi
}

# Try to stop the codes from spinning instead of sleeping
export OMPI_MCA_mpi_yield_when_idle=1
export OMPI_MCA_mpi_wait_mode=1

# ---------------------------------------------------------------------------
# Each side's program runs as a child, so that the other side can stop it
# ---------------------------------------------------------------------------
PAIR_KEY="$(printf '%s' "${PMIX_NAMESPACE:-${OMPI_MCA_ess_base_jobid:-${SLURM_JOB_ID:-local}.${SLURM_STEP_ID:-0}}}" | tr -c 'A-Za-z0-9._-' '_')_gpu${GPU_ID}"
ENGINE_PID_FILE=".mdi_bind.${PAIR_KEY}.engine.pid"
DRIVER_PID_FILE=".mdi_bind.${PAIR_KEY}.driver.pid"
CHILD_PID=""

# Pass a scancel (or mpirun's TERM) on to the program
trap '[ -n "$CHILD_PID" ] && kill -TERM "$CHILD_PID" 2>/dev/null' TERM INT

# Stop the other side's program, if it is still running: TERM, then KILL.
stop_partner() {
    local pid_file="$1" grace="$2" pid i
    [ -f "$pid_file" ] || return 0
    pid="$(cat "$pid_file" 2>/dev/null)"
    [ -n "$pid" ] || return 0
    for (( i=0; i<grace; i++ )); do
        kill -0 "$pid" 2>/dev/null || return 0
        sleep 1
    done
    echo "Stopping the partner process $pid." >&2
    kill -TERM "$pid" 2>/dev/null || return 0
    for (( i=0; i<10; i++ )); do
        kill -0 "$pid" 2>/dev/null || return 0
        sleep 1
    done
    kill -KILL "$pid" 2>/dev/null || true
}

# Read a copy of the driver's output; if LAMMPS reports an error but does not
# exit (it waits in MPI_Finalize for the engine), stop it after a grace period.
# It reads to the end whatever it finds, so that tee never writes to a closed pipe.
scan_for_errors() {
    local line stopping=""
    while IFS= read -r line || [ -n "$line" ]; do
        if [ -z "$stopping" ] && [[ $line == ERROR* ]]; then
            stopping=yes
            (
                sleep "${SEAMM_ERROR_GRACE:-30}"
                pid="$(cat "$DRIVER_PID_FILE" 2>/dev/null)"
                if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then
                    echo "LAMMPS reported an error and has not exited; stopping it." >&2
                    kill -TERM "$pid" 2>/dev/null
                    sleep 10
                    kill -KILL "$pid" 2>/dev/null || true
                fi
            ) &
        fi
    done
}

# Run the program in the background, record its PID, and wait for it. With
# "watch" as the first argument, copies of the program's stdout and stderr are
# scanned for errors by scan_for_errors.
run_child() {
    local watch=""
    if [ "$1" = watch ]; then
        watch=yes; shift
    fi
    local pid_file="$1" cpus="$2"; shift 2
    if [ -n "$watch" ]; then
        # The output goes on unchanged (tee copies bytes); the scanners read a copy.
        # Writing to a pipe rather than mpirun's terminal, the program would buffer
        # its output in blocks -- and a LAMMPS hung after an error would never
        # write its ERROR line at all -- so it is made line-buffered, as on a
        # terminal, where stdbuf exists.
        local line_buffered=()
        command -v stdbuf > /dev/null 2>&1 && line_buffered=(stdbuf -oL -eL)
        exec_bound "$cpus" ${line_buffered[@]+"${line_buffered[@]}"} "$@" \
            > >(tee >(scan_for_errors)) \
            2> >(tee >(scan_for_errors) >&2) &
    else
        exec_bound "$cpus" "$@" &
    fi
    CHILD_PID=$!
    echo "$CHILD_PID" > "$pid_file"
    local status=0
    wait "$CHILD_PID" || status=$?
    # A trapped signal interrupts wait: wait again for the program itself
    while kill -0 "$CHILD_PID" 2>/dev/null; do
        wait "$CHILD_PID" || status=$?
    done
    rm -f "$pid_file"
    return $status
}

# ---------------------------------------------------------------------------
# Rank 0 = engine (gets the GPU + monitoring)
# Rank 1 = LAMMPS driver (no GPU)
# ---------------------------------------------------------------------------
if [ "$LOCAL_RANK" -eq 0 ]; then
    # ---- Engine process ----
    CPU_BIND="$(choose_cpus engine)"

    # Leave an allocation the scheduler already made alone.  CUDA_VISIBLE_DEVICES
    # is not composable: setting it again is interpreted against the machine's
    # full set of devices, not against the allocation we were given, so
    # re-exporting an index here is how a job ends up on a GPU it does not hold.
    # When SLURM (or another scheduler) has set it, SEAMM_GPUS already counts
    # within that allocation and there is nothing to add.
    if [ -z "${CUDA_VISIBLE_DEVICES+set}" ]; then
        export CUDA_VISIBLE_DEVICES="$GPU_ID"
    fi
    MONITOR_GPU="$PHYSICAL_GPU"   # nvidia-smi indexes physically

    export OMP_NUM_THREADS=1
    export TORCH_NUM_THREADS=4
    export MKL_NUM_THREADS=4

    echo "Engine (rank $LOCAL_RANK) -> GPU $GPU_ID (physical $PHYSICAL_GPU), CPUs ${CPU_BIND:-not bound}" >&2

    # Start GPU memory/utilization monitor
    MEMORY_LOG="${SEAMM_MEMORY_LOG:-./gpu_${MONITOR_GPU}_engine.log}"
    MONITOR_PID=""
    if command -v nvidia-smi > /dev/null 2>&1; then
        nvidia-smi --query-gpu=timestamp,memory.used,memory.free,utilization.gpu \
                   --format=csv -l 5 -i "$MONITOR_GPU" > "$MEMORY_LOG" &
        MONITOR_PID=$!
        echo "GPU monitor PID = $MONITOR_PID" >&2
        trap 'kill $MONITOR_PID 2>/dev/null || true; wait $MONITOR_PID 2>/dev/null || true' EXIT
    fi

    echo "$@" > engine.cmd

    # Run the engine.  Keep its exit status: the cleanup below, and the EXIT
    # trap, kill the nvidia-smi monitor on purpose, so those kill/wait calls
    # necessarily report failure.  Left unguarded their status becomes this
    # script's, which makes mpirun treat every successful engine run as a
    # failure (rank 0 exited 1 even when the engine finished cleanly).
    rc=0
    run_child "$ENGINE_PID_FILE" "$CPU_BIND" "$@" || rc=$?
    # A failed engine leaves the driver waiting for it
    [ "$rc" -ne 0 ] && stop_partner "$DRIVER_PID_FILE" 0

    # Clean up monitor
    if [ -n "$MONITOR_PID" ]; then
        kill "$MONITOR_PID" 2>/dev/null || true
        wait "$MONITOR_PID" 2>/dev/null || true
        echo "Done!" >> "$MEMORY_LOG"
    fi
    echo "Engine finished." >&2
    exit $rc

else
    # ---- Driver process (LAMMPS) ----
    CPU_BIND="$(choose_cpus driver)"

    # LAMMPS doesn't need GPU access
    export CUDA_VISIBLE_DEVICES=""
    export OMP_NUM_THREADS=1

    echo "Driver (rank $LOCAL_RANK) -> no GPU (GPU $PHYSICAL_GPU's group), CPUs ${CPU_BIND:-not bound}" >&2

    echo "$@" > driver.cmd

    # Run LAMMPS
    rc=0
    run_child watch "$DRIVER_PID_FILE" "$CPU_BIND" "$@" || rc=$?
    # A failed driver leaves the engine waiting for it, holding the GPU. After a
    # normal finish the engine has been told to exit; give it time to.
    if [ "$rc" -ne 0 ]; then
        stop_partner "$ENGINE_PID_FILE" 0
    else
        stop_partner "$ENGINE_PID_FILE" "${SEAMM_ENGINE_GRACE:-120}"
    fi

    echo "Driver finished." >&2
    exit $rc
fi
