#!/bin/bash
#MolSSI lammps_step:gpu_bind 2.0
# gpu_bind.sh — one LAMMPS rank per GPU (e.g. Kokkos), with a memory monitor
#
# Rank i runs on GPU i of SEAMM_GPUS. CPU binding is worked out from the machine,
# never assumed:
#   * Under a batch scheduler (SLURM, PBS, LSF) the job already has its own cores
#     and GPUs, so nothing is bound.
#   * Otherwise each rank uses the CPUs next to its GPU (from `nvidia-smi topo -m`),
#     limited to the CPUs this process may use; GPUs sharing those CPUs split them.
#   * If that cannot be determined, or taskset is missing, nothing is bound.
#
# Environment:
#   SEAMM_GPUS          Comma-separated list of GPU IDs (required)
#   SEAMM_MEMORY_LOG    Path for the GPU monitoring log (optional)
#   SEAMM_BIND          "none" to never bind CPUs

set -uo pipefail

if [ -z "${SEAMM_GPUS:-}" ]; then
    echo "Error: SEAMM_GPUS is not set" >&2
    exit 1
fi

IFS=',' read -ra GPU_ARRAY <<< "$SEAMM_GPUS"
N_GPUS=${#GPU_ARRAY[@]}
LOCAL_RANK="${OMPI_COMM_WORLD_LOCAL_RANK:-${MPI_LOCALRANKID:-${PMIX_RANK:-${PMI_RANK:-${SLURM_LOCALID:-0}}}}}"

if [ "$LOCAL_RANK" -ge "$N_GPUS" ]; then
    echo "Error: local rank $LOCAL_RANK exceeds number of GPUs in SEAMM_GPUS ($N_GPUS)" >&2
    exit 1
fi
GPU_ID=${GPU_ARRAY[$LOCAL_RANK]}

# The machine's own number for that GPU (SEAMM_GPUS counts within the allocation).
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then
    PHYSICAL_GPU=$(echo "$CUDA_VISIBLE_DEVICES" | cut -d, -f$((GPU_ID + 1)))
    [ -z "$PHYSICAL_GPU" ] && PHYSICAL_GPU="$GPU_ID"
else
    PHYSICAL_GPU="$GPU_ID"
fi

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

allowed_cpus() {
    awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status 2>/dev/null
}

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

# This GPU's share of its nearby CPUs, as a taskset list, or nothing.
choose_cpus() {
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
    local chosen=("${cpus[@]:$start:$share}")
    local IFS=,
    echo "${chosen[*]}"
}

CPU_BIND="$(choose_cpus)"

# Leave an allocation the scheduler already made alone: CUDA_VISIBLE_DEVICES is
# interpreted against the machine's full set of devices, so re-exporting an
# index inside an allocation can put the rank on a GPU the job does not hold.
if [ -z "${CUDA_VISIBLE_DEVICES+set}" ]; then
    export CUDA_VISIBLE_DEVICES=$GPU_ID
fi
export OMP_NUM_THREADS=1
export TORCH_NUM_THREADS=4
export MKL_NUM_THREADS=4

echo "Rank $LOCAL_RANK -> GPU $GPU_ID (physical $PHYSICAL_GPU), CPUs ${CPU_BIND:-not bound}" >&2

# Start background memory monitor
MEMORY_LOG="${SEAMM_MEMORY_LOG:-./gpu_${PHYSICAL_GPU}_rank_${LOCAL_RANK}.log}"
MONITOR_PID=""
if command -v nvidia-smi > /dev/null 2>&1; then
    nvidia-smi --query-gpu=timestamp,memory.used,memory.free,utilization.gpu \
               --format=csv -l 5 -i "$PHYSICAL_GPU" > "$MEMORY_LOG" &
    MONITOR_PID=$!
    echo "The monitor PID = $MONITOR_PID" >&2
    trap 'kill $MONITOR_PID 2>/dev/null || true; wait $MONITOR_PID 2>/dev/null || true' EXIT
fi

# Run LAMMPS, keeping its exit status (killing the monitor necessarily "fails").
rc=0
if [ -n "$CPU_BIND" ]; then
    taskset -c "$CPU_BIND" "$@" || rc=$?
else
    "$@" || rc=$?
fi

if [ -n "$MONITOR_PID" ]; then
    kill "$MONITOR_PID" 2>/dev/null || true
    wait "$MONITOR_PID" 2>/dev/null || true
    echo "Done!" >> "$MEMORY_LOG"
fi
echo "The LAMMPS run has finished." >&2
exit $rc
