#!/bin/bash
#MolSSI lammps_step:cpu_bind 2.0
# cpu_bind.sh — bind a CPU-only LAMMPS run to SEAMM_NP cores
#
# CPU binding is worked out from the machine, never assumed:
#   * Under a batch scheduler (SLURM, PBS, LSF) the job already has its own cores,
#     so nothing is bound.
#   * Otherwise the run uses SEAMM_NP of the CPUs this process may use, preferring
#     CPUs that no GPU lists as its own (per `nvidia-smi topo -m`), so GPU jobs on
#     the same machine keep theirs.
#   * If that cannot be determined, or numactl is missing, nothing is bound.
#
# Environment:
#   SEAMM_NP       Number of CPUs to use (required)
#   SEAMM_BIND     "none" to never bind CPUs
#   SEAMM_CPUS     CPUs to use, e.g. "8-31" (overrides the choice)
#   SEAMM_DEBUG    1 to report the binding

set -uo pipefail

if [ -z "${SEAMM_NP:-}" ]; then
    echo "Error: SEAMM_NP is not set" >&2
    exit 1
fi

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

under_scheduler() {
    [ -n "${SLURM_JOB_ID:-}" ] || [ -n "${PBS_JOBID:-}" ] || [ -n "${LSB_JOBID:-}" ]
}

# The first CPU of each GPU's nearby set, which the GPU jobs would use first.
gpu_first_cpus() {
    command -v nvidia-smi > /dev/null 2>&1 || return 0
    nvidia-smi topo -m 2>/dev/null | sed 's/\x1b\[[0-9;]*m//g' | awk '
        $1 ~ /^GPU[0-9]+$/ {
            for (i = 2; i <= NF; i++)
                if ($i ~ /^[0-9]+([-,][0-9]+)+$/ || ($i ~ /^[0-9]+$/ && i > 2)) {
                    print $i; break
                }
        }' | sort -u
}

choose_cpus() {
    if [ -n "${SEAMM_CPUS:-}" ]; then echo "$SEAMM_CPUS"; return; fi
    [ "${SEAMM_BIND:-}" = none ] && return
    under_scheduler && return
    command -v numactl > /dev/null 2>&1 || return

    local allowed=() gpu=" " c preferred=() others=()
    read -ra allowed <<< "$(expand_cpus "$(awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status 2>/dev/null)")"
    [ ${#allowed[@]} -eq 0 ] && return
    [ "$SEAMM_NP" -gt ${#allowed[@]} ] && return   # more ranks than CPUs: let MPI decide
    # Keep clear of the first CPUs of each GPU's set (where GPU jobs bind first).
    local set list
    while read -r set; do
        [ -z "$set" ] && continue
        read -ra list <<< "$(expand_cpus "$set")"
        for c in "${list[@]:0:8}"; do gpu+="$c "; done
    done <<< "$(gpu_first_cpus)"
    for c in "${allowed[@]}"; do
        if [[ $gpu == *" $c "* ]]; then others+=("$c"); else preferred+=("$c"); fi
    done
    local all=("${preferred[@]}" "${others[@]}")
    local chosen=("${all[@]:0:$SEAMM_NP}")
    local IFS=,
    echo "${chosen[*]}"
}

CPU_BIND="$(choose_cpus)"

if [ "${SEAMM_DEBUG:-0}" = "1" ]; then
    echo "SEAMM_NP=$SEAMM_NP -> CPUs ${CPU_BIND:-not bound}" >&2
fi

if [ -n "$CPU_BIND" ]; then
    exec numactl --physcpubind="$CPU_BIND" "$@"
else
    exec "$@"
fi
