#!/bin/bash
#MolSSI lammps_step:mdi_monitor 1.0
# mdi_monitor.sh — Lightweight wrapper for SLURM/PBS managed environments
#
# Assumes the scheduler and MPI launcher handle GPU/CPU binding.
# This script only starts nvidia-smi monitoring for the engine rank.
#
# Environment (set by scheduler):
#   CUDA_VISIBLE_DEVICES   Set automatically by SLURM --gpu-bind or PBS
#   OMPI_COMM_WORLD_LOCAL_RANK   Set by OpenMPI
#
# Optional:
#   SEAMM_MEMORY_LOG       Custom path for GPU monitoring log
#
# Usage (SLURM example):
#   srun -n 1 mdi_monitor.sh python mace_mdi.py -mdi "..." \
#     : -n 1 mdi_monitor.sh lmp -mdi "..." -in input.dat

set -euo pipefail

LOCAL_RANK="${OMPI_COMM_WORLD_LOCAL_RANK:-${SLURM_LOCALID:-${PMI_LOCAL_RANK:-${MPI_LOCALRANKID:-0}}}}"

if [ "$LOCAL_RANK" -eq 0 ] && [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then
    # ---- Engine rank: start GPU monitor ----
    GPU_ID="${CUDA_VISIBLE_DEVICES%%,*}"  # First visible GPU
    MEMORY_LOG="${SEAMM_MEMORY_LOG:-./gpu_${GPU_ID}_engine.log}"

    echo "Engine (rank $LOCAL_RANK) -> GPU $CUDA_VISIBLE_DEVICES" >&2

    nvidia-smi --query-gpu=timestamp,memory.used,memory.free,utilization.gpu \
               --format=csv -l 5 -i "$GPU_ID" > "$MEMORY_LOG" 2>/dev/null &
    MONITOR_PID=$!
    trap "kill $MONITOR_PID 2>/dev/null; wait $MONITOR_PID 2>/dev/null" EXIT

    # Run the engine
    "$@"

    kill "$MONITOR_PID" 2>/dev/null
    wait "$MONITOR_PID" 2>/dev/null
    echo "Done!" >> "$MEMORY_LOG"
    echo "Engine finished." >&2
else
    # ---- Driver rank or no GPU: just run the command ----
    echo "Driver (rank $LOCAL_RANK)" >&2
    "$@"
    echo "Driver finished." >&2
fi
