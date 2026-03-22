#!/bin/bash
#MolSSI lammps_step:mdi_bind 1.0
# mdi_bind.sh — Resource binding for MACE MDI engine + LAMMPS driver
#
# Binds the MACE engine (rank 0) to the selected GPU and its NUMA-local CPUs,
# and the LAMMPS driver (rank 1) to adjacent CPUs with no GPU.
# Also starts an nvidia-smi monitor for the engine's GPU.
#
# Environment:
#   SEAMM_GPUS          Comma-separated list of GPU IDs (default: "0")
#   SEAMM_MEMORY_LOG    Path for GPU monitoring log (optional)
#
# Usage:
#   SEAMM_GPUS=0 mpirun --mca mpi_yield_when_idle 1 \
#     -np 1 mdi_bind.sh python mace_mdi.py -mdi "..." \
#     : -np 1 mdi_bind.sh lmp -mdi "..." -in input.dat
#
# Or for multi-GPU (simultaneous independent runs):
#   SEAMM_GPUS=0 mpirun ... -np 1 mdi_bind.sh python ... : -np 1 mdi_bind.sh lmp ...  &
#   SEAMM_GPUS=1 mpirun ... -np 1 mdi_bind.sh python ... : -np 1 mdi_bind.sh lmp ...  &

set -euo pipefail

SEAMM_GPUS="${SEAMM_GPUS:-0}"
LOCAL_RANK="${OMPI_COMM_WORLD_LOCAL_RANK:-0}"

IFS=',' read -ra GPU_ARRAY <<< "$SEAMM_GPUS"
GPU_ID="${GPU_ARRAY[0]}"  # Use first GPU for this simulation

# ---------------------------------------------------------------------------
# Map GPU index to L3 cache groups for EPYC 7763
# Each GPU gets two groups: one for the engine, one for the LAMMPS driver
# ---------------------------------------------------------------------------
declare -A GPU_TO_ENGINE_CPU
declare -A GPU_TO_DRIVER_CPU

# GPU 0: engine on L3 group 0 (cores 0-7), driver on L3 group 1 (cores 8-15)
GPU_TO_ENGINE_CPU[0]="0-7"
GPU_TO_DRIVER_CPU[0]="8-15"

# GPU 1: engine on L3 group 4 (cores 32-39), driver on L3 group 5 (cores 40-47)
GPU_TO_ENGINE_CPU[1]="32-39"
GPU_TO_DRIVER_CPU[1]="40-47"

# Try to stop the codes from spinning instead of sleeping
export OMPI_MCA_mpi_yield_when_idle=1
export OMPI_MCA_mpi_wait_mode=1

# ---------------------------------------------------------------------------
# Rank 0 = MACE engine (gets GPU + NUMA-local CPUs + monitoring)
# Rank 1 = LAMMPS driver (gets adjacent CPUs, no GPU)
# ---------------------------------------------------------------------------
if [ "$LOCAL_RANK" -eq 0 ]; then
    # ---- Engine process ----
    CPU_BIND="${GPU_TO_ENGINE_CPU[$GPU_ID]}"

    if [ -z "$CPU_BIND" ]; then
        echo "Error: no CPU binding defined for GPU $GPU_ID" >&2
        exit 1
    fi

    export CUDA_VISIBLE_DEVICES="$GPU_ID"
    export OMP_NUM_THREADS=1
    export TORCH_NUM_THREADS=4
    export MKL_NUM_THREADS=4

    echo "Engine (rank $LOCAL_RANK) -> GPU $GPU_ID, CPUs $CPU_BIND" >&2

    # Start GPU memory/utilization monitor
    MEMORY_LOG="${SEAMM_MEMORY_LOG:-./gpu_${GPU_ID}_engine.log}"
    nvidia-smi --query-gpu=timestamp,memory.used,memory.free,utilization.gpu \
               --format=csv -l 5 -i "$GPU_ID" > "$MEMORY_LOG" &
    MONITOR_PID=$!
    echo "GPU monitor PID = $MONITOR_PID" >&2
    trap "kill $MONITOR_PID 2>/dev/null; wait $MONITOR_PID 2>/dev/null" EXIT

    echo "$@" > engine.cmd
    
    # Run the engine
    taskset -c "$CPU_BIND" "$@"

    # Clean up monitor
    kill "$MONITOR_PID" 2>/dev/null
    wait "$MONITOR_PID" 2>/dev/null
    echo "Done!" >> "$MEMORY_LOG"
    echo "Engine finished." >&2

else
    # ---- Driver process (LAMMPS) ----
    CPU_BIND="${GPU_TO_DRIVER_CPU[$GPU_ID]}"

    if [ -z "$CPU_BIND" ]; then
        echo "Error: no CPU binding defined for driver with GPU $GPU_ID" >&2
        exit 1
    fi

    # LAMMPS doesn't need GPU access
    export CUDA_VISIBLE_DEVICES=""
    export OMP_NUM_THREADS=1

    echo "Driver (rank $LOCAL_RANK) -> no GPU, CPUs $CPU_BIND" >&2

    echo "$@" > driver.cmd
    
    # Run LAMMPS
    taskset -c "$CPU_BIND" "$@"

    echo "Driver finished." >&2
fi
