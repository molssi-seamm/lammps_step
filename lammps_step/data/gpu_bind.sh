#!/bin/bash
#MolSSI lammps_step:gpu_bind 1.0
if [ -z "$SEAMM_GPUS" ]; then
    echo "Error: SEAMM_GPUS is not set" >&2
    exit 1
fi

IFS=',' read -ra GPU_ARRAY <<< "$SEAMM_GPUS"
N_GPUS=${#GPU_ARRAY[@]}
LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK

if [ $LOCAL_RANK -ge $N_GPUS ]; then
    echo "Error: local rank $LOCAL_RANK exceeds number of GPUs in SEAMM_GPUS ($N_GPUS)" >&2
    exit 1
fi

# Map GPU index to L3 cache group for EPYC 7763
# GPU 0 -> L3 group 0 (cores 0-7), GPU 1 -> L3 group 4 (cores 32-39)
declare -A GPU_TO_CPU
GPU_TO_CPU[0]="0-7"
GPU_TO_CPU[1]="32-39"

GPU_ID=${GPU_ARRAY[$LOCAL_RANK]}
CPU_BIND=${GPU_TO_CPU[$GPU_ID]}

if [ -z "$CPU_BIND" ]; then
    echo "Error: no CPU binding defined for GPU $GPU_ID" >&2
    exit 1
fi

export CUDA_VISIBLE_DEVICES=$GPU_ID
export OMP_NUM_THREADS=1
export TORCH_NUM_THREADS=4
export MKL_NUM_THREADS=4

echo "Rank $LOCAL_RANK -> GPU $GPU_ID, CPUs $CPU_BIND" >&2

# Start background memory monitor, writing peak to a log file
MEMORY_LOG="${SEAMM_MEMORY_LOG:-./gpu_${GPU_ID}_rank_${LOCAL_RANK}.log}"

nvidia-smi --query-gpu=timestamp,memory.used,memory.free,utilization.gpu \
	   --format=csv -l 5 -i "$GPU_ID" > "$MEMORY_LOG" &
MONITOR_PID=$!
echo "The monitor PID = $MONITOR_PID" >&2
trap "kill $MONITOR_PID 2>/dev/null" EXIT

# Run LAMMPS
taskset -c "$CPU_BIND" "$@"

kill $MONITOR_PID
wait $MONITOR_PID
echo "Done!" >> "$MEMORY_LOG"
echo "The LAMMPS run has finished." >&2
