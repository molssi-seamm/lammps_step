#!/bin/bash
#MolSSI lammps_step:cpu_bind 1.0

if [ -z "$SEAMM_NP" ]; then
    echo "Error: SEAMM_NP is not set" >&2
    exit 1
fi

# L3 cache groups for EPYC 7763: 8 CCDs, 8 physical cores + HT siblings each
# Groups 0 and 4 are reserved for GPU jobs
CPU_RANGES=(
    "8-15,72-79"     # L3 group 1
    "16-23,80-87"    # L3 group 2
    "24-31,88-95"    # L3 group 3
    "40-47,104-111"  # L3 group 5
    "48-55,112-119"  # L3 group 6
    "56-63,120-127"  # L3 group 7
)

# Each group has 8 physical cores
PHYSICAL_CORES_PER_GROUP=8
N_GROUPS=$(( ($SEAMM_NP + $PHYSICAL_CORES_PER_GROUP - 1) / $PHYSICAL_CORES_PER_GROUP ))
N_AVAILABLE=${#CPU_RANGES[@]}

if [ $N_GROUPS -gt $N_AVAILABLE ]; then
    echo "Error: SEAMM_NP=$SEAMM_NP requires $N_GROUPS L3 groups but only $N_AVAILABLE are available for CPU jobs (max $((N_AVAILABLE * PHYSICAL_CORES_PER_GROUP)) physical cores)" >&2
    exit 1
fi

# Build the CPU bind string from as many groups as needed
CPU_BIND=${CPU_RANGES[0]}
for (( i=1; i<$N_GROUPS; i++ )); do
    CPU_BIND="$CPU_BIND,${CPU_RANGES[$i]}"
done

if [ "${SEAMM_DEBUG:-0}" = "1" ]; then
    echo "SEAMM_NP=$SEAMM_NP -> using $N_GROUPS L3 groups, CPUs $CPU_BIND" >&2
fi    echo "SEAMM_NP=$SEAMM_NP -> using $N_GROUPS L3 groups, CPUs $CPU_BIND" >&2
echo "SEAMM_NP=$SEAMM_NP -> using $N_GROUPS L3 groups, CPUs $CPU_BIND" >&2


exec numactl --physcpubind=$CPU_BIND "$@"
