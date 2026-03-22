"""
gpu_memory_parser.py

Parse nvidia-smi memory monitor logs produced by gpu_bind.sh and summarize
GPU memory usage and utilization, including multi-GPU summaries.

The command for nvidia-smi looks like this::

    # Start background memory monitor, writing peak to a log file
    MEMORY_LOG="${SEAMM_MEMORY_LOG:-./gpu_${GPU_ID}_rank_${LOCAL_RANK}.log}"

    nvidia-smi --query-gpu=timestamp,memory.used,memory.free,utilization.gpu \
               --format=csv -l 5 -i "$GPU_ID" > "$MEMORY_LOG" &
"""

import re
from pathlib import Path
from datetime import datetime


def parse_gpu_memory_log(log_path):
    """
    Parse a gpu_bind.sh memory monitor log file and return summary statistics.

    Initialization phase: samples before GPU utilization first becomes non-zero
    AND memory has stabilized (no longer growing significantly).
    Production phase: all samples after initialization.

    Args:
        log_path: Path to the nvidia-smi log file

    Returns:
        dict with statistics, or None if file cannot be parsed
    """
    log_path = Path(log_path)
    if not log_path.exists():
        return None

    # Try to extract GPU index and rank from filename
    # Expected patterns: gpu_0_rank_0.log, gpu_memory_rank0.log, etc.
    gpu_id = None
    rank = None
    name = log_path.stem
    m = re.search(r"gpu[_]?(\d+)", name, re.IGNORECASE)
    if m:
        gpu_id = int(m.group(1))
    m = re.search(r"rank[_]?(\d+)", name, re.IGNORECASE)
    if m:
        rank = int(m.group(1))

    samples = []
    with open(log_path) as f:
        for line in f:
            line = line.strip()
            # Skip header and "Done!" line
            if not line or line.startswith("timestamp") or not line[0].isdigit():
                continue
            try:
                parts = [p.strip() for p in line.split(",")]
                timestamp = datetime.strptime(parts[0], "%Y/%m/%d %H:%M:%S.%f")
                mem_used = int(parts[1].split()[0])
                mem_free = int(parts[2].split()[0])
                util = int(parts[3].split()[0])
                samples.append(
                    {
                        "timestamp": timestamp,
                        "mem_used": mem_used,
                        "mem_free": mem_free,
                        "util": util,
                    }
                )
            except (ValueError, IndexError):
                continue

    if not samples:
        return None

    # Find where initialization ends:
    # First sample where utilization > 0 AND memory has stabilized
    # Memory is "stable" when the increase from previous sample is small (<5%)
    MEM_GROWTH_THRESHOLD = 0.05  # 5% growth = still initializing

    init_end_idx = None
    for i, s in enumerate(samples):
        if s["util"] == 0:
            continue
        # Utilization is non-zero -- check if memory is still growing
        if i > 0:
            prev_mem = samples[i - 1]["mem_used"]
            if prev_mem > 0:
                growth = (s["mem_used"] - prev_mem) / prev_mem
                if growth > MEM_GROWTH_THRESHOLD:
                    continue  # still initializing -- memory growing fast
        init_end_idx = i
        break

    # Initialization phase
    init_samples = samples[:init_end_idx] if init_end_idx is not None else samples
    prod_samples = samples[init_end_idx:] if init_end_idx is not None else []

    # Initialization duration
    if len(init_samples) >= 1 and init_end_idx is not None:
        init_duration_s = (
            samples[init_end_idx]["timestamp"] - samples[0]["timestamp"]
        ).total_seconds()
    elif len(init_samples) >= 2:
        init_duration_s = (
            init_samples[-1]["timestamp"] - init_samples[0]["timestamp"]
        ).total_seconds()
    else:
        init_duration_s = None

    # Total run duration
    if len(samples) >= 2:
        total_duration_s = (
            samples[-1]["timestamp"] - samples[0]["timestamp"]
        ).total_seconds()
    else:
        total_duration_s = None

    # Production statistics
    if prod_samples:
        prod_mem_used = [s["mem_used"] for s in prod_samples]
        prod_mem_free = [s["mem_free"] for s in prod_samples]
        prod_util = [s["util"] for s in prod_samples]
        prod_stats = {
            "max_memory_used_mb": max(prod_mem_used),
            "min_memory_free_mb": min(prod_mem_free),
            "avg_gpu_utilization": sum(prod_util) / len(prod_util),
            "max_gpu_utilization": max(prod_util),
            "n_samples": len(prod_samples),
        }
    else:
        prod_stats = None

    return {
        "log_path": str(log_path),
        "gpu_id": gpu_id,
        "rank": rank,
        "init_duration_s": init_duration_s,
        "init_n_samples": len(init_samples),
        "total_duration_s": total_duration_s,
        "total_n_samples": len(samples),
        "production": prod_stats,
        "samples": samples,
    }


def print_gpu_summary(stats, label=None, fd=None):
    """Print summary for a single GPU log."""
    if stats is None:
        print("Could not parse log file.", file=fd)
        return

    # Header
    if label:
        header = label
    elif stats["gpu_id"] is not None and stats["rank"] is not None:
        header = f"GPU {stats['gpu_id']} (rank {stats['rank']})"
    else:
        header = Path(stats["log_path"]).name

    print(f"\n{header}", file=fd)
    print(f"  {'Total samples:':25s} {stats['total_n_samples']}", file=fd)

    if stats["total_duration_s"] is not None:
        print(f"  {'Total duration:':25s} {stats['total_duration_s']:.1f} s", file=fd)

    if stats["init_duration_s"] is not None:
        print(
            f"  {'Initialization:':25s} {stats['init_duration_s']:.1f} s "
            f"({stats['init_n_samples']} samples)",
            file=fd,
        )

    p = stats["production"]
    if p:
        prod_duration = (
            stats["total_duration_s"] - stats["init_duration_s"]
            if stats["total_duration_s"] and stats["init_duration_s"]
            else None
        )
        print(
            f"  {'Production phase:':25s} "
            + (f"{prod_duration:.1f} s " if prod_duration else "")
            + f"({p['n_samples']} samples)",
            file=fd,
        )
        print(f"  {'Max memory used:':25s} {p['max_memory_used_mb']:,} MiB", file=fd)
        print(f"  {'Min memory free:':25s} {p['min_memory_free_mb']:,} MiB", file=fd)
        print(
            f"  {'Avg GPU utilization:':25s} {p['avg_gpu_utilization']:.1f}%", file=fd
        )
        print(f"  {'Max GPU utilization:':25s} {p['max_gpu_utilization']}%", file=fd)
    else:
        print(
            "  No production samples found (GPU utilization never exceeded 0%)", file=fd
        )


def print_multi_gpu_summary(log_paths, fd=None):
    """
    Parse multiple GPU log files and print per-GPU and combined summary.

    Args:
        log_paths: list of paths to log files
    """
    results = []
    for path in log_paths:
        stats = parse_gpu_memory_log(path)
        if stats is None:
            print(f"WARNING: Could not parse {path}")
            continue
        results.append(stats)

    if not results:
        print("No valid log files found.", file=fd)
        return

    # Sort by rank if available
    results.sort(key=lambda x: x["rank"] if x["rank"] is not None else 0)

    print("=" * 50, file=fd)
    print(f"GPU Memory Summary ({len(results)} GPU(s))", file=fd)
    print("=" * 50, file=fd)

    # Per-GPU summaries
    for stats in results:
        print_gpu_summary(stats, fd=fd)

    # Combined summary across all GPUs
    ngpus = len(results)
    print(f"\nStatistics for {ngpus} GPUs:", file=fd)
    print("-" * 50, file=fd)

    all_prod = [r["production"] for r in results if r["production"]]
    if all_prod:
        max_mem = max(p["max_memory_used_mb"] for p in all_prod)
        min_mem_free = min(p["min_memory_free_mb"] for p in all_prod)
        avg_util = sum(p["avg_gpu_utilization"] for p in all_prod) / len(all_prod)
        max_util = max(p["max_gpu_utilization"] for p in all_prod)
        min_util = min(p["avg_gpu_utilization"] for p in all_prod)
        util_imbalance = max(p["avg_gpu_utilization"] for p in all_prod) - min_util

        # Initialization — use the longest init time (slowest GPU sets the pace)
        init_times = [
            r["init_duration_s"] for r in results if r["init_duration_s"] is not None
        ]
        max_init = max(init_times) if init_times else None

        if ngpus == 1:
            print(
                f"  {'GPU memory used:':25s} {max_mem:,} MiB "
                f"({max_mem/1024:.1f} GiB)",
                file=fd,
            )
            print(
                f"  {'Minimum GPU memory free:':25s} {min_mem_free:,} MiB "
                f"({min_mem_free/1024:.1f} GiB)",
                file=fd,
            )
            print(f"  {'GPU utilization:':25s} {avg_util:.1f}%", file=fd)
            print(f"  {'Max GPU utilization:':25s} {max_util}%", file=fd)
        else:
            print(
                f"  {'Maximum GPU memory used:':25s} {max_mem:,} MiB "
                f"({max_mem/1024:.1f} GiB)",
                file=fd,
            )
            print(
                f"  {'Minimum GPU memory free:':25s} {min_mem_free:,} MiB "
                f"({min_mem_free/1024:.1f} GiB)",
                file=fd,
            )
            print(f"  {'Avg GPU utilization:':25s} {avg_util:.1f}%", file=fd)
            print(f"  {'Max GPU utilization:':25s} {max_util}%", file=fd)
            print(
                f"  {'GPU util imbalance:':25s} {util_imbalance:.1f}% "
                + ("(good)" if util_imbalance < 5 else "(significant)"),
                file=fd,
            )
        if max_init is not None:
            if ngpus == 1:
                print(f"  {'Initialization time:':25s} {max_init:.1f} s", file=fd)
            else:
                print(f"  {'Max initialization time:':25s} {max_init:.1f} s", file=fd)

        return {
            "Number of GPUs": ngpus,
            "Maximum GPU memory used": round(max_mem / 1024, 1),
            "Minimum free GPU memory": round(min_mem_free / 1024, 1),
            "% GPU memory used": round(100 * max_mem / (max_mem + min_mem_free), 1),
            "Average GPU utilization": round(avg_util, 1),
            "Maximum GPU utilization": round(max_util, 1),
            "GPU utilization imbalance": round(util_imbalance, 1),
            "Maximum initialization time": round(max_init, 0),
        }
    else:
        print("  No production data available.", file=fd)

    print("=" * 50, file=fd)


if __name__ == "__main__":
    import sys
    import glob

    if len(sys.argv) < 2:
        print("Usage: python gpu_memory_parser.py <log_file> [log_file2 ...]")
        print("       python gpu_memory_parser.py 'gpu_*.log'")
        sys.exit(1)

    # Expand globs
    paths = []
    for pattern in sys.argv[1:]:
        expanded = glob.glob(pattern)
        if not expanded:
            paths.append(pattern)  # keep as-is, will fail gracefully
        else:
            paths.extend(expanded)

    if len(paths) == 0:
        stats = parse_gpu_memory_log(paths[0])
        print_gpu_summary(stats)
    else:
        print_multi_gpu_summary(paths)
