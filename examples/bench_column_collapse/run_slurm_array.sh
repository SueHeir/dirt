#!/usr/bin/env bash
# Execute one immutable column-collapse witness in a Slurm job array.
#
# Usage (from a clean DIRT checkout):
#   python3 examples/bench_column_collapse/sweep.py generate
#   python3 examples/bench_column_collapse/sweep.py emit-jobs
#   sbatch --array=1-33 --cpus-per-task=1 --time=08:00:00 \
#     examples/bench_column_collapse/run_slurm_array.sh
#
# This script deliberately runs exactly one `start --case` command.  It does
# not aggregate data, create runout.csv, or decide PASS/FAIL; `sweep.py graph`
# remains the fail-closed admission and experimental-law gate after every array
# task has completed successfully.
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
manifest="$repo_root/examples/bench_column_collapse/sweep/column_collapse_jobs.tsv"
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required (submit with sbatch --array=1-33)}

if ! [[ "$task_id" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID must be a positive integer, got '$task_id'" >&2
    exit 2
fi
if [[ ! -f "$manifest" ]]; then
    echo "ERROR: missing $manifest; run sweep.py generate then sweep.py emit-jobs" >&2
    exit 2
fi

# The first field is the immutable 1-based case index.  Do not select by line
# number: a manually edited, duplicate, or re-ordered manifest must fail rather
# than silently execute the wrong realization under a valid array index.
row=$(awk -F '\t' -v wanted="$task_id" '$1 == wanted { if (++n == 1) value=$0 } END { if (n == 1) print value; else exit 1 }' "$manifest") \
    || { echo "ERROR: expected exactly one manifest row for array index $task_id" >&2; exit 2; }

IFS=$'\t' read -r index aspect seed active_count source_sha protocol_sha command <<<"$row"
if [[ "$index" != "$task_id" || ! "$aspect" =~ ^[0-9]+([.][0-9]+)?$ || ! "$seed" =~ ^[0-9]+$ ]]; then
    echo "ERROR: malformed manifest row for array index $task_id" >&2
    exit 2
fi

echo "column-collapse array task=$task_id aspect=$aspect seed=$seed particles=$active_count protocol=${protocol_sha:0:12}"
cd "$repo_root"
if [[ "${COLUMN_COLLAPSE_ARRAY_DRY_RUN:-0}" == "1" ]]; then
    echo "DRY RUN: python3 examples/bench_column_collapse/sweep.py start --case $aspect,$seed"
    exit 0
fi
exec python3 examples/bench_column_collapse/sweep.py start --case "$aspect,$seed"
