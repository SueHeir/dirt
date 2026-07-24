#!/usr/bin/env bash
# Run the complete, evidence-qualified column-collapse ensemble on a host that
# has no batch scheduler.  This is intentionally a thin launcher: sweep.py
# remains the sole owner of source generation, receipts, admission, fitting,
# and the external experimental gate.
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
JOBS=${1:-1}

if ! [[ "$JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "usage: $0 [positive-worker-count]" >&2
    exit 2
fi

# Match the array launcher: do not let a non-interactive shell silently build
# a different precision/configuration from the one used to retain witnesses.
BUILD_ENV=${DIRT_BUILD_ENV:-"$HOME/projects/.build-env"}
if [[ ! -r "$BUILD_ENV" ]]; then
    echo "ERROR: missing readable build environment: $BUILD_ENV" >&2
    exit 1
fi
source "$BUILD_ENV"

cd "$ROOT"
python3 examples/bench_column_collapse/sweep.py generate
python3 examples/bench_column_collapse/sweep.py emit-jobs
python3 examples/bench_column_collapse/sweep.py validate-jobs
python3 examples/bench_column_collapse/sweep.py start --jobs "$JOBS"
python3 examples/bench_column_collapse/sweep.py graph
