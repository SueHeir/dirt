#!/usr/bin/env python3
"""Run the DIRT validation sweeps that are expected to be green in stock CI.

The full benchmark directory intentionally also contains honest known-fail and
environment-specific checks. Keep those exclusions explicit here so the CI
coverage decision is reviewable instead of hiding behind a one-off command.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

# Representative gated validation drivers for every PR. This is deliberately a
# bounded smoke suite, not one token benchmark: it covers normal/tangential
# contact, wall impact, friction/rolling/twisting, Haff cooling, thermal
# conductivity, restart determinism, damping, a no-history tangential model, and
# the peri/DEM interop gate.
SMOKE_SWEEPS = [
    "examples/bench_hertz_rebound/sweep.py",
    "examples/bench_hooke_rebound/sweep.py",
    "examples/bench_oblique_impact/sweep.py",
    "examples/bench_kharaz_oblique/sweep.py",
    "examples/bench_sliding_friction/sweep.py",
    "examples/bench_twisting_friction/sweep.py",
    "examples/bench_sphere_haff_cooling/sweep.py",
    "examples/bench_granular_conductivity/sweep.py",
    "examples/bench_restart_determinism/sweep.py",
    "examples/bench_cundall_damping/sweep.py",
    "examples/bench_nohistory_tangential/sweep.py",
    "examples/peri_dem_interop/sweep.py",
]

# Broader green no-MPI suite for scheduled/manual CI. It intentionally excludes
# documented known-fail and MPI-only checks listed below.
FULL_SWEEPS = [
    "examples/bench_angle_of_repose/sweep.py",
    "examples/bench_bond_breakage/sweep.py",
    "examples/bench_chung_ooi_impact/sweep.py",
    "examples/bench_clump_haff_cooling/sweep.py",
    "examples/bench_convergence/sweep.py",
    "examples/bench_cundall_damping/sweep.py",
    "examples/bench_dmt_sjkr_cohesion/sweep.py",
    "examples/bench_fiber_crossover/sweep.py",
    "examples/bench_granular_conductivity/sweep.py",
    "examples/bench_hertz_rebound/sweep.py",
    "examples/bench_hooke_rebound/sweep.py",
    "examples/bench_hopper_beverloo/sweep.py",
    "examples/bench_jkr_adhesion/sweep.py",
    "examples/bench_kharaz_oblique/sweep.py",
    "examples/bench_lebc_shear/sweep.py",
    "examples/bench_marshall_twisting/sweep.py",
    "examples/bench_nohistory_tangential/sweep.py",
    "examples/bench_oblique_impact/sweep.py",
    "examples/bench_plate_sinkage/sweep.py",
    "examples/bench_polydisperse_mixing/sweep.py",
    "examples/bench_restart_determinism/sweep.py",
    "examples/bench_rolling_decay/sweep.py",
    "examples/bench_sds_rolling/sweep.py",
    "examples/bench_sliding_friction/sweep.py",
    "examples/bench_sphere_haff_cooling/sweep.py",
    "examples/bench_twisting_friction/sweep.py",
    "examples/peri_dem_interop/sweep.py",
]

DOCUMENTED_EXCLUSIONS = {
    "examples/bench_column_collapse/sweep.py": (
        "documented honest FAIL in examples/VALIDATION.md; retained outside "
        "green PR CI so it cannot be mistaken for a regression"
    ),
    "examples/bench_rod_haff_cooling/sweep.py": (
        "documented honest FAIL in examples/VALIDATION.md; retained outside "
        "green PR CI so it cannot be mistaken for a regression"
    ),
    "examples/bench_mpi_decomposition/sweep.py": (
        "requires MPI/default-feature runtime, while this workflow deliberately "
        "uses the stock no-MPI quickstart configuration"
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--suite",
        choices=("smoke", "full"),
        default="smoke",
        help="validation suite to run (default: smoke)",
    )
    parser.add_argument("--list", action="store_true", help="list selected sweeps and exit")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sweeps = SMOKE_SWEEPS if args.suite == "smoke" else FULL_SWEEPS

    if args.list:
        for sweep in sweeps:
            print(sweep)
        return 0

    python = os.environ.get("BENCH_PYTHON", sys.executable)
    timeout_s = int(os.environ.get("DIRT_CI_BENCH_TIMEOUT", "1800"))
    failures: list[str] = []

    print(f"DIRT CI validation sweeps ({args.suite}):")
    for sweep in sweeps:
        print(f"  RUN {sweep}")
        try:
            proc = subprocess.run(
                [python, sweep],
                cwd=REPO_ROOT,
                timeout=timeout_s,
                check=False,
            )
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT {sweep} after {timeout_s}s")
            failures.append(f"{sweep} (timeout)")
            continue

        if proc.returncode == 0:
            print(f"  PASS {sweep}")
        else:
            print(f"  FAIL {sweep} exit={proc.returncode}")
            failures.append(f"{sweep} (exit {proc.returncode})")

    print("\nDocumented exclusions from stock no-MPI green CI:")
    for sweep, reason in DOCUMENTED_EXCLUSIONS.items():
        print(f"  SKIP {sweep}: {reason}")

    if failures:
        print("\nCI VALIDATION FAILED")
        for failure in failures:
            print(f"  {failure}")
        return 1

    print("\nCI VALIDATION PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
