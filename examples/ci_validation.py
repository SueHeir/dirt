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

Sweep = tuple[str, tuple[str, ...]]

# Representative gated validation drivers for every PR. This is deliberately a
# bounded smoke suite, not one token benchmark: it covers normal/tangential
# contact, wall impact, friction/rolling/twisting, Haff cooling, thermal
# conductivity, restart determinism, damping, a no-history tangential model, and
# the peri/DEM interop gate.
SMOKE_SWEEPS: list[Sweep] = [
    ("examples/bench_hertz_rebound/sweep.py", ()),
    ("examples/bench_hooke_rebound/sweep.py", ()),
    ("examples/bench_oblique_impact/sweep.py", ()),
    ("examples/bench_kharaz_oblique/sweep.py", ()),
    ("examples/bench_sliding_friction/sweep.py", ()),
    ("examples/bench_twisting_friction/sweep.py", ()),
    ("examples/bench_sphere_haff_cooling/sweep.py", ()),
    ("examples/bench_granular_conductivity/sweep.py", ()),
    ("examples/bench_restart_determinism/sweep.py", ()),
    ("examples/bench_cundall_damping/sweep.py", ()),
    ("examples/bench_nohistory_tangential/sweep.py", ()),
    ("examples/peri_dem_interop/sweep.py", ()),
]

# Broader green no-MPI suite for scheduled/manual CI. It intentionally excludes
# documented known-fail and MPI-only checks listed below.
FULL_SWEEPS: list[Sweep] = [
    ("examples/bench_angle_of_repose/sweep.py", ()),
    ("examples/bench_bond_breakage/sweep.py", ()),
    ("examples/bench_chung_ooi_impact/sweep.py", ()),
    ("examples/bench_clump_haff_cooling/sweep.py", ()),
    ("examples/bench_clump_inertia_sampler/sweep.py", ()),
    ("examples/bench_clump_insertion_determinism/sweep.py", ()),
    ("examples/bench_convergence/sweep.py", ()),
    ("examples/bench_cundall_damping/sweep.py", ()),
    ("examples/bench_curtis_cantilever/sweep.py", ()),
    ("examples/bench_dmt_sjkr_cohesion/sweep.py", ()),
    ("examples/bench_fiber_crossover/sweep.py", ()),
    ("examples/bench_granular_conductivity/sweep.py", ()),
    ("examples/bench_hertz_rebound/sweep.py", ()),
    ("examples/bench_hooke_wall_rebound/sweep.py", ()),
    ("examples/bench_hooke_rebound/sweep.py", ()),
    ("examples/bench_hopper_beverloo/sweep.py", ()),
    ("examples/bench_jkr_adhesion/sweep.py", ()),
    ("examples/bench_kharaz_oblique/sweep.py", ()),
    ("examples/bench_lebc_shear/sweep.py", ()),
    ("examples/bench_liquid_bridge_cohesion/sweep.py", ()),
    ("examples/bench_marshall_twisting/sweep.py", ()),
    ("examples/bench_mdr_elastoplastic_normal/sweep.py", ()),
    ("examples/bench_mindlin_rescale_tangential/sweep.py", ()),
    ("examples/bench_nohistory_tangential/sweep.py", ()),
    ("examples/bench_oblique_impact/sweep.py", ()),
    ("examples/bench_plate_sinkage/sweep.py", ()),
    ("examples/bench_polydisperse_mixing/sweep.py", ()),
    ("examples/bench_restart_determinism/sweep.py", ()),
    ("examples/bench_rod_haff_cooling/sweep.py", ()),
    ("examples/bench_rolling_decay/sweep.py", ()),
    ("examples/bench_sds_rolling/sweep.py", ()),
    ("examples/bench_sliding_friction/sweep.py", ()),
    ("examples/bench_sphere_haff_cooling/sweep.py", ()),
    ("examples/bench_twisting_friction/sweep.py", ()),
    ("examples/bench_wall_activate_by_name/sweep.py", ()),
    ("examples/bench_wall_twisting_parity/sweep.py", ()),
    ("examples/peri_dem_interop/sweep.py", ()),
    (
        "examples/SPH_glass_sphere_calibration/01_shear_rheology/sweep.py",
        ("smoke",),
    ),
    ("examples/SPH_glass_sphere_calibration/02_compressibility/sweep.py", ()),
    (
        "examples/SPH_glass_sphere_calibration/04_enduring_contact/sweep.py",
        ("smoke",),
    ),
    ("examples/SPH_glass_sphere_calibration/05_cooling_dissipation/sweep.py", ()),
    ("examples/SPH_glass_sphere_calibration/06_conductivity/sweep.py", ("smoke",)),
    ("examples/bond_cantilever/sweep.py", ()),
    ("examples/fiber_bond/sweep.py", ()),
    ("examples/fiber_bond_breakage/sweep.py", ()),
]

DOCUMENTED_EXCLUSIONS = {
    "examples/bench_column_collapse/sweep.py": (
        "documented honest FAIL in examples/VALIDATION.md; retained outside "
        "green PR CI so it cannot be mistaken for a regression"
    ),
    "examples/bench_mpi_decomposition/sweep.py": (
        "requires MPI/default-feature runtime, while this workflow deliberately "
        "uses the stock no-MPI quickstart configuration"
    ),
    "examples/bond_mpi_drift/sweep.py": (
        "requires a two-rank MPI run, while this workflow deliberately uses the "
        "stock no-MPI quickstart configuration"
    ),
    "examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py": (
        "documented honest FAIL in its README: measured glass sliding friction "
        "does not reach the target repose band for any swept rolling friction"
    ),
    "examples/SPH_glass_sphere_calibration/07_column_collapse/sweep.py": (
        "documented provisional macro validation: depends on the calibrated "
        "rolling friction from 03_angle_of_repose, which is currently unavailable"
    ),
    "examples/SPH_glass_sphere_calibration/08_cooperativity_length/sweep.py": (
        "documented exploratory calibration rig: no independent reference value "
        "or justified tolerance exists yet for A or the g∝sqrt(T) bridge"
    ),
}


def format_sweep(sweep: Sweep) -> str:
    path, extra_args = sweep
    if extra_args:
        return " ".join((path, *extra_args))
    return path


def sweep_paths(sweeps: list[Sweep]) -> set[str]:
    return {path for path, _ in sweeps}


def discover_sweeps() -> set[str]:
    return {
        path.relative_to(REPO_ROOT).as_posix()
        for path in (REPO_ROOT / "examples").rglob("sweep.py")
    }


def validate_manifest() -> bool:
    discovered = discover_sweeps()
    run_paths = sweep_paths(FULL_SWEEPS)
    excluded_paths = set(DOCUMENTED_EXCLUSIONS)
    accounted = run_paths | excluded_paths

    missing = sorted(discovered - accounted)
    stale = sorted(accounted - discovered)
    duplicate_runs = sorted(
        path
        for path in run_paths
        if sum(1 for sweep_path, _ in FULL_SWEEPS if sweep_path == path) > 1
    )

    if not missing and not stale and not duplicate_runs:
        return True

    print("CI validation manifest is inconsistent:")
    for path in missing:
        print(f"  MISSING {path}")
    for path in stale:
        print(f"  STALE {path}")
    for path in duplicate_runs:
        print(f"  DUPLICATE {path}")
    return False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--suite",
        choices=("smoke", "full"),
        default="smoke",
        help="validation suite to run (default: smoke)",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list selected sweeps, documented exclusions, and accounting",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sweeps = SMOKE_SWEEPS if args.suite == "smoke" else FULL_SWEEPS

    if not validate_manifest():
        return 2

    if args.list:
        for sweep in sweeps:
            print(format_sweep(sweep))
        print("\nDocumented exclusions from stock no-MPI green CI:")
        for sweep, reason in DOCUMENTED_EXCLUSIONS.items():
            print(f"{sweep}: {reason}")
        print(
            "\nAccounting: "
            f"{len(FULL_SWEEPS)} run + {len(DOCUMENTED_EXCLUSIONS)} excluded = "
            f"{len(discover_sweeps())} discovered sweep.py drivers"
        )
        return 0

    python = os.environ.get("BENCH_PYTHON", sys.executable)
    timeout_s = int(os.environ.get("DIRT_CI_BENCH_TIMEOUT", "1800"))
    failures: list[str] = []

    print(f"DIRT CI validation sweeps ({args.suite}):")
    for sweep in sweeps:
        sweep_path, extra_args = sweep
        print(f"  RUN {format_sweep(sweep)}")
        try:
            proc = subprocess.run(
                [python, sweep_path, *extra_args],
                cwd=REPO_ROOT,
                timeout=timeout_s,
                check=False,
            )
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT {format_sweep(sweep)} after {timeout_s}s")
            failures.append(f"{format_sweep(sweep)} (timeout)")
            continue

        if proc.returncode == 0:
            print(f"  PASS {format_sweep(sweep)}")
        else:
            print(f"  FAIL {format_sweep(sweep)} exit={proc.returncode}")
            failures.append(f"{format_sweep(sweep)} (exit {proc.returncode})")

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
