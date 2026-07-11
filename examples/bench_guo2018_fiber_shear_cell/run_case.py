#!/usr/bin/env python3
"""Run one source-parameterized Guo shear-cell case with the DIRT executable.

This is the bridge deliberately absent from the original PR: it creates the
audited topology, writes a case-specific DIRT configuration, and invokes the
solver.  It does not fabricate histories or call the acceptance validator.
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEPTH_M = 0.036


def decomposition(ranks: int) -> tuple[int, int]:
    """Return an x/z decomposition for the periodic horizontal cell axes."""
    if ranks < 1:
        raise ValueError("ranks must be positive")
    x = math.isqrt(ranks)
    while x > 1 and ranks % x:
        x -= 1
    return x, ranks // x


def materialize(pressure_pa: float, width_mm: int, output: Path, ranks: int) -> Path:
    output = output.resolve()
    if output.exists():
        shutil.rmtree(output)
    subprocess.run([sys.executable, str(HERE / "prepare.py"), "--width-mm", str(width_mm),
                    "--output", str(output)], check=True)
    subprocess.run([sys.executable, str(HERE / "prepare.py"), "--audit", str(output)], check=True)
    area = width_mm / 1000.0 * DEPTH_M
    config = (HERE / "config.toml").read_text()
    processors_x, processors_z = decomposition(ranks)
    config = config.replace("processors_x = 1", f"processors_x = {processors_x}")
    config = config.replace("processors_z = 1", f"processors_z = {processors_z}")
    config = config.replace("x_high = 0.064", f"x_high = {width_mm / 1000.0:.3f}")
    config = config.replace("target_force = 1.499904", f"target_force = {pressure_pa * area:.12g}")
    config = config.replace(
        "examples/bench_guo2018_fiber_shear_cell/generated/pack.csv", str(output / "pack.csv"))
    config = config.replace(
        "examples/bench_guo2018_fiber_shear_cell/generated/pack.bonds", str(output / "pack.bonds"))
    config = config.replace(
        'dir = "examples/bench_guo2018_fiber_shear_cell/generated"', f'dir = "{output}"')
    path = output / "case.toml"
    path.write_text(config)
    (output / "case_manifest.json").write_text(json.dumps({
        "pressure_pa": pressure_pa,
        "width_mm": width_mm,
        "ranks": ranks,
        "processors_x": processors_x,
        "processors_z": processors_z,
        "runner": "run_case.py",
        "solver_history": "cell_history.csv",
    }, indent=2) + "\n")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pressure-pa", type=float, choices=(651.0, 1735.0, 3470.0), required=True)
    parser.add_argument("--width-mm", type=int, choices=(64, 96), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=1,
                        help="MPI ranks; decomposed only over periodic x/z axes")
    parser.add_argument("--prepare-only", action="store_true", help="write and audit input but do not invoke DIRT")
    args = parser.parse_args()
    config = materialize(args.pressure_pa, args.width_mm, args.output, args.ranks)
    print(f"prepared solver case: {config}")
    if not args.prepare_only:
        subprocess.run(["cargo", "build", "--release", "--example", "bench_guo2018_fiber_shear_cell",
                        "--no-default-features", "--features", "precision-double"], cwd=ROOT, check=True)
        solver = ROOT / "target/release/examples/bench_guo2018_fiber_shear_cell"
        command = [str(solver), str(config)]
        if args.ranks > 1:
            command = [os.environ.get("MPIEXEC", "mpirun"), "-n", str(args.ranks), *command]
        subprocess.run(command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
