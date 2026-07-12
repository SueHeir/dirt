#!/usr/bin/env python3
"""Run one source-parameterized Guo shear-cell case with the DIRT executable.

This is the bridge deliberately absent from the original PR: it creates the
audited topology, writes a case-specific DIRT configuration, and invokes the
solver.  It does not fabricate histories or call the acceptance validator.
"""
import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def require_published_control_cell() -> None:
    """Refuse an input that differs from the paper's published DEM reduction."""
    sys.path.insert(0, str(HERE))
    from source_contract import require_published_control_cell as require_protocol
    require_protocol(HERE / "config.toml")


def require_primary_reference(source_pdf: Path | None) -> None:
    """Stop before creating an allegedly source-derived solver input."""
    sys.path.insert(0, str(HERE))
    from evidence_contract import verify_primary_reference
    # Do not accept a DOI, web landing page, or arbitrary file. The caller
    # supplies the local paper while the committed manifest supplies its hash.
    if source_pdf is None:
        raise RuntimeError(
            "BLOCKED: --source-pdf is required; refusing to "
            "materialize an allegedly Guo-derived geometry"
        )
    try:
        verify_primary_reference(source_pdf)
    except ValueError as error:
        raise RuntimeError(f"BLOCKED: Guo primary-method provenance failed: {error}") from error



def decomposition(ranks: int) -> tuple[int, int]:
    """Return an x/z decomposition for the periodic horizontal cell axes."""
    if ranks < 1:
        raise ValueError("ranks must be positive")
    x = math.isqrt(ranks)
    while x > 1 and ranks % x:
        x -= 1
    return x, ranks // x


def sha256(path: Path) -> str:
    """Content digest used to bind a receipt to its exact inputs."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def materialize(pressure_pa: float, width_mm: int, output: Path, ranks: int) -> Path:
    output = output.resolve()
    if output.exists():
        shutil.rmtree(output)
    subprocess.run([sys.executable, str(HERE / "prepare.py"), "--width-mm", str(width_mm),
                    "--output", str(output)], check=True)
    subprocess.run([sys.executable, str(HERE / "prepare.py"), "--audit", str(output)], check=True)
    width = width_mm / 1000.0
    # Primary source: Lx=64 mm (or 96-mm sensitivity case), Lz=36 mm.
    # This is a planar periodic cell, not the historical circular-cup area.
    area = width * 0.036
    config = (HERE / "config.toml").read_text()
    processors_x, processors_z = decomposition(ranks)
    config = config.replace("processors_x = 1", f"processors_x = {processors_x}")
    config = config.replace("processors_z = 1", f"processors_z = {processors_z}")
    config = config.replace("x_high = 0.064", f"x_high = {width:.3f}")
    config = config.replace("target_force = 1499.904", f"target_force = {pressure_pa * area:.12g}")
    config = config.replace(
        "examples/bench_guo2018_fiber_shear_cell/generated/pack.csv", str(output / "pack.csv"))
    config = config.replace(
        "examples/bench_guo2018_fiber_shear_cell/generated/pack.bonds", str(output / "pack.bonds"))
    config = config.replace(
        'dir = "examples/bench_guo2018_fiber_shear_cell/generated"', f'dir = "{output}"')
    path = output / "case.toml"
    path.write_text(config)
    source = json.loads((output / "source_parameters.json").read_text())
    (output / "case_manifest.json").write_text(json.dumps({
        "pressure_pa": pressure_pa,
        "width_mm": width_mm,
        "ranks": ranks,
        "processors_x": processors_x,
        "processors_z": processors_z, "planform_area_m2": area,
        "runner": "run_case.py",
        "solver_history": "cell_history.csv",
        "expected_global_atoms": source["n_atoms"],
        "input_sha256": {
            "pack.csv": sha256(output / "pack.csv"),
            "pack.bonds": sha256(output / "pack.bonds"),
            "case.toml": sha256(path),
        },
    }, indent=2) + "\n")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pressure-pa", type=float, choices=(651.0, 1735.0, 3470.0), required=True)
    parser.add_argument("--width-mm", type=int, choices=(64, 96), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=1,
                        help="MPI ranks; decomposed only over periodic x/z axes")
    parser.add_argument("--source-pdf", type=Path,
                        help="local primary PDF matching data/reference_provenance.json")
    parser.add_argument("--prepare-only", action="store_true", help="write and audit input but do not invoke DIRT")
    args = parser.parse_args()
    # Both preparation and execution need an authenticated method source. A
    # solver receipt cannot repair an unauditable input geometry afterwards.
    require_primary_reference(args.source_pdf)
    # Preparation itself creates an allegedly source-derived topology and
    # case.toml.  Do not let --prepare-only bypass the geometry-equivalence
    # boundary: a non-equivalent input is misleading even without a solver
    # history, and could later be mistaken for a runnable campaign case.
    require_published_control_cell()
    config = materialize(args.pressure_pa, args.width_mm, args.output, args.ranks)
    print(f"prepared solver case: {config}")
    if not args.prepare_only:
        # The default DIRT feature set includes mpi_backend.  A serial-only
        # binary launched through mpirun creates one independent one-rank
        # communicator per process, so it must never be used for `--ranks > 1`.
        # Building the same MPI-capable executable for rank one keeps the
        # campaign command unambiguous and lets the manifest describe the
        # communicator that actually ran the case.
        subprocess.run(["cargo", "build", "--release", "--example", "bench_guo2018_fiber_shear_cell"],
                       cwd=ROOT, check=True)
        solver = ROOT / "target/release/examples/bench_guo2018_fiber_shear_cell"
        command = [str(solver), str(config)]
        if args.ranks > 1:
            command = [os.environ.get("MPIEXEC", "mpirun"), "-n", str(args.ranks), *command]
        subprocess.run(command, cwd=ROOT, check=True)
        history = args.output.resolve() / "cell_history.csv"
        if not history.is_file():
            raise RuntimeError("DIRT exited without writing cell_history.csv")
        # Written only after a successful solver invocation. This binds the
        # history to the materialized input; it is not a physics result.
        (args.output.resolve() / "solver_receipt.json").write_text(json.dumps({
            "runner": "run_case.py", "command": command,
            "history_sha256": sha256(history),
            "input_sha256": json.loads((args.output.resolve() / "case_manifest.json").read_text())["input_sha256"],
        }, indent=2) + "\n")


if __name__ == "__main__":
    main()
