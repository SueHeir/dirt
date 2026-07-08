#!/usr/bin/env python3
"""Run the fiber-bond breakage validation scenarios end-to-end."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

SCENARIOS = [
    "axial_stress_constant",
    "axial_strain_constant",
    "axial_stress_weibull",
    "combined_stress",
    "combined_strain",
    "interaction_linear_stress",
    "axial_plastic_stress_constant",
]


def run(cmd: list[str], env: dict[str, str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)


def main() -> int:
    env = os.environ.copy()
    run([
        "cargo", "build", "--release", "--example", "fiber_bond",
        "--no-default-features", "--features", "precision-double",
    ], env)

    for name in SCENARIOS:
        cfg = HERE / f"{name}.toml"
        csv = HERE / name / "data" / "fiber_bond.csv"
        run([
            "cargo", "run", "--release", "--example", "fiber_bond",
            "--no-default-features", "--features", "precision-double",
            "--", str(cfg.relative_to(ROOT)),
        ], env)
        run([sys.executable, str(HERE / "validate.py"), str(csv)], env)

    print(f"VALIDATION: PASS ({len(SCENARIOS)}/{len(SCENARIOS)} scenarios)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
