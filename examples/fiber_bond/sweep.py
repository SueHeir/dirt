#!/usr/bin/env python3
"""Run the Guo/Curtis fiber-bond validation cases used by CI/regression."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
CONFIG = SCRIPT_DIR / "bending_plastic_guo.toml"
CSV = SCRIPT_DIR / "bending_plastic_guo" / "data" / "fiber_bond.csv"


def run(cmd: list[str], env: dict[str, str] | None = None) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=REPO_ROOT, env=env, check=True)


def main() -> int:
    env = os.environ.copy()
    cargo = env.get("CARGO", "cargo")
    python = env.get("BENCH_PYTHON", sys.executable)

    run([
        cargo,
        "run",
        "--release",
        "--example",
        "fiber_bond",
        "--no-default-features",
        "--features",
        "precision-mixed",
        "--",
        str(CONFIG.relative_to(REPO_ROOT)),
    ], env=env)
    run([python, str((SCRIPT_DIR / "validate.py").relative_to(REPO_ROOT)), str(CSV.relative_to(REPO_ROOT))], env=env)
    print("VALIDATION: PASS fiber_bond Guo/Curtis bending plasticity")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
