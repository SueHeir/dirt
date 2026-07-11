#!/usr/bin/env python3
"""Run and plot the BPM MPI bond-migration check."""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples" / "bond_mpi_drift"
DATA = EXAMPLE / "data" / "bond_drift.csv"
PLOT_DIR = EXAMPLE / "plots"
PLOT = PLOT_DIR / "bond_mpi_drift_counts.png"
CONFIG_MPI2 = EXAMPLE / "config_mpi2.toml"
EXPECTED_BONDS = 2
EXPECTED_MISSING = 0


def _run(cmd: list[str]) -> None:
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def start() -> None:
    DATA.unlink(missing_ok=True)
    _run(["cargo", "build", "--release", "--example", "bond_mpi_drift"])
    _run(
        [
            "mpiexec",
            "-n",
            "2",
            str(ROOT / "target" / "release" / "examples" / "bond_mpi_drift"),
            str(CONFIG_MPI2),
        ]
    )


def _load_rows() -> list[dict[str, float]]:
    if not DATA.exists():
        raise FileNotFoundError(f"{DATA} does not exist; run 'sweep.py start' first")
    with DATA.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"{DATA} has no samples")
    return [
        {
            "step": float(r["step"]),
            "time_ms": float(r["t"]) * 1.0e3,
            "bond_count": float(r["bond_count"]),
            "bond_missing": float(r["bond_missing"]),
            "nlocal_global": float(r["nlocal_global"]),
            "ranks": float(r["ranks"]),
        }
        for r in rows
    ]


def validate(rows: list[dict[str, float]]) -> bool:
    bad_count = [r for r in rows if int(r["bond_count"]) != EXPECTED_BONDS]
    bad_missing = [r for r in rows if int(r["bond_missing"]) != EXPECTED_MISSING]
    bad_atoms = [r for r in rows if int(r["nlocal_global"]) != 3]

    min_count = min(int(r["bond_count"]) for r in rows)
    max_count = max(int(r["bond_count"]) for r in rows)
    max_missing = max(int(r["bond_missing"]) for r in rows)
    samples = len(rows)

    print(
        "bond_mpi_drift: "
        f"{samples} samples, bond_count min/max={min_count}/{max_count} "
        f"(reference {EXPECTED_BONDS}), bond_missing max={max_missing} "
        f"(reference {EXPECTED_MISSING})"
    )

    ok = not bad_count and not bad_missing and not bad_atoms
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    return ok


def plot(rows: list[dict[str, float]]) -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    x = [r["time_ms"] for r in rows]
    count = [r["bond_count"] for r in rows]
    missing = [r["bond_missing"] for r in rows]

    fig, (ax_count, ax_missing) = plt.subplots(
        2, 1, figsize=(7.2, 5.4), sharex=True, constrained_layout=True
    )

    ax_count.axhline(EXPECTED_BONDS, color="black", lw=1.4, ls="--", label="reference = 2")
    ax_count.plot(x, count, "o-", color="C0", ms=3, lw=1.2, label="measured")
    ax_count.set_ylabel("bond_count")
    ax_count.set_ylim(1.2, 2.8)
    ax_count.legend(loc="upper right", ncol=2, fontsize=8)
    ax_count.grid(True, alpha=0.25)

    ax_missing.axhline(EXPECTED_MISSING, color="black", lw=1.4, ls="--", label="reference = 0")
    ax_missing.plot(x, missing, "o-", color="C1", ms=3, lw=1.2, label="measured")
    ax_missing.set_xlabel("time (ms)")
    ax_missing.set_ylabel("bond_missing")
    ax_missing.set_ylim(-0.08, 0.8)
    ax_missing.legend(loc="upper right", ncol=2, fontsize=8)
    ax_missing.grid(True, alpha=0.25)

    fig.suptitle("BPM MPI bond migration: measured vs reference")
    fig.savefig(PLOT, dpi=180)
    plt.close(fig)
    print(f"Figure -> {PLOT}")


def graph() -> bool:
    rows = _load_rows()
    ok = validate(rows)
    plot(rows)
    return ok


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd == "start":
        start()
    elif cmd == "graph":
        sys.exit(0 if graph() else 1)
    elif cmd == "all":
        start()
        sys.exit(0 if graph() else 1)
    else:
        print("Usage: sweep.py [start|graph|all]   (no arg = all)")
        sys.exit(2)


if __name__ == "__main__":
    main()
