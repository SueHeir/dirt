#!/usr/bin/env python3
"""Run and plot the BPM cantilever example.

The simulation writes ``data/cantilever.csv``.  ``graph`` reads that CSV,
compares the latest tip deflection with the Euler-Bernoulli cantilever
reference used by the README, and writes the committed validation figure.
"""

from __future__ import annotations

import csv
import math
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
CSV = HERE / "data" / "cantilever.csv"
PLOT = HERE / "plots" / "tip_deflection_vs_beam.png"

N_ATOMS = 10
N_BONDS = 9
RADIUS = 1.0e-3
DENSITY = 2500.0
GRAVITY = 9.81
YOUNGS_MODULUS = 1.0e9
SPAN = 18.0e-3
TOL = 0.05


def beam_tip_deflection() -> float:
    """Euler-Bernoulli tip deflection for the committed 10-sphere chain."""
    sphere_mass = (4.0 / 3.0) * math.pi * RADIUS**3 * DENSITY
    q = N_ATOMS * sphere_mass * GRAVITY / SPAN
    i_area = math.pi * RADIUS**4 / 4.0
    return -(q * SPAN**4) / (8.0 * YOUNGS_MODULUS * i_area)


def start() -> None:
    cmd = [
        "cargo",
        "run",
        "--release",
        "--example",
        "bond_cantilever",
        "--no-default-features",
        "--features",
        "precision-double",
        "--",
        "examples/bond_cantilever/config.toml",
    ]
    subprocess.run(cmd, cwd=REPO, check=True)


def read_rows() -> list[dict[str, float]]:
    if not CSV.exists():
        raise SystemExit(
            f"missing {CSV}; run `python3 examples/bond_cantilever/sweep.py start` first"
        )
    with CSV.open(newline="") as f:
        rows = []
        for row in csv.DictReader(f):
            rows.append(
                {
                    "step": int(row["step"]),
                    "t": float(row["t"]),
                    "tip_z": float(row["tip_z"]),
                    "max_strain": float(row["max_strain"]),
                    "bond_count": int(row["bond_count"]),
                    "bond_missing": int(row["bond_missing"]),
                }
            )
    if not rows:
        raise SystemExit(f"{CSV} has no samples")
    return rows


def graph() -> bool:
    rows = read_rows()
    ref = beam_tip_deflection()
    last = rows[-1]
    rel_err = abs(last["tip_z"] - ref) / abs(ref)
    all_bonds = all(r["bond_count"] == N_BONDS for r in rows)
    no_missing = all(r["bond_missing"] == 0 for r in rows)
    passed = rel_err <= TOL and all_bonds and no_missing

    print("=== BPM cantilever beam-theory check ===")
    print(f"  samples                 : {len(rows)}")
    print(f"  latest step             : {last['step']}")
    print(f"  measured tip_z          : {last['tip_z']:+.6e} m")
    print(f"  Euler-Bernoulli tip_z   : {ref:+.6e} m")
    print(f"  relative error          : {rel_err:.3%}")
    print(f"  tolerance               : {TOL:.0%}")
    print(f"  bonds present           : {N_BONDS if all_bonds else 'not all samples'}")
    print(f"  missing partner skips   : {max(r['bond_missing'] for r in rows)}")
    print(f"  status                  : {'PASS' if passed else 'FAIL'}")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    PLOT.parent.mkdir(parents=True, exist_ok=True)
    ts_ms = [r["t"] * 1e3 for r in rows]
    tips_um = [r["tip_z"] * 1e6 for r in rows]
    ref_um = ref * 1e6
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.plot(ts_ms, tips_um, lw=1.8, color="C0", label="DIRT tip deflection")
    ax.axhline(ref_um, color="C1", lw=1.8, ls="--", label="Euler-Bernoulli reference")
    ax.plot(
        [last["t"] * 1e3],
        [last["tip_z"] * 1e6],
        "o",
        ms=7,
        color="C0",
        mec="white",
        mew=1.0,
        label=f"latest sample ({rel_err:.2%} error)",
    )
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("tip deflection z (um)")
    ax.set_title("Bonded cantilever tip deflection vs beam theory")
    ax.legend(loc="best", framealpha=0.95)
    ax.grid(True, alpha=0.25)
    fig.savefig(PLOT, dpi=160, bbox_inches="tight")
    print(f"  plot saved              : {PLOT}")
    return passed


def main() -> int:
    action = sys.argv[1] if len(sys.argv) > 1 else "all"
    if action == "start":
        start()
        return 0
    if action == "graph":
        return 0 if graph() else 1
    if action == "all":
        start()
        return 0 if graph() else 1
    print("Usage: sweep.py [all|start|graph]", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
