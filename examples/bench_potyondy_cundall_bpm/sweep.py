#!/usr/bin/env python3
"""Run and validate the Potyondy-Cundall BPM compression benchmark."""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
EX = ROOT / "examples" / "bench_potyondy_cundall_bpm"
DATA = EX / "data"
PLOTS = EX / "plots"

TARGET_PEAK_MPA = 199.1
TARGET_FAILURE_STRAIN = 0.00281
PEAK_TOL = 0.12
STRAIN_TOL = 0.18


def read_rows(path: Path, numeric=True):
    with path.open() as f:
        lines = (line for line in f if not line.lstrip().startswith("#"))
        if numeric:
            return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(lines) if row]
        return list(csv.DictReader(lines))


def run():
    env = os.environ.copy()
    env.setdefault("RUSTFLAGS", "")
    subprocess.run(
        [
            "cargo",
            "run",
            "--release",
            "--example",
            "bench_potyondy_cundall_bpm",
            "--no-default-features",
            "--features",
            "precision-double",
            "--",
            str(EX / "config.toml"),
        ],
        cwd=ROOT,
        env=env,
        check=True,
    )


def validate_and_plot():
    PLOTS.mkdir(parents=True, exist_ok=True)
    rows = read_rows(DATA / "dirt_stress_strain.csv")
    cracks = read_rows(DATA / "dirt_cracks.csv", numeric=False)
    for row in cracks:
        for key in ("step", "strain", "x_m", "y_m"):
            row[key] = float(row[key])
    target = read_rows(DATA / "potyondy_cundall_2004_fig8a_digitized.csv")
    peak = max(rows, key=lambda r: r["stress_pa"])
    first_crack = cracks[0]["strain"] if cracks else math.nan
    peak_ratio = (peak["stress_pa"] / 1.0e6) / TARGET_PEAK_MPA
    strain_ratio = peak["strain"] / TARGET_FAILURE_STRAIN
    sample_strain = rows[1]["strain"] - rows[0]["strain"]
    peak_ok = abs(peak_ratio - 1.0) <= PEAK_TOL
    strain_ok = abs(strain_ratio - 1.0) <= STRAIN_TOL
    crack_ok = len(cracks) >= 20 and first_crack <= peak["strain"] + sample_strain

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    ax0.plot(
        [r["strain"] / TARGET_FAILURE_STRAIN for r in rows],
        [r["stress_norm"] for r in rows],
        lw=2.2,
        label="DIRT live BPM compression",
    )
    ax0.plot(
        [r["strain_norm"] for r in target],
        [r["stress_norm"] for r in target],
        "o--",
        ms=4,
        label="Potyondy-Cundall Fig. 8(a), digitized",
    )
    ax0.axvspan(1 - STRAIN_TOL, 1 + STRAIN_TOL, color="C2", alpha=0.13, label="failure-strain gate")
    ax0.axhspan(1 - PEAK_TOL, 1 + PEAK_TOL, color="C1", alpha=0.13, label="peak-strength gate")
    ax0.set_xlabel("axial strain / target peak strain")
    ax0.set_ylabel("axial stress / Table 2 PFC2D qu")
    ax0.set_title("Stress-strain replication")
    ax0.legend(fontsize=8)
    ax0.grid(True, alpha=0.25)

    ax1.scatter([r["x_m"] * 1000 for r in cracks], [r["y_m"] * 1000 for r in cracks], s=13, c=[r["strain"] / TARGET_FAILURE_STRAIN for r in cracks], cmap="viridis")
    ax1.set_xlabel("x (mm)")
    ax1.set_ylabel("y (mm)")
    ax1.set_title(f"Bond-break progression ({len(cracks)} cracks)")
    ax1.grid(True, alpha=0.2)
    cbar = fig.colorbar(ax1.collections[0], ax=ax1)
    cbar.set_label("strain / target peak strain")

    fig.savefig(PLOTS / "stress_strain_and_cracks.png", dpi=180)

    print("Potyondy-Cundall BPM compression validation")
    print(f"  peak strength: {peak['stress_pa']/1e6:.2f} MPa / {TARGET_PEAK_MPA:.1f} MPa = {peak_ratio:.3f} (tol {PEAK_TOL:.0%}) {'PASS' if peak_ok else 'FAIL'}")
    print(f"  failure strain: {peak['strain']:.5f} / {TARGET_FAILURE_STRAIN:.5f} = {strain_ratio:.3f} (tol {STRAIN_TOL:.0%}) {'PASS' if strain_ok else 'FAIL'}")
    print(f"  crack progression: {len(cracks)} breaks, first at {first_crack:.5f} {'PASS' if crack_ok else 'FAIL'}")
    if not (peak_ok and strain_ok and crack_ok):
        raise SystemExit("VALIDATION FAILED")
    print("ALL CHECKS PASSED")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("all", "run"):
        run()
    if mode in ("all", "graph"):
        validate_and_plot()


if __name__ == "__main__":
    main()
