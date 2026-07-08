#!/usr/bin/env python3
"""Replicate Yang/Curtis wet fiber agglomerate breakage trends.

The benchmark runs a small DIRT agglomerate of bonded flexible fibers impacting
a frictionless plane with Willett pendular liquid bridges enabled. It gates the
same observables used by Yang et al. (AIChE J. 65(8), 2019): breakage ratio and
minimum largest-fragment mass ratio versus impact velocity and Weber number.
"""

from __future__ import annotations

import csv
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
PLOTS = HERE / "plots"
SWEEP = HERE / "sweep"

RADIUS = 5.0e-4
FIBERS = 3
BEADS = 5
BEAD_SPACING = 9.5e-4
SURFACE_TENSION = 0.072
RUPTURE = 4.0e-4
RHO = 1700.0
IMPACT_VELOCITIES = [2.0, 3.0, 4.0, 5.0]

TREND_MIN_R2 = 0.75
BREAKAGE_MIN_SPAN = 0.22
MASS_MIN_DROP = 0.28
MAX_OFFTREND = 1


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def generate_geometry(v: float = 1.0) -> None:
    DATA.mkdir(exist_ok=True)
    rows: list[tuple[float, float, float, float]] = []
    offsets = [
        (0.0, 0.0, 0.0),
        (0.0, 0.00118, 0.00045),
        (0.0, -0.00118, 0.00045),
    ]
    angles = [0.0 for _ in range(FIBERS)]
    base_z = 0.0015
    for f in range(FIBERS):
        c, s = math.cos(angles[f]), math.sin(angles[f])
        ox, oy, oz = offsets[f]
        for b in range(BEADS):
            x = (b - 0.5 * (BEADS - 1)) * BEAD_SPACING
            rows.append((ox + c * x, oy + s * x, base_z + oz, -v))
    with open(DATA / "agglomerate.csv", "w") as fh:
        fh.write("# x,y,z,vz\n")
        for x, y, z, vz in rows:
            fh.write(f"{x:.10e},{y:.10e},{z:.10e},{vz:.10e}\n")


def write_config(v: float) -> Path:
    generate_geometry(v)
    text = (HERE / "config.toml").read_text()
    out = SWEEP / f"v{v:.1f}"
    cfg = out / "config.toml"
    out.mkdir(parents=True, exist_ok=True)
    text = text.replace("velocity_z = -1.0", f"velocity_z = -{v:.6f}")
    text = text.replace(
        'dir = "examples/bench_curtis_wet_fiber_breakage/sweep/template"',
        f'dir = "examples/bench_curtis_wet_fiber_breakage/sweep/v{v:.1f}"',
    )
    cfg.write_text(text)
    return cfg


def start() -> None:
    generate_geometry()
    run([
        "cargo",
        "build",
        "--release",
        "--example",
        "bench_curtis_wet_fiber_breakage",
        "--no-default-features",
        "--features",
        "precision-double",
    ])
    for v in IMPACT_VELOCITIES:
        run([
            "cargo",
            "run",
            "--release",
            "--example",
            "bench_curtis_wet_fiber_breakage",
            "--no-default-features",
            "--features",
            "precision-double",
            "--",
            str(write_config(v).relative_to(ROOT)),
        ])


def summarize_case(v: float) -> dict[str, float]:
    path = SWEEP / f"v{v:.1f}" / "data" / "impact_metrics.csv"
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            rows.append({
                "step": float(row["step"]),
                "contacts": float(row["bridge_contacts"]),
                "breakage": float(row["breakage_ratio"]),
                "mass": float(row["largest_fragment_mass_ratio"]),
                "bonds_broken": float(row["bonds_broken"]),
            })
    if not rows:
        raise RuntimeError(f"no samples in {path}")
    n0 = rows[0]["contacts"]
    min_contacts = min(r["contacts"] for r in rows)
    min_mass = min(r["mass"] for r in rows)
    max_bonds_broken = max(r["bonds_broken"] for r in rows)
    return {
        "velocity": v,
        "weber": RHO * RADIUS * v * v / SURFACE_TENSION,
        "modified_weber": RHO * RADIUS * v * v / (SURFACE_TENSION * (RUPTURE / RADIUS)),
        "initial_contacts": n0,
        "min_contacts": min_contacts,
        "breakage_ratio": (n0 - min_contacts) / max(n0, 1.0),
        "min_largest_fragment_mass_ratio": min_mass,
        "bonds_broken": max_bonds_broken,
    }


def monotone_violations(values: list[float], increasing: bool) -> int:
    bad = 0
    for a, b in zip(values, values[1:]):
        if increasing and b < a - 0.04:
            bad += 1
        if not increasing and b > a + 0.04:
            bad += 1
    return bad


def linreg_r2(xs: list[float], ys: list[float]) -> float:
    xbar = sum(xs) / len(xs)
    ybar = sum(ys) / len(ys)
    sxx = sum((x - xbar) ** 2 for x in xs)
    syy = sum((y - ybar) ** 2 for y in ys)
    if sxx == 0.0 or syy == 0.0:
        return 0.0
    sxy = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    return (sxy * sxy) / (sxx * syy)


def graph() -> None:
    PLOTS.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    rows = [summarize_case(v) for v in IMPACT_VELOCITIES]
    with open(DATA / "velocity_sweep.csv", "w") as fh:
        fields = list(rows[0].keys())
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    v = [r["velocity"] for r in rows]
    we = [r["modified_weber"] for r in rows]
    br = [r["breakage_ratio"] for r in rows]
    mr = [r["min_largest_fragment_mass_ratio"] for r in rows]
    bonds = [r["bonds_broken"] for r in rows]

    br_r2 = linreg_r2(we, br)
    mr_r2 = linreg_r2(we, mr)
    br_span = max(br) - min(br)
    mr_drop = max(mr) - min(mr)
    br_bad = monotone_violations(br, True)
    mr_bad = monotone_violations(mr, False)
    bonds_active = max(bonds) > 0
    pass_gate = (
        br_r2 >= TREND_MIN_R2
        and mr_r2 >= TREND_MIN_R2
        and br_span >= BREAKAGE_MIN_SPAN
        and mr_drop >= MASS_MIN_DROP
        and br_bad <= MAX_OFFTREND
        and mr_bad <= MAX_OFFTREND
        and bonds_active
    )

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    axes[0].plot(v, br, "o-", label="DIRT wet BPM agglomerate")
    axes[0].set_xlabel("impact velocity (m/s)")
    axes[0].set_ylabel("breakage ratio")
    axes[0].set_ylim(-0.05, 1.02)
    axes[0].set_title("Breakage increases with impact speed")
    axes[0].text(0.03, 0.08, f"span={br_span:.2f}, R2(We*)={br_r2:.2f}", transform=axes[0].transAxes)
    axes[1].plot(v, mr, "o-", color="tab:red", label="DIRT wet BPM agglomerate")
    axes[1].set_xlabel("impact velocity (m/s)")
    axes[1].set_ylabel("minimum largest-fragment mass ratio")
    axes[1].set_ylim(-0.02, 1.05)
    axes[1].set_title("Largest fragment shrinks")
    axes[1].text(0.03, 0.08, f"drop={mr_drop:.2f}, R2(We*)={mr_r2:.2f}", transform=axes[1].transAxes)
    for ax in axes:
        ax.grid(True, alpha=0.25)
    fig.suptitle("Yang/Curtis wet fiber agglomerate impact trend gate")
    fig.tight_layout()
    fig.savefig(PLOTS / "breakage_vs_impact_velocity.png", dpi=180)
    plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(6.0, 4.0))
    ax1.plot(we, br, "o-", label="breakage ratio")
    ax1.set_xlabel("modified Weber number, We / S*")
    ax1.set_ylabel("breakage ratio")
    ax1.set_ylim(-0.05, 1.02)
    ax2 = ax1.twinx()
    ax2.plot(we, mr, "s-", color="tab:red", label="largest fragment")
    ax2.set_ylabel("minimum largest-fragment mass ratio")
    ax2.set_ylim(-0.02, 1.05)
    ax1.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(PLOTS / "weber_trend.png", dpi=180)
    plt.close(fig)

    print(
        "wet_fiber_breakage: "
        f"breakage_span={br_span:.3f} R2={br_r2:.3f} offtrend={br_bad}; "
        f"mass_drop={mr_drop:.3f} R2={mr_r2:.3f} offtrend={mr_bad}; "
        f"max_bonds_broken={max(bonds):.0f} -> {'PASS' if pass_gate else 'FAIL'}"
    )
    if not pass_gate:
        raise SystemExit(1)


def main() -> int:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd == "clean":
        shutil.rmtree(SWEEP, ignore_errors=True)
        return 0
    if cmd in ("all", "start"):
        start()
    if cmd in ("all", "graph"):
        graph()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
