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

REFERENCE_POINTS = DATA / "yang_curtis_fig13_digitized.csv"
BREAKAGE_TOL = 0.30
MASS_TOL = 0.40


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def generate_geometry(v: float = 1.0, path: Path | None = None) -> None:
    path = path or (DATA / "agglomerate.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
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
    with open(path, "w") as fh:
        fh.write("# x,y,z,vz\n")
        for x, y, z, vz in rows:
            fh.write(f"{x:.10e},{y:.10e},{z:.10e},{vz:.10e}\n")


def write_config(v: float) -> Path:
    text = (HERE / "config.toml").read_text()
    out = SWEEP / f"v{v:.1f}"
    geometry = out / "agglomerate.csv"
    generate_geometry(v, geometry)
    cfg = out / "config.toml"
    out.mkdir(parents=True, exist_ok=True)
    text = text.replace("velocity_z = -1.0", f"velocity_z = -{v:.6f}")
    text = text.replace(
        'file = "examples/bench_curtis_wet_fiber_breakage/data/agglomerate.csv"',
        f'file = "examples/bench_curtis_wet_fiber_breakage/sweep/v{v:.1f}/agglomerate.csv"',
    )
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


def load_reference_points() -> list[dict[str, float | str]]:
    points: list[dict[str, float | str]] = []
    with open(REFERENCE_POINTS) as fh:
        for row in csv.DictReader(fh):
            points.append({
                "modified_weber": float(row["modified_weber"]),
                "breakage_ratio": float(row["breakage_ratio"]),
                "min_largest_fragment_mass_ratio": float(row["min_largest_fragment_mass_ratio"]),
                "source": row["source"],
            })
    if not points:
        raise RuntimeError(f"no reference points in {REFERENCE_POINTS}")
    return sorted(points, key=lambda r: r["modified_weber"])


def interp_reference(
    modified_weber: float, points: list[dict[str, float | str]]
) -> tuple[float, float]:
    if modified_weber < points[0]["modified_weber"] or modified_weber > points[-1]["modified_weber"]:
        raise RuntimeError(
            f"We*={modified_weber:.1f} is outside the committed Yang/Curtis Fig. 13 curve"
        )
    for lo, hi in zip(points, points[1:]):
        if lo["modified_weber"] <= modified_weber <= hi["modified_weber"]:
            span = hi["modified_weber"] - lo["modified_weber"]
            t = 0.0 if span == 0.0 else (modified_weber - lo["modified_weber"]) / span
            br = lo["breakage_ratio"] + t * (hi["breakage_ratio"] - lo["breakage_ratio"])
            mr = lo["min_largest_fragment_mass_ratio"] + t * (
                hi["min_largest_fragment_mass_ratio"] - lo["min_largest_fragment_mass_ratio"]
            )
            return br, mr
    last = points[-1]
    return last["breakage_ratio"], last["min_largest_fragment_mass_ratio"]


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
        writer = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    try:
        import matplotlib
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "matplotlib is required for graphing; run with the repo bench Python, "
            "for example: source ~/projects/.build-env && "
            "$BENCH_PYTHON examples/bench_curtis_wet_fiber_breakage/sweep.py"
        ) from exc

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ref_points = load_reference_points()
    for row in rows:
        br_ref, mr_ref = interp_reference(row["modified_weber"], ref_points)
        row["breakage_ref"] = br_ref
        row["mass_ref"] = mr_ref
        row["breakage_abs_err"] = abs(row["breakage_ratio"] - br_ref)
        row["mass_abs_err"] = abs(row["min_largest_fragment_mass_ratio"] - mr_ref)
        row["reference_pass"] = 1.0 if (
            row["breakage_abs_err"] <= BREAKAGE_TOL and row["mass_abs_err"] <= MASS_TOL
        ) else 0.0
    with open(DATA / "reference_comparison.csv", "w") as fh:
        fields = [
            "velocity",
            "modified_weber",
            "breakage_ratio",
            "breakage_ref",
            "breakage_abs_err",
            "min_largest_fragment_mass_ratio",
            "mass_ref",
            "mass_abs_err",
            "reference_pass",
        ]
        writer = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})

    v = [r["velocity"] for r in rows]
    we = [r["modified_weber"] for r in rows]
    br = [r["breakage_ratio"] for r in rows]
    mr = [r["min_largest_fragment_mass_ratio"] for r in rows]
    br_ref = [r["breakage_ref"] for r in rows]
    mr_ref = [r["mass_ref"] for r in rows]
    bonds = [r["bonds_broken"] for r in rows]

    br_r2 = linreg_r2(we, br)
    mr_r2 = linreg_r2(we, mr)
    br_bad = monotone_violations(br, True)
    mr_bad = monotone_violations(mr, False)
    bonds_active = max(bonds) > 0
    reference_passes = int(sum(r["reference_pass"] for r in rows))
    pass_gate = reference_passes == len(rows) and bonds_active and br_bad == 0 and mr_bad == 0

    def velocity_from_modified_weber(westar: float) -> float:
        hstar = RUPTURE / RADIUS
        return math.sqrt(max(westar, 0.0) * SURFACE_TENSION * hstar / (RHO * RADIUS))

    ref_we = [p["modified_weber"] for p in ref_points]
    ref_br = [p["breakage_ratio"] for p in ref_points]
    ref_mr = [p["min_largest_fragment_mass_ratio"] for p in ref_points]
    ref_v = [velocity_from_modified_weber(w) for w in ref_we]

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    axes[0].plot(ref_v, ref_br, "k--", lw=1.1, label="Yang/Curtis Fig. 13")
    axes[1].plot(ref_v, ref_mr, "k--", lw=1.1, label="Yang/Curtis Fig. 13")
    axes[0].plot(v, br, "o-", label="DIRT wet BPM agglomerate")
    axes[0].set_xlabel("impact velocity (m/s)")
    axes[0].set_ylabel("breakage ratio")
    axes[0].set_ylim(-0.05, 1.02)
    axes[0].set_title("Breakage vs Yang/Curtis Fig. 13")
    axes[1].plot(v, mr, "o-", color="tab:red", label="DIRT wet BPM agglomerate")
    axes[1].set_xlabel("impact velocity (m/s)")
    axes[1].set_ylabel("minimum largest-fragment mass ratio")
    axes[1].set_ylim(-0.02, 1.05)
    axes[1].set_title("Largest fragment vs Yang/Curtis Fig. 13")
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7)
    fig.suptitle("Yang/Curtis wet fiber agglomerate impact trend")
    fig.tight_layout()
    fig.savefig(PLOTS / "breakage_vs_impact_velocity.png", dpi=180)
    plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(6.8, 4.3))
    ax1.plot(ref_we, ref_br, "k--", lw=1.1, label="Yang/Curtis Fig. 13 breakage")
    ax1.plot(we, br, "o-", label="breakage ratio")
    ax1.set_xlabel("modified Weber number, We / S*")
    ax1.set_ylabel("breakage ratio")
    ax1.set_ylim(-0.05, 1.02)
    ax2 = ax1.twinx()
    ax2.plot(ref_we, ref_mr, "k:", lw=1.1, label="Yang/Curtis Fig. 13 mass")
    ax2.plot(we, mr, "s-", color="tab:red", label="largest fragment")
    ax2.set_ylabel("minimum largest-fragment mass ratio")
    ax2.set_ylim(-0.02, 1.05)
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="upper left", fontsize=7)
    ax2.legend(loc="upper right", fontsize=7)
    fig.tight_layout()
    fig.savefig(PLOTS / "weber_trend.png", dpi=180)
    plt.close(fig)

    print(
        "wet_fiber_breakage: "
        f"yang_curtis_reference={reference_passes}/{len(rows)}; "
        f"breakage_R2={br_r2:.3f} offtrend={br_bad}; "
        f"mass_R2={mr_r2:.3f} offtrend={mr_bad}; "
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
