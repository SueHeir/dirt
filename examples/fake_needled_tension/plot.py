#!/usr/bin/env python3
"""Plot needled tension and overlay the otherwise identical control."""

from __future__ import annotations

import csv
import argparse
from pathlib import Path

import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
mode = parser.add_mutually_exclusive_group()
mode.add_argument(
    "--slow",
    action="store_true",
    help="plot the paired one-percent, lower-rate runs",
)
mode.add_argument(
    "--fifty",
    action="store_true",
    help="plot the paired 50-percent, original-rate runs",
)
args = parser.parse_args()
output_name = (
    "output_50pct" if args.fifty else "output_slow" if args.slow else "output"
)
NEEDLED_CSV = HERE / output_name / "tension.csv"
CONTROL_CSV = HERE.parent / "fake_felted_tension" / output_name / "tension.csv"
FIGURE_PATH = HERE / output_name / "stress_strain.png"


def read_pull(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return [row for row in csv.DictReader(stream) if row["phase"] == "pull"]


needled = read_pull(NEEDLED_CSV)
if not needled:
    raise SystemExit(f"no pull samples in {NEEDLED_CSV}")

figure, (stress_axis, force_axis) = plt.subplots(1, 2, figsize=(11, 4.5))


def add_curve(
    rows: list[dict[str, str]],
    displacements: list[float],
    gauge_length: float,
    label: str,
    color: str,
    linestyle: str = "-",
    alpha: float = 1.0,
) -> None:
    strain_percent = [100.0 * displacement / gauge_length for displacement in displacements]
    stress_kpa = [1.0e-3 * float(row["engineering_stress"]) for row in rows]
    displacement_um = [1.0e6 * displacement for displacement in displacements]
    force_mn = [1.0e3 * float(row["tensile_force"]) for row in rows]
    stress_axis.plot(
        strain_percent,
        stress_kpa,
        label=label,
        color=color,
        linestyle=linestyle,
        linewidth=1.5,
        alpha=alpha,
    )
    force_axis.plot(
        displacement_um,
        force_mn,
        label=label,
        color=color,
        linestyle=linestyle,
        linewidth=1.5,
        alpha=alpha,
    )


def add_history(rows: list[dict[str, str]], label: str, color: str) -> None:
    """Plot prescribed loading continuously across any boundary-loss event."""
    # New runs report the internal, non-gripped gauge length explicitly.
    # Retain center-gap fallback so historical datasets remain plottable.
    gauge_length = float(
        rows[0].get("gauge_length")
        or float(rows[0]["gap"]) - float(rows[0]["displacement"])
    )
    initial_bonds = int(rows[0]["bond_count"])
    loss = next(
        (index for index, row in enumerate(rows) if int(row["bond_count"]) != initial_bonds),
        len(rows),
    )
    displacements = [float(row["displacement"]) for row in rows]
    if loss < len(rows):
        last_valid = rows[loss - 1]
        previous = rows[loss - 2]
        last_step = int(last_valid["step"])
        last_displacement = float(last_valid["displacement"])
        increment = (last_displacement - float(previous["displacement"])) / (
            last_step - int(previous["step"])
        )
        for index in range(loss, len(rows)):
            displacements[index] = last_displacement + increment * (
                int(rows[index]["step"]) - last_step
            )

    add_curve(rows[:loss], displacements[:loss], gauge_length, label, color)
    if loss < len(rows):
        add_curve(
            rows[loss:],
            displacements[loss:],
            gauge_length,
            f"{label} after boundary loss",
            color,
            linestyle="--",
            alpha=0.6,
        )


if CONTROL_CSV.exists():
    add_history(read_pull(CONTROL_CSV), "unneedled control", "#777777")
add_history(needled, "needled", "#8f3b76")

stress_axis.set_xlabel("Applied engineering strain, Δh/h₀ (%)")
stress_axis.set_ylabel("Through-thickness stress (kPa)")
stress_axis.set_title("Through-thickness stress–strain")
stress_axis.grid(alpha=0.25)
stress_axis.legend()

force_axis.set_xlabel("Applied grip separation increase (µm)")
force_axis.set_ylabel("Tensile force (mN)")
force_axis.set_title("Force–displacement")
force_axis.grid(alpha=0.25)
force_axis.legend()

rate_label = (
    "50% original-rate pull"
    if args.fifty
    else "1% lower-rate pull"
    if args.slow
    else "original-rate pull"
)
figure.suptitle(f"Paired fake felt: effect of needling ({rate_label})")
figure.tight_layout()
figure.savefig(FIGURE_PATH, dpi=200)
print(FIGURE_PATH)
