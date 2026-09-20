#!/usr/bin/env python3
"""Plot the fake-felted through-thickness tension history."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
CSV_PATH = HERE / "output" / "tension.csv"
FIGURE_PATH = HERE / "output" / "stress_strain.png"


with CSV_PATH.open(newline="") as stream:
    rows = list(csv.DictReader(stream))

pull = [row for row in rows if row["phase"] == "pull"]
if not pull:
    raise SystemExit(f"no pull samples in {CSV_PATH}")

# Internal fiber bonds are unbreakable in this test. A changed count therefore
# marks boundary loss, not physical damage. Split the curve at that point so
# the complete history remains visible without connecting across the jump.
initial_bond_count = int(pull[0]["bond_count"])
integrity_loss = next(
    (index for index, row in enumerate(pull) if int(row["bond_count"]) != initial_bond_count),
    len(pull),
)
valid_pull = pull[:integrity_loss]
post_loss_pull = pull[integrity_loss:]
if integrity_loss < len(pull):
    print(
        f"warning: bond count changed at step {pull[integrity_loss]['step']}; "
        "later data are plotted dashed"
    )

initial_gap = float(pull[0]["gap"]) - float(pull[0]["displacement"])
last_valid_step = int(valid_pull[-1]["step"])
last_valid_displacement = float(valid_pull[-1]["displacement"])
if len(valid_pull) >= 2:
    previous = valid_pull[-2]
    prescribed_increment_per_step = (
        last_valid_displacement - float(previous["displacement"])
    ) / (last_valid_step - int(previous["step"]))
else:
    prescribed_increment_per_step = 0.0


def applied_displacement(row: dict[str, str]) -> float:
    """Continuous prescribed separation, unaffected by deleted grip atoms."""
    step = int(row["step"])
    if step <= last_valid_step:
        return float(row["displacement"])
    return last_valid_displacement + prescribed_increment_per_step * (
        step - last_valid_step
    )

figure, (stress_axis, force_axis) = plt.subplots(1, 2, figsize=(11, 4.5))


def add_segment(rows: list[dict[str, str]], linestyle: str, label: str, alpha: float) -> None:
    displacements = [applied_displacement(row) for row in rows]
    strain_percent = [100.0 * displacement / initial_gap for displacement in displacements]
    stress_kpa = [1.0e-3 * float(row["engineering_stress"]) for row in rows]
    displacement_um = [1.0e6 * displacement for displacement in displacements]
    force_mn = [1.0e3 * float(row["tensile_force"]) for row in rows]
    stress_axis.plot(
        strain_percent,
        stress_kpa,
        color="#2468a2",
        linestyle=linestyle,
        linewidth=1.5,
        alpha=alpha,
        label=label,
    )
    force_axis.plot(
        displacement_um,
        force_mn,
        color="#c4492d",
        linestyle=linestyle,
        linewidth=1.5,
        alpha=alpha,
        label=label,
    )


add_segment(valid_pull, "-", "before boundary loss", 1.0)
if post_loss_pull:
    add_segment(post_loss_pull, "--", "after boundary loss", 0.75)

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

figure.suptitle("TANGLE fake felt: contact/friction/entanglement only")
figure.tight_layout()
figure.savefig(FIGURE_PATH, dpi=200)
print(FIGURE_PATH)
