#!/usr/bin/env python3
"""Plot hopper quiescence validation results from generated stats CSV files."""

from __future__ import annotations

import argparse
import csv
import os
import tomllib
from pathlib import Path


EXAMPLE_DIR = Path(__file__).resolve().parent
HEIGHT_TOL_M = 0.001
DISCHARGE_TOL_FRAC = 0.01


def read_runs(config: Path) -> tuple[int, int]:
    with config.open("rb") as f:
        data = tomllib.load(f)
    runs = data.get("run", [])
    if len(runs) < 2:
        raise ValueError(f"{config} must declare filling and flowing runs")
    fill_steps = int(runs[0]["steps"])
    flow_steps = int(runs[1]["steps"])
    return fill_steps, fill_steps + flow_steps


def read_stats(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as f:
        rows = []
        for row in csv.DictReader(f):
            rows.append({k: float(v) for k, v in row.items()})
    if not rows:
        raise ValueError(f"{path} is empty")
    return rows


def row_at_or_before(rows: list[dict[str, float]], step: int) -> dict[str, float]:
    eligible = [r for r in rows if int(r["step"]) <= step]
    if not eligible:
        raise ValueError(f"no stats row at or before step {step}")
    return eligible[-1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--prefix",
        default="val",
        help="config/stats prefix to plot: val or config (default: val)",
    )
    parser.add_argument(
        "--out",
        default=EXAMPLE_DIR / "plots" / "hopper_quiescence_validation.png",
        type=Path,
    )
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    variants = ("baseline", "coherence")
    rows = {}
    stages = {}
    for variant in variants:
        config = EXAMPLE_DIR / f"{args.prefix}_{variant}.toml"
        stats = EXAMPLE_DIR / f"{args.prefix}_{variant}_stats.csv"
        stages[variant] = read_runs(config)
        rows[variant] = read_stats(stats)

    fill_step = stages["baseline"][0]
    if stages["coherence"][0] != fill_step:
        raise ValueError("baseline/coherence fill stages differ; choose a matched config pair")

    base_fill = row_at_or_before(rows["baseline"], fill_step)
    coh_fill = row_at_or_before(rows["coherence"], fill_step)
    base_final = row_at_or_before(rows["baseline"], stages["baseline"][1])
    coh_final = row_at_or_before(rows["coherence"], stages["coherence"][1])
    n_ref = max(base_final["n_atoms"], coh_final["n_atoms"], 1.0)

    height_ref = base_fill["top_z"]
    height_delta = coh_fill["top_z"] - height_ref
    height_pass = abs(height_delta) <= HEIGHT_TOL_M

    def discharge_series(variant: str) -> tuple[list[float], list[float]]:
        xs = []
        ys = []
        for row in rows[variant]:
            step = int(row["step"])
            if fill_step <= step <= stages[variant][1]:
                xs.append((step - fill_step) / 1000.0)
                ys.append(row["n_discharged"] / n_ref)
        return xs, ys

    bx, by = discharge_series("baseline")
    cx, cy = discharge_series("coherence")
    discharge_pass = True
    baseline_by_step = {
        int(row["step"]): row["n_discharged"] / n_ref
        for row in rows["baseline"]
        if fill_step <= int(row["step"]) <= stages["baseline"][1]
    }
    for row in rows["coherence"]:
        step = int(row["step"])
        if step in baseline_by_step and fill_step <= step <= stages["coherence"][1]:
            diff = abs(row["n_discharged"] / n_ref - baseline_by_step[step])
            discharge_pass &= diff <= DISCHARGE_TOL_FRAC

    phases = {
        "baseline": (base_fill["elapsed_s"], base_final["elapsed_s"] - base_fill["elapsed_s"]),
        "coherence": (coh_fill["elapsed_s"], coh_final["elapsed_s"] - coh_fill["elapsed_s"]),
    }

    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 10})
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    ax = axes[0][0]
    ax.plot(bx, by, label="baseline", color="#365f91", linewidth=2)
    ax.fill_between(
        bx,
        [max(0.0, y - DISCHARGE_TOL_FRAC) for y in by],
        [min(1.05, y + DISCHARGE_TOL_FRAC) for y in by],
        color="#365f91",
        alpha=0.16,
        linewidth=0,
        label="baseline ±1% pass band",
    )
    ax.plot(cx, cy, label="coherence", color="#b45f06", linewidth=2)
    ax.set_title("Discharge curve")
    ax.set_xlabel("flow-stage step / 1000")
    ax.set_ylabel("discharged fraction")
    ax.set_ylim(-0.02, 1.05)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right")
    ax.text(
        0.03,
        0.92,
        "PASS" if discharge_pass else "FAIL",
        transform=ax.transAxes,
        weight="bold",
        color="#1b7837" if discharge_pass else "#b2182b",
    )

    ax = axes[0][1]
    ax.axhspan(
        height_ref - HEIGHT_TOL_M,
        height_ref + HEIGHT_TOL_M,
        color="#365f91",
        alpha=0.16,
        label="baseline ±1 mm pass band",
    )
    ax.axhline(height_ref, color="#365f91", linewidth=1.2, label="baseline reference")
    ax.scatter(["baseline", "coherence"], [height_ref, coh_fill["top_z"]], s=70, color=["#365f91", "#b45f06"])
    ax.set_title("Settled fill height")
    ax.set_ylabel("top_z at fill end [m]")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best")
    ax.text(
        0.03,
        0.92,
        f"{'PASS' if height_pass else 'FAIL'}  Δ={height_delta * 1000:.2f} mm",
        transform=ax.transAxes,
        weight="bold",
        color="#1b7837" if height_pass else "#b2182b",
    )

    ax = axes[1][0]
    labels = ["fill", "flow"]
    x = range(len(labels))
    width = 0.35
    ax.bar([i - width / 2 for i in x], phases["baseline"], width, label="baseline", color="#365f91")
    ax.bar([i + width / 2 for i in x], phases["coherence"], width, label="coherence", color="#b45f06")
    ax.set_xticks(list(x), labels)
    ax.set_title("Phase wall times")
    ax.set_ylabel("elapsed wall time [s]")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best")

    ax = axes[1][1]
    totals = {variant: sum(phases[variant]) for variant in variants}
    speedup = totals["baseline"] / totals["coherence"] if totals["coherence"] > 0 else 0.0
    skipped = [
        row_at_or_before(rows["baseline"], fill_step)["skipped_pairs"],
        row_at_or_before(rows["coherence"], fill_step)["skipped_pairs"],
        coh_final["skipped_pairs"],
    ]
    pairs = [
        max(row_at_or_before(rows["baseline"], fill_step)["pairs"], 1.0),
        max(row_at_or_before(rows["coherence"], fill_step)["pairs"], 1.0),
        max(coh_final["pairs"], 1.0),
    ]
    ax.bar(["base fill", "coh fill", "coh final"], [100.0 * s / p for s, p in zip(skipped, pairs)], color=["#365f91", "#b45f06", "#b45f06"])
    ax.set_title(f"Skipped pair fraction (total speedup {speedup:.2f}x)")
    ax.set_ylabel("skipped pairs / pairs [%]")
    ax.set_ylim(0, 105)
    ax.grid(True, axis="y", alpha=0.25)

    fig.suptitle(f"Hopper quiescence validation from {args.prefix}_*_stats.csv")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    print(f"Wrote {args.out}")
    print(
        "RESULT: "
        f"height {'PASS' if height_pass else 'FAIL'} "
        f"(abs delta {abs(height_delta) * 1000:.2f} mm <= {HEIGHT_TOL_M * 1000:.1f} mm); "
        f"discharge {'PASS' if discharge_pass else 'FAIL'} "
        f"(coherence within ±{DISCHARGE_TOL_FRAC:.0%} of baseline curve)"
    )
    return 0 if height_pass and discharge_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
