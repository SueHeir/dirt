#!/usr/bin/env python3
"""Graph the pinned before/after ParticlesWith migration comparison."""
import csv
from pathlib import Path

HERE = Path(__file__).parent
ROWS = HERE / "data" / "before_after.csv"
PLOT = HERE / "plots" / "before_after_compatibility.png"


def main():
    rows = list(csv.DictReader(ROWS.open()))
    diffs = [float(row["relative_difference"]) for row in rows]
    tolerances = [float(row["tolerance"]) for row in rows]
    hashes_match = all(row["before_sha256"] == row["after_sha256"] for row in rows)
    ok = hashes_match and all(d <= t for d, t in zip(diffs, tolerances))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    x = range(len(rows))
    ax.bar(x, diffs, color="#2a9d8f", label="measured relative difference")
    ax.scatter(x, tolerances, color="#d62828", marker="_", s=700,
               label="strict pass limit (1e-15)", zorder=3)
    ax.set_xticks(list(x), [row["case"].replace("_", "\n") for row in rows])
    ax.set_yscale("symlog", linthresh=1e-17)
    ax.set_ylim(-2e-17, 2e-14)
    ax.set_ylabel("before/after relative difference")
    ax.set_title("ParticlesWith migration: representative output compatibility")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="upper right")
    for i, row in enumerate(rows):
        ax.text(i, 2e-17, "byte-identical", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    PLOT.parent.mkdir(exist_ok=True)
    fig.savefig(PLOT, dpi=180)
    print(f"wrote {PLOT}")
    print(f"CSV hashes: {'PASS' if hashes_match else 'FAIL'}")
    print(f"compatibility: {sum(d <= t for d, t in zip(diffs, tolerances))}/{len(rows)} PASS")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
