#!/usr/bin/env python3
"""Run and graph a live before/after ParticlesWith compatibility comparison.

The baseline is a detached worktree at the pinned pre-migration revision; the
candidate is the checkout containing this script.  Both produce their CSVs
from scratch, so this check detects a current numerical regression rather than
merely re-plotting an earlier result.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parents[1]
ROWS = HERE / "data" / "before_after.csv"
PLOT = HERE / "plots" / "before_after_compatibility.png"
BASELINE = "41cb2f116480dba2b4d023a5df5e143408033c2c"
TOLERANCE = 1.0e-15

CASES = (
    ("contact_wall", "bench_wall_activate_by_name", "examples/bench_wall_activate_by_name/config.toml", "examples/bench_wall_activate_by_name/data/wall_activate_by_name_results.csv", "particle_fz"),
    ("bond", "bond_fiber_tensile", "examples/bond_fiber_tensile/config.toml", "examples/bond_fiber_tensile/data/fiber_tensile.csv", "stress_mid"),
    ("clump", "bench_clump_inertia_sampler", "examples/bench_clump_inertia_sampler/config.toml", "examples/bench_clump_inertia_sampler/data/inertia_sampler.csv", "max_diag_rel_err"),
)


def run_case(repo: Path, example: str, config: str, output: str) -> Path:
    path = repo / output
    path.unlink(missing_ok=True)
    subprocess.run(
        ["cargo", "run", "--release", "--example", example,
         "--no-default-features", "--features", "precision-double", "--", config],
        cwd=repo,
        check=True,
    )
    if not path.exists():
        raise RuntimeError(f"{example} did not write {path}")
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def terminal_observable(path: Path, column: str) -> float:
    rows = list(csv.DictReader(path.open(newline="")))
    if not rows:
        raise RuntimeError(f"{path} has no data rows")
    return float(rows[-1][column])


def relative_difference(before: float, after: float) -> float:
    return abs(after - before) / max(abs(before), 1.0)


def run_comparison() -> list[dict[str, str]]:
    with tempfile.TemporaryDirectory(prefix="dirt-particleswith-baseline-") as tmp:
        baseline = Path(tmp) / "baseline"
        subprocess.run(
            ["git", "worktree", "add", "--detach", str(baseline), BASELINE],
            cwd=REPO,
            check=True,
        )
        try:
            rows = []
            for name, example, config, output, column in CASES:
                before_path = run_case(baseline, example, config, output)
                after_path = run_case(REPO, example, config, output)
                before = terminal_observable(before_path, column)
                after = terminal_observable(after_path, column)
                rows.append({
                    "case": name,
                    "observable": column,
                    "before": f"{before:.17e}",
                    "after": f"{after:.17e}",
                    "relative_difference": f"{relative_difference(before, after):.17e}",
                    "tolerance": f"{TOLERANCE:.17e}",
                    "before_sha256": sha256(before_path),
                    "after_sha256": sha256(after_path),
                })
            return rows
        finally:
            subprocess.run(["git", "worktree", "remove", "--force", str(baseline)], cwd=REPO, check=True)


def write_rows(rows: list[dict[str, str]]) -> None:
    ROWS.parent.mkdir(exist_ok=True)
    with ROWS.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def graph(rows: list[dict[str, str]]) -> bool:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    diffs = [float(row["relative_difference"]) for row in rows]
    tolerances = [float(row["tolerance"]) for row in rows]
    hashes_match = all(row["before_sha256"] == row["after_sha256"] for row in rows)
    ok = hashes_match and all(d <= t for d, t in zip(diffs, tolerances))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    x = range(len(rows))
    ax.bar(x, diffs, color="#2a9d8f", label="live measured relative difference")
    ax.scatter(x, tolerances, color="#d62828", marker="_", s=700,
               label="strict pass limit (1e-15)", zorder=3)
    ax.set_xticks(list(x), [row["case"].replace("_", "\n") for row in rows])
    ax.set_yscale("symlog", linthresh=1e-17)
    ax.set_ylim(-2e-17, 2e-14)
    ax.set_ylabel("before/after relative difference")
    ax.set_title("ParticlesWith migration: live output compatibility")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="upper right")
    for i, row in enumerate(rows):
        ax.text(i, 2e-17, "byte-identical" if row["before_sha256"] == row["after_sha256"] else "different",
                ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    PLOT.parent.mkdir(exist_ok=True)
    fig.savefig(PLOT, dpi=180)
    print(f"wrote {PLOT}")
    print(f"CSV hashes: {'PASS' if hashes_match else 'FAIL'}")
    print(f"compatibility: {sum(d <= t for d, t in zip(diffs, tolerances))}/{len(rows)} PASS")
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("run", "graph"), nargs="?", default="run")
    args = parser.parse_args()
    if args.mode == "run":
        rows = run_comparison()
        write_rows(rows)
    else:
        rows = list(csv.DictReader(ROWS.open()))
    return 0 if graph(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
