#!/usr/bin/env python3
"""Exercise and graph the GRASS -> SOIL -> DIRT compatibility gate.

The two cases intentionally differ only in SOIL.  A reviewed source tuple must
pass both metadata and the non-MPI precision-double check.  The newer SOIL
tuple must be rejected and name the AtomData snapshot API drift.  Treating the
second failure as success is deliberate: it proves that the gate exposes drift
instead of resolving an unrelated locked Git revision.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
DATA = ROOT / "data"
PLOTS = ROOT / "plots"


def run(cmd: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=cwd, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, check=False)


def sha(repo: Path, revision: str) -> str:
    result = run(["git", "-C", str(repo), "rev-parse", revision], cwd=REPO)
    if result.returncode:
        raise RuntimeError(result.stdout)
    return result.stdout.strip()


def add_worktree(repo: Path, destination: Path, revision: str) -> None:
    result = run(["git", "-C", str(repo), "worktree", "add", "--detach",
                  str(destination), revision], cwd=REPO)
    if result.returncode:
        raise RuntimeError(result.stdout)


def remove_worktree(repo: Path, destination: Path) -> None:
    run(["git", "-C", str(repo), "worktree", "remove", "--force",
         str(destination)], cwd=REPO)


def record_case(name: str, grass: Path, soil: Path) -> dict[str, object]:
    result = run(["bash", "ci/ecosystem-head-check.sh", "--grass", str(grass),
                  "--soil", str(soil)], cwd=REPO)
    output = result.stdout
    return {
        "case": name,
        "dirt": sha(REPO, "HEAD"),
        "grass": sha(grass, "HEAD"),
        "soil": sha(soil, "HEAD"),
        "metadata": "PASS" if "running: cargo metadata" in output else "NOT REACHED",
        "check": "PASS" if result.returncode == 0 else "FAIL",
        "exit_code": result.returncode,
        "diagnostic": next((line.strip() for line in output.splitlines()
                            if "snapshot" in line.lower()), "(no snapshot diagnostic)"),
        "output": output,
    }


def validate(rows: list[dict[str, object]], expected: dict[str, str]) -> None:
    baseline, drift = rows
    failures = []
    if baseline["exit_code"] != 0 or baseline["metadata"] != "PASS":
        failures.append("compatible GRASS/SOIL/DIRT tuple did not pass metadata and check")
    if drift["exit_code"] == 0:
        failures.append("newer SOIL tuple unexpectedly passed; API drift was not exposed")
    diagnostic = str(drift["output"]).lower()
    for required in (expected["drift_diagnostic"].lower(), "dirt_granular", "dirt_bond"):
        if required not in diagnostic:
            failures.append(f"drift rejection did not name required diagnostic: {required}")
    if failures:
        raise SystemExit("CHECKS FAILED\n  " + "\n  ".join(failures))


def plot(rows: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    PLOTS.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(13, 5.2))
    ax.set_axis_off()
    labels = ["reviewed compatible tuple", "newer SOIL API tuple"]
    cells = []
    for label, row in zip(labels, rows):
        diagnostic = str(row["diagnostic"])
        if row["case"] == "soil-api-drift":
            # Keep the criterion itself visible in the committed figure rather
            # than relying on a truncated compiler sentence.
            diagnostic = "missing AtomData::snapshot\nin dirt_granular + dirt_bond"
        elif len(diagnostic) > 72:
            diagnostic = diagnostic[:69] + "..."
        cells.append([label, str(row["dirt"])[:12], str(row["grass"])[:12], str(row["soil"])[:12],
                      str(row["metadata"]), str(row["check"]), diagnostic])
    table = ax.table(cellText=cells,
                     colLabels=["case", "DIRT commit", "GRASS commit", "SOIL commit", "metadata", "precision-double check", "visible API diagnostic"],
                     cellLoc="left", loc="center", colWidths=[.17, .11, .11, .11, .10, .16, .24])
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 2.25)
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor("#1f4e79"); cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        elif row == 1:
            cell.set_facecolor("#d9ead3")
        elif row == 2:
            cell.set_facecolor("#f4cccc")
    ax.set_title("GRASS → SOIL → DIRT HEAD compatibility matrix\nPASS requires metadata + non-MPI precision-double check; deliberate SOIL drift must fail visibly", pad=20, weight="bold")
    fig.tight_layout()
    fig.savefig(PLOTS / "ecosystem_head_compatibility_matrix.png", dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recorded", action="store_true", help="validate and plot committed matrix data only")
    args = parser.parse_args()
    with (ROOT / "config.toml").open("rb") as handle:
        config = tomllib.load(handle)
    if args.recorded:
        rows = json.loads((DATA / "matrix.json").read_text())
    else:
        sources = config["sources"]
        grass_repo = Path(sources["grass_repo"]).expanduser()
        soil_repo = Path(sources["soil_repo"]).expanduser()
        with tempfile.TemporaryDirectory(prefix="dirt-ecosystem-matrix-") as tmp:
            tmp = Path(tmp)
            grass, compatible, drift = tmp / "grass", tmp / "soil-compatible", tmp / "soil-drift"
            try:
                add_worktree(grass_repo, grass, sources["grass_commit"])
                add_worktree(soil_repo, compatible, sources["soil_compatible_commit"])
                add_worktree(soil_repo, drift, sources["soil_drift_commit"])
                rows = [record_case("compatible", grass, compatible),
                        record_case("soil-api-drift", grass, drift)]
            finally:
                remove_worktree(grass_repo, grass)
                remove_worktree(soil_repo, compatible)
                remove_worktree(soil_repo, drift)
        DATA.mkdir(exist_ok=True)
        (DATA / "matrix.json").write_text(json.dumps(rows, indent=2) + "\n")
    validate(rows, config["expect"])
    plot(rows)
    print("VALIDATION: PASS")
    print("ALL CHECKS PASSED: compatible tuple passed; newer SOIL tuple exposed AtomData snapshot API drift.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
