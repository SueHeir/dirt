#!/usr/bin/env python3
"""Replicate Guo/Curtis flexible-fiber cantilever bending vs beam theory."""

from __future__ import annotations

import csv
import math
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CONFIG = HERE / "config.toml"
RUNS = HERE / "runs"
DATA = HERE / "data"
PLOTS = HERE / "plots"

E_BOND = 1.0e9
R_BOND = 1.0e-3
I_BOND = 0.25 * math.pi * R_BOND**4
EI = E_BOND * I_BOND
SPAN = 0.036
LOAD_NORM = [0.05, 0.10, 0.15, 0.20, 0.45]
TARGET_NORM_LOAD = 0.45
TIP_TOL = 0.03
PROFILE_TOL = 0.03
MOMENT_TOL = 0.03


def read_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def write_case_config(norm_load: float, force: float, out_dir: Path) -> Path:
    text = CONFIG.read_text()
    text = text.replace('fz = -0.1', f'fz = {-force:.12e}')
    text = text.replace(
        'dir = "examples/bench_curtis_cantilever/runs/template"',
        f'dir = "{out_dir.as_posix()}"',
    )
    path = out_dir / "config.toml"
    out_dir.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def run_case(norm_load: float) -> dict[str, float | str]:
    force = norm_load * EI / SPAN**2
    out_dir = RUNS / f"load_{norm_load:.2f}".replace(".", "p")
    if out_dir.exists():
        shutil.rmtree(out_dir)
    cfg = write_case_config(norm_load, force, out_dir)
    cmd = [
        "cargo",
        "run",
        "--release",
        "--example",
        "bench_curtis_cantilever",
        "--no-default-features",
        "--features",
        "precision-double",
        "--",
        str(cfg),
    ]
    print(f"running normalized load {norm_load:.2f} (F={force:.6e} N)")
    subprocess.run(cmd, cwd=ROOT, check=True)

    rows = read_csv(out_dir / "data" / "fiber_bond.csv")
    tail = rows[max(0, int(0.9 * len(rows))) :]
    tip_z = sum(r["right_z"] for r in tail) / len(tail)
    tip_norm = -tip_z / SPAN
    theory_norm = norm_load / 3.0
    rel_err = abs(tip_norm - theory_norm) / abs(theory_norm)
    last = rows[-1]
    return {
        "norm_load": norm_load,
        "force": force,
        "tip_norm": tip_norm,
        "theory_tip_norm": theory_norm,
        "tip_rel_err": rel_err,
        "bond_count": last["bond_count"],
        "bonds_broken": last["bonds_broken"],
        "out_dir": out_dir.as_posix(),
    }


def profile_errors(case: dict[str, float | str]) -> tuple[float, float, list[dict[str, float]]]:
    out_dir = Path(str(case["out_dir"]))
    force = float(case["force"])
    profile = read_csv(out_dir / "data" / "profile.csv")
    bonds = read_csv(out_dir / "data" / "bond_profile.csv")
    ymax = force * SPAN**3 / (3.0 * EI)
    mmax = force * SPAN

    profile_rows = []
    max_defl_err = 0.0
    for row in sorted(profile, key=lambda r: r["x0"]):
        x = row["x0"]
        xi = x / SPAN
        dem = -row["z"] / ymax
        theory = xi**2 * (3.0 - xi) / 2.0
        err = abs(dem - theory)
        max_defl_err = max(max_defl_err, err)
        profile_rows.append(
            {
                "kind": "deflection",
                "x_over_span": xi,
                "dem": dem,
                "theory": theory,
                "abs_err": err,
            }
        )

    max_moment_err = 0.0
    for row in sorted(bonds, key=lambda r: r["x0_center"]):
        xi = row["x0_center"] / SPAN
        dem = abs(row["m_bend"]) / mmax
        theory = 1.0 - xi
        err = abs(dem - theory)
        max_moment_err = max(max_moment_err, err)
        profile_rows.append(
            {
                "kind": "moment",
                "x_over_span": xi,
                "dem": dem,
                "theory": theory,
                "abs_err": err,
            }
        )
    return max_defl_err, max_moment_err, profile_rows


def write_outputs(results: list[dict[str, float | str]], profile_rows: list[dict[str, float]]) -> None:
    DATA.mkdir(exist_ok=True)
    PLOTS.mkdir(exist_ok=True)
    with (DATA / "results.csv").open("w", newline="") as f:
        fields = [
            "norm_load",
            "force",
            "tip_norm",
            "theory_tip_norm",
            "tip_rel_err",
            "bond_count",
            "bonds_broken",
        ]
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for r in results:
            w.writerow({k: r[k] for k in fields})
    with (DATA / "target_profile.csv").open("w", newline="") as f:
        fields = ["kind", "x_over_span", "dem", "theory", "abs_err"]
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(profile_rows)

    import matplotlib.pyplot as plt

    x = [float(r["norm_load"]) for r in results]
    y = [float(r["tip_norm"]) for r in results]
    yt = [float(r["theory_tip_norm"]) for r in results]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    xx = [i * max(x) / 200.0 for i in range(201)]
    yy = [v / 3.0 for v in xx]
    ax.fill_between(xx, [v * (1 - TIP_TOL) for v in yy], [v * (1 + TIP_TOL) for v in yy],
                    color="0.88", label=f"+/-{TIP_TOL:.0%} gate")
    ax.plot(xx, yy, "-", c="C1", lw=2.0, label="Euler-Bernoulli thin-beam theory")
    ax.plot(x, y, "o", ms=7, mfc="white", mec="C0", mew=1.6, label="DIRT bonded-sphere fiber")
    ax.set_xlabel("normalized load  F (L-rs)^2 / EI")
    ax.set_ylabel("normalized free-end deflection  |y0| / (L-rs)")
    ax.set_title("Flexible-fiber cantilever bending")
    ax.legend(loc="upper left", framealpha=0.95)
    fig.savefig(PLOTS / "tip_deflection_vs_load.png", dpi=160, bbox_inches="tight")

    defl = [r for r in profile_rows if r["kind"] == "deflection"]
    mom = [r for r in profile_rows if r["kind"] == "moment"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, rows, ylabel, title, tol in [
        (axes[0], defl, "normalized deflection  |y|/|ymax|", "Deflection distribution", PROFILE_TOL),
        (axes[1], mom, "normalized bending moment  |M|/Mmax", "Bending-moment distribution", MOMENT_TOL),
    ]:
        xs = [r["x_over_span"] for r in rows]
        dem = [r["dem"] for r in rows]
        th = [r["theory"] for r in rows]
        ax.fill_between(xs, [v - tol for v in th], [v + tol for v in th],
                        color="0.88", label=f"+/-{tol:.0%} absolute gate")
        ax.plot(xs, th, "-", c="C1", lw=2.0, label="beam theory")
        ax.plot(xs, dem, "o", ms=6, mfc="white", mec="C0", mew=1.4, label="DIRT")
        ax.set_xlabel("position along fiber  x / (L-rs)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlim(-0.02, 1.02)
    axes[0].legend(loc="upper left", framealpha=0.95)
    axes[1].legend(loc="upper right", framealpha=0.95)
    fig.suptitle(f"Guo/Curtis Fig. 4 analog at normalized load {TARGET_NORM_LOAD:.2f}")
    fig.savefig(PLOTS / "moment_deflection_profiles.png", dpi=160, bbox_inches="tight")


def main() -> int:
    results = [run_case(n) for n in LOAD_NORM]
    target = min(results, key=lambda r: abs(float(r["norm_load"]) - TARGET_NORM_LOAD))
    max_defl_err, max_moment_err, profile_rows = profile_errors(target)
    write_outputs(results, profile_rows)

    max_tip_err = max(float(r["tip_rel_err"]) for r in results)
    broken = sum(float(r["bonds_broken"]) for r in results)
    print("=== Guo/Curtis cantilever validation ===")
    print(f"max tip curve relative error : {max_tip_err:.3%} (gate {TIP_TOL:.0%})")
    print(f"max deflection profile error : {max_defl_err:.3%} absolute (gate {PROFILE_TOL:.0%})")
    print(f"max moment profile error     : {max_moment_err:.3%} absolute (gate {MOMENT_TOL:.0%})")
    print(f"broken bonds                 : {broken:.0f}")
    ok = (
        max_tip_err <= TIP_TOL
        and max_defl_err <= PROFILE_TOL
        and max_moment_err <= MOMENT_TOL
        and broken == 0.0
    )
    print(f"VALIDATION: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
