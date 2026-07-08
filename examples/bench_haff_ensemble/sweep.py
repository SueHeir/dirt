#!/usr/bin/env python3
"""Seeded Haff cooling-time ensemble validation for spheres, clumps, and rods."""

import csv
import importlib.util
import math
import os
import shutil
import subprocess
import sys
import tomllib
from dataclasses import dataclass

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
CONFIG = os.path.join(SCRIPT_DIR, "config.toml")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
GEN_DIR = os.path.join(DATA_DIR, "generated")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
PLOT = os.path.join(PLOT_DIR, "haff_ensemble.png")
SUMMARY_CSV = os.path.join(DATA_DIR, "haff_ensemble_summary.csv")


@dataclass(frozen=True)
class Case:
    key: str
    title: str
    example: str
    insert_table: str
    count: int
    diameter: float
    fit_after_equilibration: bool
    slope_gate: tuple[float, float]

    @property
    def path(self):
        return os.path.join(REPO_ROOT, "examples", self.example)

    @property
    def config(self):
        return os.path.join(self.path, "config.toml")

    @property
    def sweep(self):
        return os.path.join(self.path, "sweep.py")


CASES = [
    Case("sphere", "Single rough spheres", "bench_sphere_haff_cooling",
         "particles.insert", 800, 2.0 * 0.0011, False, (-2.3, -1.7)),
    Case("clump", "Multisphere clumps", "bench_clump_haff_cooling",
         "clump.insert", 500, 2.0 * 0.0011, True, (-2.3, -1.6)),
    Case("rod", "Rod clumps", "bench_rod_haff_cooling",
         "clump.insert", 500, 2.0 * 0.0017, True, (-2.3, -1.6)),
]


def load_driver_config():
    with open(CONFIG, "rb") as f:
        cfg = tomllib.load(f)
    return cfg["ensemble"]["seeds"], cfg["theory_band"]


def load_toml(path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def material_params(case):
    cfg = load_toml(case.config)
    mat = cfg["dem"]["materials"][0]
    domain = cfg["domain"]
    return {
        "E": float(mat["youngs_mod"]),
        "nu": float(mat["poisson_ratio"]),
        "e": float(mat["restitution"]),
        "mu": float(mat["friction"]),
        "density": float(cfg.get("particles", {}).get("insert", [{}])[0].get(
            "density", cfg.get("clump", {}).get("insert", [{}])[0].get("density", 2500.0))),
        "L": float(domain["x_high"]) - float(domain["x_low"]),
    }


def write_seed_config(case, seed, out_dir):
    text = open(case.config, encoding="utf-8").read().splitlines()
    out = []
    in_output = False
    in_insert = False
    inserted_seed = False
    table_header = f"[[{case.insert_table}]]"
    for line in text:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            if in_insert and not inserted_seed:
                out.append(f"seed = {seed}")
            in_output = stripped == "[output]"
            in_insert = stripped == table_header
            inserted_seed = False
            out.append(line)
            continue
        if in_output and stripped.startswith("dir ="):
            out.append(f'dir = "{out_dir}"')
            continue
        if in_insert and stripped.startswith("seed ="):
            out.append(f"seed = {seed}")
            inserted_seed = True
            continue
        out.append(line)
    if in_insert and not inserted_seed:
        out.append(f"seed = {seed}")
    os.makedirs(GEN_DIR, exist_ok=True)
    path = os.path.join(GEN_DIR, f"{case.key}_seed_{seed}.toml")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    return path


def import_sweep(case):
    spec = importlib.util.spec_from_file_location(f"{case.key}_haff_sweep", case.sweep)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_cmd(args, cwd=REPO_ROOT, stdout=None):
    return subprocess.run(args, cwd=cwd, stdout=stdout, stderr=subprocess.STDOUT, check=True)


def build_example(case):
    print(f"Building {case.example} ...", flush=True)
    run_cmd(["cargo", "build", "--release", "--example", case.example,
             "--no-default-features", "--features", "precision-double"])


def run_dirt(case, seed):
    out_dir = os.path.join(DATA_DIR, case.key, f"seed_{seed}")
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    cfg = write_seed_config(case, seed, out_dir)
    log = os.path.join(out_dir, "dirt_run.log")
    with open(log, "w", encoding="utf-8") as lf:
        proc = subprocess.run(
            ["cargo", "run", "--release", "--example", case.example,
             "--no-default-features", "--features", "precision-double", "--", cfg],
            cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT)
    csv_path = os.path.join(out_dir, "cooling.csv")
    if proc.returncode != 0 or not os.path.isfile(csv_path):
        raise RuntimeError(f"{case.key} seed {seed} DIRT run failed; see {log}")
    return csv_path


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def equilibration_time(rows):
    tr = [(r["time"], r.get("T_trans", 0.0), r.get("T_rot", 0.0)) for r in rows
          if r["time"] > 0 and r.get("T_trans", 0.0) > 0]
    if len(tr) < 8:
        return 0.0
    ratios = [ro / tt for _, tt, ro in tr]
    plateau = sorted(ratios[len(ratios) // 2:])[len(ratios[len(ratios) // 2:]) // 2]
    idx = next((i for i, q in enumerate(ratios) if q >= 0.9 * plateau), 0)
    return min(tr[idx][0], 0.3 * tr[-1][0])


def linfit(x, y):
    import numpy as np
    x = np.asarray(x)
    y = np.asarray(y)
    b, a = np.polyfit(x, y, 1)
    pred = b * x + a
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return a, b, pred, 1.0 - ss_res / ss_tot


def haff_fit(rows, t_min):
    import numpy as np
    t = np.array([r["time"] for r in rows])
    T = np.array([r["T_total"] for r in rows])
    T0r = T[T > 0][0]
    win = (t >= t_min) & (T > 1e-3 * T0r) & np.isfinite(T)
    tw, Tw = t[win], T[win]
    y = 1.0 / np.sqrt(Tw)
    a, b, pred, r2 = linfit(tw, y)
    T0 = 1.0 / a**2
    tc = a / b
    pos = tw > 0
    tp, Tp = tw[pos], Tw[pos]
    half = slice(len(tp) // 2, None)
    slope = float(np.polyfit(np.log(tp[half]), np.log(Tp[half]), 1)[0])
    return {
        "T0": T0, "tc": tc, "r2": r2, "slope": slope, "t_over_tc": tp[-1] / tc,
        "t": tw, "T": Tw, "linear_y": y, "linear_pred": pred,
    }


def kinetic_tc(case, T_start):
    p = material_params(case)
    n = case.count / p["L"]**3
    phi = n * math.pi / 6.0 * case.diameter**3
    g0 = (1.0 - phi / 2.0) / (1.0 - phi) ** 3
    omega0 = (4.0 / 3.0) * n * case.diameter**2 * g0 * math.sqrt(math.pi) * (1.0 - p["e"]**2)
    return 2.0 / (omega0 * math.sqrt(T_start)), phi, g0


def validate_case(case, seed_results, band):
    import numpy as np
    print(f"\n{case.title}")
    total = passed = 0

    def check(name, ok, detail=""):
        nonlocal total, passed
        total += 1
        passed += bool(ok)
        print(f"  {name:<34}{'PASS' if ok else 'FAIL'}   {detail}")

    for result in seed_results:
        seed = result["seed"]
        rows = result["rows"]
        T = np.array([r["T_total"] for r in rows])
        fit = result["fit"]
        check(f"seed {seed}: finite temperatures", bool(np.all(np.isfinite(T))))
        check(f"seed {seed}: non-negative T", bool(np.all(T >= 0)))
        check(f"seed {seed}: cooling", T[-1] < T[0], f"Ti={T[0]:.3e} Tf={T[-1]:.3e}")
        check(f"seed {seed}: no energy growth", float(np.max(T)) < 1.5 * T[0])
        check(f"seed {seed}: Haff linearity", fit["r2"] > 0.99, f"R^2={fit['r2']:.4f}")
        lo, hi = case.slope_gate
        check(f"seed {seed}: late slope", lo < fit["slope"] < hi,
              f"slope={fit['slope']:.3f} at t/tc={fit['t_over_tc']:.1f}")

    ratios = np.array([r["fit"]["tc"] / r["tc_theory"] for r in seed_results])
    median = float(np.median(ratios))
    lo, hi = band
    check("ensemble median tc/theory", lo <= median <= hi,
          f"median={median:.2f}, range=[{ratios.min():.2f}, {ratios.max():.2f}], band=[{lo}, {hi}]")
    print(f"  Result: {passed}/{total} checks passed")
    return passed == total


def run_lammps_reference(case, seed0_csv):
    mod = import_sweep(case)
    lammps = mod.find_lammps()
    if not lammps:
        return None
    p = material_params(case)
    for name, value in [("YOUNGS_MOD", p["E"]), ("POISSON", p["nu"]),
                        ("RESTITUTION", p["e"]), ("FRICTION", p["mu"])]:
        if hasattr(mod, name):
            setattr(mod, name, value)
    lmp_dir = os.path.join(DATA_DIR, case.key, "lammps")
    os.makedirs(lmp_dir, exist_ok=True)
    mod.DATA_DIR = lmp_dir
    mod.LMP_INPUT = os.path.join(lmp_dir, "in.lammps")
    mod.LMP_TRACE = os.path.join(lmp_dir, "haff_trace.txt")
    mod.LMP_CSV = os.path.join(lmp_dir, "lammps_cooling.csv")
    if hasattr(mod, "MOL_FILE"):
        mod.MOL_FILE = os.path.join(lmp_dir, os.path.basename(mod.MOL_FILE))
        mod.write_molecule()
    rows = load_csv(seed0_csv)
    dt = next((r["time"] / r["step"] for r in rows if r["step"] > 0), mod.dt_rayleigh_fraction())
    try:
        if case.key == "sphere":
            mod.write_lammps_input(dt)
        else:
            t_create = mod.calibrate_t_create(lammps, dt, rows[0]["T_total"])
            mod.write_lammps_input(dt, t_create)
        if os.path.exists(mod.LMP_TRACE):
            os.remove(mod.LMP_TRACE)
        log = os.path.join(lmp_dir, "lammps.log")
        proc = subprocess.run([lammps, "-in", mod.LMP_INPUT, "-log", log],
                              cwd=REPO_ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        if proc.returncode != 0 or not os.path.isfile(mod.LMP_TRACE):
            return None
        lrows = mod.parse_lammps_trace(dt)
        with open(mod.LMP_CSV, "w", newline="", encoding="utf-8") as f:
            keys = list(lrows[0].keys())
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(lrows)
        return load_csv(mod.LMP_CSV)
    except Exception as exc:
        print(f"  LAMMPS {case.key} skipped: {exc}")
        return None


def graph(all_results, bands):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(PLOT_DIR, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "figure.dpi": 150, "savefig.dpi": 150})
    fig, axes = plt.subplots(len(CASES), 3, figsize=(15, 11))
    for row, case in enumerate(CASES):
        results = all_results[case.key]
        ax = axes[row][0]
        for result in results["seeds"]:
            fit = result["fit"]
            t0 = result["t_min"]
            t = fit["t"] - t0
            T = fit["T"] / fit["T"][0]
            mask = t > 0
            ax.loglog(t[mask], T[mask], "o", ms=2.5, alpha=0.45, label=f"DIRT seed {result['seed']}")
        if results["lammps"]:
            lr = results["lammps"]
            t_min = results["seeds"][0]["t_min"]
            t = np.array([r["time"] for r in lr])
            T = np.array([r["T_total"] for r in lr])
            m = t >= t_min
            tt = t[m] - t[m][0]
            TT = T[m] / T[m][0]
            q = tt > 0
            ax.loglog(tt[q], TT[q], "s", ms=2.5, color="#ff7f0e", alpha=0.6, label="LAMMPS seed 0")
        med = sorted(results["seeds"], key=lambda r: r["fit"]["tc"])[len(results["seeds"]) // 2]
        tf = np.linspace(max(1e-12, (med["fit"]["t"][1] - med["t_min"])), med["fit"]["t"][-1] - med["t_min"], 200)
        ax.loglog(tf, 1.0 / (1.0 + tf / med["fit"]["tc"]) ** 2, "k-", lw=1.3, label="median Haff fit")
        ax.set_title(case.title)
        ax.set_xlabel("Fit-window time [s]")
        ax.set_ylabel("T / T_start")
        ax.legend(fontsize=7)

        ax = axes[row][1]
        for result in results["seeds"]:
            fit = result["fit"]
            resid = fit["linear_y"] - fit["linear_pred"]
            ax.plot(fit["t"] - result["t_min"], resid / np.mean(fit["linear_y"]), "o-", ms=2.5,
                    alpha=0.65, label=f"seed {result['seed']}")
        ax.axhline(0.0, color="black", lw=0.8)
        ax.set_xlabel("Fit-window time [s]")
        ax.set_ylabel("linearized residual / mean")
        ax.set_title(f"Haff residuals (R2 min {min(r['fit']['r2'] for r in results['seeds']):.4f})")

        ax = axes[row][2]
        tcs = np.array([r["fit"]["tc"] for r in results["seeds"]])
        theory = np.array([r["tc_theory"] for r in results["seeds"]])
        ratios = tcs / theory
        lo, hi = bands[case.key]
        ax.axhspan(lo, hi, color="#ffcc80", alpha=0.35, label="theory band")
        ax.axhline(1.0, color="black", lw=1.0, label="kinetic estimate")
        ax.plot([r["seed"] for r in results["seeds"]], ratios, "o", color="#1f77b4")
        ax.set_xlabel("Seed")
        ax.set_ylabel("fitted tc / kinetic tc")
        ax.set_ylim(0, max(hi * 1.15, float(ratios.max()) * 1.2))
        ax.set_title(f"tc distribution, median={np.median(ratios):.2f}")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(PLOT, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {PLOT}")


def main():
    seeds, bands = load_driver_config()
    os.makedirs(DATA_DIR, exist_ok=True)
    all_results = {}
    summary_rows = []
    for case in CASES:
        build_example(case)
        seed_results = []
        for seed in seeds:
            print(f"Running {case.key} seed {seed} ...", flush=True)
            csv_path = run_dirt(case, seed)
            rows = load_csv(csv_path)
            t_min = equilibration_time(rows) if case.fit_after_equilibration else 0.0
            fit = haff_fit(rows, t_min)
            T_start = fit["T"][0]
            tc_theory, phi, g0 = kinetic_tc(case, T_start)
            result = {"seed": seed, "csv": csv_path, "rows": rows, "t_min": t_min,
                      "fit": fit, "tc_theory": tc_theory, "phi": phi, "g0": g0}
            seed_results.append(result)
            summary_rows.append({
                "case": case.key, "seed": seed, "tc_fit": fit["tc"], "tc_theory": tc_theory,
                "tc_fit_over_theory": fit["tc"] / tc_theory, "r2": fit["r2"],
                "slope": fit["slope"], "t_over_tc": fit["t_over_tc"], "phi": phi, "g0": g0,
            })
            print(f"  tc_fit={fit['tc']:.4e}s tc_theory={tc_theory:.4e}s "
                  f"ratio={fit['tc']/tc_theory:.2f} R2={fit['r2']:.4f} slope={fit['slope']:.3f}")
        lammps = run_lammps_reference(case, seed_results[0]["csv"])
        all_results[case.key] = {"seeds": seed_results, "lammps": lammps}

    with open(SUMMARY_CSV, "w", newline="", encoding="utf-8") as f:
        keys = list(summary_rows[0].keys())
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(summary_rows)

    ok = True
    print("\n" + "=" * 72)
    print("Haff Ensemble Kinetic Cooling-Time Validation")
    print("=" * 72)
    for case in CASES:
        ok = validate_case(case, all_results[case.key]["seeds"], bands[case.key]) and ok
    graph(all_results, bands)
    print(f"\nSummary: {SUMMARY_CSV}")
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
