#!/usr/bin/env python3
"""
Hooke wall rebound benchmark driver.

Runs a single sphere against a real `dirt_wall` plane with
`contact_model = "hooke"` and validates the wall-particle normal contact against
the closed-form linear spring-dashpot solution. For a rigid wall the reduced
mass is the particle mass and DIRT uses the same `kn_ij` and `beta_ij` material
tables as the particle-particle Hooke path.
"""

import csv
import math
import os
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_hooke_wall_rebound"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep_results.csv")

DENSITY = 2500.0
RADIUS = 0.005
KN = 1.0e5
KT = 2.0 / 7.0 * KN
CORS = [0.3, 0.5, 0.7, 0.9, 1.0]
VELOCITIES = [0.5, 1.0, 2.0, 4.0]

COR_TOL = 0.01
TIME_TOL = 0.02
OVERLAP_TOL = 0.02
VEL_INDEP_COR_TOL = 0.005
VEL_INDEP_TC_TOL = 0.01

MASS = (4.0 / 3.0) * math.pi * RADIUS**3 * DENSITY
M_EFF = MASS
OMEGA0 = math.sqrt(KN / M_EFF)


def beta_of_e(e):
    if e >= 1.0:
        return 0.0
    le = math.log(e)
    return -le / math.sqrt(math.pi**2 + le**2)


def contact_time_theory(e):
    beta = beta_of_e(e)
    return math.pi / (OMEGA0 * math.sqrt(1.0 - beta**2))


def peak_overlap_theory(e, v):
    beta = beta_of_e(e)
    if beta == 0.0:
        return v / OMEGA0
    omega_d = OMEGA0 * math.sqrt(1.0 - beta**2)
    t_star = math.atan(math.sqrt(1.0 - beta**2) / beta) / omega_d
    return (v / omega_d) * math.exp(-beta * OMEGA0 * t_star) * math.sin(omega_d * t_star)


def dt_for():
    return (math.pi / OMEGA0) / 2500.0


def steps_for(v):
    gap = 0.007
    approach = gap / v
    total = approach + 5.0 * contact_time_theory(min(CORS))
    return int(total / dt_for()) + 3000


TOML_TEMPLATE = """\
# Auto-generated Hooke wall rebound case: v = {v} m/s, COR = {cor}
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[domain]
x_low = -0.01
x_high = 0.01
y_low = -0.01
y_high = 0.01
z_low = 0.0
z_high = 0.05
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[neighbor]
skin_fraction = 1.1
bin_size = 0.015
every = 1

[dem]
contact_model = "hooke"
limit_damping = false

[[dem.materials]]
name = "grain"
youngs_mod = 7.0e10
poisson_ratio = 0.30
restitution = {cor}
friction = 0.0
kn = {kn}
kt = {kt}

[[particles.insert]]
material = "grain"
count = 1
radius = {radius}
density = {density}
velocity_z = -{v}
region = {{ type = "block", min = [-1.0e-6, -1.0e-6, 0.012], max = [1.0e-6, 1.0e-6, 0.012002] }}

[[wall]]
type = "plane"
point_z = 0.0
normal_z = 1.0
material = "grain"
name = "floor"

[output]
dir = "{outdir}"

[run]
steps = {steps}
dt = {dt}
thermo = 100000
"""


def case_tag(v, cor):
    return f"v{v}_cor{cor}"


def case_dir(v, cor):
    return os.path.join(SWEEP_DIR, case_tag(v, cor))


def generate():
    os.makedirs(SWEEP_DIR, exist_ok=True)
    n = 0
    for cor in CORS:
        for v in VELOCITIES:
            cdir = case_dir(v, cor)
            os.makedirs(cdir, exist_ok=True)
            with open(os.path.join(cdir, "config.toml"), "w") as f:
                f.write(TOML_TEMPLATE.format(
                    v=v, cor=cor, kn=f"{KN:.6e}", kt=f"{KT:.6e}",
                    radius=RADIUS, density=DENSITY, outdir=cdir,
                    steps=steps_for(v), dt=f"{dt_for():.10e}",
                ))
            n += 1
    print(f"Generated {n} configs under {SWEEP_DIR}")


CSV_FIELDS = ["input_v", "input_cor", "v_impact", "v_rebound", "cor_measured",
              "max_force", "contact_time", "max_overlap", "dt", "radius", "density"]


def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE,
         "--no-default-features", "--features", "precision-double"],
        cwd=REPO_ROOT, check=True,
    )

    results = []
    n_total = len(CORS) * len(VELOCITIES)
    i = 0
    for cor in CORS:
        for v in VELOCITIES:
            i += 1
            cdir = case_dir(v, cor)
            config = os.path.join(cdir, "config.toml")
            print(f"  [{i:2d}/{n_total}] v={v:<4} COR={cor:<4}", end="  ", flush=True)
            log = os.path.join(cdir, "run.log")
            with open(log, "w") as lf:
                proc = subprocess.run(
                    ["cargo", "run", "--release", "--example", EXAMPLE,
                     "--no-default-features", "--features", "precision-double", "--", config],
                    cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
                )
            out_csv = os.path.join(cdir, "data", "rebound_results.csv")
            if proc.returncode == 0 and os.path.isfile(out_csv):
                with open(out_csv) as f:
                    row = next(csv.DictReader(f))
                row["input_v"], row["input_cor"] = str(v), str(cor)
                results.append(row)
                print(f"COR={float(row['cor_measured']):.4f}")
            else:
                print(f"FAILED ({log})")

    if not results:
        print("ERROR: no results collected.")
        sys.exit(1)

    with open(SWEEP_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerows(results)
    print(f"\nDIRT: {len(results)}/{n_total} cases -> {SWEEP_CSV}")


def load_rows():
    if not os.path.isfile(SWEEP_CSV):
        print(f"ERROR: {SWEEP_CSV} not found. Run 'start' first.")
        sys.exit(1)
    with open(SWEEP_CSV) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        print("ERROR: no data in results file.")
        sys.exit(1)
    return rows


def validate(rows):
    print("=" * 72)
    print("Hooke Wall Rebound Benchmark")
    print("=" * 72)
    print(f"  Sphere-wall contact: R={RADIUS} m rho={DENSITY} kg/m^3")
    print(f"  kn={KN:.3e} N/m m_eff={M_EFF:.4e} kg omega_0={OMEGA0:.4e} rad/s")
    print("  Reference: exact linear spring-dashpot wall collision.")
    print(f"  Tolerances: COR {COR_TOL*100:.0f}% contact-time {TIME_TOL*100:.0f}% overlap {OVERLAP_TOL*100:.0f}%\n")

    total = passed = 0
    by_cor_cor = {}
    by_cor_tc = {}

    for r in sorted(rows, key=lambda x: (float(x["input_cor"]), float(x["input_v"]))):
        e = float(r["input_cor"])
        v = float(r["v_impact"])
        cor_meas = float(r["cor_measured"])
        t_meas = float(r["contact_time"])
        d_meas = float(r["max_overlap"])
        by_cor_cor.setdefault(e, []).append(cor_meas)
        by_cor_tc.setdefault(e, []).append(t_meas)

        t_th = contact_time_theory(e)
        d_th = peak_overlap_theory(e, v)
        cor_err = abs(cor_meas - e)
        t_err = abs(t_meas - t_th) / t_th
        d_err = abs(d_meas - d_th) / d_th
        sc = "PASS" if cor_err <= COR_TOL else "FAIL"
        st = "PASS" if t_err <= TIME_TOL else "FAIL"
        sd = "PASS" if d_err <= OVERLAP_TOL else "FAIL"
        total += 3
        passed += (sc == "PASS") + (st == "PASS") + (sd == "PASS")

        print(f"  e={e:.2f} v={v:5.2f} m/s:")
        print(f"    COR:   {cor_meas:.5f} vs {e:.5f} (err {cor_err*100:5.2f}%) [{sc}]")
        print(f"    t_c:   {t_meas:.5e} vs {t_th:.5e} s (err {t_err*100:5.2f}%) [{st}]")
        print(f"    d_max: {d_meas:.5e} vs {d_th:.5e} m (err {d_err*100:5.2f}%) [{sd}]")

    print("\nVelocity-independence:")
    for e in sorted(by_cor_cor):
        cor_spread = max(by_cor_cor[e]) - min(by_cor_cor[e])
        tc_vals = by_cor_tc[e]
        tc_spread = (max(tc_vals) - min(tc_vals)) / (sum(tc_vals) / len(tc_vals))
        sc = "PASS" if cor_spread <= VEL_INDEP_COR_TOL else "FAIL"
        st = "PASS" if tc_spread <= VEL_INDEP_TC_TOL else "FAIL"
        total += 2
        passed += (sc == "PASS") + (st == "PASS")
        print(f"  e={e:.2f}: COR spread {cor_spread:.4f} [{sc}]  t_c spread {tc_spread*100:.2f}% [{st}]")

    expected = len(CORS) * len(VELOCITIES)
    complete = len(rows) == expected
    total += 1
    passed += complete
    print(f"\nCompleteness: {len(rows)}/{expected} cases [{'PASS' if complete else 'FAIL'}]")
    print(f"Overall: {passed}/{total} checks passed")
    print("ALL CHECKS PASSED" if passed == total else f"CHECKS FAILED: {total - passed} failed")
    return passed == total


def plot(rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as e:
        print(f"Plotting skipped: {e}")
        return

    os.makedirs(PLOT_DIR, exist_ok=True)
    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 150})
    data = {(float(r["input_cor"]), float(r["input_v"])): r for r in rows}
    cors = sorted({float(r["input_cor"]) for r in rows})
    vels = sorted({float(r["input_v"]) for r in rows})
    colors = ["#2166ac", "#b2182b", "#4d9221", "#762a83", "#bf812d"]
    markers = ["o", "s", "^", "D"]

    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    ax.fill_between([0.2, 1.05], [0.2 - COR_TOL, 1.05 - COR_TOL],
                    [0.2 + COR_TOL, 1.05 + COR_TOL], color="0.88",
                    label="+/- COR tolerance")
    for iv, v in enumerate(vels):
        pts = [(e, float(data[(e, v)]["cor_measured"])) for e in cors if (e, v) in data]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], markers[iv % len(markers)],
                color=colors[iv], markerfacecolor="none", label=f"v={v} m/s")
    ax.plot([0.2, 1.05], [0.2, 1.05], "k--", label="theory COR=e")
    ax.set_xlabel("Input restitution e")
    ax.set_ylabel("Measured COR")
    ax.set_title("Hooke wall rebound: restitution")
    ax.set_aspect("equal")
    ax.legend(fontsize=8)
    fig.savefig(os.path.join(PLOT_DIR, "cor_validation.png"), bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    e_grid = np.linspace(min(cors), 1.0, 200)
    theory = np.array([contact_time_theory(e) for e in e_grid]) * 1e6
    ax.fill_between(e_grid, theory * (1.0 - TIME_TOL), theory * (1.0 + TIME_TOL),
                    color="0.88", label="+/- 2% tolerance")
    ax.plot(e_grid, theory, "k-", label="linear theory")
    for iv, v in enumerate(vels):
        pts = [(e, float(data[(e, v)]["contact_time"]) * 1e6) for e in cors if (e, v) in data]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], markers[iv % len(markers)],
                color=colors[iv], markerfacecolor="none", label=f"v={v} m/s")
    ax.set_xlabel("Input restitution e")
    ax.set_ylabel("Contact duration [microseconds]")
    ax.set_title("Hooke wall rebound: contact time")
    ax.legend(fontsize=8)
    fig.savefig(os.path.join(PLOT_DIR, "contact_duration.png"), bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    v_grid = np.linspace(min(vels) * 0.9, max(vels) * 1.1, 200)
    for ic, e in enumerate(cors):
        theory = np.array([peak_overlap_theory(e, v) for v in v_grid]) * 1e6
        ax.fill_between(v_grid, theory * (1.0 - OVERLAP_TOL), theory * (1.0 + OVERLAP_TOL),
                        color=colors[ic % len(colors)], alpha=0.08)
        ax.plot(v_grid, theory, "-", color=colors[ic % len(colors)], linewidth=1.2)
        pts = [(v, float(data[(e, v)]["max_overlap"]) * 1e6) for v in vels if (e, v) in data]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o",
                color=colors[ic % len(colors)], markerfacecolor="none", label=f"e={e}")
    ax.set_xlabel("Impact velocity [m/s]")
    ax.set_ylabel("Peak overlap [micrometers]")
    ax.set_title("Hooke wall rebound: peak overlap")
    ax.legend(fontsize=8)
    fig.savefig(os.path.join(PLOT_DIR, "peak_overlap.png"), bbox_inches="tight")
    plt.close(fig)

    print(f"Saved plots under {PLOT_DIR}")


def graph():
    rows = load_rows()
    ok = validate(rows)
    plot(rows)
    return ok


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd == "generate":
        generate()
    elif cmd == "start":
        start()
    elif cmd == "graph":
        sys.exit(0 if graph() else 1)
    elif cmd == "all":
        generate()
        start()
        print()
        sys.exit(0 if graph() else 1)
    else:
        print("Usage: sweep.py [generate|start|graph]")
        sys.exit(2)


if __name__ == "__main__":
    main()
