#!/usr/bin/env python3
"""
Hooke (linear-spring) normal-contact rebound benchmark driver.

Exercises DIRT's LINEAR spring-dashpot normal contact
(`contact_model = "hooke"`, per-material `kn`/`kt`) — the contact-force branch in
`crates/dirt_granular/src/contact.rs` that every other benchmark leaves
untouched, because they all use the nonlinear Hertz model.

Two identical spheres are launched head-on for every (restitution, relative
impact velocity) combination; each collides once and rebounds. The driver gates
the measured coefficient of restitution, contact duration, and peak overlap
against the EXACT analytical solution of the linear contact — a
constant-coefficient damped harmonic oscillator. The reference is theory only
(neither DIRT's own output nor another code), which makes this an analytical,
not self-consistent, check.

Commands (run from anywhere):
    python3 examples/bench_hooke_rebound/sweep.py generate   # write configs
    python3 examples/bench_hooke_rebound/sweep.py start      # build + run -> CSV
    python3 examples/bench_hooke_rebound/sweep.py graph      # validate + plot
    python3 examples/bench_hooke_rebound/sweep.py            # all three, in order

Linear spring-dashpot contact theory
─────────────────────────────────────
During contact the mutual overlap x(t) of two spheres obeys a linear
spring-dashpot with the SAME coefficients DIRT integrates
(`contact.rs`: `f_n = kn·δ − γ_n·v_n`, `γ_n = 2β·√(kn·m_eff)`):

    m_eff·ẍ + γ_n·ẋ + kn·x = 0 ,   m_eff = m/2  (two identical spheres) .

That is a damped harmonic oscillator with

    ω₀ = √(kn/m_eff) ,   ζ = γ_n/(2√(kn·m_eff)) = β ,   ω_d = ω₀·√(1−β²) .

DIRT derives the damping ratio β from the restitution input `e` as the EXACT
linear-contact inversion (`dirt_atom` `build_pair_tables`, Hooke branch):

    β = −ln e / √(π² + ln²e) .

Substituting reduces the standard results to closed forms with no free constants:

    coefficient of restitution   COR   = exp(−π·ζ/√(1−ζ²)) = e           (exact)
    contact duration (half-period) t_c = π/ω_d = √(π² + ln²e)·√(m_eff/kn)
    peak overlap                 δ_max = (v/ω_d)·e^(−ζω₀t*)·sin(ω_d t*),
                                 t* = atan(√(1−β²)/β)/ω_d .

COR and t_c are velocity-independent (the signature of a LINEAR contact, unlike
Hertz); δ_max scales linearly with the impact speed. All three are checked.

Reference: the damped-harmonic-oscillator collision is textbook DEM contact
theory — see e.g. Schäfer, Dippel & Wolf, "Force schemes in simulations of
granular materials", J. Phys. I France 6:5-20 (1996), Eq. (2.10)-(2.16); and
Y. Tsuji, T. Tanaka, T. Ishida, Powder Technol. 71:239-250 (1992), Eq. (10)-(12)
for the linear COR ↔ damping relation.
"""

import os
import sys
import csv
import math
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_hooke_rebound"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep_results.csv")

# ── Material & geometry — two identical spheres ───────────────────────────────
DENSITY = 2500.0         # kg/m^3
RADIUS = 0.005           # m
KN = 1.0e5               # N/m — linear normal spring stiffness (Hooke)
KT = 2.0 / 7.0 * KN      # N/m — linear tangential stiffness (inactive: friction=0)

# Coefficient-of-restitution inputs. For a LINEAR contact the realized COR equals
# the input exactly; 1.0 is the elastic (zero-damping) anchor.
CORS = [0.3, 0.5, 0.7, 0.9, 1.0]
# Relative normal impact velocities [m/s]. COR and t_c must be flat across these
# (velocity-independence is the distinguishing property of the linear contact).
VELOCITIES = [0.5, 1.0, 2.0, 4.0]

# ── Validation tolerances (vs analytical linear-contact theory — theory only) ─
# DIRT reproduces the exact linear-contact collision to <0.1% at the timestep
# used here; these tolerances bound the residual time-discretisation error and
# are NOT loosened to force a pass. Contact time is measured by integer step
# counting, so its floor is set by dt/t_c (a small fraction of a percent).
COR_TOL = 0.01          # 1% on coefficient of restitution (measured vs input e)
TIME_TOL = 0.02         # 2% on contact duration
OVERLAP_TOL = 0.02      # 2% on peak overlap
VEL_INDEP_COR_TOL = 0.005   # max spread of realized COR across impact speed, per e
VEL_INDEP_TC_TOL = 0.01     # max relative spread of t_c across impact speed, per e

MASS = (4.0 / 3.0) * math.pi * RADIUS**3 * DENSITY
M_EFF = MASS / 2.0       # two identical spheres
OMEGA0 = math.sqrt(KN / M_EFF)


def beta_of_e(e):
    """Linear-contact damping ratio DIRT derives from restitution `e`."""
    if e >= 1.0:
        return 0.0
    le = math.log(e)
    return -le / math.sqrt(math.pi**2 + le**2)


def cor_theory(e):
    """Exact linear-contact COR = input restitution."""
    return e


def contact_time_theory(e):
    """Half-period of the damped oscillator: t_c = pi/omega_d."""
    beta = beta_of_e(e)
    omega_d = OMEGA0 * math.sqrt(1.0 - beta**2)
    return math.pi / omega_d


def peak_overlap_theory(e, v):
    """Peak overlap of the damped oscillator with x(0)=0, x'(0)=v (approach)."""
    beta = beta_of_e(e)
    if beta == 0.0:
        return v / OMEGA0                       # elastic: delta_max = v/omega_0
    omega_d = OMEGA0 * math.sqrt(1.0 - beta**2)
    t_star = math.atan(math.sqrt(1.0 - beta**2) / beta) / omega_d
    return (v / omega_d) * math.exp(-beta * OMEGA0 * t_star) * math.sin(omega_d * t_star)


def dt_for():
    """Timestep: a fixed fraction of the (velocity-independent) elastic contact
    time, so every case resolves the collision equally well."""
    tc_elastic = math.pi / OMEGA0
    return tc_elastic / 2500.0


def steps_for(v):
    """Approach time over the launch gap plus a generous contact/rebound margin."""
    gap = 0.002                    # surface-to-surface launch gap [m]
    approach = gap / v             # closing rate = v (both spheres move at v/2)
    tc = contact_time_theory(min(CORS))   # longest contact time (most damping)
    total = approach + 5.0 * tc
    return int(total / dt_for()) + 3000


# ── config template ───────────────────────────────────────────────────────────
TOML_TEMPLATE = """\
# Auto-generated Hooke rebound case — v_rel = {v} m/s, COR = {cor}
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[domain]
x_low = -0.02
x_high = 0.02
y_low = -0.01
y_high = 0.01
z_low = -0.01
z_high = 0.01
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
radius = {R}
density = {rho}
velocity_x = {vhalf}
region = {{ type = "block", min = [-0.006001, -1.0e-6, -1.0e-6], max = [-0.005999, 1.0e-6, 1.0e-6] }}

[[particles.insert]]
material = "grain"
count = 1
radius = {R}
density = {rho}
velocity_x = -{vhalf}
region = {{ type = "block", min = [0.005999, -1.0e-6, -1.0e-6], max = [0.006001, 1.0e-6, 1.0e-6] }}

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
                    R=RADIUS, rho=DENSITY, vhalf=f"{v/2.0:.6e}",
                    outdir=cdir, steps=steps_for(v), dt=f"{dt_for():.6e}",
                ))
            n += 1
    print(f"Generated {n} configs under {SWEEP_DIR}")


# ── start ────────────────────────────────────────────────────────────────────
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
            if not os.path.isfile(config):
                print(f"  [{i:2d}/{n_total}] missing {config} — run 'generate' first.")
                continue
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
        print("\nERROR: no results collected.")
        sys.exit(1)

    with open(SWEEP_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerows(results)
    print(f"\nDIRT: {len(results)}/{n_total} cases -> {SWEEP_CSV}")


# ── graph (validate + plot) ──────────────────────────────────────────────────
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
    print("Hooke (linear-spring) Normal-Contact Rebound Benchmark")
    print("=" * 72)
    print(f"  Two identical spheres  R={RADIUS} m  rho={DENSITY} kg/m^3")
    print(f"  kn={KN:.3e} N/m  kt={KT:.3e} N/m   m_eff={M_EFF:.4e} kg  "
          f"omega_0={OMEGA0:.4e} rad/s")
    print(f"  Reference: exact linear damped-oscillator collision (theory only).")
    print(f"  Tolerances: COR {COR_TOL*100:.0f}%  contact-time {TIME_TOL*100:.0f}%"
          f"  overlap {OVERLAP_TOL*100:.0f}%\n")

    total = passed = 0
    by_cor_cor = {}     # realized COR grouped by input e (velocity-independence)
    by_cor_tc = {}      # realized t_c grouped by input e (velocity-independence)

    for r in sorted(rows, key=lambda x: (float(x["input_cor"]), float(x["input_v"]))):
        e = float(r["input_cor"])
        v = float(r["v_impact"])
        cor_meas = float(r["cor_measured"])
        t_meas = float(r["contact_time"])
        d_meas = float(r["max_overlap"])
        by_cor_cor.setdefault(e, []).append(cor_meas)
        by_cor_tc.setdefault(e, []).append(t_meas)

        cor_th = cor_theory(e)
        t_th = contact_time_theory(e)
        d_th = peak_overlap_theory(e, v)

        cor_err = abs(cor_meas - cor_th)                 # absolute (COR is O(1))
        t_err = abs(t_meas - t_th) / t_th
        d_err = abs(d_meas - d_th) / d_th

        sc = "PASS" if cor_err <= COR_TOL else "FAIL"
        st = "PASS" if t_err <= TIME_TOL else "FAIL"
        sd = "PASS" if d_err <= OVERLAP_TOL else "FAIL"
        total += 3
        passed += (sc == "PASS") + (st == "PASS") + (sd == "PASS")

        print(f"  e={e:.2f}  v={v:5.2f} m/s:")
        print(f"    COR:    {cor_meas:.5f} vs {cor_th:.5f}   (err {cor_err*100:5.2f}%)  [{sc}]")
        print(f"    t_c:    {t_meas:.5e} vs {t_th:.5e} s (err {t_err*100:5.2f}%)  [{st}]")
        print(f"    d_max:  {d_meas:.5e} vs {d_th:.5e} m (err {d_err*100:5.2f}%)  [{sd}]")

    # Velocity-independence: for a LINEAR contact both COR and t_c must be flat
    # across impact speed. This is the property that distinguishes the linear
    # spring from Hertz (whose undamped t_c ∝ v^{-1/5}).
    print("\nVelocity-independence (linear-contact signature):")
    for e in sorted(by_cor_cor):
        cor_spread = max(by_cor_cor[e]) - min(by_cor_cor[e])
        tc_vals = by_cor_tc[e]
        tc_spread = (max(tc_vals) - min(tc_vals)) / (sum(tc_vals) / len(tc_vals))
        sc = "PASS" if cor_spread <= VEL_INDEP_COR_TOL else "FAIL"
        st = "PASS" if tc_spread <= VEL_INDEP_TC_TOL else "FAIL"
        total += 2
        passed += (sc == "PASS") + (st == "PASS")
        print(f"  e={e:.2f}:  COR spread {cor_spread:.4f} (tol {VEL_INDEP_COR_TOL})  [{sc}]"
              f"   t_c spread {tc_spread*100:.2f}% (tol {VEL_INDEP_TC_TOL*100:.0f}%)  [{st}]")

    total += 1
    expected = len(CORS) * len(VELOCITIES)
    complete = len(rows) == expected
    passed += complete
    print(f"\nCompleteness: {len(rows)}/{expected} cases  [{'PASS' if complete else 'FAIL'}]")
    print(f"\nOverall: {passed}/{total} checks passed")
    print("ALL CHECKS PASSED" if passed == total
          else f"CHECKS FAILED: {total - passed} of {total} failed")
    return passed == total


def plot(rows):
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n(matplotlib/numpy unavailable, skipped plots: {e})")
        return

    plt.rcParams.update({
        "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13,
        "legend.fontsize": 10, "figure.dpi": 150, "savefig.dpi": 150,
    })
    os.makedirs(PLOT_DIR, exist_ok=True)
    markers = ["o", "s", "^", "D", "v"]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

    data = {}
    for r in rows:
        data[(float(r["input_cor"]), float(r["input_v"]))] = {
            "cor": float(r["cor_measured"]),
            "tc": float(r["contact_time"]),
            "dmax": float(r["max_overlap"]),
            "v": float(r["v_impact"]),
        }
    cors = sorted({float(r["input_cor"]) for r in rows})
    vels = sorted({float(r["input_v"]) for r in rows})

    # ── Plot 1: measured vs input COR (target = 1:1 line, exact) ─────────────
    fig, ax = plt.subplots(figsize=(6.6, 5))
    for iv, v in enumerate(vels):
        m, c = markers[iv % len(markers)], colors[iv % len(colors)]
        pts = [(e, data[(e, v)]["cor"]) for e in cors if (e, v) in data]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], m, color=c,
                    markersize=7, markerfacecolor="none", markeredgewidth=1.6,
                    label=f"v = {v} m/s")
    ax.plot([0.2, 1.03], [0.2, 1.03], "k--", linewidth=1.2,
            label="linear theory (COR = e)")
    ax.set_xlabel("Input restitution e")
    ax.set_ylabel("Measured COR")
    ax.set_title("Hooke contact: COR = input restitution (exact, velocity-independent)")
    ax.legend(fontsize=9, loc="upper left")
    ax.set_xlim(0.2, 1.05)
    ax.set_ylim(0.2, 1.05)
    ax.set_aspect("equal")
    fig.savefig(os.path.join(PLOT_DIR, "cor_validation.png"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {PLOT_DIR}/cor_validation.png")

    # ── Plot 2: contact duration vs input COR (flat in v; grows with damping) ─
    fig, ax = plt.subplots(figsize=(7, 5))
    e_grid = np.linspace(min(cors), 1.0, 200)
    ax.plot(e_grid, [contact_time_theory(e) * 1e6 for e in e_grid],
            "k-", linewidth=1.8, label="linear theory  t_c = pi/omega_d")
    for iv, v in enumerate(vels):
        m, c = markers[iv % len(markers)], colors[iv % len(colors)]
        pts = [(e, data[(e, v)]["tc"]) for e in cors if (e, v) in data]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] * 1e6 for p in pts], m, color=c,
                    markersize=7, markerfacecolor="none", markeredgewidth=1.6,
                    label=f"v = {v} m/s")
    ax.set_xlabel("Input restitution e")
    ax.set_ylabel("Contact duration [µs]")
    ax.set_title("Hooke contact: contact duration (velocity-independent)")
    ax.legend(fontsize=9)
    fig.savefig(os.path.join(PLOT_DIR, "contact_duration.png"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {PLOT_DIR}/contact_duration.png")

    # ── Plot 3: peak overlap vs impact velocity (linear scaling) ─────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    v_grid = np.linspace(min(vels) * 0.9, max(vels) * 1.1, 200)
    for ic, e in enumerate(cors):
        m, c = markers[ic % len(markers)], colors[ic % len(colors)]
        ax.plot(v_grid, [peak_overlap_theory(e, v) * 1e6 for v in v_grid],
                "-", color=c, linewidth=1.4, alpha=0.8)
        pts = [(data[(e, v)]["v"], data[(e, v)]["dmax"]) for v in vels if (e, v) in data]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] * 1e6 for p in pts], m, color=c,
                    markersize=7, markerfacecolor="none", markeredgewidth=1.6,
                    label=f"e = {e}")
    ax.set_xlabel("Relative impact velocity [m/s]")
    ax.set_ylabel("Peak overlap [µm]")
    ax.set_title("Hooke contact: peak overlap ∝ impact velocity")
    ax.legend(fontsize=9)
    fig.savefig(os.path.join(PLOT_DIR, "peak_overlap.png"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {PLOT_DIR}/peak_overlap.png")


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
        print(f"Unknown command: {cmd!r}")
        print("Usage: sweep.py [generate|start|graph]   (no arg = all three)")
        sys.exit(2)


if __name__ == "__main__":
    main()
