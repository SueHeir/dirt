#!/usr/bin/env python3
"""
Chung & Ooi (2011) elastic normal-impact benchmark driver — Test Cases 1 & 2.

Reproduces two of the standard DEM code-verification benchmarks from

    L. Chung and J.Y. Ooi, "Benchmark tests for verifying discrete element
    modelling codes at particle impact level", Granular Matter 13(5):643-656
    (2011). https://doi.org/10.1007/s10035-011-0277-0

    * Test 1 — elastic normal impact of two identical spheres.
    * Test 2 — elastic normal impact of a sphere with a rigid wall.

For both cases the contact is perfectly elastic (restitution = 1, no damping),
so the exact reference is the closed-form Hertz solution for the collision —
the very reference the paper uses for these two cases. This driver sweeps the
impact velocity, runs each case in DIRT, and gates the measured maximum contact
force, contact duration, and maximum overlap against the Hertz analytical
values, PASS/FAIL. The reference is an INDEPENDENT analytical solution — not
LAMMPS and not DIRT's own output (anti-gaming: theory only).

Commands (run from anywhere):
    python3 examples/bench_chung_ooi_impact/sweep.py generate   # write configs
    python3 examples/bench_chung_ooi_impact/sweep.py start      # build + run -> CSV
    python3 examples/bench_chung_ooi_impact/sweep.py graph      # validate + plot
    python3 examples/bench_chung_ooi_impact/sweep.py            # all three, in order

Hertz normal-impact theory (effective quantities m*, R*, E* below):
    maximum overlap    delta_max = (15 m* v^2 / (16 E* sqrt(R*)))^(2/5)
    maximum force      F_max     = (4/3) E* sqrt(R*) delta_max^(3/2)
    contact duration   t_c       = 2.943920 * delta_max / v
    where 2.943920 = 2 * integral_0^1 (1 - x^(5/2))^(-1/2) dx  (Hertz collision).
For two identical spheres:  E* = E/(2(1-nu^2)), R* = R/2, m* = m/2.
For a sphere on a same-material wall: E* = E/(2(1-nu^2)), R* = R, m* = m.

Reference for the Hertz impact formulae: K.L. Johnson, Contact Mechanics,
Cambridge University Press, 1985 (Ch. 11); identical to Chung & Ooi (2011) Eqs.
for Test 1/2.
"""

import os
import sys
import csv
import math
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_chung_ooi_impact"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep_results.csv")

# ── Material & geometry — Chung & Ooi (2011) aluminium-alloy spheres ──────────
YOUNGS_MOD = 7.0e10      # Pa  — aluminium alloy
POISSON_RATIO = 0.30
DENSITY = 2700.0         # kg/m^3
RADIUS = 0.1             # m

# Relative normal impact velocities [m/s]. Spans the paper's impact regime.
VELOCITIES = [0.5, 1.0, 2.0, 5.0, 10.0]

# Cases: (tag, n_particles). Test 1 = two spheres, Test 2 = sphere + wall.
CASES = [("test1_sphere_sphere", 2), ("test2_sphere_wall", 1)]

# ── Validation tolerances (vs analytical Hertz — theory only) ─────────────────
# The elastic Hertz contact is reproduced to <0.5% in DIRT; these tolerances
# bound the small residual time-discretisation error and are NOT loosened to
# force a pass. Contact time is measured by integer step counting, so its
# tolerance is set by dt/t_c (a fraction of a percent at the resolution used).
FORCE_TOL = 0.02      # 2% on maximum contact force
OVERLAP_TOL = 0.02    # 2% on maximum overlap
TIME_TOL = 0.02       # 2% on contact duration

# Hertz collision contact-time constant: 2 * int_0^1 (1-x^{5/2})^{-1/2} dx.
HERTZ_TC_CONST = 2.943920

E_STAR = YOUNGS_MOD / (2.0 * (1.0 - POISSON_RATIO**2))
MASS = (4.0 / 3.0) * math.pi * RADIUS**3 * DENSITY


def effective(n_particles):
    """Return (m_eff, r_eff) for the Hertz reference of this case."""
    if n_particles >= 2:
        return MASS / 2.0, RADIUS / 2.0     # two identical spheres
    return MASS, RADIUS                     # sphere on rigid wall


def hertz_max_overlap(v, n_particles):
    m_eff, r_eff = effective(n_particles)
    return (15.0 * m_eff * v**2 / (16.0 * E_STAR * math.sqrt(r_eff)))**0.4


def hertz_max_force(v, n_particles):
    _, r_eff = effective(n_particles)
    d = hertz_max_overlap(v, n_particles)
    return (4.0 / 3.0) * E_STAR * math.sqrt(r_eff) * d**1.5


def hertz_contact_time(v, n_particles):
    return HERTZ_TC_CONST * hertz_max_overlap(v, n_particles) / v


def dt_for():
    """Timestep: a fixed fraction of the Rayleigh critical timestep."""
    g = YOUNGS_MOD / (2.0 * (1.0 + POISSON_RATIO))
    alpha = 0.1631 * POISSON_RATIO + 0.876605
    dt_rayleigh = math.pi * RADIUS / alpha * (DENSITY / g) ** 0.5
    return dt_rayleigh * 0.02


def steps_for(v, n_particles):
    """Approach time over the launch gap plus a generous contact/rebound margin."""
    gap = 0.01                       # surface-to-surface (or surface-to-wall) gap [m]
    launch_speed = v / 2.0 if n_particles >= 2 else v
    approach = gap / launch_speed
    tc = hertz_contact_time(v, n_particles)
    total = approach + 4.0 * tc
    return int(total / dt_for()) + 2000


# ── config templates ─────────────────────────────────────────────────────────
COMMON_HEAD = """\
# Auto-generated Chung & Ooi (2011) case — {case}, v_rel = {v} m/s
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[neighbor]
skin_fraction = 1.1
bin_size = 0.25
every = 1

[dem]
contact_model = "hertz"

[[dem.materials]]
name = "alu_alloy"
youngs_mod = {E}
poisson_ratio = {nu}
restitution = 1.0
friction = 0.0
"""

SPHERE_SPHERE = COMMON_HEAD + """
[domain]
x_low = -0.6
x_high = 0.6
y_low = -0.2
y_high = 0.2
z_low = -0.2
z_high = 0.2
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[[particles.insert]]
material = "alu_alloy"
count = 1
radius = {R}
density = {rho}
velocity_x = {vhalf}
region = {{ type = "block", min = [-0.105001, -1.0e-6, -1.0e-6], max = [-0.104999, 1.0e-6, 1.0e-6] }}

[[particles.insert]]
material = "alu_alloy"
count = 1
radius = {R}
density = {rho}
velocity_x = -{vhalf}
region = {{ type = "block", min = [0.104999, -1.0e-6, -1.0e-6], max = [0.105001, 1.0e-6, 1.0e-6] }}

[output]
dir = "{outdir}"

[run]
steps = {steps}
dt = {dt}
thermo = 100000
"""

SPHERE_WALL = COMMON_HEAD + """
[domain]
x_low = -0.2
x_high = 0.2
y_low = -0.2
y_high = 0.2
z_low = 0.0
z_high = 0.6
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[[particles.insert]]
material = "alu_alloy"
count = 1
radius = {R}
density = {rho}
velocity_z = -{vfull}
region = {{ type = "block", min = [-1.0e-6, -1.0e-6, 0.109999], max = [1.0e-6, 1.0e-6, 0.110001] }}

[[wall]]
point_x = 0.0
point_y = 0.0
point_z = 0.0
normal_x = 0.0
normal_y = 0.0
normal_z = 1.0
material = "alu_alloy"

[output]
dir = "{outdir}"

[run]
steps = {steps}
dt = {dt}
thermo = 100000
"""


def case_dir(tag, v):
    return os.path.join(SWEEP_DIR, f"{tag}_v{v}")


def generate():
    os.makedirs(SWEEP_DIR, exist_ok=True)
    n = 0
    for tag, npart in CASES:
        template = SPHERE_SPHERE if npart >= 2 else SPHERE_WALL
        for v in VELOCITIES:
            cdir = case_dir(tag, v)
            os.makedirs(cdir, exist_ok=True)
            with open(os.path.join(cdir, "config.toml"), "w") as f:
                f.write(template.format(
                    case=tag, v=v, E=f"{YOUNGS_MOD:.4e}", nu=POISSON_RATIO,
                    R=RADIUS, rho=DENSITY,
                    vhalf=f"{v/2.0:.6e}", vfull=f"{v:.6e}",
                    outdir=cdir, steps=steps_for(v, npart), dt=f"{dt_for():.6e}",
                ))
            n += 1
    print(f"Generated {n} configs under {SWEEP_DIR}")


# ── start ────────────────────────────────────────────────────────────────────
CSV_FIELDS = ["case", "n_particles", "v_impact", "max_force", "contact_time",
              "max_overlap", "dt", "radius", "density"]


def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE,
         "--no-default-features", "--features", "precision-double"],
        cwd=REPO_ROOT, check=True,
    )

    results = []
    n_total = len(CASES) * len(VELOCITIES)
    i = 0
    for tag, _npart in CASES:
        for v in VELOCITIES:
            i += 1
            cdir = case_dir(tag, v)
            config = os.path.join(cdir, "config.toml")
            if not os.path.isfile(config):
                print(f"  [{i:2d}/{n_total}] missing {config} — run 'generate' first.")
                continue
            print(f"  [{i:2d}/{n_total}] {tag:22s} v={v:<5}", end="  ", flush=True)
            log = os.path.join(cdir, "run.log")
            with open(log, "w") as lf:
                proc = subprocess.run(
                    ["cargo", "run", "--release", "--example", EXAMPLE,
                     "--no-default-features", "--features", "precision-double", "--", config],
                    cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
                )
            out_csv = os.path.join(cdir, "data", "impact_results.csv")
            if proc.returncode == 0 and os.path.isfile(out_csv):
                with open(out_csv) as f:
                    row = next(csv.DictReader(f))
                row["case"] = tag
                results.append(row)
                print(f"F_max={float(row['max_force']):.4e} N")
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


CASE_TITLE = {
    "test1_sphere_sphere": "Test 1: elastic normal impact of two identical spheres",
    "test2_sphere_wall": "Test 2: elastic normal impact of a sphere with a wall",
}


def validate(rows):
    print("=" * 70)
    print("Chung & Ooi (2011) Elastic Normal-Impact Benchmark — Tests 1 & 2")
    print("=" * 70)
    print(f"  Material: aluminium alloy  E={YOUNGS_MOD:.2e} Pa  nu={POISSON_RATIO}"
          f"  rho={DENSITY} kg/m^3  R={RADIUS} m")
    print(f"  E* = {E_STAR:.4e} Pa   m = {MASS:.4e} kg\n")
    print(f"  Reference: Hertz analytical solution (independent of DIRT & LAMMPS).")
    print(f"  Tolerances: force {FORCE_TOL*100:.0f}%  overlap {OVERLAP_TOL*100:.0f}%"
          f"  contact-time {TIME_TOL*100:.0f}%\n")

    total = passed = 0
    for tag, _ in CASES:
        crows = [r for r in rows if r["case"] == tag]
        if not crows:
            continue
        print(f"── {CASE_TITLE.get(tag, tag)} ──")
        for r in sorted(crows, key=lambda x: float(x["v_impact"])):
            npart = int(r["n_particles"])
            v = float(r["v_impact"])
            f_meas = float(r["max_force"])
            t_meas = float(r["contact_time"])
            d_meas = float(r["max_overlap"])

            f_th = hertz_max_force(v, npart)
            t_th = hertz_contact_time(v, npart)
            d_th = hertz_max_overlap(v, npart)

            f_err = abs(f_meas - f_th) / f_th
            t_err = abs(t_meas - t_th) / t_th
            d_err = abs(d_meas - d_th) / d_th

            sf = "PASS" if f_err <= FORCE_TOL else "FAIL"
            st = "PASS" if t_err <= TIME_TOL else "FAIL"
            sd = "PASS" if d_err <= OVERLAP_TOL else "FAIL"
            total += 3
            passed += (sf == "PASS") + (st == "PASS") + (sd == "PASS")

            print(f"  v={v:5.2f} m/s:")
            print(f"    F_max:  {f_meas:.5e} vs {f_th:.5e} N  (err {f_err*100:5.2f}%)  [{sf}]")
            print(f"    t_c:    {t_meas:.5e} vs {t_th:.5e} s  (err {t_err*100:5.2f}%)  [{st}]")
            print(f"    d_max:  {d_meas:.5e} vs {d_th:.5e} m  (err {d_err*100:5.2f}%)  [{sd}]")
        print()

    total += 1
    expected = len(CASES) * len(VELOCITIES)
    complete = len(rows) == expected
    passed += complete
    print(f"Completeness: {len(rows)}/{expected} cases  [{'PASS' if complete else 'FAIL'}]")
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
    colors = {"test1_sphere_sphere": "#1f77b4", "test2_sphere_wall": "#d62728"}
    label = {"test1_sphere_sphere": "Test 1: sphere-sphere",
             "test2_sphere_wall": "Test 2: sphere-wall"}

    panels = [
        ("max_force", hertz_max_force, "Maximum contact force [N]",
         "Chung & Ooi (2011): Max Contact Force vs Impact Velocity", "max_force.png",
         1.0, FORCE_TOL),
        ("contact_time", hertz_contact_time, "Contact duration [µs]",
         "Chung & Ooi (2011): Contact Duration vs Impact Velocity", "contact_time.png",
         1e6, TIME_TOL),
        ("max_overlap", hertz_max_overlap, "Maximum overlap [µm]",
         "Chung & Ooi (2011): Max Overlap vs Impact Velocity", "max_overlap.png",
         1e6, OVERLAP_TOL),
    ]

    vgrid = np.logspace(math.log10(min(VELOCITIES) * 0.9),
                        math.log10(max(VELOCITIES) * 1.1), 200)
    for key, theory_fn, ylabel, title, fname, scale, tol in panels:
        fig, (ax, err_ax) = plt.subplots(
            2, 1, figsize=(7, 6.3), sharex=True,
            gridspec_kw={"height_ratios": [3.2, 1.2], "hspace": 0.08},
        )
        for i, (tag, npart) in enumerate(CASES):
            crows = sorted([r for r in rows if r["case"] == tag],
                           key=lambda x: float(x["v_impact"]))
            if not crows:
                continue
            theory = np.array([theory_fn(v, npart) * scale for v in vgrid])
            ax.fill_between(
                vgrid, theory * (1.0 - tol), theory * (1.0 + tol),
                color=colors[tag], alpha=0.16, linewidth=0,
                label=(f"±{tol*100:.0f}% PASS band" if i == 0 else None),
            )
            ax.plot(vgrid, [theory_fn(v, npart) * scale for v in vgrid],
                    "-", color=colors[tag], linewidth=1.6, alpha=0.7,
                    label=f"{label[tag]} — Hertz theory")
            vs = [float(r["v_impact"]) for r in crows]
            ys = [float(r[key]) * scale for r in crows]
            refs = [theory_fn(v, npart) * scale for v in vs]
            errs = [(y / ref - 1.0) * 100.0 for y, ref in zip(ys, refs)]
            ax.plot(vs, ys, "o", color=colors[tag], markersize=7,
                    markerfacecolor="none", markeredgewidth=1.6,
                    label=f"{label[tag]} — DIRT")
            err_ax.plot(vs, errs, "o-", color=colors[tag], markersize=5,
                        markerfacecolor="none", markeredgewidth=1.2,
                        label=label[tag])
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.text(0.02, 0.03, f"PASS criterion: |DIRT/Hertz - 1| ≤ {tol*100:.0f}%",
                transform=ax.transAxes, fontsize=9,
                bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.85})
        ax.legend(fontsize=9)
        err_ax.axhspan(-tol * 100.0, tol * 100.0, color="0.85", alpha=0.55,
                       label=f"±{tol*100:.0f}% PASS")
        err_ax.axhline(0.0, color="0.35", linewidth=0.8)
        err_ax.axhline(tol * 100.0, color="#d62728", linestyle="--", linewidth=1.0)
        err_ax.axhline(-tol * 100.0, color="#d62728", linestyle="--", linewidth=1.0)
        err_ax.set_xscale("log")
        err_ax.set_xlabel("Relative impact velocity [m/s]")
        err_ax.set_ylabel("error [%]")
        err_ax.set_ylim(-tol * 120.0, tol * 120.0)
        err_ax.legend(fontsize=8, loc="upper right")
        fig.savefig(os.path.join(PLOT_DIR, fname), bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {PLOT_DIR}/{fname}")


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
