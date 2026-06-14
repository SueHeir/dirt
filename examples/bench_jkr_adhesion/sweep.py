#!/usr/bin/env python3
"""
JKR adhesion pull-off benchmark driver.

Brings two identical glass spheres into adhesive contact and slowly separates
them, recording the contact normal force versus separation. The peak tensile
(attractive) force is the pull-off force. The sweep varies the work of adhesion
w (the material `surface_energy`) and checks that the measured pull-off force is
linear in w with the JKR slope:

    F_pulloff = (3/2) * pi * w * R*      (JKR)

DIRT's contact model (dirt_granular::contact) implements adhesion as a constant
attractive force: with adhesion_model = "jkr" (the default) the force is
F_adh = (3/2) pi gamma R*, and with "dmt" it is F_dmt = 2 pi gamma R*, where
gamma = surface_energy = w. The pull-off force is that constant, sampled cleanly
in the gap (adhesion-only) regime. This benchmark uses JKR.

Commands (from anywhere):
    python3 examples/bench_jkr_adhesion/sweep.py generate   # write per-case configs
    python3 examples/bench_jkr_adhesion/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_jkr_adhesion/sweep.py graph      # validate + plot
    python3 examples/bench_jkr_adhesion/sweep.py            # all three, in order

LAMMPS is intentionally NOT run here: DIRT's JKR/DMT is a simplified constant-
force model with no exact LAMMPS counterpart (LAMMPS' jkr is the full
Maugis-area model with a different force-overlap law), so a code-to-code overlay
would compare different physics. Validation is against the analytical pull-off
line only. The LAMMPS hook is left as a stub for future work.

Outputs:
    sweep/<case>/config.toml   DIRT configs                  (gitignored)
    data/sweep_results.csv     measured pull-off per case     (gitignored)
    data/<case>_trace.csv      force-separation traces        (gitignored)
    plots/*.png                final figures                  (tracked)

Reference: K.L. Johnson, K. Kendall, A.D. Roberts, "Surface energy and the
contact of elastic solids", Proc. R. Soc. Lond. A 324:301-313, 1971.
"""

import os
import sys
import csv
import math
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_jkr_adhesion"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep_results.csv")

# LAMMPS binary candidates. LAMMPS is not used by this benchmark (see module
# docstring) but the list is kept for structural consistency with the others.
LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# Parameter sweep: work of adhesion w = surface_energy [J/m^2].
SURFACE_ENERGIES = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

# Adhesion model under test. DIRT default is "jkr".
#   jkr -> F_pulloff = 1.5 * pi * w * R*
#   dmt -> F_pulloff = 2.0 * pi * w * R*
ADHESION_MODEL = "jkr"
PULLOFF_COEFF = 1.5 if ADHESION_MODEL == "jkr" else 2.0

# Material / geometry.
YOUNGS_MOD = 70.0e9     # Pa — soda-lime glass
POISSON_RATIO = 0.22
RESTITUTION = 0.5
RADIUS = 0.005          # m
DENSITY = 2500.0        # kg/m^3
APPROACH_VEL = 0.002    # m/s — slow, so the ~10 nm adhesion window is sampled

# Two identical spheres -> R* = R/2.
R_EFF = RADIUS / 2.0

# Validation tolerances.
PULLOFF_REL_TOL = 0.02     # 2% per-case: measured vs analytical pull-off
LINEAR_R2_MIN = 0.999      # F_pulloff(w) must be linear in w
SLOPE_REL_TOL = 0.02       # fitted slope vs analytical 1.5*pi*R*

TOML_TEMPLATE = """\
# Auto-generated JKR pull-off config — w (surface_energy) = {w} J/m^2
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[domain]
x_low = -0.05
x_high = 0.05
y_low = -0.02
y_high = 0.02
z_low = -0.02
z_high = 0.02
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[neighbor]
skin_fraction = 1.1
bin_size = 0.02
every = 1

[dem]
contact_model = "hertz"
adhesion_model = "{adhesion_model}"

[[dem.materials]]
name = "glass"
youngs_mod = {youngs_mod}
poisson_ratio = {poisson_ratio}
restitution = {restitution}
friction = 0.0
surface_energy = {w}

# Frozen sphere (center x = 0).
[[particles.insert]]
material = "glass"
count = 1
radius = {radius}
density = {density}
region = {{ type = "block", min = [-1.0e-6, -1.0e-6, -1.0e-6], max = [1.0e-6, 1.0e-6, 1.0e-6] }}

# Free sphere, launched slowly inward (center x = 0.01105, surface gap 1.05 mm).
[[particles.insert]]
material = "glass"
count = 1
radius = {radius}
density = {density}
velocity_x = -{approach_vel}
region = {{ type = "block", min = [0.0110499, -1.0e-6, -1.0e-6], max = [0.0110501, 1.0e-6, 1.0e-6] }}

[[group]]
name = "frozen"
dynamic = false
region = {{ type = "block", min = [-0.002, -0.002, -0.002], max = [0.002, 0.002, 0.002] }}

[[freeze]]
group = "frozen"

[output]
dir = "{output_dir}"

[run]
steps = {steps}
thermo = 200000
"""


def case_tag(w):
    return f"w{w}"


def case_dir(w):
    return os.path.join(SWEEP_DIR, case_tag(w))


def steps_for(_w):
    """Step budget: travel the 1.05 mm gap at APPROACH_VEL plus a generous
    margin for the contact and separation, with the solver default timestep."""
    g = YOUNGS_MOD / (2.0 * (1.0 + POISSON_RATIO))
    alpha = 0.1631 * POISSON_RATIO + 0.876605
    dt_rayleigh = math.pi * RADIUS / alpha * (DENSITY / g) ** 0.5
    dt = dt_rayleigh * 0.15
    travel_time = 0.00105 / APPROACH_VEL
    total_time = travel_time * 1.4   # margin for contact + separation
    return int(total_time / dt) + 50000


def find_lammps():
    for b in LAMMPS_BINS:
        path = shutil.which(b)
        if path:
            return path
    return None


def analytical_pulloff(w):
    """Analytical pull-off force for the model under test."""
    return PULLOFF_COEFF * math.pi * w * R_EFF


# ── generate ───────────────────────────────────────────────────────────────
def generate():
    n = 0
    for w in SURFACE_ENERGIES:
        cdir = case_dir(w)
        os.makedirs(cdir, exist_ok=True)
        with open(os.path.join(cdir, "config.toml"), "w") as f:
            f.write(TOML_TEMPLATE.format(
                w=w,
                adhesion_model=ADHESION_MODEL,
                youngs_mod=f"{YOUNGS_MOD:.6e}",
                poisson_ratio=POISSON_RATIO,
                restitution=RESTITUTION,
                radius=RADIUS, density=DENSITY,
                approach_vel=APPROACH_VEL,
                output_dir=cdir, steps=steps_for(w),
            ))
        n += 1
    print(f"Generated {n} DIRT configs under {SWEEP_DIR}")


# ── start ──────────────────────────────────────────────────────────────────
CSV_FIELDS = ["w", "f_pulloff", "sep_at_pulloff", "r_eff", "radius", "density", "dt"]


def start():
    n_total = len(SURFACE_ENERGIES)
    os.makedirs(DATA_DIR, exist_ok=True)

    # Wipe stale results so a failed run can never re-plot old data.
    if os.path.isfile(SWEEP_CSV):
        os.remove(SWEEP_CSV)

    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE, "--no-default-features"],
        cwd=REPO_ROOT, check=True,
    )

    lammps = find_lammps()
    print("LAMMPS: not used by this benchmark (constant-force JKR/DMT has no "
          "matching LAMMPS model)." if lammps is None
          else f"LAMMPS: {lammps} found but not used by this benchmark.")

    results = []
    for i, w in enumerate(SURFACE_ENERGIES, 1):
        cdir = case_dir(w)
        config = os.path.join(cdir, "config.toml")
        if not os.path.isfile(config):
            print(f"  [{i}/{n_total}] missing {config} — run 'generate' first.")
            continue
        print(f"  [{i}/{n_total}] w={w:<5} J/m^2", end="  ", flush=True)

        log = os.path.join(cdir, "run.log")
        with open(log, "w") as lf:
            proc = subprocess.run(
                ["cargo", "run", "--release", "--example", EXAMPLE,
                 "--no-default-features", "--", config],
                cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
            )
        csv_path = os.path.join(cdir, "data", "jkr_results.csv")
        if proc.returncode == 0 and os.path.isfile(csv_path):
            with open(csv_path) as f:
                row = next(csv.DictReader(f))
            row["w"] = str(w)
            results.append(row)
            # Copy the per-case trace next to the aggregate data for plotting.
            trace_src = os.path.join(cdir, "data", "jkr_trace.csv")
            if os.path.isfile(trace_src):
                shutil.copyfile(trace_src,
                                os.path.join(DATA_DIR, f"{case_tag(w)}_trace.csv"))
            print(f"F_pulloff={float(row['f_pulloff']):.4e} N")
        else:
            print(f"FAILED ({log})")

    if not results:
        print("\nERROR: no results collected.")
        sys.exit(1)

    with open(SWEEP_CSV, "w", newline="") as f:
        w_ = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w_.writeheader()
        w_.writerows(results)
    print(f"\n{len(results)}/{n_total} cases -> {SWEEP_CSV}")


# ── graph (validate + plot) ──────────────────────────────────────────────────
def load_rows():
    if not os.path.isfile(SWEEP_CSV):
        print(f"ERROR: {SWEEP_CSV} not found.")
        print("Run the sweep first: python3 examples/bench_jkr_adhesion/sweep.py start")
        sys.exit(1)
    with open(SWEEP_CSV) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        print("ERROR: no data in sweep results file.")
        sys.exit(1)
    return rows


def linear_fit_through_origin(xs, ys):
    """Least-squares slope for y = m x (intercept fixed at 0), plus R^2."""
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    m = sxy / sxx
    ybar = sum(ys) / len(ys)
    ss_res = sum((y - m * x) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return m, r2


def validate(rows):
    print("=" * 66)
    print(f"JKR Adhesion Pull-off Benchmark Validation ({ADHESION_MODEL.upper()})")
    print("=" * 66)
    print(f"  Pull-off law: F = {PULLOFF_COEFF} * pi * w * R*")
    print(f"  R* = R/2 = {R_EFF*1e3:.3f} mm   (two equal spheres)\n")

    ws, f_meas = [], []
    all_pass = True
    for row in rows:
        w = float(row["w"])
        f = float(row["f_pulloff"])
        f_theory = analytical_pulloff(w)
        err = abs(f - f_theory) / f_theory
        status = "PASS" if err <= PULLOFF_REL_TOL else "FAIL"
        all_pass &= status == "PASS"
        ws.append(w)
        f_meas.append(f)
        print(f"  w={w:<5} J/m^2:  F={f:.4e} vs {f_theory:.4e} N  "
              f"(err={err*100:.3f}%)  [{status}]")

    # Linearity: fit F = slope * w through the origin and compare to theory.
    slope, r2 = linear_fit_through_origin(ws, f_meas)
    slope_theory = PULLOFF_COEFF * math.pi * R_EFF
    slope_err = abs(slope - slope_theory) / slope_theory
    lin_status = "PASS" if r2 >= LINEAR_R2_MIN else "FAIL"
    slope_status = "PASS" if slope_err <= SLOPE_REL_TOL else "FAIL"
    all_pass &= r2 >= LINEAR_R2_MIN and slope_err <= SLOPE_REL_TOL

    print()
    print(f"  Linear fit F = slope * w (through origin):")
    print(f"    slope (measured): {slope:.6e} N/(J/m^2)")
    print(f"    slope (theory):   {slope_theory:.6e} N/(J/m^2)  "
          f"(err={slope_err*100:.3f}%)  [{slope_status}]")
    print(f"    R^2 = {r2:.6f}  [{lin_status}]")
    print()
    print("ALL CHECKS PASSED" if all_pass else "WARNING: one or more checks FAILED")
    return all_pass


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
        "font.size": 12, "axes.labelsize": 14, "axes.titlesize": 14,
        "legend.fontsize": 10, "figure.dpi": 150, "savefig.dpi": 150,
    })
    os.makedirs(PLOT_DIR, exist_ok=True)
    save = dict(bbox_inches="tight")

    ws = [float(r["w"]) for r in rows]
    f_meas = [float(r["f_pulloff"]) for r in rows]

    # ── Plot 1: pull-off force vs work of adhesion (measured + theory line) ──
    fig, ax = plt.subplots(figsize=(7, 5))
    w_line = np.linspace(0, max(ws) * 1.05, 100)
    ax.plot(w_line, [analytical_pulloff(w) * 1e3 for w in w_line], "k-",
            linewidth=2, label=f"{ADHESION_MODEL.upper()} theory: "
                               f"F = {PULLOFF_COEFF}π w R*")
    ax.plot(ws, [f * 1e3 for f in f_meas], "o", color="#1f77b4",
            markersize=7, label="DIRT (measured)")
    ax.set_xlabel("Work of adhesion w [J/m$^2$]")
    ax.set_ylabel("Pull-off force [mN]")
    ax.set_title("Adhesion Pull-off Force vs Work of Adhesion")
    ax.legend(loc="upper left")
    ax.set_xlim(0, max(ws) * 1.05)
    ax.set_ylim(0, None)
    fig.savefig(os.path.join(PLOT_DIR, "pulloff_vs_surface_energy.png"), **save)
    plt.close(fig)
    print(f"Saved: {PLOT_DIR}/pulloff_vs_surface_energy.png")

    # ── Plot 2: force-separation curves for each case ───────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(ws)))
    for (w, c) in zip(ws, colors):
        tpath = os.path.join(DATA_DIR, f"w{w}_trace.csv")
        if not os.path.isfile(tpath):
            continue
        seps, fns = [], []
        with open(tpath) as f:
            for r in csv.DictReader(f):
                fn = float(r["f_normal"])
                if fn == 0.0:
                    continue
                seps.append(float(r["separation"]) * 1e9)  # nm
                fns.append(fn * 1e3)                        # mN
        if seps:
            ax.plot(seps, fns, ".", markersize=2, color=c, label=f"w = {w}")
    ax.axhline(0.0, color="0.6", linewidth=0.8)
    ax.set_xlabel("Surface separation [nm]   (<0 overlap, >0 gap)")
    ax.set_ylabel("Contact normal force [mN]   (<0 tensile)")
    ax.set_title("Force vs Separation (pull-off = min in the gap regime)")
    ax.legend(loc="lower right", fontsize=8)
    # Zoom near the contact/gap transition where the adhesion plateau lives.
    ax.set_xlim(-50, 50)
    fig.savefig(os.path.join(PLOT_DIR, "force_separation.png"), **save)
    plt.close(fig)
    print(f"Saved: {PLOT_DIR}/force_separation.png")


def graph():
    rows = load_rows()
    ok = validate(rows)
    plot(rows)
    return ok


# ── dispatch ─────────────────────────────────────────────────────────────────
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
