#!/usr/bin/env python3
"""
DMT / SJKR cohesion benchmark driver.

Validates two adhesion models that are physically distinct from JKR, and
exercises DIRT's adhesion-model selection (`[dem] adhesion_model` + the material
`surface_energy` vs `cohesion_energy` columns):

  * DMT arm (adhesion_model = "dmt", surface_energy = w):
        F_pulloff = 2 * pi * w * R*                          (DMT)
    DMT adds a constant attractive force with NO extended interaction range, so
    the pull-off is realized *inside* geometric overlap: the net normal force
    F_n(delta) = (4/3) E* sqrt(R*) delta^1.5 - 2 pi w R* is most tensile as the
    overlap delta -> 0, tending to exactly -2 pi w R*. The recorder captures that
    peak tension. Contrast: JKR (bench_jkr_adhesion) gives 1.5 pi w R*, so the
    DMT/JKR pull-off ratio is 4/3 -- selecting the model changes the physics.
    This arm also runs one JKR reference case to demonstrate that selection
    concretely (measured ratio ~ 4/3).

  * SJKR arm (cohesion_energy = c, no surface_energy):
        F_coh(delta) = c * pi * R* * delta                  (SJKR area law)
    Cohesion proportional to the circular contact area A = pi R* delta -- linear
    in overlap, vanishing at separation (no constant pull-off). It is isolated by
    differencing an SJKR run against a pure-Hertz baseline at matched overlap.
    Both run at restitution = 1, so the normal force is conservative (no velocity
    damping) and a pure function of delta; the shared Hertz term cancels exactly
    in the difference, leaving F_coh(delta) = c pi R* delta.

Commands (from anywhere):
    python3 examples/bench_dmt_sjkr_cohesion/sweep.py generate   # write configs
    python3 examples/bench_dmt_sjkr_cohesion/sweep.py start      # build + run -> CSV
    python3 examples/bench_dmt_sjkr_cohesion/sweep.py graph      # validate + plot
    python3 examples/bench_dmt_sjkr_cohesion/sweep.py            # all three

LAMMPS is intentionally NOT run: DIRT's DMT/SJKR are simplified constant- and
area-force models with no exact LAMMPS counterpart (LAMMPS' jkr is the full
Maugis-area model with a different force-overlap law), so a code-to-code overlay
would compare different physics. Validation is against the analytical laws only.

Outputs:
    sweep/<case>/config.toml   DIRT configs                  (gitignored)
    data/*.csv                 measured results + traces      (gitignored)
    plots/*.png                final figures                  (tracked)

References:
  B.V. Derjaguin, V.M. Muller, Yu.P. Toporov, "Effect of contact deformations on
    the adhesion of particles", J. Colloid Interface Sci. 53:314-326, 1975. (DMT)
  K.L. Johnson, K. Kendall, A.D. Roberts, "Surface energy and the contact of
    elastic solids", Proc. R. Soc. Lond. A 324:301-313, 1971. (JKR reference)
  For the "simplified JKR" (SJKR) area-proportional cohesion model as implemented
    in DEM codes, see the LIGGGHTS `cohesion sjkr` model documentation.
"""

import os
import sys
import csv
import math
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_dmt_sjkr_cohesion"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
DMT_CSV = os.path.join(DATA_DIR, "dmt_results.csv")
SJKR_CSV = os.path.join(DATA_DIR, "sjkr_results.csv")

# ── Material / geometry (shared by all cases) ────────────────────────────────
YOUNGS_MOD = 70.0e9     # Pa — soda-lime glass
POISSON_RATIO = 0.22
RADIUS = 0.005          # m
DENSITY = 2500.0        # kg/m^3
APPROACH_VEL = 0.002    # m/s — slow, so the contact is finely resolved
R_EFF = RADIUS / 2.0    # two equal spheres -> R* = R/2

# Effective Young's modulus 1/E* = 2 (1-nu^2)/E  (equal materials).
E_EFF = 1.0 / (2.0 * (1.0 - POISSON_RATIO ** 2) / YOUNGS_MOD)

# ── DMT arm: sweep the work of adhesion w = surface_energy [J/m^2] ────────────
DMT_W = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
# One JKR reference case (same w) to demonstrate model selection.
JKR_REF_W = 1.0

# ── SJKR arm: sweep the cohesion energy density c [J/m^3 = Pa] ────────────────
SJKR_C = [1.0e6, 2.0e6, 5.0e6, 1.0e7, 2.0e7]

# ── Validation tolerances ────────────────────────────────────────────────────
DMT_REL_TOL = 0.02        # 2% per-case: measured vs analytical DMT pull-off
DMT_SLOPE_TOL = 0.02      # fitted slope F(w) vs 2*pi*R*
DMT_R2_MIN = 0.999
RATIO_TOL = 0.02          # DMT/JKR pull-off ratio vs 4/3
SJKR_SLOPE_TOL = 0.05     # per-case fitted slope F_coh(delta) vs c*pi*R*
SJKR_R2_MIN = 0.995       # linearity of F_coh in delta
SJKR_C_SLOPE_TOL = 0.05   # slope-of-slopes vs pi*R* (linearity in c)

TOML_TEMPLATE = """\
# Auto-generated {arm} config
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
restitution = 1.0
friction = 0.0
{energy_line}

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
thermo = 400000
"""


def steps_for():
    """Step budget: travel the 1.05 mm gap at APPROACH_VEL, plus margin for the
    contact and separation, with the solver default timestep."""
    g = YOUNGS_MOD / (2.0 * (1.0 + POISSON_RATIO))
    alpha = 0.1631 * POISSON_RATIO + 0.876605
    dt_rayleigh = math.pi * RADIUS / alpha * (DENSITY / g) ** 0.5
    dt = dt_rayleigh * 0.15
    travel_time = 0.00105 / APPROACH_VEL
    total_time = travel_time * 1.4   # margin for contact + separation
    return int(total_time / dt) + 50000


def write_config(cdir, arm, adhesion_model, energy_line, output_dir):
    os.makedirs(cdir, exist_ok=True)
    with open(os.path.join(cdir, "config.toml"), "w") as f:
        f.write(TOML_TEMPLATE.format(
            arm=arm,
            adhesion_model=adhesion_model,
            energy_line=energy_line,
            youngs_mod=f"{YOUNGS_MOD:.6e}",
            poisson_ratio=POISSON_RATIO,
            radius=RADIUS, density=DENSITY,
            approach_vel=APPROACH_VEL,
            output_dir=output_dir, steps=steps_for(),
        ))


# ── case identity ────────────────────────────────────────────────────────────
def dmt_case(w):     return f"dmt_w{w}"
def jkr_case(w):     return f"jkr_w{w}"
def sjkr_case(c):    return f"sjkr_c{c:.3e}"
SJKR_BASELINE = "sjkr_hertz_baseline"   # c = 0, shared by all SJKR cases


# ── generate ─────────────────────────────────────────────────────────────────
def generate():
    n = 0
    for w in DMT_W:
        d = os.path.join(SWEEP_DIR, dmt_case(w))
        write_config(d, f"DMT w={w}", "dmt", f"surface_energy = {w}", d)
        n += 1
    # JKR reference (default model) at one w, to demonstrate model selection.
    d = os.path.join(SWEEP_DIR, jkr_case(JKR_REF_W))
    write_config(d, f"JKR ref w={JKR_REF_W}", "jkr", f"surface_energy = {JKR_REF_W}", d)
    n += 1
    # SJKR baseline (pure Hertz, c=0) — shared reference for the difference.
    # adhesion_model is irrelevant here (no surface_energy); use the default.
    d = os.path.join(SWEEP_DIR, SJKR_BASELINE)
    write_config(d, "SJKR Hertz baseline", "jkr", "", d)
    n += 1
    for c in SJKR_C:
        d = os.path.join(SWEEP_DIR, sjkr_case(c))
        write_config(d, f"SJKR c={c}", "jkr", f"cohesion_energy = {c}", d)
        n += 1
    print(f"Generated {n} DIRT configs under {SWEEP_DIR}")


# ── start ────────────────────────────────────────────────────────────────────
def run_case(tag, i, n_total):
    cdir = os.path.join(SWEEP_DIR, tag)
    config = os.path.join(cdir, "config.toml")
    if not os.path.isfile(config):
        print(f"  [{i}/{n_total}] missing {config} — run 'generate' first.")
        return None
    print(f"  [{i}/{n_total}] {tag:<24}", end="  ", flush=True)
    log = os.path.join(cdir, "run.log")
    with open(log, "w") as lf:
        proc = subprocess.run(
            ["cargo", "run", "--release", "--example", EXAMPLE,
             "--no-default-features", "--features", "precision-double", "--", config],
            cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
        )
    res = os.path.join(cdir, "data", "cohesion_results.csv")
    trace = os.path.join(cdir, "data", "cohesion_trace.csv")
    if proc.returncode == 0 and os.path.isfile(res):
        # Copy the trace next to the aggregate data for graphing.
        if os.path.isfile(trace):
            shutil.copyfile(trace, os.path.join(DATA_DIR, f"{tag}_trace.csv"))
        with open(res) as f:
            row = next(csv.DictReader(f))
        print(f"peak tension = {float(row['f_peak_tension']):.4e} N")
        return row
    print(f"FAILED ({log})")
    return None


def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    for p in (DMT_CSV, SJKR_CSV):
        if os.path.isfile(p):
            os.remove(p)

    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE,
         "--no-default-features", "--features", "precision-double"],
        cwd=REPO_ROOT, check=True,
    )

    tags = [dmt_case(w) for w in DMT_W] + [jkr_case(JKR_REF_W)] \
        + [SJKR_BASELINE] + [sjkr_case(c) for c in SJKR_C]
    n_total = len(tags)

    dmt_rows, jkr_row = [], None
    ok = True
    for i, tag in enumerate(tags, 1):
        row = run_case(tag, i, n_total)
        if row is None:
            ok = False
            continue
        if tag.startswith("dmt_"):
            row["w"] = tag[len("dmt_w"):]
            dmt_rows.append(row)
        elif tag.startswith("jkr_"):
            jkr_row = row

    # DMT arm summary (+ the JKR reference peak, tagged w = "jkr_ref").
    with open(DMT_CSV, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["w", "f_peak_tension", "sep_at_peak",
                                           "r_eff", "radius", "density", "dt"])
        wr.writeheader()
        for row in dmt_rows:
            wr.writerow({k: row.get(k, "") for k in wr.fieldnames})
        if jkr_row is not None:
            jkr_row["w"] = "jkr_ref"
            wr.writerow({k: jkr_row.get(k, "") for k in wr.fieldnames})

    # SJKR arm summary: just record which c-values ran (traces hold the data).
    with open(SJKR_CSV, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["c"])
        wr.writeheader()
        for c in SJKR_C:
            if os.path.isfile(os.path.join(DATA_DIR, f"{sjkr_case(c)}_trace.csv")):
                wr.writerow({"c": repr(c)})

    print(f"\nDMT -> {DMT_CSV}   SJKR -> {SJKR_CSV}")
    if not ok:
        print("ERROR: one or more cases failed.")
        sys.exit(1)


# ── trace helpers ────────────────────────────────────────────────────────────
def load_approach_curve(tag):
    """Return (deltas, forces) on the *approach* (v_normal < 0), overlap only,
    with delta = -separation monotonically increasing. With restitution = 1 the
    normal force is a pure function of delta, so approach and retreat coincide;
    the approach branch is monotonic and clean to interpolate."""
    path = os.path.join(DATA_DIR, f"{tag}_trace.csv")
    if not os.path.isfile(path):
        return None, None
    pts = []
    with open(path) as f:
        for r in csv.DictReader(f):
            sep = float(r["separation"])
            vn = float(r["v_normal"])
            fn = float(r["f_normal"])
            if sep < 0.0 and vn < 0.0:      # overlap, approaching
                pts.append((-sep, fn))       # (delta, f_normal)
    if len(pts) < 5:
        return None, None
    pts.sort()
    deltas = [p[0] for p in pts]
    forces = [p[1] for p in pts]
    return deltas, forces


def interp(xs, ys, x):
    """Linear interpolation of y(x) on a sorted xs; None outside range."""
    if x < xs[0] or x > xs[-1]:
        return None
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    if xs[hi] == xs[lo]:
        return ys[lo]
    t = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] + t * (ys[hi] - ys[lo])


def fit_through_origin(xs, ys):
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    m = sxy / sxx
    ybar = sum(ys) / len(ys)
    ss_res = sum((y - m * x) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return m, r2


# ── validate ─────────────────────────────────────────────────────────────────
def dmt_analytical(w):
    return 2.0 * math.pi * w * R_EFF


def validate_dmt():
    print("=" * 70)
    print("DMT Adhesion Pull-off  —  F_pulloff = 2*pi*w*R*")
    print("=" * 70)
    print(f"  R* = R/2 = {R_EFF*1e3:.3f} mm   E* = {E_EFF:.4e} Pa\n")
    if not os.path.isfile(DMT_CSV):
        print(f"ERROR: {DMT_CSV} not found. Run 'start' first.")
        sys.exit(1)
    with open(DMT_CSV) as f:
        rows = list(csv.DictReader(f))

    ws, f_meas, jkr_peak = [], [], None
    all_pass = True
    for row in rows:
        if row["w"] == "jkr_ref":
            jkr_peak = float(row["f_peak_tension"])
            continue
        w = float(row["w"])
        fm = float(row["f_peak_tension"])
        ft = dmt_analytical(w)
        err = abs(fm - ft) / ft
        status = "PASS" if err <= DMT_REL_TOL else "FAIL"
        all_pass &= status == "PASS"
        ws.append(w); f_meas.append(fm)
        print(f"  w={w:<5} J/m^2:  F={fm:.4e} vs {ft:.4e} N  (err={err*100:.3f}%)  [{status}]")

    slope, r2 = fit_through_origin(ws, f_meas)
    slope_theory = 2.0 * math.pi * R_EFF
    slope_err = abs(slope - slope_theory) / slope_theory
    lin_ok = r2 >= DMT_R2_MIN and slope_err <= DMT_SLOPE_TOL
    all_pass &= lin_ok
    print(f"\n  Linear fit F = slope * w (through origin):")
    print(f"    slope (meas)  : {slope:.6e}")
    print(f"    slope (theory): {slope_theory:.6e}  (err={slope_err*100:.3f}%)  "
          f"[{'PASS' if slope_err <= DMT_SLOPE_TOL else 'FAIL'}]")
    print(f"    R^2 = {r2:.6f}  [{'PASS' if r2 >= DMT_R2_MIN else 'FAIL'}]")

    # Model-selection demonstration: DMT vs JKR at the same w.
    if jkr_peak is not None:
        dmt_at_ref = dmt_analytical(JKR_REF_W)
        # measured DMT peak at w = JKR_REF_W (interpolate from the fit line).
        dmt_meas_ref = slope * JKR_REF_W
        ratio = dmt_meas_ref / jkr_peak
        ratio_err = abs(ratio - 4.0 / 3.0) / (4.0 / 3.0)
        rstatus = "PASS" if ratio_err <= RATIO_TOL else "FAIL"
        all_pass &= ratio_err <= RATIO_TOL
        print(f"\n  Adhesion-model selection (w={JKR_REF_W} J/m^2):")
        print(f"    DMT peak  = {dmt_meas_ref:.4e} N   (theory {dmt_at_ref:.4e})")
        print(f"    JKR peak  = {jkr_peak:.4e} N   (theory {1.5*math.pi*JKR_REF_W*R_EFF:.4e})")
        print(f"    ratio DMT/JKR = {ratio:.4f}  vs  4/3 = {4/3:.4f}  "
              f"(err={ratio_err*100:.3f}%)  [{rstatus}]")
    print()
    return all_pass, (ws, f_meas)


def sjkr_analytical_slope(c):
    return c * math.pi * R_EFF


def validate_sjkr():
    print("=" * 70)
    print("SJKR Cohesion  —  F_coh(delta) = c*pi*R**delta   (area law)")
    print("=" * 70)
    base_d, base_f = load_approach_curve(SJKR_BASELINE)
    if base_d is None:
        print("ERROR: missing/short SJKR Hertz baseline trace. Run 'start' first.")
        sys.exit(1)

    cs, slopes = [], []
    all_pass = True
    per_case = {}
    for c in SJKR_C:
        d_c, f_c = load_approach_curve(sjkr_case(c))
        if d_c is None:
            print(f"  c={c:.2e}:  missing/short trace  [FAIL]")
            all_pass = False
            continue
        # Common overlap range; skip the smallest 5% (float noise near delta=0).
        d_hi = min(base_d[-1], d_c[-1])
        d_lo = 0.05 * d_hi
        grid = [d_lo + (d_hi - d_lo) * k / 40.0 for k in range(41)]
        xs, ys = [], []
        for d in grid:
            fb = interp(base_d, base_f, d)
            fs = interp(d_c, f_c, d)
            if fb is None or fs is None:
                continue
            xs.append(d)
            ys.append(fb - fs)          # F_coh(delta) = Hertz - SJKR
        if len(xs) < 10:
            print(f"  c={c:.2e}:  too few matched points  [FAIL]")
            all_pass = False
            continue
        slope, r2 = fit_through_origin(xs, ys)
        slope_th = sjkr_analytical_slope(c)
        serr = abs(slope - slope_th) / slope_th
        status = "PASS" if (serr <= SJKR_SLOPE_TOL and r2 >= SJKR_R2_MIN) else "FAIL"
        all_pass &= status == "PASS"
        cs.append(c); slopes.append(slope)
        per_case[c] = (xs, ys, slope, slope_th)
        print(f"  c={c:.2e} Pa:  slope={slope:.4e} vs {slope_th:.4e} N/m  "
              f"(err={serr*100:.2f}%)  R^2={r2:.5f}  [{status}]")

    # Linearity of the fitted slope in c: slope(c) = (pi R*) c through origin.
    if len(cs) >= 2:
        m2, r2_2 = fit_through_origin(cs, slopes)
        m2_th = math.pi * R_EFF
        m2_err = abs(m2 - m2_th) / m2_th
        cstat = "PASS" if (m2_err <= SJKR_C_SLOPE_TOL and r2_2 >= DMT_R2_MIN) else "FAIL"
        all_pass &= cstat == "PASS"
        print(f"\n  Slope-of-slopes (linearity in c):")
        print(f"    d(slope)/dc (meas)  : {m2:.6e}")
        print(f"    pi*R* (theory)      : {m2_th:.6e}  (err={m2_err*100:.3f}%)  [{cstat}]")
        print(f"    R^2 = {r2_2:.6f}")
    print()
    return all_pass, per_case


# ── plot ─────────────────────────────────────────────────────────────────────
def plot(dmt_data, sjkr_data):
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(matplotlib/numpy unavailable, skipped plots: {e})")
        return
    plt.rcParams.update({"font.size": 12, "axes.labelsize": 13,
                         "axes.titlesize": 13, "legend.fontsize": 10,
                         "figure.dpi": 150, "savefig.dpi": 150})
    os.makedirs(PLOT_DIR, exist_ok=True)
    save = dict(bbox_inches="tight")

    # Plot 1: DMT pull-off vs w (with JKR theory line for contrast).
    ws, f_meas = dmt_data
    fig, ax = plt.subplots(figsize=(7, 5))
    wl = np.linspace(0, max(ws) * 1.05, 100)
    ax.plot(wl, [2 * math.pi * w * R_EFF * 1e3 for w in wl], "k-", lw=2,
            label="DMT theory: F = 2π w R*")
    ax.plot(wl, [1.5 * math.pi * w * R_EFF * 1e3 for w in wl], "k--", lw=1.3,
            label="JKR theory: F = 1.5π w R*")
    ax.plot(ws, [f * 1e3 for f in f_meas], "o", color="#d62728", ms=8,
            label="DIRT DMT (measured)")
    ax.set_xlabel("Work of adhesion w [J/m$^2$]")
    ax.set_ylabel("Pull-off force [mN]")
    ax.set_title("DMT Pull-off Force vs Work of Adhesion")
    ax.legend(loc="upper left"); ax.set_xlim(0, max(ws) * 1.05); ax.set_ylim(0, None)
    fig.savefig(os.path.join(PLOT_DIR, "dmt_pulloff_vs_w.png"), **save)
    plt.close(fig)
    print(f"Saved: {PLOT_DIR}/dmt_pulloff_vs_w.png")

    # Plot 2: SJKR isolated cohesion force vs overlap (area law).
    if sjkr_data:
        fig, ax = plt.subplots(figsize=(7, 5))
        cs = sorted(sjkr_data.keys())
        colors = plt.cm.viridis(np.linspace(0, 0.85, len(cs)))
        for c, col in zip(cs, colors):
            xs, ys, slope, slope_th = sjkr_data[c]
            ax.plot([x * 1e9 for x in xs], [y * 1e3 for y in ys], ".",
                    ms=4, color=col, label=f"c={c:.0e} Pa")
            dl = np.array([0, max(xs)])
            ax.plot(dl * 1e9, slope_th * dl * 1e3, "-", color=col, lw=1)
        ax.set_xlabel("Overlap δ [nm]")
        ax.set_ylabel("Isolated cohesion force F$_{coh}$ [mN]")
        ax.set_title("SJKR Cohesion (Hertz−SJKR difference) vs Overlap\n"
                     "points = DIRT, lines = c·π·R*·δ theory")
        ax.legend(loc="upper left", fontsize=8)
        fig.savefig(os.path.join(PLOT_DIR, "sjkr_cohesion_vs_overlap.png"), **save)
        plt.close(fig)
        print(f"Saved: {PLOT_DIR}/sjkr_cohesion_vs_overlap.png")


# ── graph / dispatch ─────────────────────────────────────────────────────────
def graph():
    dmt_ok, dmt_data = validate_dmt()
    sjkr_ok, sjkr_data = validate_sjkr()
    plot(dmt_data, sjkr_data)
    ok = dmt_ok and sjkr_ok
    print("=" * 70)
    print("ALL CHECKS PASSED" if ok else "WARNING: one or more checks FAILED")
    print("=" * 70)
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
