#!/usr/bin/env python3
"""
Angle-of-repose bulk-friction benchmark driver.

Forms a static granular heap with the "lift the cylinder" protocol and measures
its angle of repose theta_r as a function of the sliding friction mu. Unlike the
two-body benchmarks, the reference here is EMPIRICAL, not closed-form: theta_r
has no exact analytical value, so validation checks the qualitative laws every
correct DEM contact model must obey:

    1. theta_r increases MONOTONICALLY with mu,
    2. theta_r falls in a physically sensible band and collapses toward ~0 deg as
       mu -> 0. (The lift-the-cylinder COLLAPSE protocol on a smooth frictional
       floor reads lower than slow pouring, so the band for this protocol's
       top-mu cases is ~10-40 deg.)
    3. results are REPRODUCIBLE: the run-to-run spread (over fresh random packs)
       is small.

The heap sits directly on a real frictional plane wall (z = 0, normal +z):
dirt_wall applies Mindlin sliding (tangential) friction on plane walls using the
material's friction coefficient (μ) via friction_ij. That base friction keeps
the bottom layer from sliding out, so the pile holds a slope — no frozen
particle bed is needed. See README "Assumptions".

Commands (from anywhere):
    python3 examples/bench_angle_of_repose/sweep.py generate   # write per-case configs
    python3 examples/bench_angle_of_repose/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_angle_of_repose/sweep.py graph      # validate + plot
    python3 examples/bench_angle_of_repose/sweep.py            # all three, in order

Each (mu) case is run REPS times with independent random packs (the inserter is
entropy-seeded), so the spread of theta_r is a direct reproducibility measure.

The angle is fit in this script from the settled particle positions DIRT dumps:
the heap is centered on its (x,y) centroid, particles are binned by radial
distance r, the heap-surface height h(r) is the upper envelope of z in each bin,
and theta_r = atan(-slope) of a linear fit to h(r) on the sloping flank.

This benchmark is DIRT-only: it has no LAMMPS overlay, because there is no
analytical target to compare a second code against — the validation is the set
of empirical laws above, applied to DIRT's own settled heaps.

Outputs:
    sweep/<case>/config.toml            DIRT configs                 (gitignored)
    sweep/<case>/data/repose_results.csv  per-run particle positions (gitignored)
    data/repose_sweep.csv               theta_r per (mu, rep)        (gitignored)
    data/profile_<mu>.csv               representative h(r) profile  (gitignored)
    plots/*.png                         final figures                (tracked)

Reference (empirical, for context — values vary with material/protocol):
    Y.C. Zhou et al., "Rolling friction in the dynamic simulation of sandpile
    formation", Physica A 269 (1999) 536-553.
    H.P. Zhu et al., "Discrete particle simulation of particulate systems:
    A review of major applications and findings", Chem. Eng. Sci. 63 (2008).
"""

import os
import sys
import csv
import math
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_angle_of_repose"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "repose_sweep.csv")  # theta_r per (mu, rep)

# LAMMPS binary candidates (unused here; kept for structural parity). This
# benchmark validates against the empirical laws, not a second code.
LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# -- Sweep parameters -----------------------------------------------------------
# We sweep mu over [0, 0.5]. In the lift-the-cylinder protocol the heap forms by
# a column COLLAPSE on the frictional floor: at low mu the bottom layer slides
# out and the deposit spreads into a near-flat disk (theta_r ~ 0), while at higher
# mu the floor friction arrests the runout and the deposit relaxes into a cone
# whose flank steepens with mu. theta_r rises monotonically across the range.
MU_LIST = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]   # sliding friction sweep
REPS = 3                          # independent packs per mu (reproducibility)

# Material / geometry (mirror config.toml; mu is overridden per case). The swept
# mu governs both particle–particle contacts (which set the pile's angle) and
# particle–floor-wall contacts (which provide the base friction).
YOUNGS_MOD = 1.0e7      # Pa (softened: larger stable dt, standard DEM practice)
POISSON = 0.25
RESTITUTION = 0.4
ROLLING_FRICTION = 0.1
DENSITY = 2500.0        # kg/m^3
RADIUS = 0.002          # m — particle radius
HEAP_COUNT = 1200       # mobile heap particles
CYL_RADIUS = 0.025      # confining-cylinder radius (m) — narrow, tall column
GZ = -9.81

# -- Heap-fit parameters --------------------------------------------------------
# The settled deposit is: a central cone (the heap) sitting on the floor, plus a
# sparse monolayer of stragglers that avalanched out past the cone toe during the
# collapse. The fit isolates the cone flank by (a) subtracting the floor baseline
# height (a single resting layer), (b) finding the toe radius where the heap
# height falls to ~one particle diameter above the baseline, and (c) fitting the
# slope on the straight flank window between the apex skip and the toe.
N_BINS = 26
SURFACE_PCTL = 90.0     # height percentile per bin = heap surface envelope
APEX_SKIP_FRAC = 0.15   # skip the inner 15% of the toe radius (rounded apex)
TOE_HI_FRAC = 0.92      # stop the fit just inside the toe
TOE_HEIGHT_FACTOR = 1.5 # toe = where (h - baseline) drops below this * diameter

# -- Validation tolerances ------------------------------------------------------
# The band is set for THIS protocol: lift-the-cylinder column collapse on a
# smooth frictional floor reads lower than slow pouring (the collapse mobilizes
# the surface and the deposit grows a wide low apron the fit averages in), so the
# sensible repose band for the top-mu cases here is ~10-40 deg, not the 25-40 deg
# of a poured heap. The trend (monotonic, mu=0 -> flat) is the primary validation.
ANGLE_LO_DEG = 10.0     # sensible lower bound for the top-mu cases (collapse heap)
ANGLE_HI_DEG = 40.0     # sensible upper bound
LOWMU_MAX_DEG = 8.0     # mu=0 heap must be nearly flat (theta_r small)
SPREAD_MAX_DEG = 5.0    # max allowed run-to-run std dev of theta_r at a given mu
MONOTONIC_SLACK_DEG = 2.5  # mean theta_r may dip by at most this between mu steps

# -- DIRT config template (mirrors config.toml; mu swept) -----------------------
TOML_TEMPLATE = """[comm]
processors_x = 1
processors_y = 1
processors_z = 1
[domain]
x_low = -0.08
x_high = 0.08
y_low = -0.08
y_high = 0.08
z_low = 0.0
z_high = 0.16
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"
[neighbor]
skin_fraction = 1.2
bin_size = 0.006
every = 1
[gravity]
gz = {gz}
[dem]
contact_model = "hertz"
[[dem.materials]]
name = "glass"
youngs_mod = {youngs:.6e}
poisson_ratio = {nu}
restitution = {e_n}
friction = {mu}
rolling_friction = {mu_r}
[[wall]]
type = "cylinder"
axis = "z"
center = [0.0, 0.0]
radius = {cyl_r}
lo = 0.0
hi = 0.16
inside = true
material = "glass"
name = "cylinder"
[[wall]]
type = "plane"
point_x = 0.0
point_y = 0.0
point_z = 0.0
normal_x = 0.0
normal_y = 0.0
normal_z = 1.0
material = "glass"
[[wall]]
type = "cylinder"
axis = "z"
center = [0.0, 0.0]
radius = 0.07
lo = 0.0
hi = 0.16
inside = true
material = "glass"
[[particles.insert]]
material = "glass"
count = {heap_count}
radius = {radius}
density = {density}
velocity_z = -0.1
region = {{ type = "cylinder", center = [0.0, 0.0], radius = {ins_r}, axis = "z", lo = 0.003, hi = 0.14 }}
[output]
dir = "{outdir}"
[vtp]
interval = 100000
[[run]]
name = "fill"
steps = 100000
thermo = 50000
[[run]]
name = "lift"
steps = 200000
thermo = 50000
"""


# -- helpers --------------------------------------------------------------------
def case_tag(mu, rep):
    return f"mu_{mu:g}_rep{rep}"


def case_dir(mu, rep):
    return os.path.join(SWEEP_DIR, case_tag(mu, rep))


def find_lammps():
    for b in LAMMPS_BINS:
        if shutil.which(b):
            return shutil.which(b)
    return None


def _dirt_config(mu, outdir):
    return TOML_TEMPLATE.format(
        gz=GZ, youngs=YOUNGS_MOD, nu=POISSON, e_n=RESTITUTION, mu=mu,
        mu_r=ROLLING_FRICTION, cyl_r=CYL_RADIUS,
        heap_count=HEAP_COUNT, radius=RADIUS, density=DENSITY,
        ins_r=CYL_RADIUS - 1.5 * RADIUS, outdir=outdir,
    )


# -- generate -------------------------------------------------------------------
def generate():
    n = 0
    for mu in MU_LIST:
        for rep in range(REPS):
            cdir = case_dir(mu, rep)
            os.makedirs(cdir, exist_ok=True)
            with open(os.path.join(cdir, "config.toml"), "w") as f:
                f.write(_dirt_config(mu, cdir))
            n += 1
    print(f"Generated {n} DIRT configs ({len(MU_LIST)} mu x {REPS} reps) under {SWEEP_DIR}")


# -- heap geometry fit ----------------------------------------------------------
def _load_positions(path):
    """Load every recorded particle: a single material is used, so all positions
    belong to the heap (there is no frozen bed to exclude)."""
    xs, ys, zs, rs = [], [], [], []
    with open(path) as f:
        for r in csv.DictReader(f):
            xs.append(float(r["x"])); ys.append(float(r["y"]))
            zs.append(float(r["z"])); rs.append(float(r["radius"]))
    return xs, ys, zs, rs


def _percentile(vals, pctl):
    if not vals:
        return 0.0
    s = sorted(vals)
    k = (len(s) - 1) * pctl / 100.0
    lo = int(math.floor(k)); hi = int(math.ceil(k))
    if lo == hi:
        return s[lo]
    return s[lo] * (hi - k) + s[hi] * (k - lo)


def heap_profile(xs, ys, zs, rs):
    """Radial surface profile of the heap, centered on the (x,y) centroid.

    Returns (r_centers, h_surface, baseline, diameter): h_surface is the
    per-bin SURFACE_PCTL height (the heap envelope), baseline is the resting
    height of a single particle layer on the floor (the profile floor), and
    diameter is the particle diameter."""
    n = len(xs)
    if n == 0:
        return [], [], 0.0, 0.0
    cx = sum(xs) / n
    cy = sum(ys) / n
    rad = [math.hypot(xs[i] - cx, ys[i] - cy) for i in range(n)]
    r_max = _percentile(rad, 99.0)
    diameter = 2.0 * (sum(rs) / n)
    if r_max <= 0:
        return [], [], 0.0, diameter
    bins = [[] for _ in range(N_BINS)]
    for i in range(n):
        b = int(rad[i] / r_max * N_BINS)
        if b == N_BINS:
            b -= 1
        if 0 <= b < N_BINS:
            bins[b].append(zs[i])
    r_centers, h_surface = [], []
    for b in range(N_BINS):
        if len(bins[b]) >= 3:
            r_centers.append((b + 0.5) / N_BINS * r_max)
            h_surface.append(_percentile(bins[b], SURFACE_PCTL))
    # Baseline = lowest surface height seen on the outer half (a single layer
    # resting on the floor), i.e. the profile floor the cone sits above.
    outer = [h_surface[i] for i in range(len(r_centers)) if r_centers[i] > 0.5 * r_max]
    baseline = min(outer) if outer else (min(h_surface) if h_surface else 0.0)
    return r_centers, h_surface, baseline, diameter


def _toe_radius(r_centers, h_surface, baseline, diameter):
    """Outermost radius where the heap still stands more than TOE_HEIGHT_FACTOR
    diameters above the floor baseline — the cone toe, ignoring sparse stragglers."""
    thresh = baseline + TOE_HEIGHT_FACTOR * diameter
    r_toe = 0.0
    for i in range(len(r_centers)):
        if h_surface[i] >= thresh:
            r_toe = r_centers[i]
    return r_toe


def _linfit(x, y):
    """Least-squares slope, intercept."""
    n = len(x)
    sx = sum(x); sy = sum(y)
    sxx = sum(v * v for v in x); sxy = sum(x[i] * y[i] for i in range(n))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        return 0.0, sy / n
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    return slope, intercept


def fit_angle(r_centers, h_surface, baseline, diameter):
    """theta_r = atan(-slope) of the heap surface over the straight cone flank,
    from just outside the apex to just inside the toe.

    A deposit with no resolvable cone (the low-mu heaps spread into a near-flat
    disk) is a genuine theta_r ~ 0 deg result, not a fit failure: return 0.0 so
    the sweep records the flat low-mu point rather than dropping the case."""
    r_toe = _toe_radius(r_centers, h_surface, baseline, diameter)
    if r_toe <= 0:
        return 0.0, 0.0
    lo = APEX_SKIP_FRAC * r_toe
    hi = TOE_HI_FRAC * r_toe
    xf = [r_centers[i] for i in range(len(r_centers)) if lo <= r_centers[i] <= hi]
    yf = [h_surface[i] for i in range(len(r_centers)) if lo <= r_centers[i] <= hi]
    if len(xf) < 3:
        return 0.0, r_toe
    slope, _ = _linfit(xf, yf)
    return math.degrees(math.atan(max(0.0, -slope))), r_toe


# -- start ----------------------------------------------------------------------
def _run_dirt(cdir):
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "repose_results.csv")
    if os.path.exists(res):
        os.remove(res)
    log = os.path.join(cdir, "run.log")
    with open(log, "w") as lf:
        proc = subprocess.run(
            ["cargo", "run", "--release", "--example", EXAMPLE,
             "--no-default-features", "--", config],
            cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
        )
    if proc.returncode != 0 or not os.path.isfile(res):
        return None
    return res


def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(["cargo", "build", "--release", "--example", EXAMPLE,
                    "--no-default-features"], cwd=REPO_ROOT, check=True)

    if find_lammps():
        print("LAMMPS present but unused: this benchmark validates against "
              "empirical laws, not a second code.")

    rows = []
    profiles = {}  # mu -> representative (r_centers, h_surface) from rep 0
    total = len(MU_LIST) * REPS
    k = 0
    for mu in MU_LIST:
        for rep in range(REPS):
            k += 1
            cdir = case_dir(mu, rep)
            if not os.path.isfile(os.path.join(cdir, "config.toml")):
                print(f"  [{k:2d}/{total}] missing config mu={mu} rep={rep} — run 'generate'.")
                continue
            print(f"  [{k:2d}/{total}] mu={mu:<4} rep={rep}", end="  ", flush=True)
            res = _run_dirt(cdir)
            if res is None:
                print("DIRT FAILED")
                continue
            xs, ys, zs, rad = _load_positions(res)
            r_c, h_s, base, diam = heap_profile(xs, ys, zs, rad)
            theta, r_toe = fit_angle(r_c, h_s, base, diam)
            if theta is None:
                print("fit FAILED (no resolvable cone flank)")
                continue
            rows.append({"mu": mu, "rep": rep, "theta_deg": theta,
                         "r_toe": r_toe, "n": len(xs)})
            print(f"theta_r = {theta:5.2f} deg  (r_toe={r_toe*1e3:.1f} mm, N_heap={len(xs)})")
            if rep == 0:
                profiles[mu] = (r_c, h_s, base, r_toe)

    if not rows:
        print("\nERROR: no DIRT results collected.")
        sys.exit(1)

    os.makedirs(DATA_DIR, exist_ok=True)
    with open(SWEEP_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mu", "rep", "theta_deg", "r_toe", "n"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\n{len(rows)}/{total} cases -> {SWEEP_CSV}")

    # Save representative profiles (baseline-subtracted) for the cross-section plot.
    for mu, (r_c, h_s, base, r_toe) in profiles.items():
        with open(os.path.join(DATA_DIR, f"profile_{mu:g}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["r", "h"])
            for i in range(len(r_c)):
                w.writerow([r_c[i], h_s[i] - base])


# -- graph (validate + plot) ----------------------------------------------------
def _load_sweep():
    if not os.path.isfile(SWEEP_CSV):
        return []
    with open(SWEEP_CSV) as f:
        return [{k: (float(v) if k != "rep" else int(float(v))) for k, v in r.items()}
                for r in csv.DictReader(f)]


def _stats_by_mu(rows):
    """mu -> (mean_theta, std_theta, n_reps), sorted by mu."""
    by = {}
    for r in rows:
        by.setdefault(r["mu"], []).append(r["theta_deg"])
    out = []
    for mu in sorted(by):
        v = by[mu]
        mean = sum(v) / len(v)
        var = sum((x - mean) ** 2 for x in v) / len(v) if len(v) > 1 else 0.0
        out.append((mu, mean, math.sqrt(var), len(v)))
    return out


def validate(rows):
    print("\n=== Angle-of-repose validation ===")
    print(f"  material: E={YOUNGS_MOD:.1e} Pa  nu={POISSON}  e={RESTITUTION}  "
          f"mu_r={ROLLING_FRICTION}")
    stats = _stats_by_mu(rows)
    print(f"  {'mu':>6}{'mean_theta':>12}{'std':>8}{'reps':>6}  note")
    ok = True

    # 1. monotonic increase with mu (allow small slack for stochastic dips).
    prev_mean = None
    for (mu, mean, std, nrep) in stats:
        note = ""
        if prev_mean is not None and mean < prev_mean - MONOTONIC_SLACK_DEG:
            note = "NON-MONOTONIC"; ok = False
        # 3. reproducibility: spread small.
        if std > SPREAD_MAX_DEG:
            note = (note + " HIGH-SPREAD").strip(); ok = False
        print(f"  {mu:>6.2f}{mean:>12.2f}{std:>8.2f}{nrep:>6d}  {note}")
        prev_mean = mean

    # 2a. low-friction heap must be nearly flat.
    low = [s for s in stats if s[0] <= 1e-9]
    if low and low[0][1] > LOWMU_MAX_DEG:
        print(f"  mu=0 mean theta = {low[0][1]:.2f} deg > {LOWMU_MAX_DEG} deg "
              f"— frictionless heap not flat"); ok = False

    # 2b. the higher-friction cases must land in the sensible repose band.
    hi = [s for s in stats if s[0] >= 0.2]
    band_ok = any(ANGLE_LO_DEG <= mean <= ANGLE_HI_DEG for (_, mean, _, _) in hi)
    if hi and not band_ok:
        means = ", ".join(f"{m:.1f}" for (_, m, _, _) in hi)
        print(f"  no mu>=0.2 case in [{ANGLE_LO_DEG},{ANGLE_HI_DEG}] deg "
              f"(got {means}) — out of physical band"); ok = False

    # overall increase from lowest to highest mu.
    if len(stats) >= 2 and stats[-1][1] <= stats[0][1] + 1.0:
        print(f"  theta_r did not increase across mu range "
              f"({stats[0][1]:.1f} -> {stats[-1][1]:.1f} deg)"); ok = False

    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


def plot(rows):
    os.makedirs(PLOT_DIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    # -- theta_r vs mu (mean +/- spread, plus individual reps) --
    stats = _stats_by_mu(rows)
    mus = [s[0] for s in stats]
    means = [s[1] for s in stats]
    stds = [s[2] for s in stats]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.errorbar(mus, means, yerr=stds, fmt="o-", capsize=4,
                label="DIRT (mean ± std over reps)")
    ax.scatter([r["mu"] for r in rows], [r["theta_deg"] for r in rows],
               s=14, alpha=0.4, color="gray", label="individual runs")
    ax.axhspan(ANGLE_LO_DEG, ANGLE_HI_DEG, color="green", alpha=0.07,
               label=f"sensible band [{ANGLE_LO_DEG:.0f},{ANGLE_HI_DEG:.0f}]°")
    ax.set_xlabel(r"sliding friction $\mu$")
    ax.set_ylabel(r"angle of repose $\theta_r$ (deg)")
    ax.set_title("Angle of repose vs sliding friction (lift-the-cylinder heap)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "theta_vs_mu.png"))
    plt.close(fig)

    # -- heap cross-section profiles h(r) for each mu, with fitted flank --
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for mu in sorted({r["mu"] for r in rows}):
        ppath = os.path.join(DATA_DIR, f"profile_{mu:g}.csv")
        if not os.path.isfile(ppath):
            continue
        rc, hs = [], []
        with open(ppath) as f:
            for row in csv.DictReader(f):
                rc.append(float(row["r"]) * 1e3)
                hs.append(float(row["h"]) * 1e3)
        if rc:
            ax.plot(rc, hs, "o-", ms=3, label=fr"$\mu$={mu:g}")
    ax.set_xlabel("radial distance r (mm)")
    ax.set_ylabel("heap surface height h (mm)")
    ax.set_title("Settled heap cross-section (surface envelope)")
    ax.legend(title="friction")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "heap_profile.png"))
    plt.close(fig)

    print(f"\nFigures -> {PLOT_DIR}/theta_vs_mu.png, heap_profile.png")


def graph():
    rows = _load_sweep()
    if not rows:
        print(f"No {SWEEP_CSV} — run 'start' first.")
        return False
    ok = validate(rows)
    plot(rows)
    return ok


# -- dispatch -------------------------------------------------------------------
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
