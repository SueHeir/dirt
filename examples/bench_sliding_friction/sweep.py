#!/usr/bin/env python3
"""
Coulomb sliding-friction benchmark driver.

Launches a single sphere horizontally (v0, zero spin) onto a flat frictional
floor under gravity and validates the classic rigid-body result:

  - sliding phase: the center decelerates at a constant a = mu*g while the sphere
    spins up at alpha = mu*g / ((2/5) R)  (contact point slides, v > omega*R),
  - slip ends at t* = 2 v0 / (7 mu g),
  - thereafter the sphere rolls without slipping at v_final = (5/7) v0,
    independent of mu.

The floor is a real `dirt_wall` z-plane at z = 0 (normal +z). `dirt_wall` applies
Mindlin tangential (sliding) friction with a Coulomb cap on plane walls, using
the material's `friction` coefficient — so the flat wall decelerates the sliding
sphere at a = mu*g and spins it up. The floor is perfectly flat (no curvature
systematic). See README "Assumptions".

Commands (from anywhere):
    python3 examples/bench_sliding_friction/sweep.py generate   # write per-case configs
    python3 examples/bench_sliding_friction/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_sliding_friction/sweep.py graph      # validate + plot
    python3 examples/bench_sliding_friction/sweep.py            # all three, in order

If a LAMMPS binary (lmp_serial / lmp / lmp_mpi / lammps) is on PATH, each case is
ALSO run in LAMMPS's granular model — one sphere of the same radius/density
launched horizontally with NO spin onto a flat frictional `wall/gran` floor
(`hertz/material` normal + `mindlin` tangential, Coulomb cap mu) under gravity.
The sliding-phase deceleration a and rolling plateau v_final are fit the SAME way
as DIRT and overlaid on the plots. a = mu*g is identical in both codes, so this
is a near-exact cross-check. LAMMPS is OPTIONAL — without it, only DIRT runs and
the benchmark still fully validates against the exact rigid-body theory.

Outputs:
    sweep/<case>/config.toml                       DIRT configs         (gitignored)
    sweep/<case>/data/sliding_friction_results.csv raw DIRT time series (gitignored)
    sweep/<case>/in.lammps                          LAMMPS inputs        (gitignored)
    sweep/<case>/lammps_series.txt                  raw LAMMPS dump      (gitignored)
    data/sweep_summary.csv                         per-case DIRT fit     (gitignored)
    data/lammps_results.csv                        per-case LAMMPS fit   (gitignored)
    plots/*.png                                    final figures         (tracked)
"""

import os
import sys
import csv
import math
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_sliding_friction"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SUMMARY_CSV = os.path.join(DATA_DIR, "sweep_summary.csv")     # DIRT fit summary
LAMMPS_CSV = os.path.join(DATA_DIR, "lammps_results.csv")     # LAMMPS fit summary

# LAMMPS binary candidates, in preference order. LAMMPS is optional.
LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# ── Physics constants ─────────────────────────────────────────────────────────
G = 9.81             # gravity (m/s^2)
R = 0.005            # projectile radius (m)
DENSITY = 2500.0     # kg/m^3 — glass
YOUNGS_MOD = 70.0e9  # Pa
NU = 0.22            # Poisson ratio
RESTITUTION = 0.3    # low e: vertical settling mode damps fast
DT = 2.0e-6          # s

# ── Sweep: friction mu (at fixed v0) and impact speed v0 (at fixed mu) ────────
MU_LIST = [0.2, 0.3, 0.5, 0.7]   # at v0 = V0_REF
V0_LIST = [0.5, 1.0, 1.5]        # at mu = MU_REF
MU_REF = 0.5
V0_REF = 1.0


def t_star(v0, mu):
    """Slip -> roll transition time."""
    return 2.0 * v0 / (7.0 * mu * G)


def n_steps(v0, mu):
    """Run long enough to see slip end at t* plus a rolling plateau (~1.6 t*)."""
    t_end = 1.6 * t_star(v0, mu) + 0.005
    return int(math.ceil(t_end / DT))


# ── DIRT config template ──────────────────────────────────────────────────────
TOML_TEMPLATE = """[comm]
processors_x = 1
processors_y = 1
processors_z = 1
[domain]
x_low = -2.0
x_high = 2.0
y_low = -0.05
y_high = 0.05
z_low = 0.0
z_high = 0.1
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"
[neighbor]
skin_fraction = 1.1
bin_size = 0.015
every = 1
[gravity]
gx = 0.0
gy = 0.0
gz = -{g}
[dem]
contact_model = "hertz"
[[dem.materials]]
name = "glass"
youngs_mod = {youngs:.6e}
poisson_ratio = {nu}
restitution = {e_n}
friction = {mu}
[[particles.insert]]
material = "glass"
count = 1
radius = {radius}
density = {density}
velocity_x = {v0}
region = {{ type = "block", min = [-1.0e-6, -1.0e-6, {zproj_lo:.7f}], max = [1.0e-6, 1.0e-6, {zproj_hi:.7f}] }}
[[wall]]
point_x = 0.0
point_y = 0.0
point_z = 0.0
normal_x = 0.0
normal_y = 0.0
normal_z = 1.0
material = "glass"
[output]
dir = "{outdir}"
[run]
steps = {steps}
thermo = {thermo}
dt = {dt:.6e}
"""


# ── LAMMPS template ───────────────────────────────────────────────────────────
# Mirrors the DIRT setup exactly: one sphere of the same radius/density launched
# horizontally at v0 with NO spin (`velocity ... set v0 0 0`) onto a flat
# frictional granular floor (`fix wall/gran ... zplane 0.0 NULL`) under gravity.
# Same material: hertz/material <E> <e> <nu> normal + mindlin tangential with the
# Coulomb cap mu. While the contact slides the wall decelerates the center at
# a = mu*g and spins it up; once v = R*omega it rolls at v_final = (5/7) v0.
# Per-step dump of t, vx, omega_y (single atom -> reduce sum is just that atom).
LMP_TEMPLATE = """units           si
atom_style      sphere
boundary        f f f
comm_modify     vel yes
region          box block -2.0 2.0 -0.05 0.05 0.0 0.1 units box
create_box      1 box
create_atoms    1 single 0.0 0.0 {ztop:.9f} units box
group           mover id 1
set             group mover diameter {dmover:.9f} density {density}
pair_style      granular
pair_coeff      1 1 hertz/material {youngs:.6e} {e_n} {nu} &
                tangential mindlin NULL 1.0 {mu} &
                rolling none twisting none
# Flat frictional granular floor at z = 0, same contact law as the pair.
fix             floor all wall/gran granular &
                hertz/material {youngs:.6e} {e_n} {nu} &
                tangential mindlin NULL 1.0 {mu} &
                rolling none twisting none &
                zplane 0.0 NULL
# Launch: translational v0 along +x, ZERO spin (matches DIRT's omega0 = 0).
velocity        mover set {v0} 0.0 0.0 units box
set             group mover omega 0.0 0.0 0.0
fix             grav all gravity {g} vector 0 0 -1
fix             integ mover nve/sphere
timestep        {dt:.6e}
thermo          {steps}
compute         vxa mover property/atom vx
compute         wya mover property/atom omegay
compute         vxs mover reduce sum c_vxa
compute         wys mover reduce sum c_wya
fix             rec mover ave/time 1 1 {every} c_vxs c_wys file {out} mode scalar
run             {steps}
"""


def find_lammps():
    for b in LAMMPS_BINS:
        path = shutil.which(b)
        if path:
            return path
    return None


def case_tag(mu, v0):
    return f"mu_{mu:g}_v0_{v0:g}"


def case_dir(mu, v0):
    return os.path.join(SWEEP_DIR, case_tag(mu, v0))


def sweep_cases():
    """All (mu, v0) cases: a mu sweep at v0=V0_REF + a v0 sweep at mu=MU_REF."""
    cases = [(mu, V0_REF) for mu in MU_LIST]
    for v0 in V0_LIST:
        if (MU_REF, v0) not in cases:
            cases.append((MU_REF, v0))
    return cases


def _dirt_config(mu, v0, outdir):
    steps = n_steps(v0, mu)
    return TOML_TEMPLATE.format(
        g=G, youngs=YOUNGS_MOD, nu=NU, e_n=RESTITUTION, mu=mu,
        density=DENSITY,
        radius=R, v0=v0, zproj_lo=R + 1e-6, zproj_hi=R + 1.2e-6,
        outdir=outdir, steps=steps, thermo=max(1, steps // 10), dt=DT,
    )


# ── generate ─────────────────────────────────────────────────────────────────
def generate():
    n = 0
    for mu, v0 in sweep_cases():
        cdir = case_dir(mu, v0)
        os.makedirs(cdir, exist_ok=True)
        with open(os.path.join(cdir, "config.toml"), "w") as f:
            f.write(_dirt_config(mu, v0, cdir))
        n += 1
    print(f"Generated {n} DIRT sweep configs under {SWEEP_DIR}")


# ── start ────────────────────────────────────────────────────────────────────
def _run_dirt(cdir):
    """Run one prepared DIRT case; return the path to its results CSV or None."""
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "sliding_friction_results.csv")
    if os.path.exists(res):
        os.remove(res)  # never re-plot a stale run
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

    lammps = find_lammps()
    print(f"LAMMPS: {lammps}" if lammps else
          "LAMMPS: not found on PATH — running DIRT only.")

    # Wipe stale LAMMPS summary so old results can never be re-plotted.
    if os.path.exists(LAMMPS_CSV):
        os.remove(LAMMPS_CSV)

    cases = sweep_cases()
    rows, lmp_rows = [], []
    for i, (mu, v0) in enumerate(cases, 1):
        cdir = case_dir(mu, v0)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  [{i:2d}/{len(cases)}] missing config for {case_tag(mu, v0)} — run 'generate' first.")
            continue
        print(f"  [{i:2d}/{len(cases)}] {case_tag(mu, v0):16s}", end="  ", flush=True)
        res = _run_dirt(cdir)
        if res is None:
            print("DIRT FAILED")
            continue
        fit = _fit_case(res, mu, v0)
        if fit is None:
            print("FIT FAILED (no clean sliding phase)")
            continue
        rows.append(fit)
        print(f"a_fit={fit['a_fit']:.3f}  a_th={fit['a_theory']:.3f}  "
              f"v_final={fit['v_final']:.4f} (th {fit['v_final_theory']:.4f})", end="")
        if lammps:
            lf = _run_lammps(lammps, mu, v0)
            if lf is not None:
                lmp_rows.append(lf)
                print(f"   | LAMMPS a_fit={lf['a_fit']:.3f}  "
                      f"v_final={lf['v_final']:.4f}", end="")
            else:
                print("   | LAMMPS FAILED", end="")
        print()

    if not rows:
        print("\nERROR: no DIRT results collected.")
        sys.exit(1)
    _write_summary(rows)
    print(f"\nDIRT:   {len(rows)}/{len(cases)} cases -> {SUMMARY_CSV}")
    if lmp_rows:
        _write_summary(lmp_rows, LAMMPS_CSV)
        print(f"LAMMPS: {len(lmp_rows)}/{len(cases)} cases -> {LAMMPS_CSV}")


# ── analysis ──────────────────────────────────────────────────────────────────
def _load_series(path):
    """Load a per-case time series: t, vx, omega_y, radius, in_contact."""
    t, vx, w, contact = [], [], [], []
    R_proj = R
    with open(path) as f:
        for r in csv.DictReader(f):
            t.append(float(r["t"]))
            vx.append(float(r["vx"]))
            w.append(float(r["omega_y"]))
            contact.append(int(float(r["in_contact"])))
            R_proj = float(r["radius"])
    return t, vx, w, contact, R_proj


def _linfit(xs, ys):
    """Ordinary least-squares slope/intercept of ys vs xs."""
    n = len(xs)
    sx = sum(xs); sy = sum(ys)
    sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        return 0.0, 0.0
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    return slope, intercept


def _fit_series(t, vx, w, contact, R_proj, mu, v0):
    """Fit the sliding-phase deceleration and measure the rolling plateau from a
    loaded time series (t, vx, omega_y, in_contact). Code-agnostic core used by
    both the DIRT and the LAMMPS legs.

    The sliding phase is the interval (in contact) where the surface speed
    R*omega is still below the center speed |vx| (contact point slides). We fit
    vx(t) over the central 70% of that interval (trimming settling transients and
    the kink) to get the constant deceleration, and read v_final from the
    rolling plateau after the slip->roll transition.
    """
    if len(t) < 50:
        return None

    # Sliding while in contact and surface speed has not yet caught the center.
    # vx > 0 throughout (decelerating but never reversing). Slip s = vx - R*omega.
    slip = [vx[i] - R_proj * w[i] for i in range(len(t))]
    # Transition index: first in-contact sample where slip has essentially closed.
    slip_tol = 0.02 * v0
    i_trans = None
    for i in range(len(t)):
        if contact[i] and slip[i] <= slip_tol:
            i_trans = i
            break
    if i_trans is None or i_trans < 20:
        return None

    # Sliding window = in-contact samples before the transition. Trim the first
    # 15% (vertical settling) and last 5% (approach to the kink).
    lo = int(0.15 * i_trans)
    hi = int(0.95 * i_trans)
    if hi - lo < 10:
        return None
    ts = t[lo:hi]
    vs = vx[lo:hi]
    slope, intercept = _linfit(ts, vs)
    a_fit = -slope  # deceleration magnitude

    # Spin-up rate alpha over the same window (omega about -y: contact below
    # center, +x motion -> omega_y negative -> use magnitude).
    ws = [abs(w[i]) for i in range(lo, hi)]
    alpha_fit, _ = _linfit(ts, ws)

    # Rolling plateau: mean vx over the last in-contact samples after transition.
    tail = [vx[i] for i in range(len(t)) if contact[i] and i >= i_trans]
    if len(tail) < 5:
        return None
    # Use the last third of the tail to avoid the immediate post-kink ring-down.
    tail = tail[len(tail) // 3:]
    v_final = sum(tail) / len(tail)

    return {
        "mu": mu, "v0": v0,
        "a_fit": a_fit, "a_theory": mu * G,
        "alpha_fit": alpha_fit, "alpha_theory": mu * G / ((2.0 / 5.0) * R_proj),
        "v_final": v_final, "v_final_theory": (5.0 / 7.0) * v0,
        "t_star_meas": t[i_trans], "t_star_theory": t_star(v0, mu),
        "series": "",
    }


def _fit_case(path, mu, v0):
    """Load a DIRT results CSV and fit it (deceleration + rolling plateau)."""
    t, vx, w, contact, R_proj = _load_series(path)
    fit = _fit_series(t, vx, w, contact, R_proj, mu, v0)
    if fit is not None:
        fit["series"] = path
    return fit


def _write_summary(rows, path=SUMMARY_CSV):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fields = ["mu", "v0", "a_fit", "a_theory", "alpha_fit", "alpha_theory",
              "v_final", "v_final_theory", "t_star_meas", "t_star_theory", "series"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


# ── LAMMPS leg ────────────────────────────────────────────────────────────────
def _lammps_input(mu, v0, cdir):
    """Write the LAMMPS input for one (mu, v0) case; return (in_path, out_path)."""
    steps = n_steps(v0, mu)
    out = os.path.join(cdir, "lammps_series.txt")
    in_path = os.path.join(cdir, "in.lammps")
    # Seat the sphere a hair above z = R so normal contact engages immediately
    # without a penetration kick — same trick as the DIRT insert band.
    with open(in_path, "w") as f:
        f.write(LMP_TEMPLATE.format(
            ztop=R + 1.1e-6, dmover=2.0 * R, density=DENSITY,
            youngs=YOUNGS_MOD, e_n=RESTITUTION, nu=NU, mu=mu,
            v0=v0, g=G, dt=DT, steps=steps,
            every=max(1, steps // 2000), out=out,
        ))
    return in_path, out


def _load_lammps_series(path):
    """Read t, vx, |omega_y| from a LAMMPS ave/time scalar dump.

    Columns: timestep  c_vxs  c_wys.  The single atom is always in contact after
    a sub-microsecond settle, so we mark every sample in_contact; the fit trims
    the first 15% anyway. omega_y is taken in magnitude so the slip convention
    (slip = vx - R*omega) matches DIRT regardless of LAMMPS's sign.
    """
    t, vx, w, contact = [], [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            if len(p) >= 3:
                t.append(float(p[0]) * DT)
                vx.append(float(p[1]))
                w.append(abs(float(p[2])))
                contact.append(1)
    return t, vx, w, contact, R


def _run_lammps(lammps, mu, v0):
    """Run one LAMMPS case and fit it the SAME way as DIRT. Returns a fit dict."""
    cdir = case_dir(mu, v0)
    in_path, out = _lammps_input(mu, v0, cdir)
    if os.path.exists(out):
        os.remove(out)
    log = os.path.join(cdir, "lammps.log")
    with open(log, "w") as lf:
        proc = subprocess.run([lammps, "-in", in_path], cwd=cdir,
                              stdout=lf, stderr=subprocess.STDOUT)
    if proc.returncode != 0 or not os.path.isfile(out):
        return None
    t, vx, w, contact, R_proj = _load_lammps_series(out)
    fit = _fit_series(t, vx, w, contact, R_proj, mu, v0)
    if fit is not None:
        fit["series"] = out
    return fit


# ── graph (validate + plot) ──────────────────────────────────────────────────
def _load_summary(path=SUMMARY_CSV):
    if not os.path.isfile(path):
        return []
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            d = {k: (v if k == "series" else float(v)) for k, v in r.items()}
            rows.append(d)
    return rows


def compare_codes(dirt, lammps):
    """Print DIRT vs LAMMPS fitted a and v_final per case."""
    key = lambda r: (round(r["mu"], 6), round(r["v0"], 6))
    lmp = {key(r): r for r in lammps}
    print("\n=== DIRT vs LAMMPS (fitted a, rolling plateau v_final) ===")
    print(f"  {'mu':>5}{'v0':>6}{'a_DIRT':>9}{'a_LMP':>9}{'d_a':>9}"
          f"{'vf_DIRT':>10}{'vf_LMP':>10}{'d_vf':>9}")
    for r in sorted(dirt, key=lambda x: (x["v0"], x["mu"])):
        l = lmp.get(key(r))
        if not l:
            continue
        da = r["a_fit"] - l["a_fit"]
        dvf = r["v_final"] - l["v_final"]
        print(f"  {r['mu']:>5.2f}{r['v0']:>6.2f}{r['a_fit']:>9.3f}{l['a_fit']:>9.3f}"
              f"{da:>+9.3f}{r['v_final']:>10.4f}{l['v_final']:>10.4f}{dvf:>+9.4f}")


# Validation tolerances (relative).
TOL_A = 0.08         # deceleration a = mu*g
TOL_VFINAL = 0.03    # rolling plateau v_final = (5/7) v0
TOL_TSTAR = 0.10     # slip->roll transition time


def validate(rows):
    print("\n=== Sliding-friction validation ===")
    print(f"  g={G}  R={R} m  (floor: dirt_wall z-plane at z=0)")
    hdr = (f"  {'mu':>5}{'v0':>6}{'a_fit':>9}{'a=mu g':>9}{'err%':>7}"
           f"{'v_fin':>9}{'5/7 v0':>9}{'err%':>7}{'t*meas':>9}{'t*th':>9}{'err%':>7}  note")
    print(hdr)
    ok = True
    for r in sorted(rows, key=lambda x: (x["v0"], x["mu"])):
        ea = abs(r["a_fit"] - r["a_theory"]) / r["a_theory"]
        ev = abs(r["v_final"] - r["v_final_theory"]) / r["v_final_theory"]
        et = abs(r["t_star_meas"] - r["t_star_theory"]) / r["t_star_theory"]
        notes = []
        if ea > TOL_A:
            notes.append("A"); ok = False
        if ev > TOL_VFINAL:
            notes.append("VFINAL"); ok = False
        if et > TOL_TSTAR:
            notes.append("TSTAR"); ok = False
        print(f"  {r['mu']:>5.2f}{r['v0']:>6.2f}{r['a_fit']:>9.3f}{r['a_theory']:>9.3f}"
              f"{100*ea:>7.1f}{r['v_final']:>9.4f}{r['v_final_theory']:>9.4f}{100*ev:>7.1f}"
              f"{1e3*r['t_star_meas']:>8.2f}m{1e3*r['t_star_theory']:>8.2f}m{100*et:>7.1f}"
              f"  {' '.join(notes)}")
    print(f"\n  tolerances: a {100*TOL_A:.0f}%, v_final {100*TOL_VFINAL:.0f}%, t* {100*TOL_TSTAR:.0f}%")
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


def _lmp_index(lammps):
    """Map (mu, v0) -> LAMMPS fit row for quick overlay lookup."""
    return {(round(r["mu"], 6), round(r["v0"], 6)): r for r in (lammps or [])}


def plot(rows, lammps=None):
    os.makedirs(PLOT_DIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    lmp = _lmp_index(lammps)

    # ── (1) vx(t) and R*omega(t) for the mu sweep at v0 = V0_REF ──────────────
    mu_rows = sorted([r for r in rows if abs(r["v0"] - V0_REF) < 1e-9],
                     key=lambda x: x["mu"])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    cmap = plt.get_cmap("viridis")
    for k, r in enumerate(mu_rows):
        t, vx, w, contact, R_proj = _load_series(r["series"])
        color = cmap(k / max(1, len(mu_rows) - 1))
        ax1.plot([1e3 * x for x in t], vx, "-", color=color, label=f"DIRT μ={r['mu']:g}")
        ax1.plot([1e3 * x for x in t], [R_proj * abs(wi) for wi in w], "--", color=color, lw=1.0)
        # transition marker + theory plateau
        ax1.axvline(1e3 * r["t_star_meas"], color=color, lw=0.6, alpha=0.4)
        # LAMMPS overlay: vx(t) as open circles (sparse), Rω(t) thin dotted.
        lr = lmp.get((round(r["mu"], 6), round(r["v0"], 6)))
        if lr and lr.get("series") and os.path.isfile(lr["series"]):
            lt, lvx, lw, _, lR = _load_lammps_series(lr["series"])
            step = max(1, len(lt) // 30)
            ax1.plot([1e3 * x for x in lt[::step]], lvx[::step], "o", color=color,
                     mfc="none", ms=4, lw=0)
            ax1.plot([1e3 * x for x in lt], [lR * wi for wi in lw], ":", color=color, lw=0.9)
    ax1.axhline((5.0 / 7.0) * V0_REF, color="k", ls=":", lw=1.0, label="(5/7) v₀")
    if lmp:
        ax1.plot([], [], "ko", mfc="none", ms=4, label="LAMMPS vₓ")
    ax1.set_xlabel("time (ms)")
    ax1.set_ylabel("speed (m/s)")
    ax1.set_title(f"v_x (solid) and Rω (dashed) — v₀={V0_REF} m/s\n"
                  f"DIRT lines, LAMMPS open markers; kink at t*, plateau at (5/7)v₀")
    ax1.legend(fontsize=9)

    # zoom on a single mu to show the slip phase, kink, and rolling plateau clearly
    ref = next((r for r in mu_rows if abs(r["mu"] - MU_REF) < 1e-9), mu_rows[0])
    t, vx, w, contact, R_proj = _load_series(ref["series"])
    ax2.plot([1e3 * x for x in t], vx, "-", color="tab:blue", label="v_x (center)")
    ax2.plot([1e3 * x for x in t], [R_proj * abs(wi) for wi in w], "-", color="tab:orange", label="Rω (surface)")
    # theory sliding lines
    a_th = ref["a_theory"]
    ts_th = [x for x in t if x <= ref["t_star_theory"]]
    ax2.plot([1e3 * x for x in ts_th], [V0_REF - a_th * x for x in ts_th],
             "k--", lw=1.0, label="v₀−μg·t (theory)")
    alpha_th = ref["alpha_theory"]
    ax2.plot([1e3 * x for x in ts_th], [R_proj * alpha_th * x for x in ts_th],
             "k:", lw=1.0, label="R·αt (theory)")
    ax2.axvline(1e3 * ref["t_star_theory"], color="gray", lw=0.8)
    ax2.axhline((5.0 / 7.0) * V0_REF, color="green", ls=":", lw=1.0, label="(5/7) v₀")
    ax2.set_xlabel("time (ms)")
    ax2.set_ylabel("speed (m/s)")
    ax2.set_title(f"Slip → roll transition (μ={ref['mu']:g}, v₀={V0_REF} m/s)")
    ax2.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "slip_to_roll.png"))
    plt.close(fig)

    # ── (2) a_fit vs mu (should lie on a = g·mu) ──────────────────────────────
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ms = sorted(mu_rows, key=lambda x: x["mu"])
    mus = [r["mu"] for r in ms]
    ax.plot(mus, [r["a_fit"] for r in ms], "o", ms=8, color="tab:blue", label="DIRT (fit)")
    line = [G * m for m in mus]
    ax.plot(mus, line, "k--", label="a = g·μ (theory)")
    lmu = sorted([r for r in (lammps or []) if abs(r["v0"] - V0_REF) < 1e-9],
                 key=lambda x: x["mu"])
    if lmu:
        ax.plot([r["mu"] for r in lmu], [r["a_fit"] for r in lmu], "s", ms=8,
                mfc="none", color="tab:red", label="LAMMPS (fit)")
    ax.set_xlabel("friction coefficient μ")
    ax.set_ylabel("sliding deceleration a (m/s²)")
    ax.set_title("Sliding-phase deceleration vs Coulomb prediction")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "decel_vs_mu.png"))
    plt.close(fig)

    # ── (3) v_final vs v0 (should lie on (5/7) v0, independent of mu) ─────────
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    v0_rows = sorted([r for r in rows if abs(r["mu"] - MU_REF) < 1e-9],
                     key=lambda x: x["v0"])
    v0s = [r["v0"] for r in v0_rows]
    ax.plot(v0s, [r["v_final"] for r in v0_rows], "s", ms=8, color="tab:blue",
            label=f"DIRT (μ={MU_REF:g})")
    ax.plot(v0s, [(5.0 / 7.0) * v for v in v0s], "k--", label="(5/7) v₀ (theory)")
    lv0 = sorted([r for r in (lammps or []) if abs(r["mu"] - MU_REF) < 1e-9],
                 key=lambda x: x["v0"])
    if lv0:
        ax.plot([r["v0"] for r in lv0], [r["v_final"] for r in lv0], "D", ms=7,
                mfc="none", color="tab:red", label=f"LAMMPS (μ={MU_REF:g})")
    ax.set_xlabel("launch speed v₀ (m/s)")
    ax.set_ylabel("final rolling speed v_final (m/s)")
    ax.set_title("Rolling-plateau speed vs (5/7) v₀")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "vfinal_vs_v0.png"))
    plt.close(fig)

    print(f"\nFigures -> {PLOT_DIR}/slip_to_roll.png, decel_vs_mu.png, vfinal_vs_v0.png")


def graph():
    rows = _load_summary()
    if not rows:
        print(f"No {SUMMARY_CSV} — run 'start' first.")
        return False
    lammps = _load_summary(LAMMPS_CSV)
    ok = validate(rows)
    if lammps:
        compare_codes(rows, lammps)
    else:
        print("\n(no LAMMPS sweep — plotting DIRT only)")
    plot(rows, lammps)
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
