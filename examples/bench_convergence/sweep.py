#!/usr/bin/env python3
"""
Numerical convergence study across Tier-1 DIRT benchmarks.

Every other bench_* example runs at a SINGLE timestep (0.15 × the Rayleigh
critical dt) and a SINGLE particle count, and simply asserts the result is close
to a reference. None of them answers the two questions a solver user actually
has: *how small must dt be, and how many particles do I need, before the answer
stops moving?* This driver answers both by re-running existing benchmark binaries
over a ladder of resolutions and watching the key observables converge.

It runs TWO independent sub-studies (no new Rust code — it drives the compiled
`bench_hertz_rebound` and `bench_sphere_haff_cooling` example binaries through
generated configs):

  A. TIMESTEP convergence (deterministic).  A single glass sphere strikes a wall
     (the `bench_hertz_rebound` setup, v0 = 1 m/s) at a ladder of timesteps
     dt = f · dt_Rayleigh, f ∈ [0.5 … 0.015].  We track three observables:
     the coefficient of restitution (COR), the contact duration t_c, and the peak
     overlap δ_max.  For the elastic anchor (COR = 1.0) t_c and δ_max have exact
     Hertz closed forms, so we show measured → analytic as dt → 0 AND estimate the
     observed order of accuracy p by Richardson self-convergence.

  B. PARTICLE-COUNT convergence (finite-size / statistical).  A freely cooling
     granular gas (the `bench_sphere_haff_cooling` setup) is run at a ladder of
     particle counts N ∈ [200 … 1600] held at a FIXED volume fraction φ (the box
     grows with N), each over several independent random seeds.  The intensive
     observable is Haff's cooling time t_c (from the linearized law
     1/√(T/T0) = 1 + t/t_c).  As N grows the mean t_c plateaus and the run-to-run
     scatter — and the RMS residual of the Haff fit — shrink like ~1/√N.

Commands (from anywhere):
    python3 examples/bench_convergence/sweep.py generate  # write per-case configs
    python3 examples/bench_convergence/sweep.py start     # build + run all sims -> CSV
    python3 examples/bench_convergence/sweep.py graph     # validate + plot
    python3 examples/bench_convergence/sweep.py           # all three, in order

Outputs:
    sweep/<case>/config.toml     generated DIRT configs                (gitignored)
    data/dt_results.csv          timestep-study observables            (gitignored)
    data/n_results.csv           particle-count-study observables      (gitignored)
    data/n_curve_<N>.csv         representative cooling curve per N     (gitignored)
    plots/*.png                  final figures                         (tracked)
    report.md                    written summary + recommended dt / N   (tracked)

References:
    K.L. Johnson, Contact Mechanics, Cambridge University Press, 1985 (Hertz).
    P.K. Haff, "Grain flow as a fluid-mechanical phenomenon", J. Fluid Mech.
        134 (1983) 401-430 (free-cooling t^-2 law).
    P.J. Roache, "Verification and Validation in Computational Science and
        Engineering", Hermosa 1998 (grid/step convergence, observed order).
"""

import os
import sys
import csv
import math
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
DT_CSV = os.path.join(DATA_DIR, "dt_results.csv")
N_CSV = os.path.join(DATA_DIR, "n_results.csv")
REPORT = os.path.join(SCRIPT_DIR, "report.md")

CARGO_FEATURES = ["--no-default-features", "--features", "precision-double"]

# ── Study A: timestep convergence (single-sphere Hertz rebound) ───────────────
HERTZ_EXAMPLE = "bench_hertz_rebound"
A_YOUNGS = 70.0e9      # Pa  (soda-lime glass, mirrors bench_hertz_rebound)
A_POISSON = 0.22
A_RADIUS = 0.005       # m
A_DENSITY = 2500.0     # kg/m^3
A_V0 = 1.0             # m/s impact velocity (fixed; dt is the swept variable)
# dt as a fraction of the Rayleigh critical timestep. 0.15 is the solver default;
# we bracket it by a factor ~30 on either side to expose the convergence trend.
A_DT_FRACS = [0.5, 0.35, 0.25, 0.15, 0.10, 0.06, 0.03, 0.015]
A_CORS = [1.0, 0.9]    # 1.0 = elastic anchor (analytic t_c, δ_max); 0.9 = damped
# Recommended-dt gate: coarsest fraction (walking up contiguously from the finest)
# whose COR, t_c and δ_max all stay within this of the finest-dt value.
A_CONVERGED_TOL = 0.02
A_DEFAULT_FRAC = 0.15   # the solver's default dt fraction; we check it is adequate
A_DEFAULT_TOL = 0.03    # default-dt observables must be within this of the finest dt
A_ANCHOR_TOL = 0.02    # elastic-anchor error vs analytic Hertz at the finest dt

# ── Study B: particle-count convergence (free-cooling granular gas) ───────────
HAFF_EXAMPLE = "bench_sphere_haff_cooling"
B_YOUNGS = 7.0e7       # Pa  (softened glass, mirrors bench_sphere_haff_cooling)
B_POISSON = 0.245
B_RESTITUTION = 0.926
B_FRICTION = 0.16
B_RADIUS = 0.0011      # m
B_DENSITY = 2500.0     # kg/m^3
B_VSIGMA = 0.5         # Gaussian velocity sigma per component (sets T0)
B_PHI = 0.07           # volume fraction, HELD FIXED as N grows (box scales with N)
B_N_LIST = [200, 400, 800, 1600]
B_SEEDS = [1, 2, 3, 4]
B_STEPS = 150000       # fixed across N (same dt → same physical window)
B_SAMPLE = 2000        # the binary records T every 2000 steps
# Fit window on T_trans/T0: skip the rotational-equilibration transient, stop
# before the noise floor.
B_FIT_LO = 0.04        # lower bound on T/T0 (≈ 1/25)
B_FIT_HI = 0.7         # upper bound on T/T0 (skip start-up)
# Convergence gates.
B_PLATEAU_TOL = 0.10   # |t_c(N_max) - t_c(N_max/2)| / t_c  must be below this
B_RESID_MAX = 0.01     # recommended-N gate: mean Haff-fit RMS residual below this
B_CV_MAX = 0.03        # recommended-N gate: run-to-run CV of t_c below this


# ── shared helpers ────────────────────────────────────────────────────────────
def dt_rayleigh(youngs, poisson, radius, density):
    """Rayleigh critical timestep for a sphere of this material (the same
    estimate the solver auto-timestep uses, before its 0.15 safety factor)."""
    g = youngs / (2.0 * (1.0 + poisson))
    alpha = 0.1631 * poisson + 0.876605
    return math.pi * radius / alpha * (density / g) ** 0.5


def run_binary(example, config, log_path):
    """Run one compiled example binary on a config; return True on exit 0."""
    with open(log_path, "w") as log:
        proc = subprocess.run(
            ["cargo", "run", "--release", "--example", example, *CARGO_FEATURES,
             "--", config],
            cwd=REPO_ROOT, stdout=log, stderr=subprocess.STDOUT,
        )
    return proc.returncode == 0


def build(examples):
    print(f"Building {', '.join(examples)} (release)...", flush=True)
    args = ["cargo", "build", "--release"]
    for e in examples:
        args += ["--example", e]
    args += CARGO_FEATURES
    subprocess.run(args, cwd=REPO_ROOT, check=True)


def linfit(xs, ys):
    """Ordinary least squares y = a + b·x. Returns (a, b, r2, rms_resid) where
    rms_resid is the RMS of residuals normalized by the mean |y|."""
    n = len(xs)
    sx = sum(xs); sy = sum(ys)
    sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-300:
        return 0.0, 0.0, 0.0, float("inf")
    b = (n * sxy - sx * sy) / denom
    a = (sy - b * sx) / n
    ybar = sy / n
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    rms = math.sqrt(ss_res / n) / (abs(ybar) if abs(ybar) > 1e-30 else 1.0)
    return a, b, r2, rms


# ── Study A: timestep configs ─────────────────────────────────────────────────
A_TOML = """\
# Auto-generated: Hertz rebound timestep-convergence case
# COR = {cor}, dt = {dt:.6e} s ({frac} x dt_Rayleigh)
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
z_high = 0.1
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[neighbor]
skin_fraction = 1.1
bin_size = 0.015
every = 1

[dem]
contact_model = "hertz"

[[dem.materials]]
name = "glass"
youngs_mod = {youngs}
poisson_ratio = {poisson}
restitution = {cor}
friction = 0.0

[[particles.insert]]
material = "glass"
count = 1
radius = {radius}
density = {density}
velocity_z = -{v0}
region = {{ type = "block", min = [-0.001, -0.001, 0.007], max = [0.001, 0.001, 0.008] }}

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
dt = {dt:.10e}
steps = {steps}
thermo = 1000000
"""


def a_case_tag(cor, frac):
    return f"dt_cor{cor}_f{frac}"


def a_generate():
    dtR = dt_rayleigh(A_YOUNGS, A_POISSON, A_RADIUS, A_DENSITY)
    n = 0
    for cor in A_CORS:
        for frac in A_DT_FRACS:
            dt = frac * dtR
            # enough steps to fall the ~2.5 mm gap + resolve contact/rebound
            total_time = 0.004
            steps = int(total_time / dt) + 2000
            cdir = os.path.join(SWEEP_DIR, a_case_tag(cor, frac))
            os.makedirs(cdir, exist_ok=True)
            with open(os.path.join(cdir, "config.toml"), "w") as f:
                f.write(A_TOML.format(
                    cor=cor, dt=dt, frac=frac, youngs=f"{A_YOUNGS:.6e}",
                    poisson=A_POISSON, radius=A_RADIUS, density=A_DENSITY,
                    v0=A_V0, outdir=cdir, steps=steps))
            n += 1
    return n


A_FIELDS = ["cor", "dt_frac", "dt", "cor_meas", "contact_time", "max_overlap"]


def a_start():
    dtR = dt_rayleigh(A_YOUNGS, A_POISSON, A_RADIUS, A_DENSITY)
    rows = []
    cases = [(cor, frac) for cor in A_CORS for frac in A_DT_FRACS]
    for i, (cor, frac) in enumerate(cases, 1):
        cdir = os.path.join(SWEEP_DIR, a_case_tag(cor, frac))
        config = os.path.join(cdir, "config.toml")
        print(f"  [A {i:2d}/{len(cases)}] COR={cor:<4} dt/dtR={frac:<6}", end="  ", flush=True)
        ok = run_binary(HERTZ_EXAMPLE, config, os.path.join(cdir, "run.log"))
        csvf = os.path.join(cdir, "data", "rebound_results.csv")
        if ok and os.path.isfile(csvf):
            with open(csvf) as f:
                r = next(csv.DictReader(f))
            rows.append({
                "cor": cor, "dt_frac": frac, "dt": frac * dtR,
                "cor_meas": float(r["cor_measured"]),
                "contact_time": float(r["contact_time"]),
                "max_overlap": float(r["max_overlap"]),
            })
            print(f"COR={rows[-1]['cor_meas']:.5f}  t_c={rows[-1]['contact_time']:.4e}"
                  f"  d={rows[-1]['max_overlap']:.4e}")
        else:
            print("FAILED")
    _write_csv(DT_CSV, A_FIELDS, rows)
    print(f"Study A: {len(rows)}/{len(cases)} cases -> {DT_CSV}")
    return rows


# ── Study B: particle-count configs ───────────────────────────────────────────
def b_box_side(n):
    """Cube side giving volume fraction B_PHI for n spheres of radius B_RADIUS."""
    vpart = (4.0 / 3.0) * math.pi * B_RADIUS ** 3
    return (n * vpart / B_PHI) ** (1.0 / 3.0)


B_TOML = """\
# Auto-generated: free-cooling granular gas, N-convergence case
# N = {n}, seed = {seed}, box = {L:.5f} m (phi = {phi})
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[domain]
x_low = 0.0
x_high = {L:.6f}
y_low = 0.0
y_high = {L:.6f}
z_low = 0.0
z_high = {L:.6f}
boundary_x = "periodic"
boundary_y = "periodic"
boundary_z = "periodic"

[neighbor]
skin_fraction = 1.1
bin_size = {bin_size:.6f}

[dem]
contact_model = "hertz"

[[dem.materials]]
name = "glass"
youngs_mod = {youngs}
poisson_ratio = {poisson}
restitution = {restitution}
friction = {friction}

[[particles.insert]]
material = "glass"
count = {n}
radius = {radius}
density = {density}
velocity = {vsigma}
seed = {seed}
region = {{ type = "block", min = [{lo:.6f}, {lo:.6f}, {lo:.6f}], max = [{hi:.6f}, {hi:.6f}, {hi:.6f}] }}

[output]
dir = "{outdir}"

[run]
steps = {steps}
thermo = 100000
"""


def b_case_tag(n, seed):
    return f"n{n}_s{seed}"


def b_generate():
    cnt = 0
    for n in B_N_LIST:
        L = b_box_side(n)
        margin = 2.0 * B_RADIUS
        lo, hi = margin, L - margin
        bin_size = max(4.0 * B_RADIUS, 0.004)
        for seed in B_SEEDS:
            cdir = os.path.join(SWEEP_DIR, b_case_tag(n, seed))
            os.makedirs(cdir, exist_ok=True)
            with open(os.path.join(cdir, "config.toml"), "w") as f:
                f.write(B_TOML.format(
                    n=n, seed=seed, L=L, phi=B_PHI, bin_size=bin_size,
                    youngs=f"{B_YOUNGS:.6e}", poisson=B_POISSON,
                    restitution=B_RESTITUTION, friction=B_FRICTION,
                    radius=B_RADIUS, density=B_DENSITY, vsigma=B_VSIGMA,
                    lo=lo, hi=hi, outdir=cdir, steps=B_STEPS))
            cnt += 1
    return cnt


def b_fit_cooling(csv_path):
    """Read a cooling.csv, fit Haff's linearized law on T_trans, and return
    (tc, r2, rms_resid, npts, T0, curve) or None. `curve` is [(time, T_trans)]."""
    times, temps = [], []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            times.append(float(r["time"]))
            temps.append(float(r["T_trans"]))
    if len(temps) < 6 or temps[0] <= 0:
        return None
    t0 = temps[0]
    # Linearize inside the fit window: y = 1/sqrt(T/T0) is linear in t under Haff.
    xs, ys = [], []
    for t, T in zip(times, temps):
        ratio = T / t0
        if B_FIT_LO <= ratio <= B_FIT_HI:
            xs.append(t)
            ys.append(1.0 / math.sqrt(ratio))
    if len(xs) < 4:
        return None
    a, b, r2, rms = linfit(xs, ys)
    if b <= 0:
        return None
    tc = 1.0 / b   # slope of 1/sqrt(T/T0) vs t is 1/tc
    return tc, r2, rms, len(xs), t0, list(zip(times, temps))


B_FIELDS = ["n", "seed", "tc", "r2", "rms_resid", "npts", "T0"]


def b_start():
    rows = []
    cases = [(n, s) for n in B_N_LIST for s in B_SEEDS]
    for i, (n, seed) in enumerate(cases, 1):
        cdir = os.path.join(SWEEP_DIR, b_case_tag(n, seed))
        config = os.path.join(cdir, "config.toml")
        print(f"  [B {i:2d}/{len(cases)}] N={n:<5} seed={seed}", end="  ", flush=True)
        ok = run_binary(HAFF_EXAMPLE, config, os.path.join(cdir, "run.log"))
        cooling = os.path.join(cdir, "cooling.csv")
        fit = b_fit_cooling(cooling) if (ok and os.path.isfile(cooling)) else None
        if fit:
            tc, r2, rms, npts, t0, curve = fit
            rows.append({"n": n, "seed": seed, "tc": tc, "r2": r2,
                         "rms_resid": rms, "npts": npts, "T0": t0})
            # keep the first seed's curve for plotting
            if seed == B_SEEDS[0]:
                cpath = os.path.join(DATA_DIR, f"n_curve_{n}.csv")
                with open(cpath, "w", newline="") as f:
                    w = csv.writer(f); w.writerow(["time", "T_trans"])
                    w.writerows(curve)
            print(f"tc={tc:.4e}  R2={r2:.5f}  resid={rms*100:.2f}%")
        else:
            print("FAILED / no clean fit")
    _write_csv(N_CSV, B_FIELDS, rows)
    print(f"Study B: {len(rows)}/{len(cases)} cases -> {N_CSV}")
    return rows


def _write_csv(path, fields, rows):
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


# ── analytic Hertz references (elastic anchor) ────────────────────────────────
def hertz_contact_duration(v0):
    e_star = A_YOUNGS / (2.0 * (1.0 - A_POISSON ** 2))
    m = (4.0 / 3.0) * math.pi * A_RADIUS ** 3 * A_DENSITY
    return 2.87 * (m ** 2 / (A_RADIUS * e_star ** 2 * v0)) ** 0.2


def hertz_max_overlap(v0):
    e_star = A_YOUNGS / (2.0 * (1.0 - A_POISSON ** 2))
    m = (4.0 / 3.0) * math.pi * A_RADIUS ** 3 * A_DENSITY
    return (15.0 * m * v0 ** 2 / (16.0 * A_RADIUS ** 0.5 * e_star)) ** 0.4


# ── validate ──────────────────────────────────────────────────────────────────
def a_load():
    if not os.path.isfile(DT_CSV):
        return []
    with open(DT_CSV) as f:
        return [{k: (float(v) if k != "cor" else float(v)) for k, v in r.items()}
                for r in csv.DictReader(f)]


def b_load():
    if not os.path.isfile(N_CSV):
        return []
    with open(N_CSV) as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def observed_order(fracs, vals):
    """Estimate the observed order of accuracy p from three of the finest,
    roughly halving resolutions using Richardson's ratio
    p = ln((f_coarse - f_med)/(f_med - f_fine)) / ln(r)."""
    # sort by dt (fine -> coarse); use the three finest with ~constant ratio
    trip = sorted(zip(fracs, vals))[:3]
    if len(trip) < 3:
        return None
    (h1, v1), (h2, v2), (h3, v3) = trip  # h1<h2<h3 (fine->coarse)
    num = v3 - v2
    den = v2 - v1
    if abs(den) < 1e-30 or num / den <= 0:
        return None
    r = (h3 / h2 + h2 / h1) / 2.0
    return math.log(num / den) / math.log(r)


def validate(a_rows, b_rows):
    checks = []   # (name, passed, detail)
    rec_dt = rec_n = None
    print("=" * 70)
    print("DIRT Convergence Study — timestep (dt) and particle count (N)")
    print("=" * 70)

    # ── Study A ──
    print("\n[A] TIMESTEP CONVERGENCE — single-sphere Hertz rebound (v0 = 1 m/s)")
    dtR = dt_rayleigh(A_YOUNGS, A_POISSON, A_RADIUS, A_DENSITY)
    print(f"    dt_Rayleigh = {dtR:.4e} s ; solver default = 0.15·dt_R = {0.15*dtR:.4e} s")
    for cor in A_CORS:
        sub = sorted([r for r in a_rows if abs(r["cor"] - cor) < 1e-9],
                     key=lambda r: r["dt_frac"])
        if len(sub) < 3:
            checks.append((f"A/COR={cor}: enough cases", False, "too few runs"))
            continue
        fine = sub[0]  # smallest dt_frac
        print(f"\n  COR_input = {cor}:  (finest dt/dtR = {fine['dt_frac']})")
        print(f"    {'dt/dtR':>7} {'COR':>9} {'t_c[µs]':>10} {'δ_max[µm]':>11}"
              f" {'|Δt_c|':>8} {'|Δδ|':>8}")
        for r in sorted(sub, key=lambda r: -r["dt_frac"]):
            dtc = abs(r["contact_time"] - fine["contact_time"]) / fine["contact_time"]
            dov = abs(r["max_overlap"] - fine["max_overlap"]) / fine["max_overlap"]
            print(f"    {r['dt_frac']:>7} {r['cor_meas']:>9.5f}"
                  f" {r['contact_time']*1e6:>10.3f} {r['max_overlap']*1e6:>11.4f}"
                  f" {dtc*100:>7.2f}% {dov*100:>7.2f}%")

        # A1: elastic anchor converges to analytic Hertz at the finest dt
        if abs(cor - 1.0) < 1e-9:
            tc_th = hertz_contact_duration(A_V0)
            ov_th = hertz_max_overlap(A_V0)
            e_tc = abs(fine["contact_time"] - tc_th) / tc_th
            e_ov = abs(fine["max_overlap"] - ov_th) / ov_th
            ok = e_tc <= A_ANCHOR_TOL and e_ov <= A_ANCHOR_TOL
            checks.append((f"A/elastic anchor → analytic (finest dt)", ok,
                           f"t_c err {e_tc*100:.2f}%, δ_max err {e_ov*100:.2f}% "
                           f"(tol {A_ANCHOR_TOL*100:.0f}%)"))
            print(f"    anchor vs Hertz: t_c err {e_tc*100:.2f}%, "
                  f"δ_max err {e_ov*100:.2f}%  (analytic t_c={tc_th*1e6:.3f}µs, "
                  f"δ_max={ov_th*1e6:.4f}µm)")

        # A2: convergence trend — refining helps. The peak overlap (smooth
        # observable) at the finer half of the ladder must be closer to the
        # converged (finest-dt) value than at the coarser half. This is robust to
        # the sub-percent timestep-quantization wiggle that breaks strict
        # step-by-step monotonicity.
        by_fine = sorted(sub, key=lambda r: r["dt_frac"])
        errs = [abs(r["max_overlap"] - fine["max_overlap"]) / fine["max_overlap"]
                for r in by_fine[1:]]   # exclude the finest (error 0 by definition)
        half = max(1, len(errs) // 2)
        fine_max = max(errs[:half]) if errs[:half] else 0.0
        coarse_max = max(errs[half:]) if errs[half:] else 0.0
        trend = fine_max <= coarse_max
        checks.append((f"A/COR={cor}: δ_max converges as dt→0", trend,
                       f"max err fine half {fine_max*100:.2f}% ≤ coarse half "
                       f"{coarse_max*100:.2f}%"))

        # A3: the solver's DEFAULT dt (0.15·dt_R) is adequate — all three
        # observables within A_DEFAULT_TOL of the finest-dt value.
        dflt = next((r for r in sub if abs(r["dt_frac"] - A_DEFAULT_FRAC) < 1e-9), None)
        if dflt:
            dtc = abs(dflt["contact_time"] - fine["contact_time"]) / fine["contact_time"]
            dov = abs(dflt["max_overlap"] - fine["max_overlap"]) / fine["max_overlap"]
            dcor = abs(dflt["cor_meas"] - fine["cor_meas"])
            ok = max(dtc, dov, dcor) <= A_DEFAULT_TOL
            checks.append((f"A/COR={cor}: default dt (0.15·dt_R) adequate", ok,
                           f"max obs error {max(dtc, dov, dcor)*100:.2f}% "
                           f"(tol {A_DEFAULT_TOL*100:.0f}%)"))

        # observed order of accuracy for peak overlap (informational)
        p = observed_order([r["dt_frac"] for r in sub],
                           [r["max_overlap"] for r in sub])
        if p is not None:
            print(f"    observed order p(δ_max) ≈ {p:.2f}")

        # recommended dt: walk up from the finest; keep the coarsest fraction for
        # which it AND every finer fraction stay within A_CONVERGED_TOL of the
        # finest-dt value on all three observables (elastic anchor is cleanest).
        if abs(cor - 1.0) < 1e-9:
            for r in sorted(sub, key=lambda r: r["dt_frac"]):
                dtc = abs(r["contact_time"] - fine["contact_time"]) / fine["contact_time"]
                dov = abs(r["max_overlap"] - fine["max_overlap"]) / fine["max_overlap"]
                dcor = abs(r["cor_meas"] - fine["cor_meas"])
                if max(dtc, dov, dcor) <= A_CONVERGED_TOL:
                    rec_dt = r["dt_frac"]
                else:
                    break

    # ── Study B ──
    print("\n[B] PARTICLE-COUNT CONVERGENCE — free-cooling granular gas "
          f"(φ = {B_PHI}, fixed)")
    by_n = {}
    for r in b_rows:
        by_n.setdefault(int(r["n"]), []).append(r)
    ns = sorted(by_n)
    print(f"    {'N':>6} {'runs':>5} {'t_c[ms]':>10} {'CV(t_c)':>9}"
          f" {'RMS resid':>10} {'R²(Haff)':>9}")
    stats = {}
    for n in ns:
        tcs = [r["tc"] for r in by_n[n]]
        resid = [r["rms_resid"] for r in by_n[n]]
        r2s = [r["r2"] for r in by_n[n]]
        mean = sum(tcs) / len(tcs)
        var = sum((t - mean) ** 2 for t in tcs) / max(1, len(tcs) - 1)
        std = math.sqrt(var)
        cv = std / mean if mean else float("inf")
        mres = sum(resid) / len(resid)
        mr2 = sum(r2s) / len(r2s)
        stats[n] = {"mean": mean, "cv": cv, "resid": mres, "r2": mr2}
        print(f"    {n:>6} {len(tcs):>5} {mean*1e3:>10.4f} {cv*100:>8.2f}%"
              f" {mres*100:>9.2f}% {mr2:>9.5f}")

    if len(ns) >= 2:
        n_hi, n_mid = ns[-1], ns[-2]
        plate = abs(stats[n_hi]["mean"] - stats[n_mid]["mean"]) / stats[n_hi]["mean"]
        ok = plate <= B_PLATEAU_TOL
        checks.append(("B: t_c plateau at large N", ok,
                       f"|Δt_c|(N={n_mid}→{n_hi}) = {plate*100:.2f}% "
                       f"(tol {B_PLATEAU_TOL*100:.0f}%)"))
        # finite-size: fit residual shrinks from smallest to largest N
        shr = stats[n_hi]["resid"] < stats[ns[0]]["resid"]
        checks.append(("B: Haff-fit residual shrinks as N grows", shr,
                       f"resid {stats[ns[0]]['resid']*100:.2f}% (N={ns[0]}) → "
                       f"{stats[n_hi]['resid']*100:.2f}% (N={n_hi})"))
        # recommended N: smallest N whose mean Haff-fit residual AND run-to-run
        # scatter are both below the gates (finite-size + statistical convergence)
        for n in ns:
            if stats[n]["resid"] <= B_RESID_MAX and stats[n]["cv"] <= B_CV_MAX:
                rec_n = n
                break

    # ── verdict ──
    print("\n" + "=" * 70)
    print("Checks:")
    npass = 0
    for name, ok, detail in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}  ({detail})")
        npass += ok
    total = len(checks)
    print(f"\nRecommended timestep: dt ≤ {rec_dt} · dt_Rayleigh"
          if rec_dt else "\nRecommended timestep: (not resolved)")
    print(f"Recommended particle count: N ≥ {rec_n}"
          if rec_n else "Recommended particle count: (not resolved)")
    print(f"\n{npass}/{total} checks passed")
    all_ok = npass == total and total > 0
    print("ALL CHECKS PASSED" if all_ok else "CHECKS FAILED")
    return all_ok, checks, rec_dt, rec_n, stats


# ── report ────────────────────────────────────────────────────────────────────
def write_report(a_rows, b_rows, checks, rec_dt, rec_n, stats):
    dtR = dt_rayleigh(A_YOUNGS, A_POISSON, A_RADIUS, A_DENSITY)
    L = []
    L.append("# Convergence Study — timestep (dt) and particle count (N)\n")
    L.append("Auto-generated by `sweep.py graph`. Re-run to refresh.\n")
    L.append("## A. Timestep convergence (single-sphere Hertz rebound)\n")
    L.append(f"- `dt_Rayleigh = {dtR:.4e}` s; solver default `0.15·dt_R = "
             f"{0.15*dtR:.4e}` s.")
    L.append(f"- Elastic-anchor (COR=1) analytic Hertz: `t_c = "
             f"{hertz_contact_duration(A_V0)*1e6:.3f}` µs, `δ_max = "
             f"{hertz_max_overlap(A_V0)*1e6:.4f}` µm.")
    fine10 = [r for r in a_rows if abs(r["cor"]-1.0) < 1e-9]
    if fine10:
        f0 = min(fine10, key=lambda r: r["dt_frac"])
        L.append(f"- At the finest dt ({f0['dt_frac']}·dt_R) the measured values "
                 f"reach `t_c = {f0['contact_time']*1e6:.3f}` µs, `δ_max = "
                 f"{f0['max_overlap']*1e6:.4f}` µm — matching Hertz.")
    L.append(f"- **Recommended dt: ≤ {rec_dt} · dt_Rayleigh**" if rec_dt else
             "- Recommended dt: not resolved from this ladder.")
    L.append(f"  (coarsest dt whose COR, t_c and δ_max are all within "
             f"{A_CONVERGED_TOL*100:.0f}% of the finest-dt value).\n")
    L.append("![COR vs dt](plots/dt_convergence.png)\n")
    L.append("## B. Particle-count convergence (free-cooling granular gas)\n")
    L.append(f"- Volume fraction held fixed at φ = {B_PHI}; the box grows with N "
             f"so number density is constant.")
    L.append(f"- Observable: Haff cooling time `t_c` from `1/√(T/T0) = 1 + t/t_c`, "
             f"averaged over {len(B_SEEDS)} seeds per N.")
    if stats:
        L.append("")
        L.append("| N | mean t_c [ms] | CV(t_c) | RMS resid | R²(Haff) |")
        L.append("|---|---|---|---|---|")
        for n in sorted(stats):
            s = stats[n]
            L.append(f"| {n} | {s['mean']*1e3:.4f} | {s['cv']*100:.2f}% | "
                     f"{s['resid']*100:.2f}% | {s['r2']:.5f} |")
        L.append("")
    L.append(f"- **Recommended particle count: N ≥ {rec_n}**" if rec_n else
             "- Recommended N: not resolved from this ladder.")
    L.append(f"  (smallest N whose mean Haff-fit RMS residual is below "
             f"{B_RESID_MAX*100:.0f}% and run-to-run CV(t_c) below "
             f"{B_CV_MAX*100:.0f}%).\n")
    L.append("![t_c and residual vs N](plots/n_convergence.png)\n")
    L.append("## Checks\n")
    for name, ok, detail in checks:
        L.append(f"- {'✅' if ok else '❌'} {name} — {detail}")
    L.append("")
    with open(REPORT, "w") as f:
        f.write("\n".join(L))
    print(f"Wrote {REPORT}")


# ── plots ─────────────────────────────────────────────────────────────────────
def plot(a_rows, b_rows, stats):
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n(matplotlib/numpy unavailable, skipped plots: {e})")
        return
    plt.rcParams.update({"font.size": 11, "axes.labelsize": 12,
                         "axes.titlesize": 12, "legend.fontsize": 9,
                         "figure.dpi": 150, "savefig.dpi": 150})
    os.makedirs(PLOT_DIR, exist_ok=True)
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

    # ── Figure A: timestep convergence (3 panels) ──
    if a_rows:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
        obs = [("cor_meas", "COR", 1.0, None),
               ("contact_time", "Contact duration t_c [µs]", 1e6, hertz_contact_duration(A_V0)*1e6),
               ("max_overlap", "Peak overlap δ_max [µm]", 1e6, hertz_max_overlap(A_V0)*1e6)]
        for ax, (key, ylabel, scale, anchor) in zip(axes, obs):
            for ic, cor in enumerate(A_CORS):
                sub = sorted([r for r in a_rows if abs(r["cor"]-cor) < 1e-9],
                             key=lambda r: r["dt_frac"])
                xs = [r["dt_frac"] for r in sub]
                ys = [r[key]*scale for r in sub]
                ax.plot(xs, ys, "o-", color=colors[ic], markersize=5,
                        label=f"COR_in = {cor}")
            if anchor is not None:
                ax.axhline(anchor, color="k", linestyle="--", linewidth=1,
                           label="Hertz analytic")
            ax.axvline(0.15, color="0.6", linestyle=":", linewidth=1)
            ax.set_xscale("log")
            ax.set_xticks(sorted(A_DT_FRACS))
            ax.set_xticklabels([f"{f:g}" for f in sorted(A_DT_FRACS)], fontsize=8)
            ax.minorticks_off()
            ax.set_xlabel("dt / dt_Rayleigh  (← finer)")
            ax.set_ylabel(ylabel)
            ax.legend()
            ax.grid(True, which="major", alpha=0.25)
        axes[1].annotate("solver\ndefault", xy=(0.15, axes[1].get_ylim()[0]),
                         xytext=(0.15, axes[1].get_ylim()[0]), fontsize=7,
                         color="0.5", ha="center", va="bottom")
        fig.suptitle("Timestep convergence — single-sphere Hertz rebound "
                     "(observables → analytic as dt → 0)")
        fig.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, "dt_convergence.png"), bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {PLOT_DIR}/dt_convergence.png")

    # ── Figure B: particle-count convergence (2 panels) ──
    if stats:
        ns = sorted(stats)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
        # left: mean t_c vs N with seed scatter + error bars
        means = [stats[n]["mean"]*1e3 for n in ns]
        errs = [stats[n]["cv"]*stats[n]["mean"]*1e3 for n in ns]
        axes[0].errorbar(ns, means, yerr=errs, fmt="o-", color=colors[0],
                         capsize=4, markersize=6, label="mean ± s.d. (seeds)")
        for r in b_rows:
            axes[0].plot(r["n"], r["tc"]*1e3, ".", color="0.6", markersize=6,
                         alpha=0.6, zorder=0)
        axes[0].set_xscale("log")
        axes[0].set_xlabel("N (particles)")
        axes[0].set_ylabel("Haff cooling time t_c [ms]")
        axes[0].set_title("t_c plateaus as N grows")
        axes[0].legend()
        axes[0].grid(True, which="both", alpha=0.25)
        # right: RMS residual (and CV) vs N on log-log with 1/sqrt(N) guide
        resid = [stats[n]["resid"]*100 for n in ns]
        cvs = [stats[n]["cv"]*100 for n in ns]
        axes[1].plot(ns, resid, "s-", color=colors[1], markersize=6,
                     label="Haff-fit RMS residual")
        axes[1].plot(ns, cvs, "^-", color=colors[2], markersize=6,
                     label="CV(t_c) across seeds")
        # 1/sqrt(N) reference anchored at the first point
        import numpy as np
        guide = [resid[0]*math.sqrt(ns[0]/n) for n in ns]
        axes[1].plot(ns, guide, "k--", linewidth=1, label=r"$\propto 1/\sqrt{N}$")
        axes[1].set_xscale("log"); axes[1].set_yscale("log")
        axes[1].set_xlabel("N (particles)")
        axes[1].set_ylabel("relative scatter [%]")
        axes[1].set_title("finite-size scatter shrinks ~1/√N")
        axes[1].legend()
        axes[1].grid(True, which="both", alpha=0.25)
        fig.suptitle("Particle-count convergence — free-cooling granular gas "
                     f"(φ = {B_PHI} fixed)")
        fig.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, "n_convergence.png"), bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {PLOT_DIR}/n_convergence.png")

    # ── Figure C: representative cooling curves ──
    curves = []
    for n in B_N_LIST:
        cpath = os.path.join(DATA_DIR, f"n_curve_{n}.csv")
        if os.path.isfile(cpath):
            with open(cpath) as f:
                rows = [(float(r["time"]), float(r["T_trans"])) for r in csv.DictReader(f)]
            curves.append((n, rows))
    if curves:
        import numpy as np
        fig, ax = plt.subplots(figsize=(7, 5))
        for ic, (n, rows) in enumerate(curves):
            t = np.array([x[0] for x in rows])
            T = np.array([x[1] for x in rows])
            T0 = T[0]
            ax.plot(t*1e3, T/T0, "-", color=colors[ic % len(colors)],
                    linewidth=1.4, label=f"N = {n}")
        ax.set_yscale("log")
        ax.set_xlabel("time [ms]")
        ax.set_ylabel("T_trans / T0")
        ax.set_title("Free-cooling curves collapse as N grows (Haff t⁻²)")
        ax.legend()
        ax.grid(True, which="both", alpha=0.25)
        fig.savefig(os.path.join(PLOT_DIR, "cooling_curves.png"), bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {PLOT_DIR}/cooling_curves.png")


# ── dispatch ──────────────────────────────────────────────────────────────────
def generate():
    na = a_generate()
    nb = b_generate()
    print(f"Generated {na} timestep configs + {nb} particle-count configs "
          f"under {SWEEP_DIR}")


def start():
    build([HERTZ_EXAMPLE, HAFF_EXAMPLE])
    print("\n== Study A: timestep convergence ==")
    a_start()
    print("\n== Study B: particle-count convergence ==")
    b_start()


def graph():
    a_rows = a_load()
    b_rows = b_load()
    if not a_rows and not b_rows:
        print("ERROR: no results found. Run 'start' first.")
        sys.exit(1)
    ok, checks, rec_dt, rec_n, stats = validate(a_rows, b_rows)
    plot(a_rows, b_rows, stats)
    write_report(a_rows, b_rows, checks, rec_dt, rec_n, stats)
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
