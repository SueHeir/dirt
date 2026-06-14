#!/usr/bin/env python3
"""
Hertz contact rebound benchmark driver.

Drops a single glass sphere onto a rigid flat wall for every (impact velocity, COR)
combination and compares COR, contact duration, and peak overlap to Hertz theory.

Commands (from anywhere):
    python3 examples/bench_hertz_rebound/sweep.py generate   # write per-case configs
    python3 examples/bench_hertz_rebound/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_hertz_rebound/sweep.py graph      # validate + plot
    python3 examples/bench_hertz_rebound/sweep.py            # all three, in order

If a LAMMPS binary (lmp_serial / lmp / lmp_mpi / lammps) is on PATH, each case is
also run in LAMMPS with the equivalent granular Hertz model (damping tsuji) and
overlaid on the plots as open markers. Because no exact restitution->damping
formula exists for a nonlinear Hertz contact, the LAMMPS restitution input is
back-solved per nominal COR so its *measured* COR matches DIRT — isolating the
shared contact physics from each code's damping calibration. LAMMPS is optional —
without it, only DIRT runs.

Outputs:
    sweep/<case>/config.toml   DIRT configs                          (gitignored)
    sweep/<case>/in.lammps     LAMMPS inputs (calibrated restitution)(gitignored)
    data/sweep_results.csv     DIRT results                          (gitignored)
    data/lammps_results.csv    LAMMPS results                        (gitignored)
    plots/*.png                final figures                         (tracked)

Hertz contact theory (sphere on rigid wall):
    Contact duration:  t_c       = 2.87 * (m^2 / (R E*^2 v0))^(1/5)
    Peak overlap:      delta_max = (15 m v0^2 / (16 R^(1/2) E*))^(2/5)
    where E* = E / (2(1 - nu^2)) for an identical-material sphere on a wall.

Reference: K.L. Johnson, Contact Mechanics, Cambridge University Press, 1985.
"""

import os
import sys
import csv
import math
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_hertz_rebound"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep_results.csv")     # DIRT results
LAMMPS_CSV = os.path.join(DATA_DIR, "lammps_results.csv")   # LAMMPS (restitution tuned to DIRT)

# LAMMPS binary candidates, in preference order. LAMMPS is optional: if none is
# found, the LAMMPS leg is skipped and only DIRT is run/plotted.
LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# Parameter sweep.
VELOCITIES = [0.1, 0.5, 1.0, 2.0]   # impact velocity [m/s]
# COR = 1.0 is the elastic anchor: zero damping, so the sim must hit the
# (undamped) Hertz theory for contact duration and peak overlap exactly.
CORS = [0.5, 0.7, 0.9, 0.95, 1.0]   # input coefficient of restitution

# Material properties.
YOUNGS_MOD = 70.0e9    # Pa
POISSON_RATIO = 0.22
RADIUS = 0.005         # m
DENSITY = 2500.0       # kg/m^3

# Validation tolerances.
# COR tolerance is COR-dependent: the linear viscoelastic damping model (beta
# mapping derived from Hooke/linear theory) deviates when used with nonlinear
# Hertz stiffness — a well-known limitation shared by LAMMPS et al.
COR_TOL_HIGH = 0.03      # 3% for COR >= 0.7
COR_TOL_LOW = 0.12       # 12% for COR < 0.7
CONTACT_TIME_TOL = 0.10  # 10% for contact duration
# Elastic Hertz theory over-predicts overlap for dissipative contacts, so the
# overlap tolerance scales with COR.
OVERLAP_TOL_HIGH = 0.10  # COR >= 0.85 (near-elastic)
OVERLAP_TOL_MED = 0.15   # 0.6 <= COR < 0.85
OVERLAP_TOL_LOW = 0.25   # COR < 0.6 (heavy damping)

# Derived sphere-wall properties (wall: infinite mass and radius).
E_STAR = YOUNGS_MOD / (2.0 * (1.0 - POISSON_RATIO**2))
MASS = (4.0 / 3.0) * math.pi * RADIUS**3 * DENSITY
M_EFF = MASS
R_EFF = RADIUS

TOML_TEMPLATE = """\
# Auto-generated config for the Hertz rebound sweep — v0 = {v0} m/s, COR = {cor}
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
youngs_mod = {youngs_mod}
poisson_ratio = {poisson_ratio}
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
dir = "{output_dir}"

[run]
steps = {steps}
thermo = 5000
"""

# LAMMPS counterpart: same single sphere, same wall, same material.
# Mapping to DIRT's Hertz model:
#   hertz/material E e nu     -> Young's modulus, restitution, Poisson ratio
#   damping tsuji             -> derive normal damping from the restitution e
#                                (LAMMPS's viscoelastic-COR model; DIRT uses
#                                 beta = -ln(e)/sqrt(pi^2+ln^2 e) for the same job)
#   tangential ... 0.0        -> friction = 0 (clean normal rebound)
#   nve/sphere                -> translational Velocity Verlet, no gravity
# A per-step trace (time, z, vz) of the single atom is written for post-processing.
LMP_TEMPLATE = """\
# Auto-generated LAMMPS input for the Hertz rebound sweep
# v0 = {v0} m/s, nominal COR = {cor}, tsuji restitution input = {e_rest}
units           si
atom_style      sphere
atom_modify     map array
dimension       3
boundary        f f f
newton          off
comm_modify     vel yes

region          simbox block -0.01 0.01 -0.01 0.01 0.0 0.1 units box
create_box      1 simbox

create_atoms    1 single 0.0 0.0 0.0075 units box
set             group all diameter {diam}
set             group all density {density}

pair_style      granular
pair_coeff      1 1 hertz/material {E} {e_rest} {nu} tangential linear_nohistory 0.0 0.0 damping tsuji rolling none twisting none

fix             wall all wall/gran granular hertz/material {E} {e_rest} {nu} tangential linear_nohistory 0.0 0.0 damping tsuji rolling none twisting none zplane 0.0 NULL
fix             integrate all nve/sphere

velocity        all set 0.0 0.0 -{v0} units box
timestep        {dt}

variable        tnow equal time
variable        zpos equal z[1]
variable        zvel equal vz[1]
fix             rec all print 1 "${{tnow}} ${{zpos}} ${{zvel}}" file {trace} screen no title "t z vz"

thermo          5000
run             {steps}
"""


def write_lammps_input(path, v0, cor_label, e_rest, trace, steps):
    """Write a LAMMPS input. `e_rest` is the restitution fed to the tsuji damping
    model; `cor_label` is just the nominal-COR tag recorded in the comment."""
    with open(path, "w") as f:
        f.write(LMP_TEMPLATE.format(
            v0=v0, cor=cor_label, e_rest=e_rest,
            E=f"{YOUNGS_MOD:.6e}", nu=POISSON_RATIO,
            diam=2.0 * RADIUS, density=DENSITY,
            dt=f"{dt_for(v0):.10e}", steps=steps, trace=trace,
        ))


def case_tag(v0, cor):
    return f"v{v0}_cor{cor}"


def case_dir(v0, cor):
    return os.path.join(SWEEP_DIR, case_tag(v0, cor))


def dt_for(_v0):
    """Timestep: a fraction of the Rayleigh critical timestep (solver default).
    DIRT and LAMMPS share this dt so the comparison is on equal footing."""
    g = YOUNGS_MOD / (2.0 * (1.0 + POISSON_RATIO))
    alpha = 0.1631 * POISSON_RATIO + 0.876605
    dt_rayleigh = math.pi * RADIUS / alpha * (DENSITY / g) ** 0.5
    return dt_rayleigh * 0.15


def steps_for(v0):
    """Step budget: time to fall the ~2 mm gap plus a generous contact/rebound
    margin."""
    fall_time = 0.003 / v0          # ~2 mm gap, no gravity
    total_time = fall_time * 3.0    # generous margin for contact + rebound
    return int(total_time / dt_for(v0)) + 10000


def find_lammps():
    """Return the first available LAMMPS binary, or None."""
    for b in LAMMPS_BINS:
        path = shutil.which(b)
        if path:
            return path
    return None


# ── generate ───────────────────────────────────────────────────────────────
def generate():
    n = 0
    for cor in CORS:
        for v0 in VELOCITIES:
            cdir = case_dir(v0, cor)
            os.makedirs(cdir, exist_ok=True)
            steps = steps_for(v0)

            # DIRT config.
            with open(os.path.join(cdir, "config.toml"), "w") as f:
                f.write(TOML_TEMPLATE.format(
                    v0=v0, cor=cor,
                    youngs_mod=f"{YOUNGS_MOD:.1e}",
                    poisson_ratio=POISSON_RATIO,
                    radius=RADIUS, density=DENSITY,
                    output_dir=cdir, steps=steps,
                ))

            # The LAMMPS input needs the calibrated restitution e', which is
            # back-solved at run time, so it is written during 'start' — not here.
            n += 1
    print(f"Generated {n} DIRT configs under {SWEEP_DIR}")


# ── start ──────────────────────────────────────────────────────────────────
CSV_FIELDS = ["input_v0", "input_cor", "v_impact", "v_rebound", "cor_measured",
              "contact_time", "max_overlap", "dt", "radius", "density"]


def parse_rebound_trace(trace_path, radius, dt):
    """Extract rebound metrics from a per-step (t, z, vz) trace, using the same
    contact logic as the DIRT example's track_rebound. Returns a result dict or
    None if no clean contact-then-separation was captured."""
    was_contact = False
    prev_vz = 0.0
    v_impact = v_rebound = t_start = t_end = max_overlap = 0.0
    with open(trace_path) as f:
        next(f, None)  # header
        for line in f:
            parts = line.split()
            if len(parts) != 3:
                continue
            t, z, vz = (float(x) for x in parts)
            overlap = radius - z
            in_contact = overlap > 0.0
            if not was_contact and not in_contact:
                prev_vz = vz
            elif not was_contact and in_contact:
                was_contact = True
                v_impact = prev_vz
                t_start = t
                max_overlap = overlap
            elif was_contact and in_contact:
                max_overlap = max(max_overlap, overlap)
            elif was_contact and not in_contact:
                v_rebound = vz
                t_end = t
                break
    if not was_contact or t_end == 0.0 or v_impact == 0.0:
        return None
    return {
        "v_impact": abs(v_impact),
        "v_rebound": abs(v_rebound),
        "cor_measured": abs(v_rebound / v_impact),
        "contact_time": t_end - t_start,
        "max_overlap": max_overlap,
        "dt": dt,
        "radius": radius,
        "density": DENSITY,
    }


def run_lammps_rebound(lammps, in_path, log_path, trace_path, v0):
    """Run one LAMMPS input and parse its rebound trace. Returns a result dict
    (without input_* tags) or None on failure."""
    proc = subprocess.run(
        [lammps, "-in", in_path, "-log", log_path],
        cwd=REPO_ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
    )
    if proc.returncode != 0 or not os.path.isfile(trace_path):
        return None
    return parse_rebound_trace(trace_path, RADIUS, dt_for(v0))


# Calibration velocity. The viscoelastic COR is velocity-independent, so one
# velocity suffices to back-solve the restitution input for every speed.
CALIB_V = 1.0
CALIB_TOL = 1e-3       # |measured COR - target| convergence tolerance
CALIB_MAX_ITERS = 16


def calibrate_restitution(lammps, target_cor, workdir):
    """Bisection: find the LAMMPS tsuji restitution input e' whose *measured* COR
    equals `target_cor` (i.e. matches DIRT). Monotone in e'. Returns e' or None."""
    os.makedirs(workdir, exist_ok=True)
    in_path = os.path.join(workdir, "in.lammps")
    log_path = os.path.join(workdir, "lammps.log")
    trace = os.path.join(workdir, "trace.txt")
    lo, hi = 0.2, 1.0   # hi=1.0 (alpha_tsuji->0) reaches the elastic anchor
    mid = 0.5 * (lo + hi)
    for _ in range(CALIB_MAX_ITERS):
        mid = 0.5 * (lo + hi)
        write_lammps_input(in_path, CALIB_V, target_cor, mid, trace, steps_for(CALIB_V))
        cor = run_lammps_rebound(lammps, in_path, log_path, trace, CALIB_V)
        cor = cor["cor_measured"] if cor else None
        if cor is None:
            return None
        if abs(cor - target_cor) < CALIB_TOL:
            return mid
        if cor < target_cor:
            lo = mid
        else:
            hi = mid
    return mid


def run_calibrated_lammps(lammps, dirt_results):
    """For each nominal COR, back-solve e' to match DIRT's measured COR, then run
    the full velocity sweep at that e'. Returns (rows, {cor: e'})."""
    dirt = as_data(dirt_results)
    print("\nCalibrating LAMMPS restitution to match DIRT (per nominal COR):")
    calib_e = {}
    for cor in CORS:
        target = dirt.get((cor, CALIB_V), {}).get("cor_meas")
        if target is None:
            continue
        e_prime = calibrate_restitution(lammps, target, os.path.join(SWEEP_DIR, "_calib"))
        if e_prime is not None:
            calib_e[cor] = e_prime
            print(f"  COR {cor:<4} -> e'={e_prime:.4f}  (DIRT target {target:.4f})")
        else:
            print(f"  COR {cor:<4} -> calibration FAILED")

    rows = []
    for cor in CORS:
        if cor not in calib_e:
            continue
        for v0 in VELOCITIES:
            cdir = case_dir(v0, cor)
            in_path = os.path.join(cdir, "in.lammps")
            log_path = os.path.join(cdir, "lammps.log")
            trace = os.path.join(cdir, "lammps_trace.txt")
            write_lammps_input(in_path, v0, cor, calib_e[cor], trace, steps_for(v0))
            res = run_lammps_rebound(lammps, in_path, log_path, trace, v0)
            if res:
                res["input_v0"], res["input_cor"] = str(v0), str(cor)
                rows.append(res)
    return rows, calib_e


def start():
    n_total = len(CORS) * len(VELOCITIES)
    os.makedirs(DATA_DIR, exist_ok=True)

    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE, "--no-default-features"],
        cwd=REPO_ROOT, check=True,
    )

    lammps = find_lammps()
    print(f"LAMMPS: {lammps}" if lammps
          else "LAMMPS: not found on PATH — running DIRT only.")

    dirt_results = []
    i = 0
    for cor in CORS:
        for v0 in VELOCITIES:
            i += 1
            cdir = case_dir(v0, cor)
            config = os.path.join(cdir, "config.toml")
            if not os.path.isfile(config):
                print(f"  [{i:2d}/{n_total}] missing {config} — run 'generate' first.")
                continue
            print(f"  [{i:2d}/{n_total}] v0={v0:<4} COR={cor:<5}", end="  ", flush=True)

            dirt_log = os.path.join(cdir, "run.log")
            with open(dirt_log, "w") as log:
                proc = subprocess.run(
                    ["cargo", "run", "--release", "--example", EXAMPLE,
                     "--no-default-features", "--", config],
                    cwd=REPO_ROOT, stdout=log, stderr=subprocess.STDOUT,
                )
            dirt_csv = os.path.join(cdir, "data", "rebound_results.csv")
            if proc.returncode == 0 and os.path.isfile(dirt_csv):
                with open(dirt_csv) as f:
                    row = next(csv.DictReader(f))
                row["input_v0"], row["input_cor"] = str(v0), str(cor)
                dirt_results.append(row)
                print(f"DIRT COR={float(row['cor_measured']):.4f}")
            else:
                print(f"DIRT FAILED ({dirt_log})")

    if not dirt_results:
        print("\nERROR: no DIRT results collected.")
        sys.exit(1)

    _write_csv(SWEEP_CSV, dirt_results)
    print(f"\nDIRT:   {len(dirt_results)}/{n_total} cases -> {SWEEP_CSV}")
    if lammps:
        # LAMMPS restitution is back-solved per nominal COR so its measured COR
        # matches DIRT (no exact e->damping formula exists for Hertz contact).
        lammps_results, _ = run_calibrated_lammps(lammps, dirt_results)
        if lammps_results:
            _write_csv(LAMMPS_CSV, lammps_results)
            print(f"LAMMPS: {len(lammps_results)}/{n_total} cases -> {LAMMPS_CSV}")


def _write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerows(rows)


# ── graph (validate + plot) ──────────────────────────────────────────────────
def hertz_contact_duration(v0):
    return 2.87 * (M_EFF**2 / (R_EFF * E_STAR**2 * v0))**0.2


def hertz_max_overlap(v0):
    return (15.0 * M_EFF * v0**2 / (16.0 * R_EFF**0.5 * E_STAR))**0.4


# Inelastic (viscoelastic-Hertz) reference. The undamped formulas above ignore
# the damping; this integrates the same 1-DOF normal contact ODE the solver
# uses — m δ̈ = -[k_n δ + c(δ) δ̇], clamped ≥ 0 — to predict COR, contact
# duration, and peak overlap including dissipation.
#   k_n   = (4/3) E* √(R δ)              Hertz spring
#   c(δ)  = 2 β·CFAC·√(S_n m),  S_n = 2 E* √(R δ),  β = -ln e / √(π²+ln²e)
# CFAC is the solver's damping constant. NOTE the solver's constant named
# `SQRT_5_3` is actually valued √(5/6) ≈ 0.91287 (not √(5/3) ≈ 1.29099); the
# √(5/6) value is the physically correct one (it reproduces measured COR ≈ e).
CFAC_MODEL = math.sqrt(5.0 / 6.0)   # DIRT's actual damping constant (correct value)


def contact_ode(e, v0, cfac):
    """Integrate the 1-DOF viscoelastic Hertz contact; return (cor, t_c, dmax)."""
    A = 4.0 / 3.0 * E_STAR * math.sqrt(RADIUS)          # F_el = A δ^1.5
    beta = 0.0 if e >= 1.0 else -math.log(e) / math.sqrt(math.pi**2 + math.log(e)**2)
    c0 = 2.0 * beta * cfac * math.sqrt(2.0 * E_STAR * MASS) * RADIUS**0.25  # c(δ)=c0 δ^0.25
    dt = hertz_contact_duration(v0) / 6000.0

    def f(d, v):
        if d <= 0.0:
            return 0.0
        fr = A * d**1.5 + c0 * d**0.25 * v
        return fr if fr > 0.0 else 0.0

    d, v, t, dmax = 0.0, v0, 0.0, 0.0
    for _ in range(2_000_000):
        a1 = -f(d, v) / MASS
        a2 = -f(d + 0.5 * dt * v, v + 0.5 * dt * a1) / MASS
        a3 = -f(d + 0.5 * dt * (v + 0.5 * dt * a1), v + 0.5 * dt * a2) / MASS
        a4 = -f(d + dt * (v + 0.5 * dt * a2), v + dt * a3) / MASS
        dn = d + dt / 6.0 * (v + 2 * (v + 0.5 * dt * a1) + 2 * (v + 0.5 * dt * a2) + (v + dt * a3))
        vn = v + dt / 6.0 * (a1 + 2 * a2 + 2 * a3 + a4)
        t += dt
        if dn > dmax:
            dmax = dn
        if d > 0.0 and dn <= 0.0:
            frac = d / (d - dn)
            return abs(v + frac * (vn - v)) / v0, t - dt + frac * dt, dmax
        d, v = dn, vn
    return abs(v) / v0, t, dmax


def overlap_tol(cor_input):
    if cor_input >= 0.85:
        return OVERLAP_TOL_HIGH
    if cor_input >= 0.6:
        return OVERLAP_TOL_MED
    return OVERLAP_TOL_LOW


def load_rows():
    if not os.path.isfile(SWEEP_CSV):
        print(f"ERROR: {SWEEP_CSV} not found.")
        print("Run the sweep first: python3 examples/bench_hertz_rebound/sweep.py start")
        sys.exit(1)
    with open(SWEEP_CSV) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        print("ERROR: no data in sweep results file.")
        sys.exit(1)
    return rows


def load_optional(path):
    """Load a results CSV if it exists, else return []."""
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def as_data(rows):
    """Index rows by (cor, v0) -> measured quantities."""
    return {
        (float(r["input_cor"]), float(r["input_v0"])): {
            "cor_meas": float(r["cor_measured"]),
            "tc": float(r["contact_time"]),
            "delta_max": float(r["max_overlap"]),
            "v_impact": float(r["v_impact"]),
        }
        for r in rows
    }


def validate(rows):
    """Run the four checks; return True if everything passed."""
    print("=" * 65)
    print("Hertz Contact Rebound Benchmark Validation")
    print("=" * 65)
    print(f"  E* = {E_STAR:.3e} Pa")
    print(f"  m  = {MASS:.6e} kg")
    print(f"  R  = {RADIUS*1000:.1f} mm\n")

    total = passed = 0
    cor_pass = tc_pass = ov_pass = 0
    n = len(rows)

    for row in rows:
        v0 = float(row["input_v0"])
        cor_input = float(row["input_cor"])
        cor_meas = float(row["cor_measured"])
        tc_meas = float(row["contact_time"])
        delta_max_meas = float(row["max_overlap"])
        v_impact = float(row["v_impact"])

        tc_theory = hertz_contact_duration(v_impact)
        delta_max_theory = hertz_max_overlap(v_impact)

        cor_err = abs(cor_meas - cor_input) / cor_input
        cor_tol = COR_TOL_HIGH if cor_input >= 0.7 else COR_TOL_LOW
        status_cor = "PASS" if cor_err <= cor_tol else "FAIL"
        cor_pass += status_cor == "PASS"

        tc_err = abs(tc_meas - tc_theory) / tc_theory
        status_tc = "PASS" if tc_err <= CONTACT_TIME_TOL else "FAIL"
        tc_pass += status_tc == "PASS"

        ov_err = abs(delta_max_meas - delta_max_theory) / delta_max_theory
        status_ov = "PASS" if ov_err <= overlap_tol(cor_input) else "FAIL"
        ov_pass += status_ov == "PASS"

        total += 3
        passed += (status_cor == "PASS") + (status_tc == "PASS") + (status_ov == "PASS")

        print(f"v0={v0:.1f} m/s, COR_in={cor_input:.2f}:")
        print(f"  COR:     {cor_meas:.4f} vs {cor_input:.4f}  (err={cor_err*100:.2f}%)  [{status_cor}]")
        print(f"  t_c:     {tc_meas:.3e} vs {tc_theory:.3e} s  (err={tc_err*100:.1f}%)  [{status_tc}]")
        print(f"  d_max:   {delta_max_meas:.3e} vs {delta_max_theory:.3e} m  (err={ov_err*100:.1f}%)  [{status_ov}]")

    total += 1
    expected = len(VELOCITIES) * len(CORS)
    complete = n == expected
    passed += complete
    print(f"\nCompleteness: {n}/{expected} cases  [{'PASS' if complete else 'FAIL'}]")
    print()
    print(f"COR checks:          {cor_pass}/{n} passed")
    print(f"Contact time checks: {tc_pass}/{n} passed")
    print(f"Overlap checks:      {ov_pass}/{n} passed")
    print(f"\nOverall: {passed}/{total} checks passed")
    print("ALL CHECKS PASSED" if passed == total
          else f"WARNING: {total - passed} check(s) failed")
    return passed == total


def plot(dirt_rows, lammps_rows):
    """Write the three benchmark figures, overlaying DIRT (filled markers) and
    LAMMPS (open markers, restitution tuned so its COR matches DIRT). Skips
    gracefully without matplotlib."""
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except Exception as e:
        print(f"\n(matplotlib/numpy unavailable, skipped plots: {e})")
        return

    plt.rcParams.update({
        "font.size": 12, "axes.labelsize": 14, "axes.titlesize": 14,
        "legend.fontsize": 10, "figure.dpi": 150, "savefig.dpi": 150,
    })
    save = dict(bbox_inches="tight")
    markers = ["o", "s", "^", "D", "v"]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    MS = 5      # marker size
    ALPHA = 0.55  # transparency so overlapping DIRT/LAMMPS markers both show

    os.makedirs(PLOT_DIR, exist_ok=True)
    dirt = as_data(dirt_rows)
    lammps = as_data(lammps_rows)
    cors = sorted({float(r["input_cor"]) for r in dirt_rows})
    vels = sorted({float(r["input_v0"]) for r in dirt_rows})

    # Proxy legend handles distinguishing DIRT (filled) from LAMMPS (open).
    code_proxies = []
    if lammps:
        code_proxies = [
            Line2D([], [], color="gray", marker="o", linestyle="none",
                   markersize=MS, alpha=ALPHA, label="DIRT"),
            Line2D([], [], color="gray", marker="o", linestyle="none",
                   markerfacecolor="none", markersize=MS, alpha=ALPHA, label="LAMMPS"),
        ]

    def overlay_open(ax, src, xs, ys, m, c):
        pts = [(xs(k), ys(k)) for k in src]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], m, color=c,
                    markersize=MS, markerfacecolor="none", alpha=ALPHA)

    print()
    # ── Plot 1: measured vs input COR (target is the 1:1 line) ──────────────
    fig, ax = plt.subplots(figsize=(6.8, 5))
    for iv, v0 in enumerate(vels):
        m, c = markers[iv % 4], colors[iv % 4]
        keys = [(cc, v0) for cc in cors]
        # label only the DIRT series (it carries the per-velocity color legend)
        dp = [(cc, dirt[(cc, v0)]["cor_meas"]) for cc in cors if (cc, v0) in dirt]
        if dp:
            ax.plot([p[0] for p in dp], [p[1] for p in dp], m, color=c,
                    markersize=MS, alpha=ALPHA, label=f"v = {v0} m/s")
        overlay_open(ax, [k for k in keys if k in lammps],
                     lambda k: k[0], lambda k: lammps[k]["cor_meas"], m, c)
    ax.plot([0, 1], [0, 1], "k--", linewidth=1, label="Ideal (1:1)")
    # Viscoelastic-model curve (COR is velocity-independent → evaluate at v=1).
    e_in = np.linspace(0.45, 1.0, 25)
    ax.plot(e_in, [contact_ode(e, 1.0, CFAC_MODEL)[0] for e in e_in],
            "-", color="0.35", linewidth=1.4, label="viscoelastic model")
    ax.set_xlabel("Input COR")
    ax.set_ylabel("Measured COR")
    ax.set_title("Coefficient of Restitution: Input vs Measured")
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles=handles + code_proxies, loc="upper left", fontsize=8)
    ax.set_xlim(0.4, 1.03)
    ax.set_ylim(0.4, 1.03)
    ax.set_aspect("equal")
    fig.savefig(os.path.join(PLOT_DIR, "cor_validation.png"), **save)
    plt.close(fig)
    print(f"Saved: {PLOT_DIR}/cor_validation.png")

    v_theory = np.linspace(0.08, 2.5, 200)

    v_curve = np.logspace(math.log10(0.08), math.log10(2.5), 18)
    ode_idx = {"tc": 1, "delta_max": 2}

    def vs_velocity(key, ylabel, title, fname, theory_fn):
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(v_theory, [theory_fn(v) * 1e6 for v in v_theory],
                "k-", linewidth=2, label="Hertz theory (elastic)")
        idx = ode_idx[key]
        for ic, cor in enumerate(cors):
            m, c = markers[ic % 5], colors[ic % 5]
            keys = [(cor, v0) for v0 in vels]
            dp = [(dirt[(cor, v0)]["v_impact"], dirt[(cor, v0)][key])
                  for v0 in vels if (cor, v0) in dirt]
            if dp:
                ax.plot([p[0] for p in dp], [p[1] * 1e6 for p in dp], m, color=c,
                        markersize=MS, alpha=ALPHA, label=f"COR = {cor}")
            overlay_open(ax, [k for k in keys if k in lammps],
                         lambda k: lammps[k]["v_impact"],
                         lambda k: lammps[k][key] * 1e6, m, c)
            # viscoelastic-model curve (the inelastic correction the points match)
            ax.plot(v_curve, [contact_ode(cor, v, CFAC_MODEL)[idx] * 1e6 for v in v_curve],
                    "--", color=c, linewidth=1.0, alpha=0.9)
        model_proxy = [Line2D([], [], color="0.35", linestyle="--", linewidth=1.0,
                              label="viscoelastic model")]
        ax.set_xlabel("Impact velocity [m/s]")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        handles, _ = ax.get_legend_handles_labels()
        ax.legend(handles=handles + code_proxies + model_proxy, fontsize=8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        fig.savefig(os.path.join(PLOT_DIR, fname), **save)
        plt.close(fig)
        print(f"Saved: {PLOT_DIR}/{fname}")

    # ── Plot 2: contact duration vs impact velocity ─────────────────────────
    vs_velocity("tc", "Contact duration [µs]",
                "Contact Duration vs Impact Velocity",
                "contact_duration.png", hertz_contact_duration)
    # ── Plot 3: peak overlap vs impact velocity ─────────────────────────────
    vs_velocity("delta_max", "Peak overlap [µm]",
                "Peak Overlap vs Impact Velocity",
                "peak_overlap.png", hertz_max_overlap)


def compare_codes(dirt_rows, lammps_rows):
    """Print a per-case DIRT-vs-LAMMPS measured-COR comparison."""
    dirt = as_data(dirt_rows)
    lammps = as_data(lammps_rows)
    print("\n" + "=" * 50)
    print("Measured COR: DIRT vs LAMMPS (restitution tuned to DIRT)")
    print("=" * 50)
    print(f"  {'v0':>5}{'COR':>6} | {'DIRT':>8}{'LAMMPS':>9} | {'diff':>8}")
    for cor in sorted({float(r['input_cor']) for r in dirt_rows}):
        for v0 in sorted({float(r['input_v0']) for r in dirt_rows}):
            k = (cor, v0)
            if k not in dirt or k not in lammps:
                continue
            d, l = dirt[k]["cor_meas"], lammps[k]["cor_meas"]
            print(f"  {v0:>5}{cor:>6} | {d:>8.4f}{l:>9.4f} | {l - d:>+8.4f}")


def graph():
    dirt_rows = load_rows()
    lammps_rows = load_optional(LAMMPS_CSV)
    ok = validate(dirt_rows)
    if lammps_rows:
        compare_codes(dirt_rows, lammps_rows)
    else:
        print(f"\n(no {os.path.basename(LAMMPS_CSV)} — plotting DIRT only)")
    plot(dirt_rows, lammps_rows)
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
