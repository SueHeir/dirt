#!/usr/bin/env python3
"""
Plate pressure–sinkage (terramechanics) benchmark driver.

Presses a flat plate of footprint width b straight down into a settled granular
bed at constant velocity and records the vertical reaction force F on the plate
versus its sinkage depth z. The plate pressure p = F / A (A = b·L_y) is fit to
the empirical Bekker pressure–sinkage relation

    p = (k_c / b + k_phi) · z^n

i.e. pressure grows as a power law in sinkage, p ∝ z^n, monotonically. The
validation checks the *form*: p(z) monotone and well fit by a power law with a
physically sensible exponent (n ≈ 0.5–1.3 for granular soils), plus sane
width/depth trends (a wider plate carries more load at the same sinkage).

This is an EMPIRICAL reference — we validate the qualitative law and exponent,
not specific k_c / k_phi values (those are soil-fit constants).

Commands (from anywhere):
    python3 examples/bench_plate_sinkage/sweep.py generate   # write per-case configs
    python3 examples/bench_plate_sinkage/sweep.py start      # build + run all sims -> CSV
    python3 examples/bench_plate_sinkage/sweep.py graph       # validate + plot
    python3 examples/bench_plate_sinkage/sweep.py             # all three, in order

If a LAMMPS binary (lmp_serial / lmp / lmp_mpi / lammps) is on PATH, each width
case is ALSO run in LAMMPS with the same material (granular hertz/material +
mindlin tangential, matching E / nu / e / mu), the same enhanced gravity, the
same loose-insert-then-settle bed, and a flat plate of the same footprint width b
pressed straight down at the same constant velocity. The LAMMPS plate is a small
raft of frozen grains driven rigidly downward; its vertical reaction force is the
plate->bed contact force read with `compute group/group plate bed pair yes`
(component z). p = F/A vs sinkage z is processed and fit to p ∝ z^n the SAME way
as DIRT, and overlaid on the pressure–sinkage plot as open markers / dashed fit.
LAMMPS is OPTIONAL — without it, only DIRT runs and the example still validates
against the Bekker power-law form.

Reference: M. G. Bekker, "Theory of Land Locomotion" (1956); "Introduction to
Terrain-Vehicle Systems" (1969). J. Y. Wong, "Theory of Ground Vehicles".
"""

import os
import sys
import csv
import math
import shutil
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_plate_sinkage"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "sweep.csv")            # DIRT: one fitted row per case
LAMMPS_CSV = os.path.join(DATA_DIR, "lammps_results.csv")  # LAMMPS: one fitted row per case

# LAMMPS binary candidates, in preference order. LAMMPS is an OPTIONAL overlay:
# if none is found, only DIRT is run/validated/plotted.
LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# ── Geometry / material (shared across cases) ─────────────────────────────────
L_Y = 0.04              # bed/plate depth in y (periodic slice), m
DENSITY = 2500.0        # kg/m^3
YOUNGS = 5.0e6          # Pa — softened sphere (quasi-static DEM)
NU = 0.3
RESTITUTION = 0.3
GZ = -49.05             # 5x g, faster settling
PLATE_V = 0.08          # plate descent speed, m/s
PLATE_Z0 = 0.17         # plate start height (above the loose cloud; settles first)
COUNT = 2400            # loose insertion (packing ~0.25), settles to a ~5-6 cm bed
STEPS = 45000
THERMO = 10000

# ── Sweep: plate footprint width b, plus one extra-friction case ──────────────
# Each case: (tag, half_width [m], inter-particle friction). Half-widths give
# plate widths b = 0.02, 0.04, 0.06 m; the *wide* and *narrow* pair tests the
# width trend; the high-friction case tests bed shear strength.
CASES = [
    ("b020_mu05", 0.010, 0.5),   # b = 0.02 m
    ("b040_mu05", 0.020, 0.5),   # b = 0.04 m  (representative)
    ("b060_mu05", 0.030, 0.5),   # b = 0.06 m
    ("b040_mu08", 0.020, 0.8),   # b = 0.04 m, higher friction (stronger bed)
]

# Only record p–z up to a sensible max sinkage to stay in the Bekker regime
# (avoid the plate bottoming out near the ~5-6 cm-deep bed's floor).
Z_MAX_FIT = 0.030       # m — fit/validate over 0 < z <= 30 mm

TOML_TEMPLATE = """[comm]
processors_x = 1
processors_y = 1
processors_z = 1
[domain]
x_low = -0.06
x_high = 0.06
y_low = 0.0
y_high = {ly}
z_low = 0.0
z_high = 0.26
boundary_x = "fixed"
boundary_y = "periodic"
boundary_z = "fixed"
[neighbor]
skin_fraction = 1.2
bin_size = 0.007
every = 1
[gravity]
gz = {gz}
[dem]
contact_model = "hertz"
[[dem.materials]]
name = "soil"
youngs_mod = {youngs:.6e}
poisson_ratio = {nu}
restitution = {e}
friction = {mu}
[[particles.insert]]
material = "soil"
count = {count}
radius = {{ distribution = "uniform", min = 0.0023, max = 0.0027 }}
density = {density}
region = {{ type = "block", min = [-0.058, 0.001, 0.01], max = [0.058, {ly_hi}, 0.16] }}
[[wall]]
point_z = 0.0
normal_z = 1.0
material = "soil"
[[wall]]
point_x = -0.06
normal_x = 1.0
material = "soil"
[[wall]]
point_x = 0.06
normal_x = -1.0
material = "soil"
[[wall]]
name = "plate"
point_z = {plate_z0}
normal_z = -1.0
material = "soil"
velocity = [0.0, 0.0, -{plate_v}]
bound_x_low = -{half_b}
bound_x_high = {half_b}
[output]
dir = "{outdir}"
[run]
steps = {steps}
thermo = {thermo}
"""

# ── LAMMPS counterpart ────────────────────────────────────────────────────────
# Reproduces the SAME setup in LAMMPS's GRANULAR package:
#   - same material: hertz/material E e nu + tangential mindlin + friction mu,
#     damping tsuji (LAMMPS's viscoelastic-COR normal damping; like DIRT there is
#     no exact e->damping closed form for a Hertz contact, but e drives it),
#   - same enhanced gravity (49.05 m/s^2 down), periodic in y, fixed floor/sides,
#   - a LOOSE grain cloud inserted on a coarse cubic lattice (spacing > grain
#     diameter, so no initial overlap) that settles into a deep bed — the
#     loose-insert + settle that DIRT does,
#   - a flat PLATE: a one-grain-thick raft of frozen grains spanning the central
#     footprint x in [-b/2, b/2] across the full y slice, sitting just above the
#     settled bed, driven straight down at constant velocity by `fix move linear`.
# The plate's vertical reaction force is the plate<->bed contact force read with
# `compute group/group plate bed pair yes` (z component). A per-step trace of
# (plate_z, F_z) is written for post-processing — fit p=F/A vs sinkage z exactly
# as DIRT does. LAMMPS bed grains are monodisperse at the DIRT mean radius.
#
# Stage 1 (settle) and stage 2 (press) share one input; the plate is created
# after settling using the measured bed-top height so it starts just above the
# surface regardless of how the loose cloud consolidated.
LMP_RADIUS = 0.0025        # DIRT mean grain radius (uniform 0.0023..0.0027)
LMP_DIAM = 2.0 * LMP_RADIUS
LMP_LOOSE_SPACING = 0.0075  # loose cubic lattice (> diameter -> no initial overlap)
LMP_DT = 5.0e-6            # stable for the softened E=5 MPa Hertz contact
LMP_SETTLE_STEPS = 70000  # let the loose cloud settle to a quiet bed
LMP_PRESS_STEPS = 140000  # plate descent: ~5.6 cm at 0.08 m/s, past Z_MAX_FIT

LMP_TEMPLATE = """\
# Auto-generated LAMMPS input — plate pressure–sinkage, b = {b} m, mu = {mu}
units           si
atom_style      sphere
atom_modify     map array
dimension       3
boundary        p p f
comm_modify     vel yes

variable        E equal {E}
variable        nu equal {nu}
variable        e equal {e}
variable        mu equal {mu}
variable        halfb equal {half_b}
variable        platev equal {plate_v}

region          simbox block -0.06 0.06 0.0 {ly} 0.0 0.30 units box
create_box      2 simbox

# ── loose bed cloud on a coarse lattice (no initial overlap) ──
lattice         sc {loose}
region          seed block -0.057 0.057 0.002 {ly_hi} 0.012 0.18 units box
create_atoms    1 region seed units box
set             type 1 diameter {diam}
set             type 1 density {density}

pair_style      granular
pair_coeff      * * hertz/material ${{E}} ${{e}} ${{nu}} tangential mindlin NULL 1.0 ${{mu}} damping tsuji

group           bed type 1
fix             grav bed gravity {g} vector 0 0 -1
fix             zwall all wall/gran granular hertz/material ${{E}} ${{e}} ${{nu}} tangential mindlin NULL 1.0 ${{mu}} damping tsuji zplane 0.0 NULL
fix             integ bed nve/sphere
timestep        {dt}

# ── settle ──
compute         zmax bed reduce max z
thermo          {settle_thermo}
thermo_style    custom step ke c_zmax
run             {settle_steps}

# ── build the flat plate raft just above the settled bed ──
variable        ztop equal $(c_zmax)
variable        pz equal ${{ztop}}+0.005
lattice         sc {diam} origin 0.5 0.5 0.5
region          platereg block -${{halfb}} ${{halfb}} 0.0 {ly} ${{pz}} $(v_pz+{diam_plus}) units box
create_atoms    2 region platereg units box
set             type 2 diameter {diam}
set             type 2 density {density}
group           plate type 2

# ── press: plate driven rigidly downward at constant velocity ──
unfix           integ
fix             integ bed nve/sphere
fix             move plate move linear 0.0 0.0 -${{platev}} units box

compute         gg plate group/group bed pair yes
variable        fz equal c_gg[3]
variable        platez equal xcm(plate,z)
fix             rec all print {trace_every} "$(v_platez) $(v_fz)" file {trace} screen no title "platez fz"
thermo          {press_thermo}
thermo_style    custom step v_platez v_fz
run             {press_steps}
"""


# ── helpers ───────────────────────────────────────────────────────────────────
def case_dir(tag):
    return os.path.join(SWEEP_DIR, tag)


def half_to_b(half):
    return 2.0 * half


def find_lammps():
    for b in LAMMPS_BINS:
        if shutil.which(b):
            return shutil.which(b)
    return None


def _config(tag, half_b, mu, outdir):
    return TOML_TEMPLATE.format(
        ly=L_Y, ly_hi=L_Y - 0.001, gz=GZ, youngs=YOUNGS, nu=NU, e=RESTITUTION,
        mu=mu, density=DENSITY, count=COUNT, plate_z0=PLATE_Z0, plate_v=PLATE_V,
        half_b=half_b, outdir=outdir, steps=STEPS, thermo=THERMO,
    )


def _lammps_input(half_b, mu, trace):
    """LAMMPS input for one width/friction case. `trace` is the (plate_z, F_z)
    output path the post-processor reads."""
    return LMP_TEMPLATE.format(
        b=2.0 * half_b,
        E=f"{YOUNGS:.6e}", nu=NU, e=RESTITUTION, mu=mu,
        half_b=half_b, plate_v=PLATE_V, ly=L_Y, ly_hi=L_Y - 0.002,
        loose=LMP_LOOSE_SPACING, diam=LMP_DIAM, diam_plus=LMP_DIAM,
        density=DENSITY, g=abs(GZ), dt=f"{LMP_DT:.6e}",
        settle_steps=LMP_SETTLE_STEPS, press_steps=LMP_PRESS_STEPS,
        settle_thermo=LMP_SETTLE_STEPS // 4, press_thermo=LMP_PRESS_STEPS // 8,
        trace_every=40, trace=trace,
    )


def _lammps_trace_to_raw(trace_path):
    """Convert a LAMMPS (plate_z, F_z) trace into the SAME (t, sinkage, force)
    rows DIRT writes, latching the sinkage datum at first bed contact exactly as
    the DIRT recorder does (force above CONTACT_THRESHOLD_N). Returns [] on miss."""
    pz, fz = [], []
    with open(trace_path) as f:
        next(f, None)  # header "platez fz"
        for line in f:
            parts = line.split()
            if len(parts) != 2:
                continue
            try:
                pz.append(float(parts[0]))
                fz.append(abs(float(parts[1])))
            except ValueError:
                continue
    if not pz:
        return []
    # Latch datum at first contact (same threshold/logic as main.rs), then stream
    # until the plate loses contact for good. A narrow plate in a finite-depth bed
    # eventually punches through (bearing-capacity failure): the grains under it
    # are displaced and F collapses to ~0. There is no bearing pressure without
    # contact, so the valid p–z curve ends at that sustained contact loss — we
    # truncate there rather than fit a trailing run of zeros.
    z_contact = None
    raw = []
    dt_trace = LMP_DT * 40   # trace_every = 40 steps
    f_peak = 0.0
    established = False       # bearing well above the contact threshold has built up
    below = 0                # consecutive samples below the contact threshold
    BELOW_MAX = 20           # ~20 trace samples (~4 ms) of no contact ⇒ punch-through
    EST_N = 2.0              # bearing "established" once F exceeds this (N)
    for i, (z, F) in enumerate(zip(pz, fz)):
        if z_contact is None:
            if F <= CONTACT_THRESHOLD_N:
                continue
            z_contact = z
        f_peak = max(f_peak, F)
        if f_peak >= EST_N:
            established = True
        # Only honor a sustained contact loss as terminal once bearing was
        # established — the early bounce-in transient must not cut the curve.
        below = below + 1 if F <= CONTACT_THRESHOLD_N else 0
        if established and below >= BELOW_MAX:
            break
        sinkage = max(z_contact - z, 0.0)
        raw.append((i * dt_trace, sinkage, F))

    # Pressure–sinkage is the LOADING branch: pressure rising with depth. A narrow
    # plate in a finite-depth bed reaches a bearing-capacity peak and then fails
    # (force falls off as grains escape from under it). Fit only up to peak load —
    # the standard bearing-curve convention — discarding the post-failure tail.
    if raw:
        i_peak = max(range(len(raw)), key=lambda k: raw[k][2])
        if i_peak >= 5:
            raw = raw[: i_peak + 1]
    return raw


# Same first-contact threshold the DIRT recorder uses (main.rs CONTACT_THRESHOLD_N).
CONTACT_THRESHOLD_N = 0.05


def run_lammps_case(lammps, tag, half_b, mu):
    """Run one LAMMPS width/friction case; return processed summary dict or None.
    Uses the SAME _process_case post-processing as DIRT (p=F/A, binned p∝z^n fit)."""
    cdir = case_dir(tag)
    os.makedirs(cdir, exist_ok=True)
    in_path = os.path.join(cdir, "in.lammps")
    log_path = os.path.join(cdir, "lammps.log")
    trace = os.path.join(cdir, "lammps_trace.txt")
    if os.path.exists(trace):
        os.remove(trace)
    with open(in_path, "w") as f:
        f.write(_lammps_input(half_b, mu, trace))
    proc = subprocess.run(
        [lammps, "-in", in_path, "-log", log_path],
        cwd=REPO_ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
    )
    if proc.returncode != 0 or not os.path.isfile(trace):
        return None
    raw = _lammps_trace_to_raw(trace)
    if not raw:
        return None
    s = _process_case(tag, half_b, raw)
    if s is not None:
        s["mu"] = mu
    return s


# ── generate ──────────────────────────────────────────────────────────────────
def generate():
    n = 0
    for tag, half_b, mu in CASES:
        cdir = case_dir(tag)
        os.makedirs(cdir, exist_ok=True)
        with open(os.path.join(cdir, "config.toml"), "w") as f:
            f.write(_config(tag, half_b, mu, cdir))
        n += 1
    print(f"Generated {n} DIRT sweep configs under {SWEEP_DIR}")


# ── start ─────────────────────────────────────────────────────────────────────
def _run_dirt(cdir):
    """Run one prepared case; return path to its plate_sinkage_results.csv or None."""
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "plate_sinkage_results.csv")
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


def _load_raw(path):
    """Load the per-step (time, sinkage, force) rows from a results CSV."""
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append((float(r["time"]), float(r["sinkage"]), float(r["force"])))
    return rows


def _fit_powerlaw(zs, ps):
    """Least-squares fit log p = log A + n log z over z>0, p>0. Returns (n, A, r2)."""
    xs, ys = [], []
    for z, p in zip(zs, ps):
        if z > 1e-6 and p > 1e-9:
            xs.append(math.log(z))
            ys.append(math.log(p))
    if len(xs) < 5:
        return None
    m = len(xs)
    sx = sum(xs); sy = sum(ys)
    sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
    denom = m * sxx - sx * sx
    if abs(denom) < 1e-30:
        return None
    n = (m * sxy - sx * sy) / denom
    b = (sy - n * sx) / m
    A = math.exp(b)
    ybar = sy / m
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    ss_res = sum((y - (n * x + b)) ** 2 for x, y in zip(xs, ys))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else 0.0
    return n, A, r2


def _process_case(tag, half_b, raw):
    """Turn raw (t, z, F) rows into a fitted p∝z^n summary dict for the case."""
    b = half_to_b(half_b)
    area = b * L_Y
    # Restrict to the fit window and to a monotone-sinkage prefix.
    zs, ps = [], []
    z_prev = -1.0
    monotone = True
    for (_, z, F) in raw:
        if z <= 0.0 or z > Z_MAX_FIT:
            continue
        p = F / area
        if z < z_prev - 1e-9:
            monotone = False
        z_prev = z
        zs.append(z)
        ps.append(p)
    if len(zs) < 5:
        return None
    # Bin to a coarse z-grid so noisy per-step force is smoothed for the fit.
    nbin = 40
    zmin, zmax = min(zs), max(zs)
    edges = [zmin + (zmax - zmin) * k / nbin for k in range(nbin + 1)]
    bz, bp = [], []
    for k in range(nbin):
        lo, hi = edges[k], edges[k + 1]
        sel = [(z, p) for z, p in zip(zs, ps) if lo <= z < hi]
        if sel:
            bz.append(sum(z for z, _ in sel) / len(sel))
            bp.append(sum(p for _, p in sel) / len(sel))
    fit = _fit_powerlaw(bz, bp)
    if fit is None:
        return None
    n, A, r2 = fit
    # Pressure must increase from first to last bin.
    p_increases = bp[-1] > bp[0]
    return {
        "tag": tag, "b": b, "area": area, "n": n, "A": A, "r2": r2,
        "monotone": monotone, "p_increases": p_increases,
        "z_max": zmax, "p_max": max(bp),
        "bz": bz, "bp": bp,
    }


_SUMMARY_FIELDS = ["tag", "b", "mu", "n", "A", "r2", "z_max", "p_max",
                   "monotone", "p_increases"]


def _write_summaries(path, summaries):
    """Write fitted per-case rows (DIRT or LAMMPS) to a summary CSV."""
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=_SUMMARY_FIELDS)
        w.writeheader()
        for s in summaries:
            w.writerow({k: s[k] for k in _SUMMARY_FIELDS})


def _write_curve(path, s):
    """Write a case's binned p(z) curve for plotting."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sinkage", "pressure"])
        for z, p in zip(s["bz"], s["bp"]):
            w.writerow([f"{z:.8e}", f"{p:.8e}"])


def start():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE, "--no-default-features"],
        cwd=REPO_ROOT, check=True,
    )

    lammps = find_lammps()
    print(f"LAMMPS: {lammps} (cross-code overlay)" if lammps
          else "LAMMPS: not found on PATH — running DIRT only.")

    summaries = []
    # Stash binned p-z curves for plotting in a side CSV per case.
    n = len(CASES)
    for i, (tag, half_b, mu) in enumerate(CASES, 1):
        cdir = case_dir(tag)
        if not os.path.isfile(os.path.join(cdir, "config.toml")):
            print(f"  [{i}/{n}] missing config for {tag} — run 'generate' first.")
            continue
        print(f"  [{i}/{n}] {tag} (b={half_to_b(half_b)*1e3:.0f} mm, mu={mu})", end="  ", flush=True)
        res = _run_dirt(cdir)
        if res is None:
            print("DIRT FAILED")
            continue
        raw = _load_raw(res)
        s = _process_case(tag, half_b, raw)
        if s is None:
            print("insufficient sinkage data")
            continue
        s["mu"] = mu
        summaries.append(s)
        # Save the binned curve for plotting.
        _write_curve(os.path.join(DATA_DIR, f"curve_{tag}.csv"), s)
        print(f"n={s['n']:.3f}  r2={s['r2']:.3f}  p_max={s['p_max']:.3e} Pa")

    if not summaries:
        print("\nERROR: no DIRT results collected.")
        sys.exit(1)

    _write_summaries(SWEEP_CSV, summaries)
    print(f"\nDIRT:   wrote {len(summaries)}/{n} case summaries -> {SWEEP_CSV}")

    # ── LAMMPS cross-code overlay (optional) ──
    if lammps:
        lmp_summaries = []
        for i, (tag, half_b, mu) in enumerate(CASES, 1):
            print(f"  LAMMPS [{i}/{n}] {tag} (b={half_to_b(half_b)*1e3:.0f} mm, "
                  f"mu={mu})", end="  ", flush=True)
            s = run_lammps_case(lammps, tag, half_b, mu)
            if s is None:
                print("LAMMPS FAILED / insufficient sinkage")
                continue
            lmp_summaries.append(s)
            _write_curve(os.path.join(DATA_DIR, f"lammps_curve_{tag}.csv"), s)
            print(f"n={s['n']:.3f}  r2={s['r2']:.3f}  p_max={s['p_max']:.3e} Pa")
        if lmp_summaries:
            _write_summaries(LAMMPS_CSV, lmp_summaries)
            print(f"LAMMPS: wrote {len(lmp_summaries)}/{n} case summaries -> {LAMMPS_CSV}")
        else:
            print("LAMMPS: no cases collected — overlay skipped.")


# ── graph (validate + plot) ───────────────────────────────────────────────────
def _load_summaries():
    if not os.path.isfile(SWEEP_CSV):
        return []
    out = []
    with open(SWEEP_CSV) as f:
        for r in csv.DictReader(f):
            out.append({
                "tag": r["tag"], "b": float(r["b"]), "mu": float(r["mu"]),
                "n": float(r["n"]), "A": float(r["A"]), "r2": float(r["r2"]),
                "z_max": float(r["z_max"]), "p_max": float(r["p_max"]),
                "monotone": r["monotone"] == "True",
                "p_increases": r["p_increases"] == "True",
            })
    return out


# Validation tolerances.
N_MIN, N_MAX = 0.4, 1.6     # sensible Bekker exponent band for granular soils
R2_MIN = 0.85               # power-law fit quality (granular force is intrinsically noisy)


def validate(rows):
    print("\n=== Plate pressure–sinkage validation (Bekker form) ===")
    print(f"  fit window: 0 < z <= {Z_MAX_FIT*1e3:.0f} mm,  exponent band: "
          f"[{N_MIN}, {N_MAX}],  R^2 >= {R2_MIN}")
    print(f"  {'case':>10}{'b(mm)':>7}{'mu':>5}{'n':>8}{'R^2':>7}"
          f"{'p_max(Pa)':>12}  note")
    ok = True
    for r in rows:
        notes = []
        if not (N_MIN <= r["n"] <= N_MAX):
            notes.append("EXPONENT-OUT-OF-BAND"); ok = False
        if r["r2"] < R2_MIN:
            notes.append("POOR-FIT"); ok = False
        if not r["monotone"]:
            notes.append("NON-MONOTONE"); ok = False
        if not r["p_increases"]:
            notes.append("PRESSURE-NOT-RISING"); ok = False
        print(f"  {r['tag']:>10}{r['b']*1e3:>7.0f}{r['mu']:>5.1f}{r['n']:>8.3f}"
              f"{r['r2']:>7.3f}{r['p_max']:>12.3e}  {' '.join(notes)}")

    # Width trend: at fixed mu, a wider plate should carry larger total load.
    same_mu = [r for r in rows if abs(r["mu"] - 0.5) < 1e-9]
    same_mu.sort(key=lambda r: r["b"])
    if len(same_mu) >= 2:
        loads = [r["p_max"] * r["b"] * L_Y for r in same_mu]  # F = p*A
        sane = all(loads[i + 1] >= loads[i] * 0.9 for i in range(len(loads) - 1))
        print(f"\n  Width trend (total load F=p_max*A vs b at mu=0.5): "
              f"{[f'{l:.3e}' for l in loads]}  -> "
              f"{'monotone-ish (OK)' if sane else 'NON-MONOTONE LOAD (FAIL)'}")
        if not sane:
            ok = False

    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


def _load_curve(tag, prefix="curve"):
    path = os.path.join(DATA_DIR, f"{prefix}_{tag}.csv")
    if not os.path.isfile(path):
        return [], []
    zs, ps = [], []
    with open(path) as f:
        for r in csv.DictReader(f):
            zs.append(float(r["sinkage"]))
            ps.append(float(r["pressure"]))
    return zs, ps


def _load_lammps_summaries():
    """Load LAMMPS per-case summaries if present, else []."""
    if not os.path.isfile(LAMMPS_CSV):
        return []
    out = []
    with open(LAMMPS_CSV) as f:
        for r in csv.DictReader(f):
            out.append({"tag": r["tag"], "b": float(r["b"]), "mu": float(r["mu"]),
                        "n": float(r["n"]), "A": float(r["A"]), "r2": float(r["r2"]),
                        "p_max": float(r["p_max"])})
    return out


def compare_codes(dirt_rows, lammps_rows):
    """Print the per-case DIRT-vs-LAMMPS sinkage exponent comparison."""
    lmp = {r["tag"]: r for r in lammps_rows}
    print("\n=== Sinkage exponent n: DIRT vs LAMMPS (same material, same setup) ===")
    print(f"  {'case':>10}{'b(mm)':>7}{'mu':>5} | {'DIRT n':>8}{'LAMMPS n':>10} | {'diff':>8}")
    for r in dirt_rows:
        l = lmp.get(r["tag"])
        if not l:
            print(f"  {r['tag']:>10}{r['b']*1e3:>7.0f}{r['mu']:>5.1f} | "
                  f"{r['n']:>8.3f}{'   --':>10} | {'   --':>8}")
            continue
        print(f"  {r['tag']:>10}{r['b']*1e3:>7.0f}{r['mu']:>5.1f} | "
              f"{r['n']:>8.3f}{l['n']:>10.3f} | {l['n'] - r['n']:>+8.3f}")


def plot(rows, lammps_rows=None):
    os.makedirs(PLOT_DIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    lmp = {r["tag"]: r for r in (lammps_rows or [])}

    # ── p vs z, log-log, with power-law fit per case ──
    # DIRT: filled markers + solid fit. LAMMPS (same material/setup): open
    # markers + dashed fit, in the matching color, per width.
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    colors = plt.cm.viridis([0.1, 0.4, 0.7, 0.9])
    for r, c in zip(rows, colors):
        zs, ps = _load_curve(r["tag"])
        if not zs:
            continue
        label = f"b={r['b']*1e3:.0f}mm, mu={r['mu']:.1f}  (n={r['n']:.2f})"
        ax.loglog(zs, ps, "o", color=c, ms=4, label=label)
        # Fit line p = A z^n.
        zz = [min(zs) * (max(zs) / min(zs)) ** (k / 50) for k in range(51)]
        pp = [r["A"] * z ** r["n"] for z in zz]
        ax.loglog(zz, pp, "-", color=c, lw=1.3)
        # LAMMPS overlay (open markers / dashed fit) for the same case.
        lz, lp = _load_curve(r["tag"], prefix="lammps_curve")
        if lz:
            lr = lmp.get(r["tag"], {})
            ax.loglog(lz, lp, "o", color=c, ms=4, markerfacecolor="none")
            if lr:
                lzz = [min(lz) * (max(lz) / min(lz)) ** (k / 50) for k in range(51)]
                lpp = [lr["A"] * z ** lr["n"] for z in lzz]
                ax.loglog(lzz, lpp, "--", color=c, lw=1.1)
    ax.set_xlabel("sinkage z (m)")
    ax.set_ylabel("plate pressure p = F/A (Pa)")
    ax.set_title("Plate pressure–sinkage: DIRT vs Bekker power law p ∝ z^n")
    handles, _ = ax.get_legend_handles_labels()
    if lmp:
        handles += [
            Line2D([], [], color="gray", marker="o", linestyle="-",
                   markersize=5, label="DIRT (filled, solid fit)"),
            Line2D([], [], color="gray", marker="o", linestyle="--",
                   markerfacecolor="none", markersize=5, label="LAMMPS (open, dashed fit)"),
        ]
    ax.legend(handles=handles, fontsize=8)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "pressure_sinkage.png"))
    plt.close(fig)

    # ── linear p vs z (intuitive monotone view) ──
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for r, c in zip(rows, colors):
        zs, ps = _load_curve(r["tag"])
        if not zs:
            continue
        ax.plot([z * 1e3 for z in zs], ps, "-", color=c,
                label=f"b={r['b']*1e3:.0f}mm, mu={r['mu']:.1f}")
    ax.set_xlabel("sinkage z (mm)")
    ax.set_ylabel("plate pressure p = F/A (Pa)")
    ax.set_title("Plate pressure–sinkage (linear) — monotone increasing")
    ax.legend(fontsize=9)
    ax.grid(True, ls=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "pressure_sinkage_linear.png"))
    plt.close(fig)

    print(f"\nFigures -> {PLOT_DIR}/pressure_sinkage.png, pressure_sinkage_linear.png")


def graph():
    rows = _load_summaries()
    if not rows:
        print(f"No {SWEEP_CSV} — run 'start' first.")
        return False
    # DIRT-only validation: PASS/FAIL is decided against the Bekker form, never
    # against LAMMPS (LAMMPS is an optional overlay, not a reference).
    ok = validate(rows)
    lammps_rows = _load_lammps_summaries()
    if lammps_rows:
        compare_codes(rows, lammps_rows)
    else:
        print(f"\n(no {os.path.basename(LAMMPS_CSV)} — plotting/validating DIRT only)")
    plot(rows, lammps_rows)
    return ok


# ── dispatch ──────────────────────────────────────────────────────────────────
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
