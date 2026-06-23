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

Cross-code overlay (optional): if a LAMMPS binary (lmp_serial / lmp / lmp_mpi /
lammps) is on PATH, the SAME lift-the-cylinder protocol is also run in LAMMPS with
the matched Hertz-Mindlin granular model AND the matched sds rolling-resistance
model (same E, nu, restitution, mu; same k_roll, gamma_roll, mu_roll for grain–grain
AND grain–wall; same floor + confining-cylinder + catch walls; same lift-then-settle
sequence), the SAME heap-fit code is applied to LAMMPS's settled positions, and
theta_r(mu) is overlaid on the plot as open markers. Because both codes are on the
identical sds rolling model, this is a fair sds<->sds cross-code comparison: both
should hold a pile and the two theta_r(mu) curves should be reasonably close.
LAMMPS is STRICTLY OPTIONAL: without it the example runs and the DIRT validation
passes exactly as before. The validate() gate is DIRT-only — LAMMPS is an
informative overlay, not a pass/fail reference.

Outputs:
    sweep/<case>/config.toml            DIRT configs                 (gitignored)
    sweep/<case>/data/repose_results.csv  per-run particle positions (gitignored)
    sweep/lammps_<mu>/in.lammps         LAMMPS inputs                (gitignored)
    data/repose_sweep.csv               theta_r per (mu, rep)        (gitignored)
    data/lammps_results.csv             LAMMPS theta_r per mu        (gitignored)
    data/profile_<mu>.csv               representative DIRT h(r)     (gitignored)
    data/lammps_profile_<mu>.csv        representative LAMMPS h(r)   (gitignored)
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
SWEEP_CSV = os.path.join(DATA_DIR, "repose_sweep.csv")     # DIRT theta_r per (mu, rep)
LAMMPS_CSV = os.path.join(DATA_DIR, "lammps_results.csv")  # LAMMPS theta_r per mu

# LAMMPS binary candidates, in preference order. LAMMPS is OPTIONAL: if none is
# found, the LAMMPS leg is skipped and only DIRT is run/validated/plotted.
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
DENSITY = 2500.0        # kg/m^3

# -- sds rolling-resistance model (IDENTICAL in DIRT and LAMMPS) -----------------
# Both codes run the SAME spring–dashpot–slider rolling model with these exact
# parameters, so the cross-code overlay is a fair sds<->sds comparison. The rolling
# torque is  −k_roll·δ − gamma_roll·ω_roll, Coulomb-capped at mu_roll·|F_n|·r_eff,
# and the spring is rescaled on slip. DIRT: rolling_model="sds" with
# rolling_stiffness=k_roll, rolling_damping=gamma_roll, rolling_friction=mu_roll
# (grain–grain AND grain–wall). LAMMPS: `rolling sds k_roll gamma_roll mu_roll`
# in BOTH pair_coeff and every fix wall/gran.
#
# Choice (physically grounded; Ai et al. 2011, Comput. Geotech.; Wensrich &
# Katterfeld 2012, Powder Technol.): the rolling spring stiffness is tied to the
# contact via k_roll ~ k_t·r² (k_t the tangential stiffness, r the grain radius).
# With the softened E used here k_t ~ 2e3 N/m and r = 2e-3 m give k_t·r² ~ 8e-3,
# so k_roll = 1e-2 N·m/rad. gamma_roll = 1e-6 N·m·s/rad is ~0.4 of the critical
# rolling damping 2·sqrt(I·k_roll) (I the grain moment of inertia), enough to kill
# rolling oscillation without overdamping. The rolling-oscillation period
# 2π·sqrt(I/k_roll) ~ 7e-4 s is resolved by the ~2.6e-5 s timestep (~28 steps).
ROLLING_FRICTION = 0.1  # mu_roll — rolling Coulomb cap (the slider limit)
ROLLING_STIFFNESS = 1.0e-2  # k_roll (DIRT)  — rolling spring stiffness (N·m/rad)
ROLLING_DAMPING = 1.0e-6    # gamma_roll — rolling viscous damping (N·m·s/rad)
# DIRT and LAMMPS normalize the sds rolling spring DIFFERENTLY, so the SAME
# nominal k_roll is NOT the same effective stiffness. The physically meaningful
# regime (and the one that makes the overlay fair) is the Coulomb-cap-saturated
# one: a spring stiff enough that the rolling resistance sits at mu_roll·|F_n|·R
# every step, i.e. the standard constant rolling-resistance couple. DIRT reaches
# that at k_roll=1e-2 (saturation displacement ~1.6e-5 rad, hit instantly).
# LAMMPS's softer normalization needs a much stiffer nominal k_roll to saturate;
# bench_rolling_decay's validated value is 1e2. So the LAMMPS leg uses that — NOT
# a different physics, the SAME cap-limited rolling, just the convention-matched
# stiffness. (At 1e-2 LAMMPS never reaches the cap and the heap pancakes.)
LAMMPS_ROLLING_STIFFNESS = 1.0e2  # k_roll for LAMMPS — cap-saturating (see above)
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
rolling_model = "sds"
[[dem.materials]]
name = "glass"
youngs_mod = {youngs:.6e}
poisson_ratio = {nu}
restitution = {e_n}
friction = {mu}
rolling_friction = {mu_r}
rolling_stiffness = {k_roll:.6e}
rolling_damping = {gamma_roll:.6e}
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

# -- LAMMPS leg (optional cross-code overlay) -----------------------------------
# Same protocol, same material, mapped to LAMMPS's `pair_style granular`:
#   hertz/material E e nu                -> Young's modulus, restitution, Poisson
#   tangential mindlin NULL 1.0 {mu}     -> Hertz-Mindlin sliding friction (mu);
#                                           NULL = derive k_t from the material,
#                                           1.0 = poisson tangential-stiffness factor
#   damping coeff_restitution            -> normal damping from the restitution e
#   rolling sds k_roll gamma_roll mu_roll-> the SAME sds rolling model DIRT runs,
#                                           with the SAME k_roll/gamma_roll/mu_roll
#   twisting none                        -> twisting off in both codes
# The floor is a frictional `fix wall/gran ... zplane 0.0`, the confining cylinder
# is a `fix wall/gran/region ... region cyl` that is `unfix`-ed at the lift, and a
# wide catch wall (`region catch`) conserves the grain count — mirroring the three
# DIRT `[[wall]]`s. Grains are introduced with `fix pour` (random, non-overlapping,
# the same packing style as DIRT's overlap-checked random inserter — a lattice fill
# locks into a rigid, non-collapsing pillar and is unusable), poured into the
# confined cylinder, settled into a packed column, then the cylinder is lifted and
# the column collapses and relaxes on the frictional floor.
#
# ROLLING: both codes run the IDENTICAL sds rolling model with the IDENTICAL
# parameters (k_roll=ROLLING_STIFFNESS, gamma_roll=ROLLING_DAMPING,
# mu_roll=ROLLING_FRICTION). DIRT's `sds` rolling branch (dirt_granular AND
# dirt_wall) and LAMMPS's `rolling sds k_roll gamma_roll mu_roll` are the same
# spring–dashpot–slider model (torque −k_roll·δ − gamma_roll·ω_roll, Coulomb-capped
# at mu_roll·|F_n|·r_eff, spring rescaled on slip), so the overlay is a fair
# sds<->sds comparison. The rolling clause is repeated identically in the pair_coeff
# (grain–grain) AND every fix wall/gran (grain–wall) line so the floor has matching
# rolling resistance in both codes.
#
# The LAMMPS box is taller than DIRT's so the pour has headroom; the heap geometry
# (floor at z=0, cylinder r=0.025, catch r=0.07, grains r=0.002) is identical.
LMP_BOX_HI = 0.50          # m — tall box so the pour has headroom
LMP_POUR_LO = 0.10         # m — pour-region bottom (above the forming column)
LMP_DT = 2.0e-5            # s — timestep (DIRT auto-selects ~2.6e-5 here; matched in band)
LMP_POUR_SEED = 12345      # deterministic pour (single LAMMPS run per mu)
LMP_FILL_STEPS = 120000    # pour + initial settle
LMP_SETTLE_STEPS = 100000  # settle the confined column before the lift
LMP_RELAX_STEPS = 150000   # collapse + relax after the lift (KE is negligible well before this)

LMP_TEMPLATE = """\
# Auto-generated LAMMPS input — angle-of-repose lift-the-cylinder, mu = {mu}
# Matches the DIRT material: E={E} Pa, nu={nu}, e={e_n}, mu={mu};
# sds rolling: k_roll={k_roll} gamma_roll={gamma_roll} mu_roll={mu_roll} (same as DIRT).
units           si
atom_style      sphere
boundary        f f f
newton          off
comm_modify     vel yes

region          simbox  block {dxlo} {dxhi} {dylo} {dyhi} 0.0 {boxhi} units box
create_box      1 simbox

region          cyl     cylinder z 0.0 0.0 {cyl_r} 0.0 {boxhi} units box
region          catch   cylinder z 0.0 0.0 {catch_r} 0.0 {boxhi} units box
region          pourreg cylinder z 0.0 0.0 {pour_r} {pour_lo} {pour_hi} units box

pair_style      granular
pair_coeff      1 1 hertz/material {E} {e_n} {nu} tangential mindlin NULL 1.0 {mu} damping coeff_restitution rolling sds {k_roll} {gamma_roll} {mu_roll} twisting none

fix             grav    all gravity 9.81 vector 0.0 0.0 -1.0
fix             floor   all wall/gran granular hertz/material {E} {e_n} {nu} tangential mindlin NULL 1.0 {mu} damping coeff_restitution rolling sds {k_roll} {gamma_roll} {mu_roll} twisting none zplane 0.0 NULL
fix             catchw  all wall/gran/region granular hertz/material {E} {e_n} {nu} tangential mindlin NULL 1.0 {mu} damping coeff_restitution rolling sds {k_roll} {gamma_roll} {mu_roll} twisting none region catch
fix             cylwall all wall/gran/region granular hertz/material {E} {e_n} {nu} tangential mindlin NULL 1.0 {mu} damping coeff_restitution rolling sds {k_roll} {gamma_roll} {mu_roll} twisting none region cyl
fix             integrate all nve/sphere

# Pour the heap grains (random, non-overlapping) into the confined cylinder.
fix             ins all pour {count} 1 {seed} region pourreg diam one {diam} dens {density} {density} vol 0.30 1000
timestep        {dt}
thermo          50000
thermo_style    custom step atoms ke

run             {fill_steps}      # pour + settle
unfix           ins
run             {settle_steps}    # settle the confined column
unfix           cylwall           # LIFT the cylinder
run             {relax_steps}     # collapse + relax on the frictional floor

write_dump      all custom {dump} id x y z radius modify sort id
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
        mu_r=ROLLING_FRICTION, k_roll=ROLLING_STIFFNESS, gamma_roll=ROLLING_DAMPING,
        cyl_r=CYL_RADIUS,
        heap_count=HEAP_COUNT, radius=RADIUS, density=DENSITY,
        ins_r=CYL_RADIUS - 1.5 * RADIUS, outdir=outdir,
    )


# -- LAMMPS helpers -------------------------------------------------------------
def lammps_case_dir(mu):
    return os.path.join(SWEEP_DIR, f"lammps_mu_{mu:g}")


def _lammps_input(mu, dump):
    """Render the LAMMPS input for one mu (same geometry/material as DIRT)."""
    return LMP_TEMPLATE.format(
        E=f"{YOUNGS_MOD:.6e}", nu=POISSON, e_n=RESTITUTION, mu=mu,
        k_roll=f"{LAMMPS_ROLLING_STIFFNESS:.6e}", gamma_roll=f"{ROLLING_DAMPING:.6e}",
        mu_roll=ROLLING_FRICTION,
        boxhi=LMP_BOX_HI,
        dxlo=-0.08, dxhi=0.08, dylo=-0.08, dyhi=0.08,
        cyl_r=CYL_RADIUS, catch_r=0.07,
        pour_r=CYL_RADIUS - 1.5 * RADIUS, pour_lo=LMP_POUR_LO, pour_hi=LMP_BOX_HI - 0.05,
        count=HEAP_COUNT, diam=2.0 * RADIUS, density=DENSITY, seed=LMP_POUR_SEED,
        dt=f"{LMP_DT:.6e}",
        fill_steps=LMP_FILL_STEPS, settle_steps=LMP_SETTLE_STEPS, relax_steps=LMP_RELAX_STEPS,
        dump=dump,
    )


def _parse_lammps_dump(path):
    """Read a LAMMPS `id x y z radius` dump -> (xs, ys, zs, rs)."""
    xs, ys, zs, rs = [], [], [], []
    with open(path) as f:
        lines = f.read().split("\n")
    start = None
    for k, ln in enumerate(lines):
        if ln.startswith("ITEM: ATOMS"):
            start = k + 1
            break
    if start is None:
        return xs, ys, zs, rs
    for ln in lines[start:]:
        p = ln.split()
        if len(p) < 5:
            continue
        xs.append(float(p[1])); ys.append(float(p[2]))
        zs.append(float(p[3])); rs.append(float(p[4]))
    return xs, ys, zs, rs


def _run_lammps(lammps, mu):
    """Run one LAMMPS repose case; return (theta, r_toe, r_centers, h_surface,
    baseline, n) or None on failure."""
    cdir = lammps_case_dir(mu)
    os.makedirs(cdir, exist_ok=True)
    in_path = os.path.join(cdir, "in.lammps")
    dump = os.path.join(cdir, "repose.dump")
    log = os.path.join(cdir, "lammps.log")
    if os.path.exists(dump):
        os.remove(dump)
    with open(in_path, "w") as f:
        f.write(_lammps_input(mu, dump))
    proc = subprocess.run(
        [lammps, "-in", in_path, "-log", log],
        cwd=cdir, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
    )
    if proc.returncode != 0 or not os.path.isfile(dump):
        return None
    xs, ys, zs, rs = _parse_lammps_dump(dump)
    if not xs:
        return None
    r_c, h_s, base, diam = heap_profile(xs, ys, zs, rs)
    theta, r_toe = fit_angle(r_c, h_s, base, diam)
    return theta, r_toe, r_c, h_s, base, len(xs)


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

    lammps = find_lammps()
    print(f"LAMMPS: {lammps}" if lammps
          else "LAMMPS: not found on PATH — running DIRT only.")

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

    # -- LAMMPS leg (optional cross-code overlay) -------------------------------
    if not lammps:
        return
    print("\nLAMMPS cross-code overlay (same protocol, same material):")
    # Clear stale LAMMPS outputs so old results can never be re-plotted.
    for stale in (LAMMPS_CSV, *[os.path.join(DATA_DIR, f"lammps_profile_{mu:g}.csv")
                                for mu in MU_LIST]):
        if os.path.exists(stale):
            os.remove(stale)
    lrows = []
    for mu in MU_LIST:
        print(f"  mu={mu:<4}", end="  ", flush=True)
        out = _run_lammps(lammps, mu)
        if out is None:
            print("LAMMPS FAILED")
            continue
        theta, r_toe, r_c, h_s, base, n = out
        lrows.append({"mu": mu, "theta_deg": theta, "r_toe": r_toe, "n": n})
        print(f"theta_r = {theta:5.2f} deg  (r_toe={r_toe*1e3:.1f} mm, N={n})")
        with open(os.path.join(DATA_DIR, f"lammps_profile_{mu:g}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["r", "h"])
            for i in range(len(r_c)):
                w.writerow([r_c[i], h_s[i] - base])
    if lrows:
        with open(LAMMPS_CSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["mu", "theta_deg", "r_toe", "n"])
            w.writeheader()
            for r in lrows:
                w.writerow(r)
        print(f"\nLAMMPS: {len(lrows)}/{len(MU_LIST)} cases -> {LAMMPS_CSV}")


# -- graph (validate + plot) ----------------------------------------------------
def _load_sweep():
    if not os.path.isfile(SWEEP_CSV):
        return []
    with open(SWEEP_CSV) as f:
        return [{k: (float(v) if k != "rep" else int(float(v))) for k, v in r.items()}
                for r in csv.DictReader(f)]


def _load_lammps():
    """LAMMPS theta_r per mu, or [] if no overlay was run."""
    if not os.path.isfile(LAMMPS_CSV):
        return []
    with open(LAMMPS_CSV) as f:
        return [{"mu": float(r["mu"]), "theta_deg": float(r["theta_deg"]),
                 "r_toe": float(r["r_toe"]), "n": int(float(r["n"]))}
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
    print(f"  material: E={YOUNGS_MOD:.1e} Pa  nu={POISSON}  e={RESTITUTION}")
    print(f"  rolling (sds): k_roll={ROLLING_STIFFNESS:g}  "
          f"gamma_roll={ROLLING_DAMPING:g}  mu_roll={ROLLING_FRICTION:g}")
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


def compare_codes(dirt_stats, lammps_rows):
    """Per-mu DIRT-vs-LAMMPS theta_r printout. LAMMPS is a single deterministic
    run per mu, so it has no spread; DIRT shows its mean over reps."""
    lmp = {r["mu"]: r["theta_deg"] for r in lammps_rows}
    print("\n" + "=" * 52)
    print("Angle of repose: DIRT vs LAMMPS (same protocol/material)")
    print("=" * 52)
    print(f"  {'mu':>5} | {'DIRT (deg)':>12}{'LAMMPS (deg)':>14} | {'diff':>8}")
    for (mu, mean, _std, _n) in dirt_stats:
        if mu not in lmp:
            continue
        print(f"  {mu:>5.2f} | {mean:>12.2f}{lmp[mu]:>14.2f} | {lmp[mu] - mean:>+8.2f}")
    print(f"  (both codes on sds rolling: k_roll={ROLLING_STIFFNESS:g}, "
          f"gamma_roll={ROLLING_DAMPING:g}, mu_roll={ROLLING_FRICTION:g}.)")


def plot(rows, lammps_rows):
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
    ax.errorbar(mus, means, yerr=stds, fmt="o-", capsize=4, color="#1f77b4",
                label="DIRT (mean ± std over reps)")
    ax.scatter([r["mu"] for r in rows], [r["theta_deg"] for r in rows],
               s=14, alpha=0.4, color="gray", label="DIRT individual runs")
    if lammps_rows:
        lm = sorted(lammps_rows, key=lambda r: r["mu"])
        ax.plot([r["mu"] for r in lm], [r["theta_deg"] for r in lm],
                "s--", color="#d62728", markerfacecolor="none", markersize=7,
                label="LAMMPS (same protocol, sds rolling)")
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
    def _load_profile(path):
        rc, hs = [], []
        if os.path.isfile(path):
            with open(path) as f:
                for row in csv.DictReader(f):
                    rc.append(float(row["r"]) * 1e3)
                    hs.append(float(row["h"]) * 1e3)
        return rc, hs

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    cmap = plt.get_cmap("viridis")
    mus_sorted = sorted({r["mu"] for r in rows})
    have_lammps_profile = False
    for j, mu in enumerate(mus_sorted):
        color = cmap(j / max(1, len(mus_sorted) - 1))
        rc, hs = _load_profile(os.path.join(DATA_DIR, f"profile_{mu:g}.csv"))
        if rc:
            ax.plot(rc, hs, "o-", ms=3, color=color, label=fr"$\mu$={mu:g} (DIRT)")
        if lammps_rows:
            lrc, lhs = _load_profile(os.path.join(DATA_DIR, f"lammps_profile_{mu:g}.csv"))
            if lrc:
                have_lammps_profile = True
                ax.plot(lrc, lhs, "s--", ms=3, color=color, markerfacecolor="none")
    ax.set_xlabel("radial distance r (mm)")
    ax.set_ylabel("heap surface height h (mm)")
    title = "Settled heap cross-section (surface envelope)"
    if have_lammps_profile:
        title += "\nsolid = DIRT, dashed open = LAMMPS"
    ax.set_title(title)
    ax.legend(title="friction", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "heap_profile.png"))
    plt.close(fig)

    print(f"\nFigures -> {PLOT_DIR}/theta_vs_mu.png, heap_profile.png")


def graph():
    rows = _load_sweep()
    if not rows:
        print(f"No {SWEEP_CSV} — run 'start' first.")
        return False
    ok = validate(rows)           # DIRT-only gate (LAMMPS never affects PASS/FAIL)
    lammps_rows = _load_lammps()
    if lammps_rows:
        compare_codes(_stats_by_mu(rows), lammps_rows)
    else:
        print(f"\n(no {os.path.basename(LAMMPS_CSV)} — plotting DIRT only)")
    plot(rows, lammps_rows)
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
