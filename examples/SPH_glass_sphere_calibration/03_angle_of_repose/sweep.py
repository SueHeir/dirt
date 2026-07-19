#!/usr/bin/env python3
"""
SPH glass-sphere calibration — angle-of-repose ROLLING-friction (mu_r) gate.

This is the CALIBRATION GATE for the glass-sphere suite. It forms a static
granular heap with the proven "lift the cylinder" protocol and measures its
angle of repose theta_r as a function of the ROLLING friction mu_r, at a FIXED
sliding friction mu_p = 0.16 (the measured glass value), to pin the rolling
friction that reproduces the measured glass repose band [22,26] deg.

The deposit is formed at e=0.4, not the canonical e=0.926. This harness does not
claim that a value fitted under that different formation protocol transfers to
canonical glass; a protocol-matched external validation is required for that.

Validation (the gate):

    1. theta_r increases MONOTONICALLY with mu_r,
    2. the formation-protocol result is compared with the measured glass band
       [22,26] deg, but cannot by itself close canonical glass calibration: the
       formation restitution differs from the canonical material,
    3. results are REPRODUCIBLE: the run-to-run spread (over independent packs)
       is small but NONZERO. The per-case seeds are derived from a recorded base
       seed or read back from a seed manifest, so a campaign can be regenerated.

The heap sits directly on a real frictional plane wall (z = 0, normal +z):
dirt_wall applies Mindlin sliding (tangential) friction on plane walls using the
material's friction coefficient (μ) via friction_ij. That base friction keeps
the bottom layer from sliding out, so the pile holds a slope — no frozen
particle bed is needed. See README "Assumptions".

Commands (from anywhere):
    python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py generate --base-seed 20260706
    python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py generate --seed-manifest data/seed_manifest.csv
    python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py seed-check
    python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py estimator-check
    python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py start
    python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py graph
    python3 examples/SPH_glass_sphere_calibration/03_angle_of_repose/sweep.py   # all

Each (mu_r) case is run REPS times with independent random packs. The inserter is
deterministic given its seed; `generate` assigns a distinct seed to every
(mu_r, rep) case from a base seed, records those seeds in a manifest, and can read
that manifest back to reproduce the exact same configs.

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
import argparse
import hashlib
import json
import math
import subprocess
import tempfile
import reference_audit

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# This example lives THREE levels under the repo root:
#   examples/SPH_glass_sphere_calibration/03_angle_of_repose/ -> repo root
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
EXAMPLE = "sphcal_angle_of_repose"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
REPLAY_DIR = os.path.join(SWEEP_DIR, "zero_history_replay")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
SWEEP_CSV = os.path.join(DATA_DIR, "repose_sweep.csv")     # DIRT theta_r per (mu, rep)
LAMMPS_CSV = os.path.join(DATA_DIR, "lammps_results.csv")  # LAMMPS theta_r per mu
SEED_MANIFEST = os.path.join(DATA_DIR, "seed_manifest.csv")
CAMPAIGN_LEDGER = os.path.join(DATA_DIR, "campaign_ledger.json")
PROTOCOL_REFERENCE = os.path.join(SCRIPT_DIR, "external_records", "protocol_matched_glass.json")

# Cross-code evidence is deliberately opt-in (`external` below): it is too slow
# for the bounded local sweep, but must not be silently disabled.  A matching
# LAMMPS executable is an adversarial implementation check, not an experimental
# reference and therefore never participates in the glass-band PASS decision.
LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# -- Sweep parameters -----------------------------------------------------------
# CALIBRATION GATE: we sweep ROLLING friction mu_r at a FIXED sliding friction
# mu_p = 0.16 (the measured glass value). The STATIC angle of repose is governed
# by FRICTION, not restitution, so the heap is formed with the proven dissipative
# restitution e=0.4. This is a formation-specific measurement, not a transfer
# calibration for canonical e=0.926 glass. In the lift-the-cylinder protocol the
# heap forms by a column
# COLLAPSE on the frictional floor: with sliding friction fixed, raising mu_r
# arrests the surface grains' rolling and the deposit relaxes into a steeper
# cone, so theta_r rises monotonically with mu_r.
MU_P = 0.16                                  # FIXED sliding (Coulomb) friction (measured glass)
MU_R_LIST = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]  # rolling friction sweep (the swept variable)
REPS = 2                          # independent packs per mu_r (distinct recorded seeds)
DEFAULT_BASE_SEED = 20260706       # deterministic campaign seed used by generate unless overridden

# Measured glass angle-of-repose band the calibration targets.
GLASS_BAND_LO_DEG = 22.0
GLASS_BAND_HI_DEG = 26.0

# Material / geometry (mirror config.toml; mu_r is overridden per case). The fixed
# sliding mu_p governs both particle–particle and particle–floor-wall sliding; the
# swept mu_r is the rolling Coulomb cap that sets the pile's static angle.
YOUNGS_MOD = 1.0e7      # Pa (softened: larger stable dt, standard DEM practice)
POISSON = 0.25
RESTITUTION = 0.4       # FORMATION AID — dissipative so the heap settles (see above)
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
# LAMMPS SDS is a *pseudo-force* law: its history is a length and its coefficients
# have units N/m and N s/m (not N m/rad and N m s/rad).  The original campaign
# values were torque-law coefficients.  After DIRT was corrected to LAMMPS's
# pseudo-force form, retaining those bare numbers weakened the spring and dashpot
# by R².  Convert the declared torque-law values at R=2 mm instead of silently
# changing the intended resistance: k_F = k_tau/R² = 0.01/(0.002²)=2500 N/m and
# gamma_F = gamma_tau/R² = 1e-6/(0.002²)=0.25 N s/m.  This is a dimensional
# correction, not a fit to the repose target.  k_F is also of the same order as
# the softened Hertz-Mindlin tangential stiffness (~2e3 N/m).
#
# mu_roll (the force Coulomb cap / slider limit) remains the swept variable.
ROLLING_STIFFNESS = 2.5e3   # k_roll pseudo-force stiffness (N/m)
ROLLING_DAMPING = 2.5e-1    # gamma_roll pseudo-force damping (N s/m)
# `GranSubModRollingSDS` in the checked LAMMPS source uses the same
# length-valued history and force law as DIRT: F_r = -k*xi - gamma*v_r with
# v_r = -R*omega_rel x n (the documented LAMMPS sign convention).
# Therefore a comparison must use the *same numerical coefficient*.  Replacing
# it with a cap-saturating value is a different model parameter, not a mapping.
LAMMPS_ROLLING_STIFFNESS = ROLLING_STIFFNESS
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
# The formation-study checks demand a monotonic, reproducible response and report
# whether it reaches the measured glass band.  They are intentionally insufficient
# to approve a canonical-material edit because the formation restitution differs.
# At mu_r=0 (no rolling resistance) the surface grains roll freely and
# the collapse deposit is a shallow cone, so the low-mu_r case reads small.
GLASS_BAND_LO_DEG = GLASS_BAND_LO_DEG  # measured glass band lower bound (22 deg)
GLASS_BAND_HI_DEG = GLASS_BAND_HI_DEG  # measured glass band upper bound (26 deg)
SPREAD_MAX_DEG = 5.0    # max allowed run-to-run std dev of theta_r at a given mu_r
MONOTONIC_SLACK_DEG = 2.5  # mean theta_r may dip by at most this between mu_r steps
REST_MAX_SPEED_M_S = 2.0e-3  # must match main.rs; common pre/post-release rest bound

# A source record must match this frozen, pre-existing protocol.  It is not
# derived from solver output and it cannot change the retained 22--26° band.
EXTERNAL_PROTOCOL = {
    "sphere_radius_m": RADIUS,
    "density_kg_m3": DENSITY,
    "restitution": RESTITUTION,
    "sliding_friction": MU_P,
    "container": "lifted cylinder, radius 0.025 m, frictional plane floor",
    "deposition": "settled 1200 monodisperse spheres in cylinder",
    "release": "remove confining cylinder after fill-rest witness",
    "angle_definition": "linear fit to settled radial upper-envelope flank",
}

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
# LAMMPS documents displacement-history `mindlin_rescale` as a historical law
# that double-rescales unloading. Its force-history counterpart is declared
# explicitly in both DIRT and the independent replay.
tangential_model = "mindlin_rescale/force"
rolling_model = "sds"
# Keep the standard DIRT/LAMMPS `limit_damping` contact boundary explicit: a
# separating, non-cohesive contact cannot acquire a tensile normal force solely
# from viscous damping. The LAMMPS replay below carries the same keyword for
# both particle and wall contacts.
limit_damping = true
[[dem.materials]]
name = "glass"
youngs_mod = {youngs:.6e}
poisson_ratio = {nu}
restitution = {e_n}
friction = {mu_p}
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
# Per-case insertion seed. The inserter RNG is deterministic given a seed, so a
# distinct seed per rep yields an INDEPENDENT random pack while the seed manifest
# makes the campaign exactly regenerable.
seed = {seed}
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
#   tangential mindlin_rescale/force NULL 1.0 {mu}
#                                        -> unloading-rescaled Hertz-Mindlin
#                                           sliding friction (mu);
#                                           NULL = derive k_t from the material,
#                                           1.0 = poisson tangential-stiffness factor
#   damping coeff_restitution limit_damping -> normal damping from restitution e,
#                                               with DIRT's declared repulsive-only
#                                               normal-force cutoff
#   rolling sds k_roll gamma_roll mu_roll-> the SAME sds rolling model DIRT runs,
#                                           with the SAME k_roll/gamma_roll/mu_roll
#   twisting none                        -> twisting off in both codes
# The floor is a frictional `fix wall/gran ... zplane 0.0`, the confining cylinder
# is a `fix wall/gran/region ... region cyl` that is `unfix`-ed at the lift, and a
# wide catch wall (`region catch`) conserves the grain count — mirroring the three
# DIRT `[[wall]]`s. LAMMPS does not create a separate packing: it reads DIRT's
# solver-authored settled pre-lift snapshot, then removes the same cylinder and
# relaxes the same positions on the frictional floor. This deliberately narrows
# the independent comparison to the post-lift constitutive/integration path.
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
LMP_DT = 2.0e-5            # s — timestep (DIRT auto-selects ~2.6e-5 here; matched in band)
# The independent replay must use the same post-lift stopping observable as
# DIRT. A fixed duration cannot establish that two integrators reached comparable
# states, so it is unsuitable for a contact-model comparison.
LMP_REST_SPEED_M_S = 1.0e-2
LMP_RELAX_POLL_STEPS = 2000
LMP_RELAX_MAX_STEPS = 150000

LMP_REPLAY_TEMPLATE = """\
# Auto-generated LAMMPS input — replay DIRT's settled pre-lift state, mu = {mu}
# Matches the DIRT material: E={E} Pa, nu={nu}, e={e_n}, mu={mu};
# sds rolling: k_roll={k_roll} gamma_roll={gamma_roll} mu_roll={mu_roll} (same as DIRT).
units           si
atom_style      sphere
boundary        f f f
newton          off
read_data       {datafile}
comm_modify     vel yes

region          cyl     cylinder z 0.0 0.0 {cyl_r} 0.0 {boxhi} units box
region          catch   cylinder z 0.0 0.0 {catch_r} 0.0 {boxhi} units box

pair_style      granular
pair_coeff      1 1 hertz/material {E} {e_n} {nu} tangential mindlin_rescale/force NULL 1.0 {mu} damping coeff_restitution rolling sds {k_roll} {gamma_roll} {mu_roll} twisting none limit_damping

fix             grav    all gravity 9.81 vector 0.0 0.0 -1.0
fix             floor   all wall/gran granular hertz/material {E} {e_n} {nu} tangential mindlin_rescale/force NULL 1.0 {mu} damping coeff_restitution rolling sds {k_roll} {gamma_roll} {mu_roll} twisting none limit_damping zplane 0.0 NULL
fix             catchw  all wall/gran/region granular hertz/material {E} {e_n} {nu} tangential mindlin_rescale/force NULL 1.0 {mu} damping coeff_restitution rolling sds {k_roll} {gamma_roll} {mu_roll} twisting none limit_damping region catch
fix             cylwall all wall/gran/region granular hertz/material {E} {e_n} {nu} tangential mindlin_rescale/force NULL 1.0 {mu} damping coeff_restitution rolling sds {k_roll} {gamma_roll} {mu_roll} twisting none limit_damping region cyl
fix             integrate all nve/sphere

# Physical qualification, not a parity tolerance: emit a geometry only after
# LAMMPS independently reaches DIRT's declared global heap-rest speed.
variable        atom_speed atom sqrt(vx*vx+vy*vy+vz*vz)
compute         vmax all reduce max v_atom_speed
timestep        {dt}
thermo          50000
thermo_style    custom step atoms ke c_vmax

unfix           cylwall           # LIFT the cylinder
variable        rest_loop loop {rest_polls}
label           relax_until_rest
run             {rest_poll_steps}
variable        rest_speed equal c_vmax
if              "${{rest_speed}} < {rest_speed_limit}" then "jump SELF replay_settled"
next            rest_loop
jump            SELF relax_until_rest
print           "ERROR: LAMMPS replay did not meet the DIRT heap-rest speed criterion"
quit            2
label           replay_settled
print           "LAMMPS_REPLAY_REST step=$(step) vmax=$(v_rest_speed)"

write_dump      all custom {dump} id x y z radius vx vy vz modify sort id
"""


# -- helpers --------------------------------------------------------------------
def case_tag(mu_r, rep):
    return f"mur_{mu_r:g}_rep{rep}"


def case_dir(mu_r, rep):
    return os.path.join(SWEEP_DIR, case_tag(mu_r, rep))


def _case_dir_in(root, mu_r, rep):
    return os.path.join(root, case_tag(mu_r, rep))


def _dirt_config(mu_r, outdir, seed):
    return TOML_TEMPLATE.format(
        gz=GZ, youngs=YOUNGS_MOD, nu=POISSON, e_n=RESTITUTION, mu_p=MU_P,
        mu_r=mu_r, k_roll=ROLLING_STIFFNESS, gamma_roll=ROLLING_DAMPING,
        cyl_r=CYL_RADIUS,
        heap_count=HEAP_COUNT, radius=RADIUS, density=DENSITY,
        ins_r=CYL_RADIUS - 1.5 * RADIUS, outdir=outdir, seed=seed,
    )


def _replay_config(mu_r, outdir, seed):
    """A separate zero-history release for the LAMMPS implementation check.

    The calibration campaign retains formation histories. LAMMPS `read_data`
    cannot receive DIRT's contact-history records, so this diagnostic clears
    DIRT pair and wall histories at lift instead of pretending they transfer.
    """
    return _dirt_config(mu_r, outdir, seed) + (
        "\n[repose_replay]\nzero_contact_history_at_lift = true\n")


def replay_generate(argv=None):
    """Generate one isolated zero-history DIRT release for cross-code diagnosis."""
    parser = argparse.ArgumentParser(prog="sweep.py replay-generate")
    parser.add_argument("--mu-r", type=float, default=MU_R_LIST[-1], choices=MU_R_LIST)
    parser.add_argument("--rep", type=int, default=0, choices=range(REPS))
    args = parser.parse_args(argv or [])
    seed = _derive_seed_manifest(DEFAULT_BASE_SEED)[_case_key(args.mu_r, args.rep)]
    out = os.path.join(REPLAY_DIR, case_tag(args.mu_r, args.rep))
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "config.toml"), "w", encoding="utf-8") as handle:
        handle.write(_replay_config(args.mu_r, out, seed))
    print("Generated zero-history replay config -> " + os.path.join(out, "config.toml"))
    return out


def _lammps_binary():
    """Return a runnable independent granular solver, if one is installed."""
    for candidate in LAMMPS_BINS:
        try:
            probe = subprocess.run([candidate, "-h"], stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL, timeout=20)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
        if probe.returncode == 0:
            return candidate
    return None


def _load_lammps_dump(path, expected_particles):
    """Read exactly the final complete `id x y z radius vx vy vz` LAMMPS frame.

    A dump may contain several frames when a run is resumed or a user changes
    the output cadence.  Pooling them would fabricate a taller, denser heap.
    Keep the last complete frame only, require its declared atom count, and
    reject non-finite coordinates before applying the shared heap fit.
    """
    frames = []
    with open(path, encoding="utf-8") as handle:
        lines = iter(handle)
        for line in lines:
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            try:
                timestep = int(next(lines).strip())
                if not next(lines).startswith("ITEM: NUMBER OF ATOMS"):
                    raise ValueError("missing NUMBER OF ATOMS section")
                count = int(next(lines).strip())
                if not next(lines).startswith("ITEM: BOX BOUNDS"):
                    raise ValueError("missing BOX BOUNDS section")
                next(lines); next(lines); next(lines)
                atom_header = next(lines)
            except StopIteration as exc:
                raise ValueError("truncated LAMMPS dump frame") from exc
            if not atom_header.startswith("ITEM: ATOMS"):
                raise ValueError("missing ATOMS section")
            index = {name: i for i, name in enumerate(atom_header.split()[2:])}
            required = ("x", "y", "z", "radius", "vx", "vy", "vz")
            if any(name not in index for name in required):
                raise ValueError("LAMMPS dump lacks required x/y/z/radius columns")
            rows = []
            try:
                for _ in range(count):
                    values = next(lines).split()
                    point = tuple(float(values[index[name]]) for name in required)
                    if not all(math.isfinite(value) for value in point):
                        raise ValueError("non-finite coordinate in LAMMPS dump")
                    rows.append(point)
            except (StopIteration, IndexError, ValueError) as exc:
                raise ValueError("malformed LAMMPS atom frame") from exc
            frames.append((timestep, rows))
    if not frames:
        raise ValueError("LAMMPS dump has no complete frame")
    timestep, rows = frames[-1]
    if len(rows) != expected_particles:
        raise ValueError(
            f"LAMMPS final frame has {len(rows)} particles; expected {expected_particles}")
    return tuple(zip(*rows)), timestep, len(frames)


def _write_lammps_replay_data(source_csv, destination):
    """Translate DIRT's solver-authored pre-lift snapshot to LAMMPS data.

    Transfer position, radius, linear velocity, and angular velocity.  A
    fill-rest threshold bounds motion but does not make it identically zero;
    discarding it silently changes the collapse initial-value problem.  LAMMPS
    `atom_style sphere` accepts angular velocity in the six-field Velocities
    section, so the adapter can preserve the complete DIRT-authored state.
    """
    required = ("x", "y", "z", "radius", "vx", "vy", "vz",
                "omega_x", "omega_y", "omega_z")
    with open(source_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or any(key not in reader.fieldnames for key in required):
            raise ValueError("pre-lift snapshot lacks complete x/y/z/radius/velocity/omega state")
        try:
            rows = [tuple(float(row[key]) for key in required) for row in reader]
        except (TypeError, ValueError) as exc:
            raise ValueError("pre-lift snapshot has non-numeric state") from exc
    if len(rows) != HEAP_COUNT:
        raise ValueError(f"pre-lift snapshot has {len(rows)} particles; expected {HEAP_COUNT}")
    if not all(math.isfinite(value) for row in rows for value in row):
        raise ValueError("pre-lift snapshot has non-finite state")
    if any(row[3] <= 0.0 for row in rows):
        raise ValueError("pre-lift snapshot has non-positive radius")
    with open(destination, "w", encoding="utf-8") as f:
        f.write("DIRT angle-of-repose pre-lift replay\n\n")
        f.write(f"{len(rows)} atoms\n1 atom types\n\n")
        f.write("-0.08 0.08 xlo xhi\n-0.08 0.08 ylo yhi\n0.0 0.16 zlo zhi\n\n")
        f.write("Atoms # sphere\n\n")
        for i, (x, y, z, radius, _vx, _vy, _vz, _wx, _wy, _wz) in enumerate(rows, 1):
            f.write(f"{i} 1 {2.0 * radius:.9e} {DENSITY:.9e} {x:.9e} {y:.9e} {z:.9e}\n")
        f.write("\nVelocities\n\n")
        for i, (_x, _y, _z, _radius, vx, vy, vz, wx, wy, wz) in enumerate(rows, 1):
            f.write(f"{i} {vx:.9e} {vy:.9e} {vz:.9e} {wx:.9e} {wy:.9e} {wz:.9e}\n")
    return len(rows)


def _lammps_replay_state_error(source_csv, datafile):
    """Independently check complete mechanical-state transfer to LAMMPS data.

    The receipt hashes both files, but hashes alone cannot expose a systematic
    exporter omission when both artifacts were freshly generated.  This parser
    checks the solver-authored CSV against LAMMPS's distinct data syntax before
    the independent calculation starts.
    """
    required = ("x", "y", "z", "radius", "vx", "vy", "vz",
                "omega_x", "omega_y", "omega_z")
    with open(source_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or any(key not in reader.fieldnames for key in required):
            raise ValueError("pre-lift CSV lacks complete mechanical state")
        source = [tuple(float(row[key]) for key in required) for row in reader]

    atoms, velocities, section = {}, {}, None
    with open(datafile, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line == "Atoms # sphere":
                section = "atoms"
                continue
            if line == "Velocities":
                section = "velocities"
                continue
            if not line:
                continue
            fields = line.split()
            if section == "atoms" and len(fields) == 7:
                ident = int(fields[0])
                # id type diameter density x y z
                atoms[ident] = (float(fields[4]), float(fields[5]), float(fields[6]),
                                0.5 * float(fields[2]))
            elif section == "velocities" and len(fields) == 7:
                velocities[int(fields[0])] = tuple(float(value) for value in fields[1:])
    if len(source) != len(atoms) or len(atoms) != len(velocities):
        raise ValueError("LAMMPS data lacks one complete state per source particle")
    max_error = 0.0
    for ident, row in enumerate(source, 1):
        if ident not in atoms or ident not in velocities:
            raise ValueError(f"LAMMPS data is missing particle id {ident}")
        x, y, z, radius = atoms[ident]
        got = (x, y, z, radius, *velocities[ident])
        max_error = max(max_error, *(abs(a - b) for a, b in zip(row, got)))
    return max_error


def external(argv=None):
    """Run a LAMMPS replay from DIRT's exact qualified pre-lift geometry.

    It is intentionally a *comparison receipt*, not an acceptance criterion:
    the source positions come from DIRT's qualified pre-lift snapshot, so the
    receipt can test the post-lift contact/integration path without conflating
    different inserters. It remains a software comparison, not an experiment.
    No independently justified DIRT--LAMMPS agreement criterion has been
    declared for this protocol, so a completed replay is a diagnostic receipt,
    never a passing cross-code validation.
    """
    parser = argparse.ArgumentParser(prog="sweep.py external")
    parser.add_argument("--mu-r", type=float, default=MU_R_LIST[-1], choices=MU_R_LIST)
    parser.add_argument("--rep", type=int, default=0, choices=range(REPS))
    args = parser.parse_args(argv or [])
    binary = _lammps_binary()
    if binary is None:
        print("SKIPPED: no LAMMPS executable found; no external comparison claimed.")
        return False
    cdir = os.path.join(REPLAY_DIR, case_tag(args.mu_r, args.rep))
    prelift = os.path.join(cdir, "data", "repose_prelift.csv")
    final = os.path.join(cdir, "data", "repose_results.csv")
    qualification = os.path.join(cdir, "data", "repose_qualification.json")
    # A pre-lift CSV is not independent evidence merely because it is named
    # ``prelift``.  Require the same solver-written final-rest witness used by
    # the DIRT campaign, then bind *both* snapshots into this receipt.  This
    # prevents a hand-edited or stale source state from being replayed as a
    # supposedly qualified cross-code comparison.
    if not all(os.path.isfile(path) for path in (prelift, final, qualification)):
        print("FAIL: external replay requires DIRT pre-lift, final, and qualification artifacts; run the corresponding DIRT case first.")
        return False
    try:
        with open(qualification, encoding="utf-8") as handle:
            witness = json.load(handle)
        witness_ok = (
            witness["schema"] == 2
            and witness["history_at_lift"] == "cleared"
            and witness["fill_step"] == witness["lift_step"] > 0
            and witness["heap_step"] >= witness["lift_step"] + 2000
            and 0.0 <= float(witness["fill_vmax_m_s"]) < 2e-3
            and 0.0 <= float(witness["heap_vmax_m_s"]) < 1e-2
            and int(witness["particle_count"]) == HEAP_COUNT
        )
        dirt_x, dirt_y, dirt_z, dirt_r = _load_positions(final)
        if len(dirt_x) != HEAP_COUNT:
            witness_ok = False
    except (OSError, ValueError, KeyError, TypeError):
        witness_ok = False
    if not witness_ok:
        print("FAIL: external replay requires the declared zero-history DIRT release plus population/rest qualification; run replay-generate and the generated DIRT case first.")
        return False
    out = os.path.join(SWEEP_DIR, f"lammps_replay_mur_{args.mu_r:g}_rep{args.rep}")
    os.makedirs(out, exist_ok=True)
    dump = os.path.join(out, "repose.dump")
    datafile = os.path.join(out, "prelift.data")
    try:
        replay_count = _write_lammps_replay_data(prelift, datafile)
        state_error = _lammps_replay_state_error(prelift, datafile)
    except ValueError as exc:
        print(f"FAIL: ineligible DIRT pre-lift snapshot: {exc}")
        return False
    # The adapter writes decimal values with nine digits after the decimal
    # exponent.  This is a serialization bound, not a DIRT--LAMMPS physics
    # tolerance: any larger discrepancy means the replay changed its initial
    # value problem and must not be described as an independent comparison.
    if state_error > 1.0e-9:
        print(f"FAIL: LAMMPS replay state transfer error {state_error:.3e} exceeds 1e-9")
        return False
    source = LMP_REPLAY_TEMPLATE.format(
        mu=args.mu_r, E=YOUNGS_MOD, nu=POISSON, e_n=RESTITUTION,
        k_roll=ROLLING_STIFFNESS, gamma_roll=ROLLING_DAMPING, mu_roll=args.mu_r,
        boxhi=LMP_BOX_HI, cyl_r=CYL_RADIUS, catch_r=0.07, datafile=datafile,
        dt=LMP_DT, rest_poll_steps=LMP_RELAX_POLL_STEPS,
        rest_polls=LMP_RELAX_MAX_STEPS // LMP_RELAX_POLL_STEPS,
        rest_speed_limit=LMP_REST_SPEED_M_S, dump=dump,
    )
    input_path = os.path.join(out, "in.lammps")
    with open(input_path, "w", encoding="utf-8") as handle:
        handle.write(source)
    run = subprocess.run([binary, "-in", input_path], cwd=out,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with open(os.path.join(out, "run.log"), "w", encoding="utf-8") as handle:
        handle.write(run.stdout)
    if run.returncode != 0 or not os.path.isfile(dump):
        print(f"FAIL: LAMMPS sentinel failed (exit={run.returncode}); see {out}/run.log")
        return False
    try:
        (xs, ys, zs, radii, vxs, vys, vzs), timestep, frames = _load_lammps_dump(
            dump, HEAP_COUNT)
    except ValueError as exc:
        print(f"FAIL: LAMMPS sentinel produced ineligible geometry: {exc}")
        return False
    lammps_vmax = max(math.sqrt(vx * vx + vy * vy + vz * vz)
                      for vx, vy, vz in zip(vxs, vys, vzs))
    if timestep < LMP_RELAX_POLL_STEPS or not lammps_vmax < LMP_REST_SPEED_M_S:
        print("FAIL: LAMMPS replay final dump does not meet the declared heap-rest qualification.")
        return False
    rc, hs, baseline, diameter = heap_profile(xs, ys, zs, radii)
    theta, toe = fit_angle(rc, hs, baseline, diameter)
    dirt_rc, dirt_hs, dirt_baseline, dirt_diameter = heap_profile(
        dirt_x, dirt_y, dirt_z, dirt_r)
    dirt_theta, dirt_toe = fit_angle(
        dirt_rc, dirt_hs, dirt_baseline, dirt_diameter)
    receipt = {"solver": binary, "mu_r": args.mu_r, "rep": args.rep, "particles": len(xs),
               "formation_protocol": "LAMMPS replay from DIRT qualified pre-lift particle state with both solvers explicitly released from empty contact histories",
               "prelift_particles": replay_count,
               "state_transfer_max_abs_error": state_error,
               "state_transfer_serialization_bound": 1.0e-9,
               "prelift_sha256": _sha256_file(prelift),
               "dirt_final_sha256": _sha256_file(final),
               "dirt_qualification_sha256": _sha256_file(qualification),
               "dirt_theta_deg": dirt_theta,
               "dirt_r_toe_m": dirt_toe,
               "lammps_data_sha256": _sha256_file(datafile),
               "final_timestep": timestep, "complete_frames": frames,
               "lammps_heap_vmax_m_s": lammps_vmax,
               "lammps_heap_rest_limit_m_s": LMP_REST_SPEED_M_S,
               "theta_deg": theta, "r_toe_m": toe,
               "comparison_status": "not_assessed",
               "comparison_reason": (
                   "No predeclared, independently justified DIRT/LAMMPS "
                   "agreement criterion exists for this protocol"),
               "input_sha256": _sha256_file(input_path),
               "dump_sha256": _sha256_file(dump)}
    path = os.path.join(DATA_DIR, "lammps_sentinel.json")
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    lammps_value = "unresolved flank" if theta is None else f"theta={theta:.4f} deg"
    dirt_value = "unresolved flank" if dirt_theta is None else f"theta={dirt_theta:.4f} deg"
    print(f"Cross-code replay: DIRT {dirt_value}; LAMMPS {lammps_value}; receipt -> {path}")
    print("RECEIPT WRITTEN; CROSS-CODE VALIDATION: NOT ASSESSED")
    print("No predeclared independent agreement criterion exists, so this diagnostic cannot pass.")
    return False


# -- generate -------------------------------------------------------------------
def _fmt_mu(mu_r):
    return f"{mu_r:g}"


def _case_key(mu_r, rep):
    return (_fmt_mu(mu_r), int(rep))


def _seed_from_base(base_seed, mu_r, rep):
    """Stable 63-bit positive seed for one case, independent across reps."""
    payload = f"sphcal_angle_of_repose|base={int(base_seed)}|mu_r={_fmt_mu(mu_r)}|rep={int(rep)}"
    digest = hashlib.blake2b(payload.encode("ascii"), digest_size=8).digest()
    seed = int.from_bytes(digest, "big") & ((1 << 63) - 1)
    return seed or 1


def _expected_cases():
    return [(mu_r, rep) for mu_r in MU_R_LIST for rep in range(REPS)]


def _derive_seed_manifest(base_seed):
    return {
        _case_key(mu_r, rep): _seed_from_base(base_seed, mu_r, rep)
        for (mu_r, rep) in _expected_cases()
    }


def _load_seed_manifest(path):
    seeds = {}
    base_seeds = set()
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        required = {"mu_r", "rep", "seed"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            mu_s = _fmt_mu(float(row["mu_r"]))
            rep = int(row["rep"])
            seed = int(row["seed"])
            if seed <= 0 or seed >= (1 << 63):
                raise ValueError(f"{path} has out-of-range seed for mu_r={mu_s} rep={rep}: {seed}")
            key = (mu_s, rep)
            if key in seeds:
                raise ValueError(f"{path} has duplicate seed row for mu_r={mu_s} rep={rep}")
            seeds[key] = seed
            base = row.get("base_seed", "").strip()
            if base:
                base_seeds.add(int(base))
    expected = {_case_key(mu_r, rep) for (mu_r, rep) in _expected_cases()}
    got = set(seeds)
    missing = sorted(expected - got)
    extra = sorted(got - expected)
    if missing or extra:
        msg = []
        if missing:
            msg.append("missing " + ", ".join(f"mu_r={m} rep={r}" for (m, r) in missing))
        if extra:
            msg.append("extra " + ", ".join(f"mu_r={m} rep={r}" for (m, r) in extra))
        raise ValueError(f"{path} does not match this sweep: {'; '.join(msg)}")
    if len(base_seeds) > 1:
        raise ValueError(f"{path} has inconsistent base_seed values: {sorted(base_seeds)}")
    base_seed = next(iter(base_seeds)) if base_seeds else None
    return seeds, base_seed


def _write_seed_manifest(path, seeds, base_seed=None):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["case", "mu_r", "rep", "seed", "base_seed"])
        w.writeheader()
        for (mu_r, rep) in _expected_cases():
            w.writerow({
                "case": case_tag(mu_r, rep),
                "mu_r": _fmt_mu(mu_r),
                "rep": rep,
                "seed": seeds[_case_key(mu_r, rep)],
                "base_seed": "" if base_seed is None else int(base_seed),
            })


def _write_configs(root, seeds):
    n = 0
    for mu_r in MU_R_LIST:
        for rep in range(REPS):
            cdir = _case_dir_in(root, mu_r, rep)
            os.makedirs(cdir, exist_ok=True)
            seed = seeds[_case_key(mu_r, rep)]
            with open(os.path.join(cdir, "config.toml"), "w") as f:
                f.write(_dirt_config(mu_r, cdir, seed))
            n += 1
    return n


def _parse_generate_args(argv):
    p = argparse.ArgumentParser(
        prog="sweep.py generate",
        description="Generate reproducible DIRT configs and a per-case seed manifest.",
    )
    p.add_argument("--base-seed", type=int, default=DEFAULT_BASE_SEED,
                   help=f"base seed used to derive per-case seeds (default: {DEFAULT_BASE_SEED})")
    p.add_argument("--seed-manifest", default=None,
                   help="CSV manifest to read exact per-case seeds from")
    p.add_argument("--write-seed-manifest", default=SEED_MANIFEST,
                   help=f"where to record the generated/read seed table (default: {SEED_MANIFEST})")
    return p.parse_args(argv)


def generate(argv=None):
    args = _parse_generate_args(argv or [])
    if args.seed_manifest:
        seeds, base_seed = _load_seed_manifest(args.seed_manifest)
        manifest_source = f"manifest {args.seed_manifest}"
    else:
        seeds = _derive_seed_manifest(args.base_seed)
        manifest_source = f"base seed {args.base_seed}"
        base_seed = args.base_seed

    n = _write_configs(SWEEP_DIR, seeds)
    _write_seed_manifest(args.write_seed_manifest, seeds, base_seed=base_seed)
    print(f"Generated {n} DIRT configs ({len(MU_R_LIST)} mu_r x {REPS} reps, "
          f"fixed sliding mu_p={MU_P}, seeds from {manifest_source}) under {SWEEP_DIR}")
    print(f"Seed manifest -> {args.write_seed_manifest}")


# -- seed reproducibility evidence ---------------------------------------------
def _config_hashes(root, seeds):
    rows = []
    for mu_r, rep in _expected_cases():
        rel = os.path.join(case_tag(mu_r, rep), "config.toml")
        path = os.path.join(root, rel)
        with open(path, "rb") as f:
            blob = f.read()
        rows.append({
            "case": case_tag(mu_r, rep),
            "mu_r": _fmt_mu(mu_r),
            "rep": rep,
            "seed": seeds[_case_key(mu_r, rep)],
            "relpath": rel,
            "sha256": hashlib.sha256(blob).hexdigest(),
        })
    return rows


def _same_hashes(a, b):
    return [x["sha256"] == y["sha256"] for x, y in zip(a, b)]


def _same_seeds(a, b):
    return [int(x["seed"]) == int(y["seed"]) for x, y in zip(a, b)]


def _write_seed_check_csv(path, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fields = ["case", "mu_r", "rep", "base_seed", "seed",
              "same_base_config_match", "manifest_config_match",
              "changed_base_config_match", "changed_base_seed_match"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def _plot_seed_check(rows, out_png):
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    # Presentation must not mask the deterministic reproducibility assertion in
    # minimal Python environments that intentionally omit matplotlib.
    try:
        import matplotlib
    except ModuleNotFoundError:
        print(f"  plot skipped (matplotlib unavailable): {out_png}")
        return False
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 10})

    labels = [
        "same base seed\n(config hash)",
        "manifest replay\n(config hash)",
        "changed base seed\n(config hash)",
        "changed base seed\n(per-case seed)",
    ]
    values = [
        sum(int(r["same_base_config_match"]) for r in rows),
        sum(int(r["manifest_config_match"]) for r in rows),
        sum(int(r["changed_base_config_match"]) for r in rows),
        sum(int(r["changed_base_seed_match"]) for r in rows),
    ]
    expected = [len(rows), len(rows), 0, 0]
    colors = ["#2ca02c" if v == e else "#d62728" for v, e in zip(values, expected)]

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.bar(range(len(labels)), values, color=colors, width=0.65)
    for i, (v, e) in enumerate(zip(values, expected)):
        ax.text(i, v + 0.25, f"{v}/{len(rows)}\nexpected {e}",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, len(rows) + 2.0)
    ax.set_ylabel("matching cases out of %d" % len(rows))
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_title("Seed-manifest generation reproducibility")
    ax.axhline(len(rows), color="0.75", lw=0.8, zorder=0)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    return True


def _parse_seed_check_args(argv):
    p = argparse.ArgumentParser(
        prog="sweep.py seed-check",
        description="Check and plot seed-manifest config reproducibility.",
    )
    p.add_argument("--base-seed", type=int, default=DEFAULT_BASE_SEED,
                   help=f"reference base seed (default: {DEFAULT_BASE_SEED})")
    p.add_argument("--changed-base-seed", type=int, default=DEFAULT_BASE_SEED + 1,
                   help="different base seed used to prove a new campaign changes")
    p.add_argument("--csv", default=os.path.join(DATA_DIR, "seed_reproducibility.csv"),
                   help="where to write the check table")
    p.add_argument("--plot", default=os.path.join(PLOT_DIR, "seed_reproducibility.png"),
                   help="where to write the committed evidence figure")
    return p.parse_args(argv)


def seed_check(argv=None):
    args = _parse_seed_check_args(argv or [])
    ref_seeds = _derive_seed_manifest(args.base_seed)
    repeat_seeds = _derive_seed_manifest(args.base_seed)
    changed_seeds = _derive_seed_manifest(args.changed_base_seed)

    with tempfile.TemporaryDirectory(prefix="sphcal_seedcheck_") as tmp:
        # Byte identity is defined for regenerating the same checkout paths.
        # Keep the output root fixed so the comparison isolates seed selection.
        root = os.path.join(tmp, "checkout", "sweep")
        manifest = os.path.join(tmp, "seed_manifest.csv")
        _write_seed_manifest(manifest, ref_seeds, base_seed=args.base_seed)
        replay_seeds, _ = _load_seed_manifest(manifest)

        _write_configs(root, ref_seeds)
        ref = _config_hashes(root, ref_seeds)
        _write_configs(root, repeat_seeds)
        repeat = _config_hashes(root, repeat_seeds)
        _write_configs(root, replay_seeds)
        replay = _config_hashes(root, replay_seeds)
        _write_configs(root, changed_seeds)
        changed = _config_hashes(root, changed_seeds)

    same_base_config = _same_hashes(ref, repeat)
    manifest_config = _same_hashes(ref, replay)
    changed_base_config = _same_hashes(ref, changed)
    changed_base_seed = _same_seeds(ref, changed)

    rows = []
    for i, r in enumerate(ref):
        rows.append({
            "case": r["case"],
            "mu_r": r["mu_r"],
            "rep": r["rep"],
            "base_seed": args.base_seed,
            "seed": r["seed"],
            "same_base_config_match": int(same_base_config[i]),
            "manifest_config_match": int(manifest_config[i]),
            "changed_base_config_match": int(changed_base_config[i]),
            "changed_base_seed_match": int(changed_base_seed[i]),
        })

    _write_seed_check_csv(args.csv, rows)
    plotted = _plot_seed_check(rows, args.plot)

    n = len(rows)
    ok = (
        sum(same_base_config) == n
        and sum(manifest_config) == n
        and sum(changed_base_config) == 0
        and sum(changed_base_seed) == 0
        and len({r["seed"] for r in ref}) == n
    )
    print("\n=== Seed-manifest reproducibility check ===")
    print(f"  cases: {n} ({len(MU_R_LIST)} mu_r x {REPS} reps)")
    print(f"  same base seed config matches: {sum(same_base_config)}/{n} (expected {n})")
    print(f"  manifest replay config matches: {sum(manifest_config)}/{n} (expected {n})")
    print(f"  changed base seed config matches: {sum(changed_base_config)}/{n} (expected 0)")
    print(f"  changed base seed per-case seed matches: {sum(changed_base_seed)}/{n} (expected 0)")
    print(f"  distinct seeds in reference campaign: {len({r['seed'] for r in ref})}/{n}")
    print(f"  table -> {args.csv}")
    if plotted:
        print(f"  figure -> {args.plot}")
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


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

    An angle is only an observable when there is a cone flank to fit.  A
    monolayer/pancake has no such flank: treating its missing geometry as 0 deg
    manufactures a numerical observation and can hide a protocol failure.
    Return None for that explicitly unmeasurable state; it is a hard failure of
    the unchanged calibration gate, never a value that can close it."""
    r_toe = _toe_radius(r_centers, h_surface, baseline, diameter)
    if r_toe <= 0:
        return None, 0.0
    lo = APEX_SKIP_FRAC * r_toe
    hi = TOE_HI_FRAC * r_toe
    xf = [r_centers[i] for i in range(len(r_centers)) if lo <= r_centers[i] <= hi]
    yf = [h_surface[i] for i in range(len(r_centers)) if lo <= r_centers[i] <= hi]
    if len(xf) < 3:
        return None, r_toe
    slope, _ = _linfit(xf, yf)
    return math.degrees(math.atan(max(0.0, -slope))), r_toe


def estimator_check():
    """Check the estimator against an analytic cone and a no-flank control.

    This is deliberately independent of solver output: a 24-degree conical
    surface has a known slope, while a monolayer has no angle-of-repose flank.
    The check guards the analysis layer only; it is not a substitute for the
    particle simulation, replicate-spread, monotonicity, or glass-band gates.
    """
    target_deg = 24.0
    target_slope = math.tan(math.radians(target_deg))
    diameter = 0.004
    toe = 0.030
    # A regular radial profile avoids binning/contact-model assumptions and
    # exercises precisely the regression's fitting and flank-selection logic.
    radii = [toe * i / 40.0 for i in range(1, 41)]
    surface = [diameter + target_slope * max(0.0, toe - r) for r in radii]
    recovered, recovered_toe = fit_angle(radii, surface, diameter, diameter)
    pancake, pancake_toe = fit_angle(radii, [diameter] * len(radii), diameter, diameter)
    cone_ok = recovered is not None and abs(recovered - target_deg) <= 1e-10 and recovered_toe > 0.0
    pancake_ok = pancake is None and pancake_toe == 0.0
    print("=== Angle estimator analytical regression ===")
    print(f"  analytic cone: target={target_deg:.6f} deg, recovered="
          f"{('unresolved' if recovered is None else f'{recovered:.6f} deg')}")
    print(f"  monolayer negative control: {'unresolved flank' if pancake is None else f'{pancake:.6f} deg'}")
    print("RESULT:", "PASS" if cone_ok and pancake_ok else "FAIL")
    return cone_ok and pancake_ok


def lammps_mapping_check():
    """Guard the declared DIRT/LAMMPS tangential-law correspondence.

    This is a configuration audit, not a numerical agreement criterion: the
    replay remains an adversarial diagnostic and may disagree.  It catches a
    constitutive mismatch before a replay is interpreted at all.
    """
    mapped_contacts = [
        line.strip() for line in LMP_REPLAY_TEMPLATE.splitlines()
        if "tangential " in line and (line.lstrip().startswith("pair_coeff") or line.lstrip().startswith("fix"))
    ]
    expected = 4
    rescaled = sum("tangential mindlin_rescale/force" in line for line in mapped_contacts)
    plain = sum("tangential mindlin NULL" in line for line in mapped_contacts)
    ok = len(mapped_contacts) == expected and rescaled == expected and plain == 0
    print("=== LAMMPS tangential-law mapping audit ===")
    print(f"  mapped grain/wall contacts: {len(mapped_contacts)}/{expected}")
    print(f"  mindlin_rescale/force declarations: {rescaled}/{expected}")
    print(f"  plain-mindlin declarations: {plain} (expected 0)")
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


# -- start ----------------------------------------------------------------------
def _run_dirt(cdir):
    config = os.path.join(cdir, "config.toml")
    res = os.path.join(cdir, "data", "repose_results.csv")
    witness = os.path.join(cdir, "data", "repose_qualification.json")
    if os.path.exists(res):
        os.remove(res)
    if os.path.exists(witness):
        os.remove(witness)
    log = os.path.join(cdir, "run.log")
    with open(log, "w") as lf:
        proc = subprocess.run(
            ["cargo", "run", "--release", "--example", EXAMPLE,
             "--no-default-features", "--features", "precision-double", "--", config],
            cwd=REPO_ROOT, stdout=lf, stderr=subprocess.STDOUT,
        )
    # A CSV alone is not enough: it could have been emitted after an automatic
    # stage transition rather than after the declared two rest criteria. Require
    # the solver-written event witness; never infer qualification from stdout.
    try:
        with open(witness, encoding="utf-8") as f:
            q = json.load(f)
        # Schema 3 binds whether the contact histories were retained (ordinary
        # campaign) or cleared (the separate LAMMPS diagnostic). The campaign
        # must retain its physical formation history.
        qualified = (
            q["schema"] == 3
            and q["history_at_lift"] == "retained"
            and q["fill_step"] == q["lift_step"] > 0
            and q["heap_step"] >= q["lift_step"] + 2000
            and 0.0 <= float(q["fill_vmax_m_s"]) < REST_MAX_SPEED_M_S
            and 0.0 <= float(q["heap_vmax_m_s"]) < REST_MAX_SPEED_M_S
            and int(q["fill_rest_samples"]) >= 10
            and int(q["heap_rest_samples"]) >= 10
            and int(q["particle_count"]) == HEAP_COUNT
        )
    except (OSError, ValueError, KeyError, TypeError):
        qualified = False
    if proc.returncode != 0 or not os.path.isfile(res) or not qualified:
        if proc.returncode == 0 and not qualified:
            print("unqualified snapshot (missing declared rest events)", end="; ")
        return None
    return res


def _sha256_file(path):
    """Digest an emitted artifact without loading a full trajectory in memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_seed_from_config(path):
    """Read the single inserter seed without adding a TOML dependency."""
    with open(path, encoding="utf-8") as f:
        for line in f:
            key, sep, value = line.partition("=")
            if sep and key.strip() == "seed":
                return int(value.strip().split("#", 1)[0])
    raise ValueError(f"no insertion seed in {path}")


def _case_ledger_entry(mu_r, rep, cdir, results_path):
    """Bind one fitted result to its exact input and solver witness.

    A CSV and a human-readable run log can be copied independently. This local
    ledger is written only after the solver witness has qualified; `graph`
    re-hashes all three artifacts before accepting the campaign. It is provenance
    protection, not an experimental reference or a physics gate.
    """
    config = os.path.join(cdir, "config.toml")
    witness = os.path.join(cdir, "data", "repose_qualification.json")
    prelift = os.path.join(cdir, "data", "repose_prelift.csv")
    return {
        "mu_r": mu_r,
        "rep": rep,
        "seed": _load_seed_from_config(config),
        "config_sha256": _sha256_file(config),
        "results_sha256": _sha256_file(results_path),
        "qualification_sha256": _sha256_file(witness),
        "prelift_sha256": _sha256_file(prelift),
    }


def _write_campaign_ledger(entries):
    with open(CAMPAIGN_LEDGER, "w", encoding="utf-8") as f:
        json.dump({"schema": 1, "cases": entries}, f, indent=2, sort_keys=True)
        f.write("\n")


def _campaign_ledger_ok():
    """Reject stale, mixed, or edited artifacts before any numerical claim."""
    try:
        with open(CAMPAIGN_LEDGER, encoding="utf-8") as f:
            ledger = json.load(f)
        if ledger.get("schema") != 1 or not isinstance(ledger.get("cases"), list):
            raise ValueError("unsupported ledger schema")
        expected = {(mu, rep) for mu, rep in _expected_cases()}
        seen = set()
        for entry in ledger["cases"]:
            mu_r, rep = float(entry["mu_r"]), int(entry["rep"])
            key = (mu_r, rep)
            if key not in expected or key in seen:
                raise ValueError(f"invalid ledger case {key}")
            cdir = case_dir(mu_r, rep)
            artifact_paths = {
                "config_sha256": os.path.join(cdir, "config.toml"),
                "results_sha256": os.path.join(cdir, "data", "repose_results.csv"),
                "qualification_sha256": os.path.join(cdir, "data", "repose_qualification.json"),
                "prelift_sha256": os.path.join(cdir, "data", "repose_prelift.csv"),
            }
            if any(entry[name] != _sha256_file(path)
                   for name, path in artifact_paths.items()):
                raise ValueError(f"digest mismatch for {case_tag(mu_r, rep)}")
            # Keep the path separate from its ledger field name: the prior
            # spelling made it look as though a SHA-256 string was parsed as a
            # TOML path, obscuring that seed replay is independently checked.
            config_path = artifact_paths["config_sha256"]
            if entry["seed"] != _load_seed_from_config(config_path):
                raise ValueError(f"seed mismatch for {case_tag(mu_r, rep)}")
            seen.add(key)
        if seen != expected:
            raise ValueError("ledger is incomplete")
    except (OSError, ValueError, KeyError, TypeError):
        return False
    return True


def start():
    """Run the declared campaign and return whether every case qualified.

    This return value is deliberately about execution provenance only.  A
    complete set of qualified snapshots is still not a calibration pass:
    ``graph`` owns the independent flank, spread, monotonicity, band, and
    external-reference checks.  Keeping the two outcomes separate prevents a
    failed solver run from being reported as a successful campaign merely
    because the driver reached the end of its loop.
    """
    # The independent control is informative only: a real granular heap can
    # recruit a multi-contact force network, so it must not replace solver or
    # external-reference evidence and it does not gate the campaign.
    subprocess.run(
        [sys.executable, os.path.join(SCRIPT_DIR, "physical_feasibility.py")],
        cwd=SCRIPT_DIR,
        check=True,
    )
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Building {EXAMPLE} (release)...", flush=True)
    subprocess.run(["cargo", "build", "--release", "--example", EXAMPLE,
                    "--no-default-features", "--features", "precision-double"], cwd=REPO_ROOT, check=True)

    print("LAMMPS: opt-in external sentinel; running the DIRT formation study only.")

    rows = []
    ledger_entries = []
    profiles = {}  # mu_r -> representative (r_centers, h_surface) from rep 0
    total = len(MU_R_LIST) * REPS
    k = 0
    for mu_r in MU_R_LIST:
        for rep in range(REPS):
            k += 1
            cdir = case_dir(mu_r, rep)
            if not os.path.isfile(os.path.join(cdir, "config.toml")):
                print(f"  [{k:2d}/{total}] missing config mu_r={mu_r} rep={rep} — run 'generate'.")
                continue
            print(f"  [{k:2d}/{total}] mu_r={mu_r:<5} rep={rep}", end="  ", flush=True)
            res = _run_dirt(cdir)
            if res is None:
                print("DIRT FAILED")
                continue
            xs, ys, zs, rad = _load_positions(res)
            r_c, h_s, base, diam = heap_profile(xs, ys, zs, rad)
            theta, r_toe = fit_angle(r_c, h_s, base, diam)
            status = "resolved" if theta is not None else "unresolved_flank"
            # Preserve the case and its qualified raw witness in the record;
            # graph() then fails closed on the state instead of silently
            # dropping it or replacing it with a self-serving zero angle.
            rows.append({"mu_r": mu_r, "rep": rep, "theta_deg": theta,
                         "r_toe": r_toe, "n": len(xs), "status": status})
            ledger_entries.append(_case_ledger_entry(mu_r, rep, cdir, res))
            if theta is None:
                print(f"UNRESOLVED FLANK  (r_toe={r_toe*1e3:.1f} mm, N_heap={len(xs)})")
            else:
                print(f"theta_r = {theta:5.2f} deg  (r_toe={r_toe*1e3:.1f} mm, N_heap={len(xs)})")
            if rep == 0:
                profiles[mu_r] = (r_c, h_s, base, r_toe)

    if not rows:
        print("\nERROR: no DIRT results collected.")
        return False

    os.makedirs(DATA_DIR, exist_ok=True)
    with open(SWEEP_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mu_r", "rep", "theta_deg", "r_toe", "n", "status"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    # A partial campaign has no ledger and therefore cannot be graphed as a
    # qualified calibration result.
    if len(ledger_entries) == total:
        _write_campaign_ledger(ledger_entries)
    print(f"\n{len(rows)}/{total} cases -> {SWEEP_CSV}")

    # Save representative profiles (baseline-subtracted) for the cross-section plot.
    for mu_r, (r_c, h_s, base, r_toe) in profiles.items():
        with open(os.path.join(DATA_DIR, f"profile_{mu_r:g}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["r", "h"])
            for i in range(len(r_c)):
                w.writerow([r_c[i], h_s[i] - base])

    complete = len(ledger_entries) == total
    if not complete:
        print("RESULT: FAIL (one or more cases lacked a solver-qualified snapshot)")
    return complete


# -- graph (validate + plot) ----------------------------------------------------
def _load_sweep():
    if not os.path.isfile(SWEEP_CSV):
        return []
    with open(SWEEP_CSV) as f:
        rows = []
        for r in csv.DictReader(f):
            try:
                r["mu_r"] = float(r["mu_r"])
                r["rep"] = int(float(r["rep"]))
                r["r_toe"] = float(r["r_toe"])
                r["n"] = int(float(r["n"]))
                r["theta_deg"] = None if not r["theta_deg"] else float(r["theta_deg"])
            except (KeyError, TypeError, ValueError):
                return []
            rows.append(r)
        return rows


def _stats_by_mu_r(rows):
    """mu_r -> (mean_theta, std_theta, n_reps), sorted by mu_r."""
    by = {}
    for r in rows:
        by.setdefault(r["mu_r"], []).append(r["theta_deg"])
    out = []
    for mu_r in sorted(by):
        v = by[mu_r]
        mean = sum(v) / len(v)
        var = sum((x - mean) ** 2 for x in v) / len(v) if len(v) > 1 else 0.0
        out.append((mu_r, mean, math.sqrt(var), len(v)))
    return out


def _campaign_shape_ok(rows):
    """Reject partial, duplicated, or malformed data before any calibration gate.

    A graph of a convenient subset of runs must never be able to close the
    calibration.  The required surface is deliberately defined next to the
    sweep constants: every declared (mu_r, rep) pair must occur once, and every
    reported observable must be finite and each realization must retain the
    declared 1,200 mobile grains.  A truncated deposit can look deceptively
    steep or flat depending on which grains escaped, so its fitted angle is not
    a measurement of the specified experiment.  This is a provenance/integrity check;
    it does not relax any of the physical acceptance gates below.
    """
    expected = {(mu, rep) for mu in MU_R_LIST for rep in range(REPS)}
    seen = set()
    for row in rows:
        try:
            key = (float(row["mu_r"]), int(row["rep"]))
            theta = row["theta_deg"]
            n = int(float(row["n"]))
        except (KeyError, TypeError, ValueError):
            print("  malformed campaign row")
            return False
        if row.get("status") == "unresolved_flank":
            print(f"  no resolved cone flank: {key}; angle-of-repose is unmeasurable")
            return False
        if (key not in expected or key in seen or row.get("status") != "resolved"
                or theta is None or not math.isfinite(theta) or n != HEAP_COUNT):
            print(f"  inadmissible campaign row: {key}")
            return False
        seen.add(key)
    missing = expected - seen
    if missing:
        print("  incomplete campaign; missing " + ", ".join(
            f"(mu_r={mu:g}, rep={rep})" for mu, rep in sorted(missing)))
        return False
    return True


def _protocol_reference_ok():
    """Audit a *candidate* record without granting calibration closure.

    A repository JSON transcription plus a Crossref metadata response is useful
    provenance evidence, but cannot independently validate the transcription,
    protocol equivalence, or material transfer.  This helper must therefore
    never promote a numerical sweep to a calibrated material value.
    """
    if not os.path.isfile(PROTOCOL_REFERENCE):
        print("  no protocol-matched external glass record: " + PROTOCOL_REFERENCE)
        return False
    try:
        # Do not treat a locally authored JSON file as external evidence.  A
        # closure must at least survive an independent bibliographic lookup at
        # validation time; a network failure fails closed rather than silently
        # downgrading this requirement to a local self-attestation.
        record = reference_audit.audit_record(PROTOCOL_REFERENCE)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print("  external primary-record audit failed: " + str(error))
        return False
    if not reference_audit.protocol_matches(record, EXTERNAL_PROTOCOL):
        print("  external record does not match the frozen glass/lift/measurement protocol")
        return False
    if record["band_deg"] != [GLASS_BAND_LO_DEG, GLASS_BAND_HI_DEG]:
        print("  external record band differs from the existing 22--26 degree gate")
        return False
    return True


def validate(rows):
    """Evaluate the frozen, solver-local formation-study criteria only.

    This function deliberately does not decide whether a material has been
    calibrated.  Completeness, monotonicity, replicate spread, and the frozen
    22--26 degree target are statements about a DIRT campaign; experimental
    comparability and an SPH transfer remain separate scientific questions.
    Keeping the two results separate makes a successful local run visible
    without allowing it to certify itself as a glass-material calibration.
    """
    print("\n=== Angle-of-repose mu_r formation study ===")
    print(f"  material: E={YOUNGS_MOD:.1e} Pa  nu={POISSON}  e={RESTITUTION} (formation aid)")
    print(f"  FIXED sliding friction: mu_p={MU_P}  (measured glass)")
    print(f"  rolling (sds, swept): mu_r in {MU_R_LIST}  "
          f"k_roll={ROLLING_STIFFNESS:g}  gamma_roll={ROLLING_DAMPING:g}")
    print(f"  target glass band: [{GLASS_BAND_LO_DEG:.0f},{GLASS_BAND_HI_DEG:.0f}] deg")
    if not _campaign_shape_ok(rows):
        print("RESULT: FAIL (campaign integrity)")
        return False
    stats = _stats_by_mu_r(rows)
    print(f"  {'mu_r':>6}{'mean_theta':>12}{'std':>8}{'reps':>6}  note")
    formation_ok = True

    # 1. monotonic increase with mu_r (allow small slack for stochastic dips).
    prev_mean = None
    for (mu_r, mean, std, nrep) in stats:
        note = ""
        if prev_mean is not None and mean < prev_mean - MONOTONIC_SLACK_DEG:
            note = "NON-MONOTONIC"; formation_ok = False
        # 3. reproducibility: spread small (but nonzero — reps must differ).
        if std > SPREAD_MAX_DEG:
            note = (note + " HIGH-SPREAD").strip(); formation_ok = False
        print(f"  {mu_r:>6.2f}{mean:>12.2f}{std:>8.2f}{nrep:>6d}  {note}")
        prev_mean = mean

    # The fixed external band is evaluated only after the ordinary campaign
    # gates; it is never inferred from this solver's output.
    in_band = [(mu_r, mean) for (mu_r, mean, _, _) in stats
               if GLASS_BAND_LO_DEG <= mean <= GLASS_BAND_HI_DEG]
    if not in_band:
        means = ", ".join(f"{m:.1f}" for (_, m, _, _) in stats)
        print(f"  no mu_r lands theta_r in [{GLASS_BAND_LO_DEG:.0f},"
              f"{GLASS_BAND_HI_DEG:.0f}] deg (got {means}) — no calibration closure")
        formation_ok = False

    # overall increase from lowest to highest mu_r.
    if len(stats) >= 2 and stats[-1][1] <= stats[0][1] + 1.0:
        print(f"  theta_r did not increase across mu_r range "
              f"({stats[0][1]:.1f} -> {stats[-1][1]:.1f} deg)")
        formation_ok = False

    if in_band:
        band_mid = 0.5 * (GLASS_BAND_LO_DEG + GLASS_BAND_HI_DEG)
        best = min(in_band, key=lambda t: abs(t[1] - band_mid))
        print(f"\n  in-band formation observation: mu_r = {best[0]:g} -> theta_r = {best[1]:.2f} deg "
              f"(historical target [{GLASS_BAND_LO_DEG:.0f},{GLASS_BAND_HI_DEG:.0f}] deg; not calibrated)")

    print("FORMATION STUDY:", "PASS" if formation_ok else "FAIL")
    return formation_ok


def plot(rows):
    # An unresolved flank is not a zero-degree measurement.  Still render the
    # campaign figure: it shows each missing observable below the physical axis
    # and the required 22--26 degree band, so an incomplete fit cannot disappear
    # from the visual record or be mistaken for a passed calibration.
    os.makedirs(PLOT_DIR, exist_ok=True)
    try:
        import matplotlib
    except ModuleNotFoundError:
        print(f"Figures skipped (matplotlib unavailable): {PLOT_DIR}")
        return False
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 150, "font.size": 11})

    # -- theta_r vs mu_r, plus explicit unresolved-flank markers --
    resolved = [r for r in rows if r.get("status") == "resolved" and r["theta_deg"] is not None]
    unresolved = [r for r in rows if r not in resolved]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    if resolved:
        ax.scatter([r["mu_r"] for r in resolved], [r["theta_deg"] for r in resolved],
                   s=42, alpha=0.85, color="#1f77b4", label="resolved DIRT fit")
    if unresolved:
        ax.scatter([r["mu_r"] for r in unresolved], [-2.0] * len(unresolved),
                   marker="x", s=48, linewidths=1.5, color="#d62728",
                   label="unresolved flank (not 0°)")
    ax.axhspan(GLASS_BAND_LO_DEG, GLASS_BAND_HI_DEG, color="green", alpha=0.10,
               label=f"measured glass band [{GLASS_BAND_LO_DEG:.0f},"
                     f"{GLASS_BAND_HI_DEG:.0f}]°")
    ax.set_xlabel(r"rolling friction $\mu_r$  (sliding $\mu_p$=%.2f fixed)" % MU_P)
    ax.set_ylabel(r"angle of repose $\theta_r$ (deg)")
    ax.set_ylim(-4.0, max(GLASS_BAND_HI_DEG + 3.0,
                          max((r["theta_deg"] for r in resolved), default=0.0) + 3.0))
    ax.axhline(0.0, color="0.55", lw=0.8)
    ax.text(0.305, -2.0, "× = no measurable cone", va="center", ha="right", fontsize=8)
    ax.set_title("Angle of repose vs rolling friction — FAIL if any flank is unresolved")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "theta_vs_mu.png"))
    plt.close(fig)

    # -- heap cross-section profiles h(r) for each mu_r --
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
    murs_sorted = sorted({r["mu_r"] for r in rows})
    for j, mu_r in enumerate(murs_sorted):
        color = cmap(j / max(1, len(murs_sorted) - 1))
        rc, hs = _load_profile(os.path.join(DATA_DIR, f"profile_{mu_r:g}.csv"))
        if rc:
            ax.plot(rc, hs, "o-", ms=3, color=color, label=fr"$\mu_r$={mu_r:g}")
    ax.set_xlabel("radial distance r (mm)")
    ax.set_ylabel("heap surface height h (mm)")
    ax.set_title("Settled heap cross-section (surface envelope)")
    ax.legend(title="rolling friction", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "heap_profile.png"))
    plt.close(fig)

    print(f"\nFigures -> {PLOT_DIR}/theta_vs_mu.png, heap_profile.png")
    return True


def graph():
    rows = _load_sweep()
    if not rows:
        print(f"No {SWEEP_CSV} — run 'start' first.")
        return False
    if not _campaign_ledger_ok():
        print(f"No complete, hash-matched solver ledger at {CAMPAIGN_LEDGER} — run 'generate' then 'start' again.")
        return False
    formation_ok = validate(rows)
    plot(rows)
    if not formation_ok:
        print("CALIBRATION: FAIL (formation-study criteria failed)")
        return False

    # A candidate record can be audited automatically for bibliographic
    # identity and declared protocol fields.  That is valuable adversarial
    # provenance checking, but it cannot establish the human scientific review
    # of source extraction, protocol equivalence, or a DIRT-to-SPH transfer.
    # Do not turn this local check into a numerical pass criterion.
    _protocol_reference_ok()
    print("CALIBRATION: WITHHELD")
    print("  A passing DIRT formation study is not an experimental validation or an SPH transfer.")
    print("  No independently reviewed protocol-matched primary record and cross-substrate")
    print("  validation are committed; no mu_r is pinned.")
    return False


def reference_audit_command(argv):
    """Verify every committed reference record's bibliography with Crossref."""
    if argv:
        print("Usage: sweep.py reference-audit")
        return False
    records_dir = os.path.join(SCRIPT_DIR, "external_records")
    paths = sorted(
        os.path.join(records_dir, name) for name in os.listdir(records_dir)
        if name.endswith(".json")) if os.path.isdir(records_dir) else []
    if not paths:
        print("FAIL: no external reference records found")
        return False
    ok = True
    for path in paths:
        try:
            record = reference_audit.audit_record(path)
            print(f"VERIFIED BIBLIOGRAPHY: {os.path.basename(path)} ({record['doi']}) "
                  f"comparability={record['comparability']}")
        except (OSError, ValueError, json.JSONDecodeError) as error:
            print(f"FAIL: {os.path.basename(path)}: {error}")
            ok = False
    if ok:
        print("REFERENCE AUDIT: PASS (bibliographic identity only; no protocol match or calibration implied)")
    else:
        print("REFERENCE AUDIT: FAIL")
    return ok


# -- dispatch -------------------------------------------------------------------
def main():
    args = sys.argv[1:]
    if args and not args[0].startswith("-"):
        cmd = args[0]
        rest = args[1:]
    else:
        cmd = "all"
        rest = args
    if cmd == "generate":
        generate(rest)
    elif cmd == "seed-check":
        sys.exit(0 if seed_check(rest) else 1)
    elif cmd == "estimator-check":
        if rest:
            print("Usage: sweep.py estimator-check")
            sys.exit(2)
        sys.exit(0 if estimator_check() else 1)
    elif cmd == "lammps-mapping-check":
        if rest:
            print("Usage: sweep.py lammps-mapping-check")
            sys.exit(2)
        sys.exit(0 if lammps_mapping_check() else 1)
    elif cmd == "reference-audit":
        sys.exit(0 if reference_audit_command(rest) else 1)
    elif cmd == "start":
        if rest:
            print("Usage: sweep.py start")
            sys.exit(2)
        sys.exit(0 if start() else 1)
    elif cmd == "replay-generate":
        replay_generate(rest)
    elif cmd == "external":
        sys.exit(0 if external(rest) else 1)
    elif cmd == "graph":
        if rest:
            print("Usage: sweep.py graph")
            sys.exit(2)
        sys.exit(0 if graph() else 1)
    elif cmd == "all":
        generate(rest)
        start()
        print()
        sys.exit(0 if graph() else 1)
    else:
        print(f"Unknown command: {cmd!r}")
        print("Usage: sweep.py [generate|seed-check|estimator-check|lammps-mapping-check|start|external|graph]   (no arg = all three; start prints independent Coulomb control)")
        sys.exit(2)


if __name__ == "__main__":
    main()
