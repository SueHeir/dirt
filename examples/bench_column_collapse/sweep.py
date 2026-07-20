#!/usr/bin/env python3
"""
Granular column-collapse benchmark driver.

Releases a quasi-2D rectangular column of grains (released width L_i, height H_i) on a
flat floor for a range of aspect ratios a = H/L0, then extracts the final runout
L_f from each settled deposit and checks the dimensionless runout against the
experimental planar aspect-ratio scaling law of Lajeunesse et al. (2004):

    (L_f - L0)/L0 ~ 1.2 * a          (a <~ 2-3, linear regime)
    (L_f - L0)/L0 ~ 1.6 * a^(2/3)    (a >~ 3,   power-law regime)

Commands (from anywhere):
    python3 examples/bench_column_collapse/sweep.py generate   # write per-case configs
    python3 examples/bench_column_collapse/sweep.py start      # build + resume qualified sims -> CSV
    python3 examples/bench_column_collapse/sweep.py start --rerun # rerun every DIRT witness
    python3 examples/bench_column_collapse/sweep.py start --case 2,0 # one immutable witness
    python3 examples/bench_column_collapse/sweep.py graph       # extract L_f, validate + plot
    python3 examples/bench_column_collapse/sweep.py            # all three, in order

The aspect ratio is swept by changing the particle count at fixed scheduled gate
width L0. Each realization is fitted using its measured released H_i/L_i. Each aspect is run at several insertion seeds and the
runout is AVERAGED over seeds. Each DIRT run dumps the rest-state deposit
(x, y, z, radius) to data/<case>/column_collapse_results.csv; this script reads
those, computes L_f as the far edge of the deposit toe on a sub-diameter grid
(see measure_column), and fits the runout exponent in each regime.

Measurement quality (not the tolerance) was hardened to remove three possible
fit artifacts: diameter-scale runout quantization (now a sub-diameter toe
metric), single-fabric scatter (three deterministic fabrics), and a coarse
aspect sweep (11 points).  The ±0.25 exponent tolerance is unchanged.  Results
from the superseded 8d × 3d protocol are deliberately not reported by this
driver: they are neither evidence for nor a numerical prediction of the
current 32d × 10d rough-base campaign.

If a LAMMPS binary (lmp_serial / lmp / lmp_mpi / lammps) is on PATH, each aspect
ratio is ALSO run in LAMMPS with the equivalent granular model (pair_style
granular hertz/material ... tangential mindlin ... damping tsuji, same E/nu/e/mu,
gravity, and frictional floor + back + side + removable-gate walls via
fix wall/gran). LAMMPS's final deposit is parsed into the SAME (x, y, z, radius)
form and runout is extracted with the SAME measure_column() the DIRT leg uses, so
the two codes are compared on equal footing and overlaid (open markers) on
plots/runout_scaling.png. LAMMPS is optional: with no binary present, only DIRT
runs and the validation (DIRT-vs-theory) is unchanged.

Outputs:
    sweep/<case>/config.toml            DIRT configs                  (gitignored)
    sweep/<case>/in.lammps              LAMMPS inputs                 (gitignored)
    data/<case>/column_collapse_*.csv   per-case DIRT deposits        (gitignored)
    data/runout.csv                     L0, H, a, L_f per DIRT case   (gitignored)
    data/lammps_results.csv             L0, H, a, L_f per LAMMPS case (gitignored)
    plots/*.png                         final figures                 (tracked)
"""

import os
import sys
import csv
import math
import shutil
import subprocess
import contextlib
import fcntl
import hashlib
import json
import concurrent.futures
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
EXAMPLE = "bench_column_collapse"

SWEEP_DIR = os.path.join(SCRIPT_DIR, "sweep")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PLOT_DIR = os.path.join(SCRIPT_DIR, "plots")
RUNOUT_CSV = os.path.join(DATA_DIR, "runout.csv")          # DIRT runout per aspect
LAMMPS_CSV = os.path.join(DATA_DIR, "lammps_results.csv")  # LAMMPS runout per aspect
LAMMPS_RECEIPT_NAME = "lammps_receipt.json"
MANIFEST_NAME = "column_collapse_protocol.json"
CASE_RECEIPT_NAME = "column_collapse_case_receipt.json"
CAMPAIGN_LOCK_NAME = ".column_collapse_campaign.lock"
CASE_LOCK_NAME = ".column_collapse_case.lock"
JOB_MANIFEST_NAME = "column_collapse_jobs.tsv"

# LAMMPS binary candidates, in preference order. LAMMPS is optional: if none is
# found, the LAMMPS leg is skipped and only DIRT is run/plotted.
LAMMPS_BINS = ["lmp_serial", "lmp", "lmp_mpi", "lammps"]

# ── Geometry / material (shared by every case) ───────────────────────────────
RADIUS = 0.0015            # m (d = 3 mm; Lajeunesse used ~1–3 mm glass beads)
DENSITY = 2500.0           # kg/m^3 (glass)
# The previous 16d x 6d section left a toe only a few grains wide/deep.  Its
# fitted exponent changed when the front definition changed, which is direct
# evidence that it was not in the continuum regime of the experiments.  Keep
# the material, aspect schedule and acceptance band fixed, but resolve the
# initial cross-section at 32d x 10d.  At fixed aspect this is 6.7x as many
# active grains and makes the measured toe a bulk deposit rather than a
# two-to-three-grain feature.
L0 = 0.096                 # initial column width [m] (= 32 diameters)
W = 0.030                  # slab width in y [m] (= 10 diameters, quasi-3D)
# Before the named gate is removed, no active grain may pass its physical face.
# This admits one radius for contact with the plane plus a small numerical margin;
# it is a boundary-condition witness, not a fitted acceptance tolerance.
GATE_RELEASE_WIDTH_MAX = L0 + 1.05 * RADIUS
# The bead bed must cover every place the deposit may come to rest.  Stopping it
# at L0 would change the basal boundary from rough grains to the smooth safety
# plane immediately after release, precisely where runout is determined.
BASE_X_HIGH = 0.60          # m; the downstream fixed-domain boundary
# Canonical glass-bead (ballotini) material — measured properties, shared across
# all DIRT calibrations (shear/cooling/conduction/collapse). E softened from the
# real ~65 GPa (rigid-grain limit; keeps dt tractable). e and μ_p are measured
# glass–glass values (Wu et al. 2019, Meas. of restitution & friction for glass beads).
YOUNGS_MOD = 7.0e7         # Pa (softened from ~65 GPa real glass)
POISSON = 0.245
RESTITUTION = 0.926        # measured glass–glass COR
FRICTION = 0.16            # measured glass–glass sliding friction
# The exact source/base coordinates are non-overlapping (their minimum centre
# distance is one diameter), but the 32d x 10d source creates thousands of
# simultaneous Hertz contacts on its first loaded step.  A 1 us step resolves
# that loaded network.  Keep the same step after gate removal until a separate
# timestep-convergence study has established that a coarser released step
# preserves the measured runout.  A successful trajectory at a larger step is
# only an execution check; it is not convergence evidence for this experiment.
SETTLE_DT = 1.0e-6
COLLAPSE_DT = 1.0e-6
SETTLE_STEPS = 800000      # 0.8 s overdamped preparation; released law unchanged
COLLAPSE_STEPS = 4000000   # 4 s at COLLAPSE_DT: enough for a high-e glass bed to arrest
# The inserted fabric is deliberately dense and supports a full 32d x 10d
# column from its first gravity step.  At 10 um/step the former limiter allowed
# a newly loaded contact to advance several Hertz overlaps before Cundall
# damping could dissipate its insertion transient; the reference a=0.5 source
# then aborted during settling with the force kernel's >500-large-overlap
# safeguard.  Limit the *settling-only* displacement to one micrometre.  The
# 0.8 s duration gives that overdamped preparation time to dissipate the load
# wave before the independently checked release frame.  Both this limiter and
# Cundall damping are removed before gate release, so this is a stable
# preparation integrator, not a change to the released contact law.
PREPARATION_MAX_DISPLACEMENT = 1.0e-6

PACKING = 0.60             # settled solid fraction used to size the particle count
# The source is a non-overlapping close packing, not a dilute airborne cloud.
# Its packing fraction is fixed by fcc geometry; PACKING remains the requested
# released-column fraction used to translate aspect ratio to particle count.
INSERT_PACKING = math.pi / (3.0 * math.sqrt(2.0))
# A perfect close packing has a mechanically special contact network.  It is a
# useful non-overlap construction, but it is not an independently prepared
# granular fabric: merely changing ABC registries leaves every grain exactly at
# contact.  Start each realization just above close packing and apply a small,
# deterministic in-plane perturbation.  The 8% dilation gives this prepared
# source its intended roughly-0.60 solid fraction.  Its scheduled source height
# is compensated below, because its measured *released* height is the physical
# control variable and the release witness remains the acceptance test.
SOURCE_DILATION = 1.08
SOURCE_JITTER = 0.005      # fraction of a diameter, per horizontal coordinate
BASE_Z = 2.0 * RADIUS
BASE_SELECT_Z = 2.5 * RADIUS
# The removable gate starts at the *top surface* of the frozen rough base.  A
# plane wall is otherwise infinite, so after the containment-side repair it
# would also push on the downstream base beads that deliberately extend beyond
# the gate.  Those beads are support, not retained-column particles.  Bounding
# the gate above the base makes the two boundary conditions disjoint while
# retaining a full-height barrier for every mobile grain.
GATE_Z_LOW = 2.0 * RADIUS

# Aspect ratios to sweep. Spans both regimes (linear a<~2-3, power-law a>~3) so a
# regime change is resolvable. The sweep is deliberately FINE — extra points in
# BOTH regimes (7 in the linear regime, 5 in the power regime) so each least-
# squares exponent is fit from many points instead of a coarse handful, which was
# a dominant source of fit noise in the earlier 6-point sweep.
ASPECTS = [0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

# Insertion RNG seeds. Every aspect ratio is run once per seed and the runout is
# AVERAGED over seeds before the exponent is fit. Packing randomness in a quasi-2D
# (3-diameter-deep) column produced ~±20-25% run-to-run scatter in the runout at a
# single seed — the other dominant source of fit noise — so seed-averaging is
# essential to a stable exponent. (The DIRT inserter is seeded by the `seed` field
# on [[particles.insert]]; LAMMPS, if present, is run once at its own fixed seed.)
SEEDS = [0, 1, 2]

# ── Runout metric ────────────────────────────────────────────────────────────
# Sub-diameter deposit-toe front (see measure_column). This keeps the ORIGINAL
# physical definition — the runout is the far edge of the deposit where it is at
# least ~one diameter tall (a grain centre >= one diameter above the floor, i.e.
# >= 2 grain layers, which rejects the single-grain leading smear and lone rolled
# fliers) — but REMOVES that definition's diameter-scale quantization. The earlier
# metric binned x at one particle DIAMETER and reported the bin edge, quantizing
# L_f to 3 mm steps (~0.125 in normalized runout) so whole aspect ratios shared a
# runout bin. FINE_BINS subdivides each diameter into sub-cells (sub-diameter
# resolution); GAP_TOL_D is how much clear vacuum (in diameters) may separate
# neighbouring toe grains before the contiguous deposit is judged to have ended.
FINE_BINS = 8              # x-cells per particle diameter (sub-diameter resolution)
GAP_TOL_D = 1.0            # a clear gap > 1 grain diameter ends the contiguous deposit
TOE_MIN_HEIGHT_D = 1.0     # deposit toe = grain centre >= this many diameters up (>=2 layers)

# External experimental reference. This quasi-2D, planar-gate benchmark must
# not silently borrow the high-aspect *axisymmetric* Lube et al. law (whose
# exponent is different). The two exponents and crossover below are the planar
# relation reported by Lajeunesse, Mangeney-Castelnau & Vilotte, Phys. Fluids 16,
# 2371–2381 (2004), doi:10.1063/1.1736611. Keep this source separate from the
# numerical fit: DIRT output is compared to it; it does not define itself.
EMPIRICAL_REFERENCE = {
    "authors": "Lajeunesse, Mangeney-Castelnau & Vilotte",
    "journal": "Physics of Fluids 16 (2004) 2371–2381",
    "doi": "10.1063/1.1736611",
    "geometry": "planar granular mass on a horizontal plane",
}

# Validation tolerances on the externally sourced fitted exponents. UNCHANGED —
# the measurement is improved, the pass band is not touched.
EXP_TOL = 0.25             # |fitted exponent - target| pass band
LINEAR_TARGET = 1.0        # (L_f-L_i)/L_i ~ a^1   for a <~ 2-3
POWER_TARGET = 2.0 / 3.0   # (L_f-L_i)/L_i ~ a^2/3 for a >~ 3
REGIME_SPLIT = 3.0         # aspect ratio dividing the two regimes
# Initialization-fidelity admission check, not a fitted-exponent tolerance.
ASPECT_REL_TOL = 0.02

# Filled lazily after the source-coordinate functions are defined.  A common
# population is retained across seeds so seed averaging does not alter mass.
_SOURCE_POPULATION_CACHE = {}


def protocol_fingerprint():
    """Stable identity of the physical/measurement contract behind a campaign.

    ``runout.csv`` is deliberately retained locally so graphing does not rerun a
    costly ensemble.  It must therefore carry enough immutable context to reject
    a CSV produced by a different base, material, geometry, seed plan, toe
    estimator, or acceptance band.  This is not a fitted quantity and is not
    included in any numerical result; it is an evidence provenance guard.
    """
    contract = {
        "geometry": [RADIUS, DENSITY, L0, W, GATE_RELEASE_WIDTH_MAX, PACKING, INSERT_PACKING,
                     SOURCE_DILATION, SOURCE_JITTER,
                     BASE_Z, BASE_SELECT_Z, GATE_Z_LOW],
                     "material": [YOUNGS_MOD, POISSON, RESTITUTION, FRICTION, SETTLE_DT, COLLAPSE_DT,
                                  PREPARATION_MAX_DISPLACEMENT],
        "schedule": [SETTLE_STEPS, COLLAPSE_STEPS, ASPECTS, SEEDS],
        # Source preparation is part of the physical protocol: a campaign made
        # with an earlier translated crystal cannot be relabelled as this
        # fabric ensemble merely because its summary rows look compatible.
        "initialization": ["deterministic-stacking-disordered-fcc-source-v3-geometry-qualified"],
        "boundary": [BASE_X_HIGH, rough_base_positions(),
                     "frozen_close_packed_bead_layer_full_runout"],
        "measurement": [FINE_BINS, GAP_TOL_D, TOE_MIN_HEIGHT_D],
        "validation": [EXP_TOL, LINEAR_TARGET, POWER_TARGET, REGIME_SPLIT,
                       ASPECT_REL_TOL, RELEASE_FROUDE_MAX,
                       PREPARATION_WINDOW_SAMPLES, PREPARATION_SAMPLE_INTERVAL,
                       REST_FROUDE_MAX, ARREST_WINDOW_SAMPLES,
                       ARREST_SAMPLE_INTERVAL],
    }
    encoded = json.dumps(contract, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def n_particles(aspect):
    """Exact source population for a scheduled release aspect.

    The old bulk-volume formula was appropriate only for an unbounded packing.
    This benchmark instead inserts a finite, phase-shifted triangular source.
    Its edge losses and registry shifts mean that the formula prepared a column
    substantially taller than the requested aspect *before it could settle*.
    Choose the largest common population whose actual generated source remains
    below the scheduled height for every seed.  The analysis still uses the
    measured post-settlement height; this function only prevents the source
    itself from silently relabelling the aspect schedule.
    """
    if aspect not in _SOURCE_POPULATION_CACHE:
        # The source is deliberately dilated to make a non-overlapping,
        # approximately 0.60-solid-fraction preparation.  Its artificial fcc
        # clearance closes during the damped settling stage, so initializing it
        # at the released target height would systematically make the measured
        # release short by the known dilation.  Compensate the *preparation*
        # height only; the independent raw release witness still has to meet the
        # unchanged 2% H_i/L_i admission before any fit is allowed.
        target = aspect * L0 * SOURCE_DILATION
        # The infinite-packing count is a safe upper bound for this deliberately
        # dilated finite source.  Search the actual coordinate generator, not a
        # second approximate packing model.
        bulk = max(1, int(math.ceil(
            PACKING * L0 * W * target / ((4.0 / 3.0) * math.pi * RADIUS**3)
        )))
        lo, hi = 1, bulk
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if max(source_preparation_height(active_column_positions(mid, aspect, seed))
                   for seed in SEEDS) <= target + 1e-12:
                lo = mid
            else:
                hi = mid - 1
        _SOURCE_POPULATION_CACHE[aspect] = lo
    return _SOURCE_POPULATION_CACHE[aspect]


def sites_per_source_layer():
    """Number of admissible triangular-lattice sites in one source layer."""
    d = 2.0 * RADIUS
    spacing = SOURCE_DILATION * d
    dy = math.sqrt(3.0) * spacing / 2.0
    # Registry shifts and alternating row offsets can only reduce capacity by
    # one site.  Use the minimum across all registries as a conservative bound.
    capacities = []
    for px, py in ((0.0, 0.0), (0.5 * spacing, dy / 3.0),
                   (spacing, 2.0 * dy / 3.0)):
        n = 0
        iy = 0
        while RADIUS + py + iy * dy <= W - RADIUS + 1e-12:
            x0 = RADIUS + px + (0.5 * spacing if iy % 2 else 0.0)
            ix = 0
            while x0 + ix * spacing <= L0 - RADIUS + 1e-12:
                n += 1
                ix += 1
            iy += 1
        capacities.append(n)
    return min(capacities)


def loose_insert_top(count, aspect):
    """Conservative top bound of the fcc source, measured from z=0.

    A simple-cubic source is only pi/6 dense, so it cannot represent the
    requested 0.60 solid fraction without making the nominal column too tall.
    The fcc stack has layer spacing sqrt(2/3)d and packing pi/(3 sqrt(2)).
    """
    layers = int(math.ceil(count / sites_per_source_layer()))
    return BASE_Z + RADIUS + layers * math.sqrt(2.0 / 3.0) * (2.0 * RADIUS)


def rough_base_positions():
    """One glued hexagonal bead layer: the rough granular substrate.

    A square grid has the right nearest-neighbour spacing in its cardinal
    directions but is not close packed: its area fraction is pi/4, rather than
    pi/(2 sqrt(3)).  That is a different basal topography and should not be
    described as the close-packed rough bed used in granular-collapse work.
    Alternate row registries give each mobile grain a genuine three-point rough
    support while retaining a minimum centre separation of one diameter.
    """
    d = 2.0 * RADIUS
    dy = math.sqrt(3.0) * d / 2.0
    points = []
    iy = 0
    while True:
        y = RADIUS + iy * dy
        if y > W - RADIUS + 1e-12:
            break
        x0 = RADIUS + (0.5 * d if iy % 2 else 0.0)
        ix = 0
        while True:
            x = x0 + ix * d
            if x > BASE_X_HIGH - RADIUS + 1e-12:
                break
            points.append((x, y, RADIUS))
            ix += 1
        iy += 1
    return points


def write_rough_base():
    os.makedirs(SWEEP_DIR, exist_ok=True)
    path = os.path.join(SWEEP_DIR, "rough_base.csv")
    with open(path, "w", newline="") as f:
        csv.writer(f).writerows(rough_base_positions())
    return path


def stacking_registries(seed, layers):
    """Return a reproducible, stacking-disordered close-packed sequence.

    The former ``(layer + seed) % 3`` rule made the three advertised seeds mere
    translations of one perfect ABC crystal.  Such translations cannot sample
    preparation variability and therefore are not an ensemble.  Here every
    source starts in the same supported registry, while the subsequent close
    packed layers choose one of the two registries different from the preceding
    layer using a hash of the seed and layer.  Adjacent layers remain exactly
    tangent (and never overlap), but each seed has a distinct stacking-fault
    realization--a standard, controlled way to decrystallize a close packing
    without changing grain size, count, material, walls, or acceptance band.
    """
    if layers <= 0:
        return []
    registries = [0]
    for layer in range(1, layers):
        choices = [p for p in range(3) if p != registries[-1]]
        digest = hashlib.sha256(f"column-collapse-v2:{seed}:{layer}".encode()).digest()
        registries.append(choices[digest[0] & 1])
    return registries


def _source_jitter(seed, layer, iy, ix):
    """Return a deterministic bounded in-plane perturbation for one source site.

    Hashing rather than consuming a process-global RNG makes a source replayable
    from its seed and site label, independent of iteration order.  The dilation
    below leaves more than twice this maximum approach margin between nearest
    sites, and the exact spatial audit remains the final authority.
    """
    raw = hashlib.sha256(
        f"column-collapse-source-v3:{seed}:{layer}:{iy}:{ix}".encode()
    ).digest()
    scale = SOURCE_JITTER * 2.0 * RADIUS
    return ((raw[0] / 255.0 - 0.5) * 2.0 * scale,
            (raw[1] / 255.0 - 0.5) * 2.0 * scale)


def active_column_positions(count, aspect, seed):
    """Return one exact-count prepared granular source for release.

    The prior 15-by-5 simple grid was non-overlapping but only pi/6 dense.  A
    fresh nominal-a=0.5 replay therefore released at H/L0=0.817, invalidating
    the requested aspect-ratio sweep.  This fcc construction has a strict
    minimum centre separation d while its compactness is sufficient for the
    target 0.60-volume-fraction column.  Seeds select distinct ABC registries
    *and* bounded in-plane source perturbations, rather than rigid translations
    or an exactly crystalline contact network.  The source is slightly dilated
    so the perturbations cannot create overlap; settling still occurs before
    the gate is removed.
    """
    d = 2.0 * RADIUS
    spacing = SOURCE_DILATION * d
    dy = math.sqrt(3.0) * spacing / 2.0
    dz = math.sqrt(2.0 / 3.0) * spacing
    phases = ((0.0, 0.0), (0.5 * spacing, dy / 3.0),
              (spacing, 2.0 * dy / 3.0))
    # The fixed-width envelope can make shifted layers hold fewer sites.  Keep a
    # source sequence long enough even when a preparation change (such as the
    # deliberate dilation above) reduces per-layer capacity; ``count`` is a
    # conservative finite upper bound because every accepted layer has >=1 site.
    registries = stacking_registries(seed, count)
    points = []
    layer = 0
    while len(points) < count:
        px, py = phases[registries[layer]]
        z = BASE_Z + RADIUS + layer * dz
        iy = 0
        while True:
            y = RADIUS + py + iy * dy
            if y > W - RADIUS + 1e-12:
                break
            ix = 0
            while True:
                # Triangular in-layer registry: without this half-diameter row
                # offset, neighbouring rows are only sqrt(3)/2 d apart.
                x = RADIUS + px + ix * spacing + (0.5 * spacing if iy % 2 else 0.0)
                if x > L0 - RADIUS + 1e-12:
                    break
                jx, jy = _source_jitter(seed, layer, iy, ix)
                # Keep every source centre inside the same fixed L0 × W
                # envelope; randomization must not change a boundary condition.
                x = min(L0 - RADIUS, max(RADIUS, x + jx))
                y = min(W - RADIUS, max(RADIUS, y + jy))
                points.append((x, y, z))
                if len(points) == count:
                    return points
                ix += 1
            iy += 1
        layer += 1
    return points


def source_min_separation(points):
    """Return the exact nearest centre spacing using diameter-sized cells."""
    d = 2.0 * RADIUS
    cells = {}
    best = float("inf")
    for point in points:
        key = tuple(math.floor(v / d) for v in point)
        for ix in range(key[0] - 1, key[0] + 2):
            for iy in range(key[1] - 1, key[1] + 2):
                for iz in range(key[2] - 1, key[2] + 2):
                    for other in cells.get((ix, iy, iz), ()):
                        best = min(best, math.dist(point, other))
        cells.setdefault(key, []).append(point)
    return best


def source_preparation_height(points):
    """Physical source height measured with particle envelopes above the bed."""
    if not points:
        raise ValueError("empty active-column source")
    return max(z + RADIUS for _, _, z in points) - BASE_Z


def audit_active_source(points, count, aspect):
    """Fail before simulation if a seed source is not an admissible packing."""
    if len(points) != count:
        raise ValueError(f"source population {len(points)} does not equal {count}")
    minimum = source_min_separation(points)
    if minimum < 2.0 * RADIUS * (1.0 - 1.0e-10):
        raise ValueError(f"source has overlapping grains: minimum separation {minimum}")
    # A finite source changes height in whole layers.  One fcc layer is the
    # unavoidable discretisation uncertainty, not an acceptance tolerance.
    # This validates only the compensated *source* geometry.  It cannot prove
    # how the bed settles; ``checked_release_dimensions`` remains the dynamic
    # evidence gate for the scheduled physical aspect.
    target = aspect * L0 * SOURCE_DILATION
    height_error = abs(source_preparation_height(points) - target)
    layer = math.sqrt(2.0 / 3.0) * SOURCE_DILATION * (2.0 * RADIUS)
    if height_error > layer + 1e-12:
        raise ValueError(
            f"source height {source_preparation_height(points) / L0:.6f} L0 "
            f"does not represent scheduled aspect {aspect:.6f}"
        )


def write_active_column(path, count, aspect, seed):
    points = active_column_positions(count, aspect, seed)
    audit_active_source(points, count, aspect)
    with open(path, "w", newline="") as f:
        csv.writer(f).writerows(points)
    return path


def source_digest(points):
    """Digest coordinates in the representation consumed by both solvers.

    The seed label is not evidence that a case used that seed: a stale or copied
    ``active_column.csv`` can retain a plausible population, width, and height.
    Bind the campaign to the *canonical coordinates* instead.  This is kept
    separate from the protocol fingerprint, which identifies the model and
    estimator but deliberately does not enumerate tens of thousands of points.
    """
    payload = "".join(
        f"{x:.10e},{y:.10e},{z:.10e}\n" for x, y, z in points
    ).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def source_file_digest(path, expected_count):
    """Read a generated source and return its canonical-coordinate digest."""
    points = active_column_from_file(path, expected_count)
    return source_digest(points)


def protocol_manifest():
    """Recompute the canonical, case-level preparation witness.

    This intentionally derives every expected source from code rather than
    trusting files left under ``sweep/``.  It is cheap relative to a 33-case DEM
    campaign and makes a copied seed, edited source, or stale generated tree
    inadmissible before runs are launched or figures are regenerated.
    """
    cases = []
    for aspect in ASPECTS:
        count = n_particles(aspect)
        for seed in SEEDS:
            points = active_column_positions(count, aspect, seed)
            audit_active_source(points, count, aspect)
            cases.append({"nominal_aspect": aspect, "seed": seed,
                          "active_count": count,
                          "active_source_sha256": source_digest(points)})
    base = rough_base_positions()
    return {"schema": 1, "protocol_sha256": protocol_fingerprint(),
            "rough_base_count": len(base),
            "rough_base_sha256": source_digest(base), "cases": cases}


def checked_protocol_manifest(write=False):
    """Require every on-disk source file to match the current protocol.

    ``write=True`` records the self-describing witness beside campaign output;
    ``False`` is used by analysis to reject a source tree that no longer matches
    the one that a valid run would have consumed.
    """
    manifest = protocol_manifest()
    expected_base = manifest["rough_base_sha256"]
    base_path = os.path.join(SWEEP_DIR, "rough_base.csv")
    if not os.path.isfile(base_path):
        raise ValueError("missing rough-base source")
    if source_file_digest(base_path, manifest["rough_base_count"]) != expected_base:
        raise ValueError("rough-base source does not match protocol")
    for case in manifest["cases"]:
        path = os.path.join(case_dir_seed(case["nominal_aspect"], case["seed"]),
                            "active_column.csv")
        if not os.path.isfile(path):
            raise ValueError(f"missing active source a={case['nominal_aspect']} seed={case['seed']}")
        if source_file_digest(path, case["active_count"]) != case["active_source_sha256"]:
            raise ValueError(f"active source does not match protocol a={case['nominal_aspect']} seed={case['seed']}")
    if write:
        os.makedirs(DATA_DIR, exist_ok=True)
        target = os.path.join(DATA_DIR, MANIFEST_NAME)
        fd, temporary = tempfile.mkstemp(prefix=".column_collapse_protocol.",
                                         dir=DATA_DIR, text=True)
        with os.fdopen(fd, "w") as f:
            json.dump(manifest, f, sort_keys=True, indent=2)
            f.write("\n")
        os.replace(temporary, target)
    return manifest


def case_tag(aspect):
    return f"a{aspect:g}".replace(".", "p")


def case_dir(aspect):
    """Seed-0 case directory (also used by the LAMMPS leg and the profile plot)."""
    return os.path.join(SWEEP_DIR, case_tag(aspect))


def case_dir_seed(aspect, seed):
    """Per-(aspect, seed) directory. Seed 0 reuses case_dir() so LAMMPS and the
    representative-deposit plot keep pointing at a stable location."""
    if seed == 0:
        return case_dir(aspect)
    return os.path.join(SWEEP_DIR, f"{case_tag(aspect)}_s{seed}")


@contextlib.contextmanager
def exclusive_campaign_lock(operation):
    """Exclude writers while replacing the shared generated source fabric.

    This lock deliberately protects only generation and the short manifest
    snapshot that follows it.  A full campaign contains 33 independent dynamic
    witnesses, so holding a single campaign lock during simulation would make
    scheduler-dispatched, non-overlapping cases needlessly conflict.
    """
    os.makedirs(SWEEP_DIR, exist_ok=True)
    path = os.path.join(SWEEP_DIR, CAMPAIGN_LOCK_NAME)
    with open(path, "a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError(
                "another column-collapse source generation owns this campaign; "
                "wait for it to finish before replacing the source fabric"
            ) from exc
        lock.seek(0)
        lock.truncate()
        lock.write(f"pid={os.getpid()} operation={operation}\n")
        lock.flush()
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


@contextlib.contextmanager
def shared_campaign_lock(operation):
    """Keep source generation out while one or more campaigns consume it."""
    os.makedirs(SWEEP_DIR, exist_ok=True)
    path = os.path.join(SWEEP_DIR, CAMPAIGN_LOCK_NAME)
    with open(path, "a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError(
                "column-collapse source generation is in progress; retry after it completes"
            ) from exc
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


@contextlib.contextmanager
def exclusive_case_lock(aspect, seed, operation):
    """Lease one deterministic witness directory for its complete write.

    The release, pre-release-rest, deposit, arrest window, and content receipt form one atomic
    evidence unit.  A non-blocking per-case lease rejects only a duplicate
    ``(aspect, seed)`` launch; other members of the declared 11x3 ensemble may
    run concurrently in the same campaign directory.
    """
    cdir = case_dir_seed(aspect, seed)
    os.makedirs(cdir, exist_ok=True)
    path = os.path.join(cdir, CASE_LOCK_NAME)
    with open(path, "a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError(
                f"column-collapse witness a={aspect} seed={seed} is already being {operation}"
            ) from exc
        lock.seek(0)
        lock.truncate()
        lock.write(f"pid={os.getpid()} operation={operation}\n")
        lock.flush()
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


@contextlib.contextmanager
def shared_case_lock(aspect, seed, operation):
    """Admit a stable completed witness without blocking an active writer."""
    cdir = case_dir_seed(aspect, seed)
    path = os.path.join(cdir, CASE_LOCK_NAME)
    if not os.path.isdir(cdir):
        raise ValueError(f"missing witness directory a={aspect} seed={seed}")
    with open(path, "a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError(
                f"column-collapse witness a={aspect} seed={seed} is being written; "
                f"cannot {operation} a mixed snapshot"
            ) from exc
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def data_case_dir(aspect):
    return os.path.join(DATA_DIR, case_tag(aspect))


# ── DIRT config template ─────────────────────────────────────────────────────
TOML_TEMPLATE = """\
# Auto-generated column-collapse config — aspect a = {aspect}, N = {count}
[comm]
processors_x = 1
processors_y = 1
processors_z = 1

[domain]
x_low = -0.01
x_high = 0.60
y_low = -0.003
y_high = {y_high}
z_low = 0.0
z_high = {z_high}
boundary_x = "fixed"
boundary_y = "fixed"
boundary_z = "fixed"

[neighbor]
skin_fraction = 1.1
bin_size = 0.005
every = 1

[gravity]
gx = 0.0
gy = 0.0
gz = -9.81

[dem]
contact_model = "hertz"

[[dem.materials]]
# Same glass contact law as the released grains.  This distinct type exists
# solely to make the frozen rough boundary's membership survive atom sorting;
# it is not a second material model.  It must remain material type 0 because
# the dynamic rough-base group below uses that type identity.
name = "rough_glass"
youngs_mod = {youngs:.6e}
poisson_ratio = {poisson}
restitution = {restitution}
friction = {friction}

[[dem.materials]]
name = "glass"
youngs_mod = {youngs:.6e}
poisson_ratio = {poisson}
restitution = {restitution}
friction = {friction}

[[particles.insert]]
source = "file"
file = "{rough_base}"
format = "csv"
material = "rough_glass"
radius = {radius}
density = {density}
columns = {{ x = 0, y = 1, z = 2 }}

[[particles.insert]]
source = "file"
file = "{active_column}"
format = "csv"
material = "glass"
radius = {radius}
density = {density}
columns = {{ x = 0, y = 1, z = 2 }}

[[group]]
name = "rough_base"
# A static spatial mask is invalid after atom sorting: its index membership
# follows slots rather than the frozen grains.  Rebuild by the dedicated,
# physically identical rough-base type each step so the support stays fixed.
# DIRT material type IDs are zero-based; ``rough_glass`` is declared first.
type = [0]
dynamic = true

[[freeze]]
group = "rough_base"

# Numerical preparation only.  A dense, deliberately non-overlapping source
# carries a load wave when gravity is first applied; Cundall damping and the
# conservative displacement cap dissipate that artificial insertion transient.
# `main.rs` removes both fixes before the gate is opened, so neither participates
# in the measured collapse dynamics.
[[cundall]]
group = "all"
gamma_l = 0.8
gamma_a = 0.8

[[nve_limit]]
group = "all"
max_displacement = {preparation_max_displacement}

[[wall]]
type = "plane"
point_z = 0.0
normal_z = 1.0
material = "glass"
name = "floor"

[[wall]]
type = "plane"
point_x = 0.0
normal_x = 1.0
material = "glass"
name = "back"

[[wall]]
type = "plane"
point_y = 0.0
normal_y = 1.0
material = "glass"
name = "side_lo"

[[wall]]
type = "plane"
point_y = {w}
normal_y = -1.0
material = "glass"
name = "side_hi"

[[wall]]
type = "plane"
point_x = {l0}
normal_x = -1.0
material = "glass"
name = "gate"
# This finite gate begins at the rough-base surface.  It must not contact the
# frozen downstream support, which intentionally spans beyond x=L0.
bound_z_low = {gate_z_low}

[output]
dir = "{output_dir}"

[vtp]
interval = 1000000

[[run]]
name = "settle"
steps = {settle_steps}
thermo = 20000
dt = {settle_dt}

[[run]]
name = "collapse"
steps = {collapse_steps}
thermo = 20000
dt = {collapse_dt}
"""


def generate():
    with exclusive_campaign_lock("generate"):
        os.makedirs(SWEEP_DIR, exist_ok=True)
        rough_base = write_rough_base()
        n_cfg = 0
        for a in ASPECTS:
            n = n_particles(a)
            insert_top = loose_insert_top(n, a)
            z_high = max(0.2, insert_top + 0.05)
            for s in SEEDS:
                cdir = case_dir_seed(a, s)
                os.makedirs(cdir, exist_ok=True)
                active_column = write_active_column(
                    os.path.join(cdir, "active_column.csv"), n, a, s
                )
                with open(os.path.join(cdir, "config.toml"), "w") as f:
                    f.write(TOML_TEMPLATE.format(
                    aspect=a, count=n, seed=s,
                    youngs=YOUNGS_MOD, poisson=POISSON,
                    restitution=RESTITUTION, friction=FRICTION,
                    radius=RADIUS, density=DENSITY,
                    l0=f"{L0:.4f}",
                    w=f"{W:.4f}", y_high=f"{W + 0.003:.4f}",
                    # Config-runner file paths are resolved from the repository
                    # root (the same cwd used by ``start``).  Keeping them
                    # repository-relative makes a generated campaign portable:
                    # a receipt can be rerun from a clean checkout instead of
                    # silently retaining the workstation that generated it.
                    rough_base=os.path.relpath(rough_base, REPO_ROOT), active_z_low=f"{BASE_Z + RADIUS:.4f}",
                    base_select_z=f"{BASE_SELECT_Z:.4f}",
                    gate_z_low=f"{GATE_Z_LOW:.4f}",
                    insert_top=f"{insert_top:.4f}", z_high=f"{z_high:.4f}",
                    active_column=os.path.relpath(active_column, REPO_ROOT),
                    output_dir=os.path.relpath(cdir, REPO_ROOT),
                    settle_dt=f"{SETTLE_DT:.3e}", collapse_dt=f"{COLLAPSE_DT:.3e}",
                    preparation_max_displacement=f"{PREPARATION_MAX_DISPLACEMENT:.3e}",
                    settle_steps=SETTLE_STEPS, collapse_steps=COLLAPSE_STEPS,
                    ))
                n_cfg += 1
        # Generation itself is a physics preflight: retain a digestible witness of
        # the exact base and all 33 source fabrics, but do not mistake it for
        # dynamic evidence.  ``start`` rewrites it immediately before execution.
        checked_protocol_manifest(write=True)
    print(f"Generated {n_cfg} configs ({len(ASPECTS)} aspects x {len(SEEDS)} seeds) "
          f"under {SWEEP_DIR}")


def emit_jobs():
    """Write an auditable, scheduler-neutral map for the 33 raw witnesses.

    A full continuum-resolution campaign is intentionally too large to hide
    behind one interactive invocation. This manifest is not a result and does
    not create ``runout.csv``: it gives a batch system exactly one immutable
    witness per row, bound to the canonical source digest that ``graph`` later
    requires. A failed or missing row remains a fail-closed incomplete ensemble.
    """
    with shared_campaign_lock("emit-jobs"):
        manifest = checked_protocol_manifest(write=True)
        by_case = {(float(row["nominal_aspect"]), int(row["seed"])): row
                   for row in manifest["cases"]}
        path = os.path.join(SWEEP_DIR, JOB_MANIFEST_NAME)
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=("index", "aspect", "seed", "active_count",
                            "active_source_sha256", "protocol_sha256", "command"),
                delimiter="\t",
            )
            writer.writeheader()
            for index, (aspect, seed) in enumerate(
                    ((a, s) for a in ASPECTS for s in SEEDS), start=1):
                witness = by_case[(float(aspect), int(seed))]
                writer.writerow({
                    "index": index,
                    "aspect": f"{aspect:g}",
                    "seed": seed,
                    "active_count": witness["active_count"],
                    "active_source_sha256": witness["active_source_sha256"],
                    "protocol_sha256": manifest["protocol_sha256"],
                    "command": (
                        "python3 examples/bench_column_collapse/sweep.py "
                        f"start --case {aspect:g},{seed}"
                    ),
                })
    print(f"Wrote 33 immutable DIRT witness jobs -> {path}")
    return path


def campaign_status():
    """Report admission of the declared raw ensemble without fitting it.

    This is deliberately a provenance/status operation, not a cheaper form of
    ``graph``.  In particular, an all-admitted ensemble is reported only as
    ``READY_FOR_GRAPH``: the independent observer and the experimental
    exponent gate must still run before either PASS or FAIL can be claimed.
    Conversely, one missing, malformed, moving, or stale witness leaves the
    physical result *incomplete*, rather than inviting a partial fit.
    """
    result = {
        "declared_cases": len(ASPECTS) * len(SEEDS),
        "protocol_sha256": protocol_fingerprint(),
        "state": "UNPREPARED",
        "admitted_cases": [],
        "inadmissible_cases": [],
    }
    try:
        with shared_campaign_lock("inspect campaign status"):
            manifest = checked_protocol_manifest(write=False)
            result["protocol_sha256"] = manifest["protocol_sha256"]
            for aspect in ASPECTS:
                for seed in SEEDS:
                    reason = _case_evidence_error(aspect, seed)
                    case = {"aspect": aspect, "seed": seed}
                    if reason is None:
                        result["admitted_cases"].append(case)
                    else:
                        case["reason"] = reason
                        result["inadmissible_cases"].append(case)
    except ValueError as exc:
        result["reason"] = str(exc)
        return result

    result["state"] = (
        "READY_FOR_GRAPH" if not result["inadmissible_cases"] else "INCOMPLETE"
    )
    return result


def print_campaign_status():
    """Print a compact, scheduler-safe ledger and return its machine state."""
    status = campaign_status()
    admitted = len(status["admitted_cases"])
    total = status["declared_cases"]
    print(f"Column-collapse campaign: {status['state']} ({admitted}/{total} admitted)")
    print(f"  protocol: {status['protocol_sha256']}")
    if status["state"] == "UNPREPARED":
        print(f"  source admission failed: {status['reason']}")
    elif status["state"] == "INCOMPLETE":
        for case in status["inadmissible_cases"]:
            print(f"  a={case['aspect']:g}, seed={case['seed']}: {case['reason']}")
        print("  No runout, exponent, PASS, or numerical-failure claim is available.")
    else:
        print("  Raw witnesses are admitted; run 'graph' for independent observation and")
        print("  the unchanged experimental exponent gate. This is not a PASS or FAIL.")
    return status["state"]


# ── LAMMPS leg (optional cross-code overlay) ─────────────────────────────────
def find_lammps():
    """Return the first available LAMMPS binary on PATH, or None."""
    for b in LAMMPS_BINS:
        path = shutil.which(b)
        if path:
            return path
    return None


# LAMMPS counterpart of the DIRT column collapse. Same material, same geometry,
# same two-stage protocol (settle against a gate, remove the gate, collapse):
#   pair_style granular hertz/material E e nu  -> Young's modulus, restitution,
#       Poisson ratio (E/nu/e identical to DIRT's [dem.materials]).
#   tangential mindlin NULL {damp} {mu}        -> Mindlin tangential spring with
#       k_t = 8 G* sqrt(R* delta) derived from the normal contact (NULL) — exactly
#       DIRT's k_t — Coulomb friction coefficient {mu}.
#   damping tsuji                              -> viscoelastic normal+tangential
#       damping from the restitution e (DIRT uses the same e-> damping mapping).
#   fix grav ... 9.81 vector 0 0 -1            -> gravity g_z = -9.81, matching
#       [gravity] in the DIRT config.
#   fix wall/gran ... zplane/xplane/yplane     -> frictional floor, back wall, and
#       the quasi-2D side walls — same granular model (incl. friction) as the pair
#       style, the LAMMPS analogue of DIRT's frictional dirt_wall planes.
#   fix gate ... xplane NULL {L0}; unfix gate  -> removable gate at x = L0, present
#       during 'settle', unfix-ed at the start of 'collapse' (mirrors
#       Walls::deactivate_by_name on the first collapse step in DIRT).
# LAMMPS receives the exact active-column coordinate file consumed by DIRT,
# including the deterministic per-seed phase relative to the same frozen rough
# base. The final deposit is dumped as (id, x, y, z, radius); runout is then
# extracted with the SAME measure_column().
LMP_TEMPLATE = """\
# Auto-generated LAMMPS input for the column-collapse sweep — aspect a = {aspect}
units           si
atom_style      sphere
dimension       3
boundary        f f f
newton          off
comm_modify     vel yes

region          simbox block {x_low} {x_high} {y_low} {y_high} 0.0 {z_high} units box
create_box      1 simbox

{active_atoms}
{base_atoms}
set             group all diameter {diam}
set             group all density {density}

region          rough_base_region block INF INF INF INF 0.0 {base_select_z} units box
group           rough_base region rough_base_region
group           mobile subtract all rough_base

pair_style      granular
pair_coeff      1 1 hertz/material {E} {e} {nu} tangential mindlin NULL {tdamp} {mu} damping tsuji rolling none twisting none

# The glued layer is an immobile granular boundary.  ``fix freeze`` makes
# pair contacts treat it as infinite-mass, but it must also be excluded from
# body, wall, and integration fixes: applying a wall force after freeze can
# reintroduce a force on the supposed fixed particles.  In particular, the
# floor contacts the tangent base layer at startup.  The mobile group is the
# exact complement, so this changes no force acting on a released grain.
fix             grav mobile gravity {g} vector 0 0 -1
fix             base_freeze rough_base freeze
fix             floor mobile wall/gran granular hertz/material {E} {e} {nu} tangential mindlin NULL {tdamp} {mu} damping tsuji rolling none twisting none zplane 0.0 NULL
fix             back mobile wall/gran granular hertz/material {E} {e} {nu} tangential mindlin NULL {tdamp} {mu} damping tsuji rolling none twisting none xplane 0.0 NULL
fix             sides mobile wall/gran granular hertz/material {E} {e} {nu} tangential mindlin NULL {tdamp} {mu} damping tsuji rolling none twisting none yplane 0.0 {W}
fix             gate mobile wall/gran granular hertz/material {E} {e} {nu} tangential mindlin NULL {tdamp} {mu} damping tsuji rolling none twisting none xplane NULL {L0}
# These two fixes are a quasi-static source-preparation aid only.  They are
# removed together before gate release, mirroring DIRT's stage-local Cundall
# damping and displacement cap below.
fix             settle_damp mobile damping/cundall 0.8 0.8
fix             integrate mobile nve/limit {preparation_max_displacement}
# Match DIRT's sustained-quiescence evidence.  Preparation and released
# dynamics get separate files: a low-speed final deposit cannot certify the
# still-gated state that physically defines a gate-release experiment.
variable        speed atom sqrt(vx*vx+vy*vy+vz*vz)
compute         max_speed mobile reduce max v_speed
fix             preparation mobile ave/time {preparation_interval} 1 {preparation_interval} c_max_speed file {preparation_file} mode scalar

thermo_modify   lost warn flush yes
timestep        {settle_dt}
thermo          {thermo}

# Stage 1: settle the loose column against the gate.
run             {settle_steps}
write_dump      all custom {release_dump} id x y z radius modify sort id

# Stage 2: remove the gate; the column collapses and spreads to rest.
unfix           gate
unfix           settle_damp
unfix           integrate
unfix           preparation
timestep        {collapse_dt}
fix             integrate mobile nve/sphere
fix             arrest mobile ave/time {arrest_interval} 1 {arrest_interval} c_max_speed file {arrest_file} mode scalar
run             {collapse_steps}

write_dump      all custom {dump} id x y z radius vx vy vz modify sort id
"""


def lammps_dump_path(aspect, seed, stage):
    return os.path.join(case_dir_seed(aspect, seed), f"lammps_{stage}.txt")


def lammps_receipt_path(aspect, seed):
    """Receipt for one external-code witness, kept beside its raw LAMMPS files."""
    return os.path.join(case_dir_seed(aspect, seed), LAMMPS_RECEIPT_NAME)


def lammps_binary_identity(binary):
    """Content identity, rather than a PATH spelling, of the independent solver."""
    resolved = os.path.realpath(binary)
    if not os.path.isfile(resolved):
        raise ValueError("LAMMPS executable is not a regular file")
    return {"path": resolved, "sha256": sha256_file(resolved)}


def lammps_case_receipt(aspect, seed, binary):
    """Bind a LAMMPS witness to its rendered input, solver and raw outputs.

    A complete dump alone does not establish which LAMMPS input produced it.
    This receipt is deliberately provenance only: it neither supplies a target
    nor participates in the DIRT-versus-experiment exponent verdict.
    """
    cdir = case_dir_seed(aspect, seed)
    paths = {
        "input": os.path.join(cdir, "in.lammps"),
        "release_dump": lammps_dump_path(aspect, seed, "release"),
        "deposit_dump": lammps_dump_path(aspect, seed, "deposit"),
        "release_csv": os.path.join(cdir, "lammps_release.csv"),
        "deposit_csv": os.path.join(cdir, "lammps_deposit.csv"),
        "preparation": os.path.join(cdir, "lammps_preparation.txt"),
        "arrest": os.path.join(cdir, "lammps_arrest.txt"),
    }
    missing = [name for name, path in paths.items() if not os.path.isfile(path)]
    if missing:
        raise ValueError("missing LAMMPS receipt input(s): " + ", ".join(missing))
    return {
        "schema": 1,
        "nominal_aspect": aspect,
        "seed": seed,
        "protocol_sha256": protocol_fingerprint(),
        "lammps": lammps_binary_identity(binary),
        "input_sha256": sha256_file(paths.pop("input")),
        "witness_sha256": {name: sha256_file(path) for name, path in paths.items()},
    }


def write_lammps_case_receipt(aspect, seed, binary):
    with open(lammps_receipt_path(aspect, seed), "w") as f:
        json.dump(lammps_case_receipt(aspect, seed, binary), f, sort_keys=True, indent=2)
        f.write("\n")


def checked_lammps_case_receipt(aspect, seed):
    path = lammps_receipt_path(aspect, seed)
    if not os.path.isfile(path):
        raise ValueError("missing LAMMPS per-case receipt")
    try:
        with open(path) as f:
            recorded = json.load(f)
        binary = recorded["lammps"]["path"]
        current = lammps_case_receipt(aspect, seed, binary)
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError("malformed LAMMPS per-case receipt") from exc
    if recorded != current:
        raise ValueError("LAMMPS receipt disagrees with executable, input, or raw witnesses")


def lammps_create_atoms(points):
    """Render one LAMMPS ``create_atoms single`` command per source coordinate.

    DIRT reads these same coordinates from ``active_column.csv``.  Rendering
    them directly is deliberately more verbose than LAMMPS's random inserter:
    it makes code-to-code agreement a comparison of the contact/integration
    implementations rather than of two different initial packings.  The
    preflight count check below makes a truncated or stale source file fail
    before LAMMPS is launched.
    """
    return "\n".join(
        f"create_atoms    1 single {x:.10e} {y:.10e} {z:.10e} units box"
        for x, y, z in points
    )


def active_column_from_file(path, expected_count):
    """Read and validate the exact active-column source shared with DIRT."""
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    if len(rows) != expected_count:
        raise ValueError(
            f"active-column source has {len(rows)} coordinates, expected {expected_count}"
        )
    points = []
    for row in rows:
        if len(row) != 3:
            raise ValueError("active-column source row is not x,y,z")
        try:
            point = tuple(float(v) for v in row)
        except ValueError as exc:
            raise ValueError("active-column source contains non-numeric coordinate") from exc
        if not all(math.isfinite(v) for v in point):
            raise ValueError("active-column source contains non-finite coordinate")
        points.append(point)
    return points


def write_lammps_input(path, aspect, seed):
    """Write a LAMMPS case from the exact DIRT active-column source.

    Both solvers receive the identical active grains and frozen rough base.  A
    seeded ensemble is still retained because each source file is a distinct
    deterministic phase relative to the base; LAMMPS must not substitute a
    different random realization for that phase.
    """
    n = n_particles(aspect)
    source_path = os.path.join(case_dir_seed(aspect, seed), "active_column.csv")
    active_points = active_column_from_file(source_path, n)
    z_high = max(z for _, _, z in active_points) + RADIUS + 0.05
    # The active source is serialized with ten significant decimal places.  The
    # glued bed must use the same representation: at d=3 mm, the former fixed
    # decimal format could round an exactly tangent source/base pair inward.
    # That turns an intended zero-overlap preparation into a large Hertz force
    # before either solver advances.  This is an initialization fidelity repair,
    # not a material or tolerance adjustment.
    base_atoms = "\n".join(
        f"create_atoms    1 single {x:.10e} {y:.10e} {z:.10e} units box"
        for x, y, z in rough_base_positions()
    )
    with open(path, "w") as f:
        f.write(LMP_TEMPLATE.format(
            aspect=aspect, count=n,
            x_low=-0.01, x_high=0.60, y_low=-0.003, y_high=W + 0.003,
            z_high=f"{z_high:.4f}", active_atoms=lammps_create_atoms(active_points),
            base_atoms=base_atoms, base_select_z=f"{BASE_SELECT_Z:.4f}",
            diam=2.0 * RADIUS, density=DENSITY,
            E=f"{YOUNGS_MOD:.6e}", e=RESTITUTION, nu=POISSON,
            tdamp=1.0, mu=FRICTION, g=9.81,
            W=W, L0=L0, settle_dt=f"{SETTLE_DT:.3e}",
            collapse_dt=f"{COLLAPSE_DT:.3e}", thermo=40000,
            preparation_max_displacement=f"{PREPARATION_MAX_DISPLACEMENT:.3e}",
            settle_steps=SETTLE_STEPS, collapse_steps=COLLAPSE_STEPS,
            release_dump=lammps_dump_path(aspect, seed, "release"),
            dump=lammps_dump_path(aspect, seed, "deposit"),
            preparation_file=os.path.join(case_dir_seed(aspect, seed), "lammps_preparation.txt"),
            preparation_interval=PREPARATION_SAMPLE_INTERVAL,
            arrest_file=os.path.join(case_dir_seed(aspect, seed), "lammps_arrest.txt"),
            arrest_interval=ARREST_SAMPLE_INTERVAL,
        ))


def lammps_dump_to_csv(dump_path, csv_path):
    """Convert a LAMMPS 'id x y z radius' dump to the same x,y,z,radius CSV that
    the DIRT recorder writes, so measure_column() can read it unchanged."""
    with open(dump_path) as f:
        lines = f.readlines()
    # Find the 'ITEM: ATOMS' header; columns follow it.
    start = None
    cols = []
    for i, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            cols = line.split()[2:]
            start = i + 1
            break
    if start is None:
        return False
    ix, iy, iz, ir = (cols.index(c) for c in ("x", "y", "z", "radius"))
    with open(csv_path, "w", newline="") as out:
        out.write("x,y,z,radius\n")
        for line in lines[start:]:
            p = line.split()
            if len(p) < len(cols):
                continue
            out.write(f"{p[ix]},{p[iy]},{p[iz]},{p[ir]}\n")
    return True


def lammps_max_speed(dump_path):
    """Read the terminal speed from a LAMMPS custom dump, if present."""
    with open(dump_path) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            cols = line.split()[2:]
            if not {"vx", "vy", "vz"}.issubset(cols):
                raise ValueError("LAMMPS terminal dump has no velocity columns")
            iv = [cols.index(k) for k in ("vx", "vy", "vz")]
            return max((math.sqrt(sum(float(p[j]) ** 2 for j in iv))
                        for p in (line.split() for line in lines[i + 1:])
                        if len(p) >= len(cols)), default=0.0)
    raise ValueError("LAMMPS terminal dump has no atom section")


def lammps_arrest_window(path):
    """Read the final sustained-speed window written by LAMMPS ``fix ave/time``."""
    rows = []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if len(fields) != 2:
                raise ValueError("malformed LAMMPS arrest-window row")
            try:
                step, speed = int(fields[0]), float(fields[1])
            except ValueError as exc:
                raise ValueError("non-numeric LAMMPS arrest-window row") from exc
            rows.append((step, speed))
    if len(rows) < ARREST_WINDOW_SAMPLES:
        raise ValueError("too few LAMMPS arrest-window samples")
    window = rows[-ARREST_WINDOW_SAMPLES:]
    if (any(step <= 0 for step, _ in window)
            or any(b[0] - a[0] != ARREST_SAMPLE_INTERVAL for a, b in zip(window, window[1:]))
            or not all(math.isfinite(speed) and speed >= 0.0 for _, speed in window)):
        raise ValueError("invalid LAMMPS arrest-window values")
    return [speed for _, speed in window]


def lammps_preparation_window(path):
    """Read the still-gated LAMMPS rest witness ending at gate removal.

    This deliberately mirrors ``checked_preparation_window`` for DIRT.  The
    LAMMPS overlay is optional, but when it is rendered it must begin from a
    demonstrably quiet, still-gated source rather than use final arrest as a
    proxy for a valid release.
    """
    rows = []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            if len(fields) != 2:
                raise ValueError("malformed LAMMPS preparation-window row")
            try:
                step, speed = int(fields[0]), float(fields[1])
            except ValueError as exc:
                raise ValueError("non-numeric LAMMPS preparation-window row") from exc
            rows.append((step, speed))
    if len(rows) < PREPARATION_WINDOW_SAMPLES:
        raise ValueError("too few LAMMPS preparation-window samples")
    window = rows[-PREPARATION_WINDOW_SAMPLES:]
    if (window[-1][0] != SETTLE_STEPS
            or any(step <= 0 for step, _ in window)
            or any(b[0] - a[0] != PREPARATION_SAMPLE_INTERVAL for a, b in zip(window, window[1:]))
            or not all(math.isfinite(speed) and speed >= 0.0 for _, speed in window)):
        raise ValueError("invalid LAMMPS preparation-window values")
    return [speed for _, speed in window]


def run_lammps_sweep(lammps):
    """Run the same aspect×seed campaign in LAMMPS or return no overlay.

    A one-seed/final-only overlay is not an independent comparison to DIRT's
    three-seed arrested ensemble.  This deliberately fails closed unless every
    realization has exact release/final populations and the same Froude gate.
    """
    rows = []
    failures = []
    for i, a in enumerate(ASPECTS, 1):
        lfs, hs = [], []
        for seed in SEEDS:
            cdir = case_dir_seed(a, seed)
            os.makedirs(cdir, exist_ok=True)
            in_path = os.path.join(cdir, "in.lammps")
            log_path = os.path.join(cdir, "lammps.log")
            release_dump = lammps_dump_path(a, seed, "release")
            deposit_dump = lammps_dump_path(a, seed, "deposit")
            release_csv = os.path.join(cdir, "lammps_release.csv")
            deposit_csv = os.path.join(cdir, "lammps_deposit.csv")
            preparation = os.path.join(cdir, "lammps_preparation.txt")
            arrest = os.path.join(cdir, "lammps_arrest.txt")
            for stale in (release_dump, deposit_dump, release_csv, deposit_csv, preparation, arrest,
                          lammps_receipt_path(a, seed)):
                if os.path.isfile(stale): os.remove(stale)
            write_lammps_input(in_path, a, seed)
            print(f"  [LAMMPS {i}/{len(ASPECTS)}] a={a:<4} seed={seed} N={n_particles(a)}", flush=True)
            proc = subprocess.run([lammps, "-in", in_path, "-log", log_path], cwd=REPO_ROOT,
                                  stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            if (proc.returncode != 0 or not all(os.path.isfile(p) for p in (release_dump, deposit_dump, preparation, arrest))
                    or not lammps_dump_to_csv(release_dump, release_csv)
                    or not lammps_dump_to_csv(deposit_dump, deposit_csv)):
                failures.append(f"a={a} seed={seed}: no parseable LAMMPS snapshots")
                continue
            expected = total_particles(a)
            if csv_particle_count(release_csv) != expected or csv_particle_count(deposit_csv) != expected:
                failures.append(f"a={a} seed={seed}: LAMMPS population not {expected} at release/final")
                continue
            try:
                # External-code admission is deliberately independent of the
                # DIRT receipt: establish from LAMMPS's own release dump that
                # its frozen granular support was neither lost nor moved.
                checked_lammps_release_support(release_csv)
                release_geometry(release_csv)
                preparation_speeds = lammps_preparation_window(preparation)
                vmax = lammps_max_speed(deposit_dump)
                arrest_speeds = lammps_arrest_window(arrest)
            except ValueError as exc:
                failures.append(f"a={a} seed={seed}: {exc}")
                continue
            if max(preparation_speeds) / math.sqrt(9.81 * 2.0 * RADIUS) > RELEASE_FROUDE_MAX:
                failures.append(f"a={a} seed={seed}: LAMMPS still-gated preparation is not at rest")
                continue
            if vmax / math.sqrt(9.81 * 2.0 * RADIUS) > REST_FROUDE_MAX:
                failures.append(f"a={a} seed={seed}: LAMMPS terminal state is not arrested")
                continue
            if max(arrest_speeds) / math.sqrt(9.81 * 2.0 * RADIUS) > REST_FROUDE_MAX:
                failures.append(f"a={a} seed={seed}: LAMMPS final sustained state is not arrested")
                continue
            # Record the external executable and the fully rendered input only
            # after every raw witness has passed its physical admission checks.
            write_lammps_case_receipt(a, seed, lammps)
            hs.append(release_height(release_csv)); _, lf = measure_column(deposit_csv); lfs.append(lf)
        if len(lfs) != len(SEEDS):
            continue
        rn = [(v - L0) / L0 for v in lfs]
        rows.append({"nominal_aspect": a, "aspect": sum(hs) / len(hs) / L0,
                     "L0": L0, "H": sum(hs) / len(hs), "L_f": sum(lfs) / len(lfs),
                     "runout_norm": sum(rn) / len(rn),
                     "runout_std": (sum((v - sum(rn) / len(rn)) ** 2 for v in rn) / len(rn)) ** 0.5,
                     "n_seeds": len(lfs), "protocol_sha256": protocol_fingerprint()})
    if failures:
        print("LAMMPS: " + "; ".join(failures))
    return rows if len(rows) == len(ASPECTS) else []


# ── start ────────────────────────────────────────────────────────────────────
REST_FROUDE_MAX = 0.05
# A gate-release experiment must start from the same quiescent state required
# at the end of the run. This is an additional admission witness, not a new
# fit tolerance.
RELEASE_FROUDE_MAX = REST_FROUDE_MAX
# A final snapshot may coincide with a transient low-speed phase.  Require four
# equally spaced collapse witnesses instead; this is a stricter arrest check,
# not a relaxation of the existing Froude threshold.
ARREST_WINDOW_SAMPLES = 4
# Preserve the 0.1 s physical sampling interval under the released 1 us step.
ARREST_SAMPLE_INTERVAL = 100_000
PREPARATION_WINDOW_SAMPLES = 4
PREPARATION_SAMPLE_INTERVAL = 100_000
# LAMMPS's default custom-dump coordinate format is six decimal places.  This
# is still 30,000 times smaller than a bead diameter (3 mm), while allowing the
# documented text round-trip (worst-case 0.5 micrometre in displayed units).
# It is an input/output comparison tolerance, never a physics acceptance band.
RELEASE_COORDINATE_TOL = 1.0e-7


def csv_particle_count(path):
    with open(path, newline="") as f:
        return sum(1 for _ in csv.DictReader(f))


def _coordinate_key(point):
    """Stable text-round-trip key for an immobile source coordinate."""
    return tuple(int(round(value / RELEASE_COORDINATE_TOL)) for value in point)


def checked_lammps_release_support(path):
    """Require the independent solver to retain the declared rough support.

    ``fix freeze`` is intended to make the emitted LAMMPS rough-base particles
    immobile.  A release/final population check alone cannot establish that it
    did: a missing support particle can be replaced by a mobile particle while
    preserving the total count.  Compare the raw release dump with the
    canonical base coordinates before using it as an external-code witness.

    This checks a boundary-condition invariant, not DIRT's runout or either
    experimental exponent.  It therefore cannot make the DIRT-vs-experiment
    gate pass; it only prevents a non-equivalent LAMMPS boundary from being
    displayed as a cross-code comparison.
    """
    expected = {_coordinate_key(point) for point in rough_base_positions()}
    if len(expected) != len(rough_base_positions()):
        raise ValueError("rough-base source has duplicate coordinates")
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    try:
        observed = {_coordinate_key((float(row["x"]), float(row["y"]),
                                     float(row["z"]))) for row in rows}
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("malformed LAMMPS release-state record") from exc
    missing = expected - observed
    if missing:
        raise ValueError(
            f"LAMMPS release does not retain {len(missing)} frozen rough-base coordinate(s)"
        )


def sha256_file(path):
    """Return the content digest of one campaign input or raw witness."""
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dirt_case_paths(a, seed):
    """Return raw witnesses for one DIRT realization, including its run trace.

    The log is not a fitted observable or an additional numerical criterion.
    It records the executable's stage/gate/final-population trace, so an
    independently dispatched receipt cannot be detached from the run that
    produced its CSV witnesses.
    """
    case_dir = case_dir_seed(a, seed)
    case_data = os.path.join(case_dir, "data")
    return {
        "deposit": os.path.join(case_data, "column_collapse_results.csv"),
        "release": os.path.join(case_data, "column_collapse_release.csv"),
        "terminal": os.path.join(case_data, "column_collapse_final_state.csv"),
        "arrest": os.path.join(case_data, "column_collapse_arrest.csv"),
        "preparation": os.path.join(case_data, "column_collapse_preparation.csv"),
        "run_log": os.path.join(case_dir, "run.log"),
    }


def case_receipt_path(a, seed):
    return os.path.join(case_dir_seed(a, seed), "data", CASE_RECEIPT_NAME)


def dirt_binary_identity(binary):
    """Content identity of the DIRT recorder used for one raw witness.

    The configuration and recorder source identify the intended protocol, but
    neither one establishes which executable actually wrote a trajectory.  A
    content digest binds the receipt to the executable that was launched,
    rather than trusting a target-path spelling that a later build can replace.
    """
    resolved = os.path.realpath(binary)
    if not os.path.isfile(resolved):
        raise ValueError("DIRT executable is not a regular file")
    return {"path": resolved, "sha256": sha256_file(resolved)}


def case_receipt(a, seed, binary):
    """Make an immutable-content receipt for a completed DIRT witness.

    A complete set of CSV files is necessary but not sufficient evidence: a
    resumed campaign must also establish which generated source, executable
    recorder, and configuration produced them.  The receipt deliberately hashes
    raw files rather than fit results, so it cannot tune or manufacture an
    exponent.  It detects ordinary stale/mixed artifacts; it is not a signature
    scheme and makes no claim to defend against a writer able to replace both
    data and receipt.
    """
    cdir = case_dir_seed(a, seed)
    paths = dirt_case_paths(a, seed)
    required = {
        "config": os.path.join(cdir, "config.toml"),
        "active_source": os.path.join(cdir, "active_column.csv"),
        "rough_base": os.path.join(SWEEP_DIR, "rough_base.csv"),
        "recorder_source": os.path.join(SCRIPT_DIR, "main.rs"),
        **paths,
    }
    missing = [name for name, path in required.items() if not os.path.isfile(path)]
    if missing:
        raise ValueError("missing receipt input(s): " + ", ".join(missing))
    return {
        "schema": 1,
        "nominal_aspect": a,
        "seed": seed,
        "protocol_sha256": protocol_fingerprint(),
        "expected_particle_count": total_particles(a),
        "dirt": dirt_binary_identity(binary),
        "input_sha256": {
            "config": sha256_file(required["config"]),
            "active_source": sha256_file(required["active_source"]),
            "rough_base": sha256_file(required["rough_base"]),
            "recorder_source": sha256_file(required["recorder_source"]),
        },
        "witness_sha256": {name: sha256_file(path) for name, path in paths.items()},
    }


def write_case_receipt(a, seed, binary):
    receipt = case_receipt(a, seed, binary)
    path = case_receipt_path(a, seed)
    with open(path, "w") as f:
        json.dump(receipt, f, sort_keys=True, indent=2)
        f.write("\n")


def checked_case_receipt(a, seed):
    """Reject a witness whose raw artifacts no longer match its run receipt."""
    path = case_receipt_path(a, seed)
    if not os.path.isfile(path):
        raise ValueError("missing per-case content receipt")
    try:
        with open(path) as f:
            recorded = json.load(f)
        binary = recorded["dirt"]["path"]
        current = case_receipt(a, seed, binary)
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError("malformed per-case content receipt") from exc
    if recorded != current:
        raise ValueError("per-case receipt disagrees with executable, inputs, or raw witnesses")


def total_particles(aspect):
    return n_particles(aspect) + len(rough_base_positions())


def release_geometry(path):
    """Measured active-bed dimensions before gate removal.

    ``L0`` is a physical control variable in both normalized runout and aspect
    ratio.  Therefore the recorded release fabric must actually span it; a
    count-only check cannot detect an empty gap at the gate.  The rough base is
    excluded by its fixed selection height.  We use particle envelopes, not
    centres, because those are the physical extents of a granular column.
    """
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError("empty release-state record")
    try:
        active = [(float(r["x"]), float(r["z"]), float(r["radius"]))
                  for r in rows if float(r["z"]) > BASE_SELECT_Z]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("malformed release-state record") from exc
    if not active:
        raise ValueError("release state has no active grains above rough base")
    top = max(z + radius for _, z, radius in active)
    h = top - BASE_Z
    if not math.isfinite(h) or h <= 0.0:
        raise ValueError("non-positive release height")
    left = min(x - radius for x, _, radius in active)
    right = max(x + radius for x, _, radius in active)
    width = right - left
    if not math.isfinite(width) or width < 0.95 * L0:
        raise ValueError(f"release width {width / L0:.3f} L0 is below 0.95 L0")
    if right > GATE_RELEASE_WIDTH_MAX:
        raise ValueError(
            f"release crossed the still-active gate: right envelope {right / L0:.3f} L0 "
            f"exceeds {GATE_RELEASE_WIDTH_MAX / L0:.3f} L0"
        )
    return h, width, left, right


def release_height(path):
    """Compatibility wrapper for callers that only need the measured height."""
    return release_geometry(path)[0]


def checked_release_dimensions(height, width, nominal_aspect):
    """Admit only a release that preserves the scheduled physical geometry.

    Lajeunesse et al.'s control parameter is the *released* H_i / L_i, and its
    normalized runout is (L_f - L_i) / L_i.  A gate location is a useful
    construction control, but is not a substitute for L_i after the source has
    settled.  In particular, a column that compacts away from the gate can pass
    a height-only check while being fitted at the wrong aspect and normalization.
    """
    if not math.isfinite(width) or width <= 0.0:
        raise ValueError("non-positive measured release width")
    width_error = abs(width - L0) / L0
    if width_error > ASPECT_REL_TOL:
        raise ValueError(
            f"release width {width:.6f} differs from scheduled {L0:.6f} "
            f"by {width_error:.2%} (limit {ASPECT_REL_TOL:.0%})"
        )
    actual = height / width
    if not math.isfinite(actual) or actual <= 0.0:
        raise ValueError("non-positive measured release aspect")
    relative_error = abs(actual - nominal_aspect) / nominal_aspect
    if relative_error > ASPECT_REL_TOL:
        raise ValueError(
            f"release aspect {actual:.6f} differs from scheduled {nominal_aspect:.6f} "
            f"by {relative_error:.2%} (limit {ASPECT_REL_TOL:.0%})"
        )
    return actual


def checked_release_aspect(height, nominal_aspect):
    """Legacy height-only helper retained for focused compatibility tests.

    Dynamic evidence must call :func:`checked_release_dimensions`; it is the
    only admission path used by campaign execution and graphing.
    """
    return checked_release_dimensions(height, L0, nominal_aspect)


def checked_final_state(path, expected_count):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 1 or set(rows[0]) != {"particle_count", "max_speed_m_s"}:
        raise ValueError("malformed terminal-state record")
    try:
        count = int(rows[0]["particle_count"])
        vmax = float(rows[0]["max_speed_m_s"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("non-numeric terminal-state record") from exc
    if count != expected_count or not math.isfinite(vmax) or vmax < 0.0:
        raise ValueError(f"invalid terminal state {count}/{expected_count}, vmax={vmax}")
    return vmax


def checked_arrest_window(path):
    """Return a sustained terminal-speed witness from the executable output.

    The recorder emits one row every fixed number of collapse steps.  Keeping
    only the last four samples makes the check insensitive to the early dynamic
    stage while refusing a deposit that merely happens to be slow at its final
    output instant.
    """
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if len(rows) < ARREST_WINDOW_SAMPLES:
        raise ValueError("too few arrest-window samples")
    window = rows[-ARREST_WINDOW_SAMPLES:]
    try:
        steps = [int(row["collapse_step"]) for row in window]
        speeds = [float(row["max_speed_m_s"]) for row in window]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("malformed arrest-window record") from exc
    if (any(step <= 0 for step in steps)
            or any(b - a != ARREST_SAMPLE_INTERVAL for a, b in zip(steps, steps[1:]))
            or not all(math.isfinite(v) and v >= 0.0 for v in speeds)):
        raise ValueError("invalid arrest-window values")
    return speeds


def checked_preparation_window(path, expected_count):
    """Require a sustained quiet tail while the gate is still present."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if len(rows) < PREPARATION_WINDOW_SAMPLES:
        raise ValueError("too few preparation-rest samples")
    window = rows[-PREPARATION_WINDOW_SAMPLES:]
    try:
        steps = [int(row["settle_step"]) for row in window]
        counts = [int(row["particle_count"]) for row in window]
        speeds = [float(row["max_speed_m_s"]) for row in window]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("malformed preparation-rest record") from exc
    if (steps[-1] != SETTLE_STEPS
            or any(b - a != PREPARATION_SAMPLE_INTERVAL for a, b in zip(steps, steps[1:]))
            or any(count != expected_count for count in counts)
            or not all(math.isfinite(v) and v >= 0.0 for v in speeds)):
        raise ValueError("invalid preparation-rest values")
    froude = max(speeds) / math.sqrt(9.81 * 2.0 * RADIUS)
    if froude > RELEASE_FROUDE_MAX:
        raise ValueError(f"preparation-rest Fr={froude:.6g} > {RELEASE_FROUDE_MAX}")
    return froude


def derive_dirt_ensemble():
    """Derive every fit input from the 11 x 3 executable witnesses.

    ``runout.csv`` is a convenience artifact, never primary evidence.  Graphing
    must therefore re-read every release, pre-release-rest, final, and terminal-state witness and
    reproduce each seed average before it is allowed to fit or write a figure.
    This makes a partial campaign and an edited summary equally inadmissible.
    """
    manifest_path = os.path.join(DATA_DIR, MANIFEST_NAME)
    if not os.path.isfile(manifest_path):
        raise ValueError("missing protocol manifest")
    try:
        with open(manifest_path) as f:
            recorded_manifest = json.load(f)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError("malformed protocol manifest") from exc
    current_manifest = checked_protocol_manifest(write=False)
    if recorded_manifest != current_manifest:
        raise ValueError("protocol manifest disagrees with current generated sources")

    rows = []
    failures = []
    for a in ASPECTS:
        lfs, hs, widths, rights, aspects = [], [], [], [], []
        for s in SEEDS:
            try:
                with shared_case_lock(a, s, "derive the ensemble"):
                    paths = dirt_case_paths(a, s)
                    deposit, release = paths["deposit"], paths["release"]
                    terminal, arrest, preparation = paths["terminal"], paths["arrest"], paths["preparation"]
                    expected = total_particles(a)
                    if not all(os.path.isfile(p) for p in (deposit, release, terminal, arrest, preparation)):
                        raise ValueError("missing release, pre-release, final, terminal, or arrest evidence")
                    checked_case_receipt(a, s)
                    release_n = csv_particle_count(release)
                    deposit_n = csv_particle_count(deposit)
                    vmax = checked_final_state(terminal, expected)
                    arrest_speeds = checked_arrest_window(arrest)
                    checked_preparation_window(preparation, expected)
                    h, width, _, right = release_geometry(release)
                    actual_aspect = checked_release_dimensions(h, width, a)
                    _, lf = measure_column(deposit)
            except (OSError, ValueError, csv.Error) as exc:
                failures.append(f"a={a} seed={s}: malformed witness ({exc})")
                continue
            if release_n != expected or deposit_n != expected:
                failures.append(f"a={a} seed={s}: population is not {expected} at release/final")
                continue
            froude = vmax / math.sqrt(9.81 * 2.0 * RADIUS)
            if froude > REST_FROUDE_MAX:
                failures.append(f"a={a} seed={s}: terminal Fr={froude:.6g} > {REST_FROUDE_MAX}")
                continue
            arrest_froude = max(arrest_speeds) / math.sqrt(9.81 * 2.0 * RADIUS)
            if arrest_froude > REST_FROUDE_MAX:
                failures.append(
                    f"a={a} seed={s}: final {ARREST_WINDOW_SAMPLES}-sample Fr={arrest_froude:.6g} "
                    f"> {REST_FROUDE_MAX}"
                )
                continue
            if not all(math.isfinite(v) for v in (h, lf)):
                failures.append(f"a={a} seed={s}: non-finite measured geometry")
                continue
            hs.append(h)
            lfs.append(lf)
            # Keep the per-realization measured geometry: averaging an H/L_i
            # surrogate would reintroduce the very substitution this gate
            # prevents.
            widths.append(width)
            rights.append(right)
            aspects.append(actual_aspect)
        if len(lfs) != len(SEEDS):
            continue
        rn = [(v - initial_right) / width
              for v, initial_right, width in zip(lfs, rights, widths)]
        rn_mean = sum(rn) / len(rn)
        rows.append({"nominal_aspect": a,
                     "L0": sum(widths) / len(widths), "H": sum(hs) / len(hs), "L_f": sum(lfs) / len(lfs),
                     "release_front": sum(rights) / len(rights),
                     "aspect": sum(aspects) / len(aspects),
                     "runout_norm": rn_mean,
                     "runout_std": (sum((v - rn_mean) ** 2 for v in rn) / len(rn)) ** 0.5,
                     "n_seeds": len(lfs), "protocol_sha256": protocol_fingerprint()})
    if failures or len(rows) != len(ASPECTS):
        detail = "; ".join(failures) if failures else "missing scheduled aspect"
        raise ValueError(f"incomplete or non-arrested 11x3 ensemble: {detail}")
    return rows


def _case_evidence_error(a, seed):
    """Return the reason a cached case cannot be reused, or ``None``.

    This is deliberately the same physical admission test used before fitting:
    a resumable campaign must not turn a merely present CSV into a completed
    realization.  In particular, source provenance is checked by the manifest
    at ``start``/``derive_dirt_ensemble`` and this per-case check verifies the
    released and terminal populations, measured release geometry, pre-release
    rest, deposit readability, and sustained Froude arrest.
    """
    paths = dirt_case_paths(a, seed)
    deposit, release = paths["deposit"], paths["release"]
    terminal, arrest, preparation = paths["terminal"], paths["arrest"], paths["preparation"]
    if not all(os.path.isfile(p) for p in (deposit, release, terminal, arrest, preparation)):
        return "missing release, pre-release, final, terminal, or arrest evidence"
    expected = total_particles(a)
    try:
        checked_case_receipt(a, seed)
        if csv_particle_count(release) != expected or csv_particle_count(deposit) != expected:
            return f"population is not {expected} at release/final"
        vmax = checked_final_state(terminal, expected)
        arrest_speeds = checked_arrest_window(arrest)
        checked_preparation_window(preparation, expected)
        height, width, _, _ = release_geometry(release)
        checked_release_dimensions(height, width, a)
        _, lf = measure_column(deposit)
        if not math.isfinite(lf):
            return "non-finite deposit measurement"
    except (OSError, ValueError, csv.Error) as exc:
        return f"malformed witness ({exc})"
    froude_scale = math.sqrt(9.81 * 2.0 * RADIUS)
    if vmax / froude_scale > REST_FROUDE_MAX:
        return f"terminal Fr={vmax / froude_scale:.6g} > {REST_FROUDE_MAX}"
    if max(arrest_speeds) / froude_scale > REST_FROUDE_MAX:
        return (f"final {ARREST_WINDOW_SAMPLES}-sample Fr="
                f"{max(arrest_speeds) / froude_scale:.6g} > {REST_FROUDE_MAX}")
    return None


def _clear_case_evidence(a, seed):
    """Remove every dynamic witness before a case is rerun.

    The preparation-rest record is as much a part of one realization as its
    release and terminal records.  In particular, retaining it across a forced
    rerun would let a newly generated trajectory borrow a quiet tail from an
    earlier preparation.  A receipt hashes files, not their causal provenance,
    so the correct boundary is to remove the entire witness set before the
    executable starts.
    """
    cdir = case_dir_seed(a, seed)
    for name in ("column_collapse_results.csv", "column_collapse_release.csv",
                 "column_collapse_final_state.csv", "column_collapse_arrest.csv",
                 "column_collapse_preparation.csv",
                 CASE_RECEIPT_NAME):
        stale = os.path.join(cdir, "data", name)
        if os.path.isfile(stale):
            os.remove(stale)
    log = os.path.join(cdir, "run.log")
    if os.path.isfile(log):
        os.remove(log)


def _run_dirt_case(a, seed, binary, env):
    """Run one independent witness after invalidating only that case's evidence."""
    with exclusive_case_lock(a, seed, "written"):
        cdir = case_dir_seed(a, seed)
        config = os.path.join(cdir, "config.toml")
        if not os.path.isfile(config):
            raise FileNotFoundError(f"missing {config} — run 'generate' first")
        _clear_case_evidence(a, seed)
        with open(os.path.join(cdir, "run.log"), "w") as log:
            subprocess.run([binary, config], cwd=REPO_ROOT, stdout=log,
                           stderr=subprocess.STDOUT, env=env, check=True)
        # Write only after every raw witness exists.  A later graph/reuse operation
        # re-hashes this receipt instead of trusting filenames or the aggregate CSV.
        write_case_receipt(a, seed, binary)
    return a, seed


def _start(jobs=1, rerun=False, selected_cases=None):
    """Launch either the full campaign or an explicit subset of its witnesses.

    A full-scale 11 x 3 collapse campaign is deliberately expensive.  The
    individual aspect/seed simulations have no data dependency, so an HPC or
    batch scheduler must be able to run one named witness without pretending
    that the partial result is a fitted ensemble.  A subset launch therefore
    writes only its raw executable evidence; it never writes ``runout.csv`` or
    declares PASS.  ``graph`` remains the sole path to aggregate all 33 raw
    witnesses and enforce the unchanged exponent gates.
    """
    if jobs < 1:
        raise ValueError("jobs must be at least one")
    os.makedirs(DATA_DIR, exist_ok=True)
    try:
        # Recompute immediately before invalidating old witnesses: a campaign
        # can only be labelled as this ensemble if every exact source file is
        # canonical at launch, not merely when ``generate`` once ran.
        checked_protocol_manifest(write=True)
    except ValueError as exc:
        print(f"ERROR: {exc}; refusing to launch non-canonical ensemble.")
        sys.exit(1)
    all_cases = [(a, s) for a in ASPECTS for s in SEEDS]
    if selected_cases is None:
        selected_cases = all_cases
    else:
        selected_cases = list(selected_cases)
        invalid = [(a, s) for a, s in selected_cases if a not in ASPECTS or s not in SEEDS]
        if invalid:
            raise ValueError(f"case(s) outside the declared 11x3 campaign: {invalid}")
        if len(set(selected_cases)) != len(selected_cases):
            raise ValueError("duplicate --case selection")
    print(f"Building {EXAMPLE} (release)...", flush=True)
    env = dict(os.environ)
    # macOS: ensure system libffi is found if the workspace needs it.
    subprocess.run(
        ["cargo", "build", "--release", "--example", EXAMPLE, "--no-default-features", "--features", "precision-double"],
        cwd=REPO_ROOT, check=True, env=env,
    )

    binary = os.path.join(REPO_ROOT, "target", "release", "examples", EXAMPLE)
    if not os.path.isfile(binary):
        raise RuntimeError(f"release binary missing after build: {binary}")
    reusable, cases = [], []
    for a, seed in selected_cases:
        # ``--rerun`` is an explicit request for a new dynamic witness.  It
        # must therefore bypass reuse *by making the cached case ineligible*,
        # not by assigning the successful ``None`` admission result.  The old
        # inversion silently reported a successful forced rerun after launching
        # zero simulations, which is especially unsafe for this expensive
        # validation campaign: users could believe raw evidence had been
        # refreshed when it was merely retained.
        reason = "forced rerun" if rerun else _case_evidence_error(a, seed)
        if reason is None:
            reusable.append((a, seed))
        else:
            cases.append((a, seed))
            if not rerun and os.path.isdir(case_dir_seed(a, seed)):
                print(f"  rerun a={a:<4} seed={seed}: {reason}", flush=True)
    if reusable:
        print(f"Reusing {len(reusable)} independently admitted DIRT case(s).", flush=True)
    print(f"Running {len(cases)} selected DIRT case(s) with {jobs} worker(s)...", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(_run_dirt_case, a, s, binary, env) for a, s in cases]
        for k, future in enumerate(concurrent.futures.as_completed(futures), 1):
            a, s = future.result()
            print(f"  [{k}/{len(cases)}] complete a={a:<4} seed={s}", flush=True)

    # A partial worker must return success when its own witness is qualified;
    # otherwise a batch scheduler cannot distinguish a sound individual run
    # from the intentionally incomplete global campaign.  It is still not a
    # validation result: only an unfiltered launch may write the aggregate CSV.
    if set(selected_cases) != set(all_cases):
        print("Subset complete: raw witnesses recorded only; run 'graph' after all 33 cases.")
        return

    try:
        rows = derive_dirt_ensemble()
    except ValueError as exc:
        print(f"\nERROR: {exc}; refusing to fit.")
        sys.exit(1)
    os.makedirs(DATA_DIR, exist_ok=True)
    _write_runout(RUNOUT_CSV, rows)
    print(f"\nDIRT:   wrote {len(rows)} seed-averaged runout rows "
          f"({len(SEEDS)} seeds/aspect; protocol {protocol_fingerprint()[:12]}) -> {RUNOUT_CSV}")

    # LAMMPS leg — optional cross-code overlay. Skipped entirely with no binary.
    lammps = find_lammps()
    if lammps:
        print(f"LAMMPS: {lammps} — running cross-code overlay.")
        if os.path.isfile(LAMMPS_CSV):
            os.remove(LAMMPS_CSV)
        lrows = run_lammps_sweep(lammps)
        if len(lrows) == len(ASPECTS):
            _write_runout(LAMMPS_CSV, lrows)
            print(f"LAMMPS: wrote {len(lrows)} runout rows -> {LAMMPS_CSV}")
        else:
            print("LAMMPS: incomplete independent campaign — refusing overlay.")
    else:
        print("LAMMPS: not found on PATH — running DIRT only.")


def start(jobs=1, rerun=False, selected_cases=None):
    """Run independent witnesses; duplicate case launches fail closed."""
    # Shared source admission lets multiple non-overlapping cases proceed while
    # preventing ``generate`` from replacing active/base coordinates underneath
    # an executable.  It is compatible with other shared campaign readers.
    with shared_campaign_lock("start"):
        return _start(jobs=jobs, rerun=rerun, selected_cases=selected_cases)


def _write_runout(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["nominal_aspect", "aspect", "L0", "H", "L_f", "runout_norm",
                        "release_front", "runout_std", "n_seeds", "protocol_sha256"],
            restval="", extrasaction="ignore",
        )
        w.writeheader()
        w.writerows(rows)


# ── deposit analysis ─────────────────────────────────────────────────────────
def measure_column(deposit_path):
    """Return (H_initial_estimate, L_f) from a settled deposit.

    H is estimated from the particle count and footprint. L_f is the far edge of
    the CONTIGUOUS deposit TOE, measured on a SUB-DIAMETER grid:

      * A grain contributes to the deposit toe only where its centre sits at least
        TOE_MIN_HEIGHT_D diameters above the floor (>= 2 grain layers). This is the
        original definition — "furthest x where the local deposit height exceeds
        one particle diameter" — and it rejects the single-grain leading smear and
        lone rolled/saltated grains (a grain resting on the floor has centre z = r
        < d, so it never counts).
      * Each qualifying grain is projected onto x as an occupied interval [x-r, x+r]
        on a fine grid of FINE_BINS cells per diameter, so the runout front is
        resolved to d/FINE_BINS rather than quantized to a full diameter (the old
        metric binned at one diameter and reported the bin edge, so whole aspect
        ratios landed in the same runout bin — a fit artifact, not physics).
      * Walking outward from the back wall, the toe is judged to end at the first
        CLEAR GAP wider than GAP_TOL_D diameters; the far edge of the last connected
        occupied cell is L_f.

    This is a pure measurement change: it preserves the original physical runout
    definition and only removes its diameter-scale quantization. No tolerance,
    target, or fit window is altered, and it is applied identically to the DIRT and
    LAMMPS deposits.
    """
    xs, zs, rs = [], [], []
    with open(deposit_path) as f:
        for r in csv.DictReader(f):
            z = float(r["z"])
            if z <= BASE_SELECT_Z:
                continue
            xs.append(float(r["x"]))
            zs.append(z)
            rs.append(float(r["radius"]))
    if not xs:
        return 0.0, L0

    r_mean = sum(rs) / len(rs)
    d = 2.0 * r_mean                   # mean particle diameter
    # Initial column height from solids volume in the L0 x W footprint.
    n = len(xs)
    vol = n * (4.0 / 3.0) * math.pi * r_mean ** 3
    h_init = vol / (PACKING * L0 * W)

    # Sub-diameter occupancy of the deposit TOE along x: only grains whose centre
    # is >= TOE_MIN_HEIGHT_D diameters above the floor (>= 2 layers) qualify.
    z_min = TOE_MIN_HEIGHT_D * d
    dx = d / FINE_BINS
    x_min = min(x - r for x, r in zip(xs, rs))
    x_max = max(x + r for x, r in zip(xs, rs))
    nb = int((x_max - x_min) / dx) + 2
    occ = [False] * nb
    for x, z, r in zip(xs, zs, rs):
        if z < z_min:
            continue
        b0 = max(0, int((x - r - x_min) / dx))
        b1 = min(nb - 1, int((x + r - x_min) / dx))
        for b in range(b0, b1 + 1):
            occ[b] = True

    # Far edge of the contiguous toe: walk outward, stop at the first vacuum run
    # longer than GAP_TOL_D diameters.
    gap_cells = max(1, int(round(GAP_TOL_D * d / dx)))
    front = 0
    run_empty = 0
    for b in range(nb):
        if occ[b]:
            front = b
            run_empty = 0
        else:
            run_empty += 1
            if run_empty > gap_cells:
                break
    lf = max(L0, x_min + (front + 1) * dx)
    return h_init, lf


# ── graph (validate + plot) ──────────────────────────────────────────────────
def load_runout():
    if not os.path.isfile(RUNOUT_CSV):
        print(f"ERROR: {RUNOUT_CSV} not found.")
        print("Run the sweep first: python3 examples/bench_column_collapse/sweep.py start")
        sys.exit(1)
    with open(RUNOUT_CSV) as f:
        rows = list(csv.DictReader(f))
    expected = {float(a) for a in ASPECTS}
    seen = set()
    required = {"nominal_aspect", "aspect", "L0", "H", "L_f", "runout_norm",
                "runout_std", "n_seeds", "protocol_sha256"}
    fingerprint = protocol_fingerprint()
    if len(rows) != len(ASPECTS):
        print("ERROR: runout CSV does not contain one row per scheduled aspect.")
        sys.exit(1)
    for row in rows:
        if not required.issubset(row):
            print("ERROR: runout CSV has an incomplete schema.")
            sys.exit(1)
        try:
            nominal = float(row["nominal_aspect"])
            aspect = float(row["aspect"])
            values = [float(row[k]) for k in ("L0", "H", "L_f", "runout_norm", "runout_std")]
            seeds = int(row["n_seeds"])
        except (TypeError, ValueError):
            print("ERROR: runout CSV contains non-numeric validation evidence.")
            sys.exit(1)
        if (nominal not in expected or nominal in seen or aspect <= 0.0
                or seeds != len(SEEDS) or not all(math.isfinite(v) for v in values)
                or values[1] <= 0.0 or row["protocol_sha256"] != fingerprint):
            print("ERROR: runout CSV is incomplete, duplicated, or inadmissible.")
            sys.exit(1)
        seen.add(nominal)
    if seen != expected:
        print("ERROR: runout CSV aspect set differs from the configured sweep.")
        sys.exit(1)
    # The summary is intentionally not accepted on its own.  Reconstruct the
    # full ensemble from all executable witnesses now, then require the cached
    # CSV to agree field-for-field.  Thus `graph` cannot turn a partial campaign
    # (or a manually edited CSV) into a PASS or an apparently fresh figure.
    try:
        witnessed = derive_dirt_ensemble()
    except ValueError as exc:
        print(f"ERROR: {exc}; refusing to graph summary-only evidence.")
        sys.exit(1)
    by_aspect = {float(r["nominal_aspect"]): r for r in rows}
    numeric = ("aspect", "L0", "H", "L_f", "runout_norm", "runout_std")
    for expected_row in witnessed:
        summary = by_aspect[float(expected_row["nominal_aspect"])]
        if (int(summary["n_seeds"]) != expected_row["n_seeds"]
                or summary["protocol_sha256"] != expected_row["protocol_sha256"]
                or any(not math.isclose(float(summary[key]), expected_row[key],
                                        rel_tol=0.0, abs_tol=1.0e-12)
                       for key in numeric)):
            print("ERROR: runout CSV disagrees with its complete raw ensemble.")
            sys.exit(1)
    return witnessed


def fit_loglog(pairs):
    """Least-squares slope of log(y) vs log(x). Returns (exponent, prefactor)."""
    pts = [(a, y) for a, y in pairs if a > 0 and y > 0]
    if len(pts) < 2:
        return float("nan"), float("nan")
    lx = [math.log(a) for a, _ in pts]
    ly = [math.log(y) for _, y in pts]
    n = len(pts)
    sx, sy = sum(lx), sum(ly)
    sxx = sum(v * v for v in lx)
    sxy = sum(a * b for a, b in zip(lx, ly))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        return float("nan"), float("nan")
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    return slope, math.exp(intercept)


def validate(rows):
    print("=" * 66)
    print("Granular Column-Collapse Runout Validation")
    print("=" * 66)
    print(f"  L0 = {L0*1000:.1f} mm, slab W = {W*1000:.1f} mm, d = {2*RADIUS*1000:.1f} mm")
    print(f"  E = {YOUNGS_MOD:.1e} Pa, e = {RESTITUTION}, mu = {FRICTION}")
    print(f"  External reference: {EMPIRICAL_REFERENCE['authors']}, "
          f"{EMPIRICAL_REFERENCE['journal']}, doi:{EMPIRICAL_REFERENCE['doi']}")
    print(f"  Reference geometry: {EMPIRICAL_REFERENCE['geometry']}\n")
    print(f"  {'a_nom':>5} {'a_rel':>6} {'H[mm]':>8} {'L_f[mm]':>9} {'(Lf-L0)/L0':>12} "
          f"{'seed_sd':>8} {'seeds':>6}")

    pairs = []
    for r in rows:
        a = float(r["aspect"])
        h = float(r["H"])
        lf = float(r["L_f"])
        rn = float(r["runout_norm"])
        sd = float(r["runout_std"]) if r.get("runout_std") not in (None, "") else float("nan")
        ns = r.get("n_seeds", "")
        pairs.append((a, rn))
        print(f"  {float(r['nominal_aspect']):>5.2f} {a:>6.3f} {h*1000:>8.2f} {lf*1000:>9.2f} {rn:>12.3f} "
              f"{sd:>8.3f} {str(ns):>6}")

    low = [(a, rn) for a, rn in pairs if a <= REGIME_SPLIT]
    high = [(a, rn) for a, rn in pairs if a >= REGIME_SPLIT]
    e_low, _ = fit_loglog(low)
    e_high, _ = fit_loglog(high)

    low_ok = abs(e_low - LINEAR_TARGET) <= EXP_TOL
    high_ok = abs(e_high - POWER_TARGET) <= EXP_TOL

    print()
    print(f"  Linear regime (a <= {REGIME_SPLIT}): fitted exponent = {e_low:.3f} "
          f"(target {LINEAR_TARGET:.2f})  [{'PASS' if low_ok else 'FAIL'}]")
    print(f"  Power regime  (a >= {REGIME_SPLIT}): fitted exponent = {e_high:.3f} "
          f"(target {POWER_TARGET:.2f})  [{'PASS' if high_ok else 'FAIL'}]")

    ok = low_ok and high_ok
    if not ok:
        print()
        print("  NOTE: this complete current campaign does NOT validate to tolerance.")
        print("  The unchanged target and ±0.25 band remain visible above.  Historical")
        print("  8d × 3d outputs are not used to explain, soften, or substitute for this")
        print("  32d × 10d rough-base result; see README/VALIDATION.md for scope.")
    print("\nALL CHECKS PASSED" if ok else "VALIDATION FAILED (see note above)")
    return ok


def compare_codes(dirt_rows, lammps_rows):
    """Print a per-aspect DIRT-vs-LAMMPS normalized-runout comparison and the
    fitted exponents for both codes."""
    dirt = {float(r["nominal_aspect"]): float(r["runout_norm"]) for r in dirt_rows}
    lammps = {float(r["nominal_aspect"]): float(r["runout_norm"]) for r in lammps_rows}
    print("\n" + "=" * 58)
    print("Normalized runout (L_f-L_i)/L_i: DIRT vs LAMMPS")
    print("=" * 58)
    print(f"  {'a':>5} | {'DIRT':>8} {'LAMMPS':>8} | {'diff':>8}")
    for a in sorted(set(dirt) & set(lammps)):
        d, l = dirt[a], lammps[a]
        print(f"  {a:>5.2f} | {d:>8.3f} {l:>8.3f} | {l - d:>+8.3f}")

    def fits(data):
        pairs = [(a, data[a]) for a in sorted(data)]
        low = [(a, v) for a, v in pairs if a <= REGIME_SPLIT]
        high = [(a, v) for a, v in pairs if a >= REGIME_SPLIT]
        return fit_loglog(low)[0], fit_loglog(high)[0]

    dl, dh = fits(dirt)
    ll, lh = fits(lammps)
    print("\n  Fitted exponents:        linear (a<=3)   power (a>=3)")
    print(f"    DIRT   :               {dl:>10.3f}    {dh:>10.3f}")
    print(f"    LAMMPS :               {ll:>10.3f}    {lh:>10.3f}")
    print(f"    targets:               {LINEAR_TARGET:>10.2f}    {POWER_TARGET:>10.3f}")


def plot(rows, lammps_rows=None):
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker
    except Exception as e:
        print(f"\n(matplotlib/numpy unavailable, skipped plots: {e})")
        return

    os.makedirs(PLOT_DIR, exist_ok=True)
    plt.rcParams.update({"font.size": 12, "figure.dpi": 150, "savefig.dpi": 150})

    a = np.array([float(r["aspect"]) for r in rows])
    rn = np.array([float(r["runout_norm"]) for r in rows])

    # ── Plot 1: normalized runout vs aspect ratio (log-log) with scaling lines.
    fig, ax = plt.subplots(figsize=(7, 5.2))
    ax.plot(a, rn, "o", color="#1f77b4", markersize=7, label="DIRT")
    if lammps_rows:
        la = np.array([float(r["aspect"]) for r in lammps_rows])
        lrn = np.array([float(r["runout_norm"]) for r in lammps_rows])
        ax.plot(la, lrn, "s", color="#d62728", markersize=8,
                markerfacecolor="none", markeredgewidth=1.6, label="LAMMPS")
    aa = np.logspace(math.log10(0.4), math.log10(6.0), 100)
    ax.plot(aa, 1.2 * aa, "k--", linewidth=1.4, label=r"$1.2\,a$ (linear, $a\lesssim3$)")
    ax.plot(aa, 1.6 * aa ** (2.0 / 3.0), "k:", linewidth=1.4,
            label=r"$1.6\,a^{2/3}$ (power, $a\gtrsim3$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.4, 6.0)
    ax.xaxis.set_major_locator(mticker.FixedLocator([0.5, 1.0, 2.0, 3.0, 4.0, 5.0]))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%g"))
    ax.xaxis.set_minor_locator(mticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_xlabel("Aspect ratio  a = H / L0")
    ax.set_ylabel(r"Normalized runout  $(L_f - L_0)/L_0$")
    ax.set_title("Column-Collapse Runout vs Aspect Ratio")
    low, _ = fit_loglog([(float(r["aspect"]), float(r["runout_norm"]))
                         for r in rows if float(r["aspect"]) <= REGIME_SPLIT])
    high, _ = fit_loglog([(float(r["aspect"]), float(r["runout_norm"]))
                          for r in rows if float(r["aspect"]) >= REGIME_SPLIT])
    criterion = (f"exponent gate: low {low:.3f} vs 1.000 ± {EXP_TOL:.2f}; "
                 f"high {high:.3f} vs 0.667 ± {EXP_TOL:.2f}")
    ax.text(0.02, 0.02, criterion, transform=ax.transAxes, fontsize=8,
            va="bottom", ha="left",
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.5"})
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    fig.savefig(os.path.join(PLOT_DIR, "runout_scaling.png"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {PLOT_DIR}/runout_scaling.png")

    # ── Plot 2: deposit-profile snapshot for the representative a = 2 case.
    target = min(rows, key=lambda r: abs(float(r["aspect"]) - 2.0))
    a_t = float(target["aspect"])
    deposit = os.path.join(case_dir(a_t), "data", "column_collapse_results.csv")
    if os.path.isfile(deposit):
        xs, zs = [], []
        with open(deposit) as f:
            for r in csv.DictReader(f):
                xs.append(float(r["x"]) * 1000)
                zs.append(float(r["z"]) * 1000)
        fig, ax = plt.subplots(figsize=(9, 3.2))
        ax.scatter(xs, zs, s=6, color="#ff7f0e")
        ax.axvline(L0 * 1000, color="0.5", linestyle="--", linewidth=1,
                   label=r"$L_0$")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("z [mm]")
        ax.set_title(f"Deposit profile (a = {a_t:g})")
        ax.set_aspect("equal")
        ax.legend(fontsize=9)
        fig.savefig(os.path.join(PLOT_DIR, "deposit_profile.png"), bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {PLOT_DIR}/deposit_profile.png")


def load_optional(path):
    """Load a runout CSV if it exists, else return []."""
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def derive_lammps_ensemble():
    """Reconstruct the optional cross-code overlay from raw LAMMPS witnesses.

    The LAMMPS CSV is a convenience cache, not independent evidence.  In
    particular, accepting it directly would make it possible to show a stale,
    partial, or hand-edited comparison beside a freshly re-derived DIRT result.
    This has no bearing on the experimental PASS/FAIL gate (which remains DIRT
    versus the published exponents); it only makes the optional external-code
    comparison auditable to the same standard as the DIRT campaign.
    """
    try:
        checked_protocol_manifest(write=False)
    except ValueError as exc:
        raise ValueError(f"LAMMPS source provenance is inadmissible: {exc}") from exc

    rows, failures = [], []
    froude_scale = math.sqrt(9.81 * 2.0 * RADIUS)
    for a in ASPECTS:
        lfs, hs, widths, rights, aspects = [], [], [], [], []
        for seed in SEEDS:
            cdir = case_dir_seed(a, seed)
            release_dump = lammps_dump_path(a, seed, "release")
            deposit_dump = lammps_dump_path(a, seed, "deposit")
            release = os.path.join(cdir, "lammps_release.csv")
            deposit = os.path.join(cdir, "lammps_deposit.csv")
            preparation = os.path.join(cdir, "lammps_preparation.txt")
            arrest = os.path.join(cdir, "lammps_arrest.txt")
            required = (release_dump, deposit_dump, release, deposit, preparation, arrest)
            if not all(os.path.isfile(p) for p in required):
                failures.append(f"a={a} seed={seed}: missing LAMMPS raw witness")
                continue
            try:
                checked_lammps_case_receipt(a, seed)
                expected = total_particles(a)
                if csv_particle_count(release) != expected or csv_particle_count(deposit) != expected:
                    raise ValueError(f"population is not {expected} at release/final")
                preparation_speeds = lammps_preparation_window(preparation)
                vmax = lammps_max_speed(deposit_dump)
                arrest_speeds = lammps_arrest_window(arrest)
                h, width, _, right = release_geometry(release)
                actual_aspect = checked_release_dimensions(h, width, a)
                _, lf = measure_column(deposit)
                if not all(math.isfinite(v) for v in (vmax, h, lf)):
                    raise ValueError("non-finite LAMMPS witness")
                if vmax / froude_scale > REST_FROUDE_MAX:
                    raise ValueError("terminal LAMMPS state is not arrested")
                if max(preparation_speeds) / froude_scale > RELEASE_FROUDE_MAX:
                    raise ValueError("still-gated LAMMPS preparation is not at rest")
                if max(arrest_speeds) / froude_scale > REST_FROUDE_MAX:
                    raise ValueError("final sustained LAMMPS state is not arrested")
            except (OSError, ValueError, csv.Error) as exc:
                failures.append(f"a={a} seed={seed}: {exc}")
                continue
            hs.append(h)
            lfs.append(lf)
            widths.append(width)
            rights.append(right)
            aspects.append(actual_aspect)
        if len(lfs) == len(SEEDS):
            rn = [(value - initial_right) / width
                  for value, initial_right, width in zip(lfs, rights, widths)]
            rn_mean = sum(rn) / len(rn)
            rows.append({"nominal_aspect": a, "aspect": sum(aspects) / len(aspects),
                         "L0": sum(widths) / len(widths), "H": sum(hs) / len(hs), "L_f": sum(lfs) / len(lfs),
                         "release_front": sum(rights) / len(rights),
                         "runout_norm": rn_mean,
                         "runout_std": (sum((value - rn_mean) ** 2 for value in rn) / len(rn)) ** 0.5,
                         "n_seeds": len(lfs), "protocol_sha256": protocol_fingerprint()})
    if failures or len(rows) != len(ASPECTS):
        detail = "; ".join(failures) if failures else "missing scheduled LAMMPS aspect"
        raise ValueError(f"incomplete LAMMPS 11x3 ensemble: {detail}")
    return rows


def load_verified_lammps():
    """Return an externally reproducible overlay, or no overlay at all."""
    cached = load_optional(LAMMPS_CSV)
    if not cached:
        return []
    try:
        witnessed = derive_lammps_ensemble()
    except ValueError as exc:
        print(f"LAMMPS overlay withheld: {exc}")
        return []
    by_aspect = {float(row["nominal_aspect"]): row for row in cached}
    numeric = ("aspect", "L0", "H", "L_f", "runout_norm", "runout_std")
    if len(by_aspect) != len(ASPECTS):
        print("LAMMPS overlay withheld: cached summary has an incomplete aspect set")
        return []
    for row in witnessed:
        cached_row = by_aspect.get(float(row["nominal_aspect"]))
        if cached_row is None:
            print("LAMMPS overlay withheld: cached summary omits a witnessed aspect")
            return []
        try:
            same = (int(cached_row["n_seeds"]) == row["n_seeds"]
                    and cached_row["protocol_sha256"] == row["protocol_sha256"]
                    and all(math.isclose(float(cached_row[key]), row[key], rel_tol=0.0, abs_tol=1.0e-12)
                            for key in numeric))
        except (KeyError, TypeError, ValueError):
            same = False
        if not same:
            print("LAMMPS overlay withheld: cached summary disagrees with raw witnesses")
            return []
    return witnessed


def graph():
    rows = load_runout()
    lammps_rows = load_verified_lammps()
    ok = validate(rows)            # DIRT-vs-theory only; LAMMPS never gates PASS.
    # Re-measure the raw witnesses in a separate program.  It neither imports
    # this driver nor accepts this aggregate CSV, so a benchmark PASS needs
    # corroboration rather than a self-referential check of one estimator.
    observer = os.path.join(SCRIPT_DIR, "independent_observer.py")
    observed = subprocess.run([sys.executable, observer], cwd=REPO_ROOT).returncode == 0
    if not observed:
        print("INDEPENDENT OBSERVER FAILED — refusing a benchmark PASS.")
    ok = ok and observed
    if lammps_rows:
        compare_codes(rows, lammps_rows)
    else:
        print(f"\n(no {os.path.basename(LAMMPS_CSV)} — plotting DIRT only)")
    plot(rows, lammps_rows)
    return ok


# ── dispatch ─────────────────────────────────────────────────────────────────
def main():
    args = sys.argv[1:]
    cmd = args[0] if args else "all"
    jobs = 1
    if "--jobs" in args:
        try:
            jobs = int(args[args.index("--jobs") + 1])
        except (IndexError, ValueError):
            print("ERROR: --jobs requires a positive integer")
            sys.exit(2)
    if jobs < 1:
        print("ERROR: --jobs requires a positive integer")
        sys.exit(2)
    selected_cases = None
    if "--case" in args:
        selected_cases = []
        case_positions = [i for i, arg in enumerate(args) if arg == "--case"]
        for i in case_positions:
            try:
                a_text, seed_text = args[i + 1].split(",", 1)
                a, seed = float(a_text), int(seed_text)
            except (IndexError, ValueError):
                print("ERROR: --case requires ASPECT,SEED (for example --case 2,0)")
                sys.exit(2)
            selected_cases.append((a, seed))
    if cmd == "generate":
        if selected_cases is not None:
            print("ERROR: --case is valid only with start")
            sys.exit(2)
        generate()
    elif cmd == "emit-jobs":
        if selected_cases is not None:
            print("ERROR: --case is valid only with start")
            sys.exit(2)
        emit_jobs()
    elif cmd == "status":
        if selected_cases is not None:
            print("ERROR: --case is valid only with start")
            sys.exit(2)
        # An incomplete or unprepared campaign is an honest, expected state;
        # make it visible to a scheduler without mislabelling it as a failed
        # physical validation.  Only malformed CLI use is an error here.
        print_campaign_status()
    elif cmd == "start":
        try:
            start(jobs, rerun="--rerun" in args, selected_cases=selected_cases)
        except ValueError as exc:
            print(f"ERROR: {exc}")
            sys.exit(2)
    elif cmd == "graph":
        if selected_cases is not None:
            print("ERROR: --case is valid only with start")
            sys.exit(2)
        sys.exit(0 if graph() else 1)
    elif cmd == "all":
        generate()
        start(jobs, rerun="--rerun" in args)
        print()
        sys.exit(0 if graph() else 1)
    else:
        print(f"Unknown command: {cmd!r}")
        print("Usage: sweep.py [generate|emit-jobs|status|start|graph] [--jobs N] [--rerun]   (no arg = all three)")
        sys.exit(2)


if __name__ == "__main__":
    main()
