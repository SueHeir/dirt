# Clumps (Multisphere)

A **clump** is a rigid body built from several overlapping spheres — the
multisphere recipe for non-spherical particles. Each sub-sphere takes part in
contact detection, but the forces it feels are aggregated to the parent body,
which integrates rigid-body dynamics; the sub-sphere positions and velocities
are then derived back from the body state. Clumps live in the `dirt_clump`
crate; add the clump plugin to enable them.

## Clumps vs. bonds

Reach for a **clump** when you want a *permanently rigid* non-spherical
particle — a dimer, a rod, an angular grain — that never deforms. Reach for
[bonds](bonds.md) when the composite must *flex, yield, or break* (fibers,
agglomerates, cohesive solids). A clump is rigid by construction; a bonded
cluster is elastic and breakable.

## Architecture

- **No phantom parent atom.** Body data lives in a separate body store, not in
  the atom arrays. Each sub-sphere atom references its body through a `body_id`.
- **Full inertia tensor.** Built analytically (parallel-axis theorem) for
  non-overlapping spheres, or by Monte Carlo for overlapping ones, then
  diagonalized to principal moments plus an axes quaternion.
- **Euler-equation integration.** Torque is mapped into the principal frame,
  the Euler angular acceleration is computed there, and `ω` is half-kicked.

## Contact-exclusion: `same_body`

Sub-spheres of the *same* body must never push on each other — they are held
rigid, so an intra-body "contact" would be a spurious internal force. The crate
does **not** silently filter the neighbor list; instead it exposes a `same_body`
predicate, and every contact-force consumer is responsible for calling it to
skip same-body pairs:

```rust
// inside a pair loop over (i, j) from the neighbor list:
if same_body(&clump_data, i, j) {
    continue; // sub-spheres of one rigid body do not interact
}
```

`same_body` returns `true` only when both atoms carry a non-zero `body_id` and
those ids match. The granular contact kernels already honor this contract.

> **Every force plugin must call `same_body`.** The crate does not filter the
> neighbor list, so the exclusion is the *force plugin's* responsibility. A
> custom force plugin that iterates the neighbor list and forgets this check
> will let the sub-spheres of a single rigid body push on one another — the body
> self-repels and explodes apart on the first step.

## Rigid-body integration

Each body is advanced with a velocity-Verlet scheme following LIGGGHTS's
`FixMultisphere`. Two design choices are worth calling out:

- **Angular momentum is the integrated state.** The body stores `angmom`
  (space-frame angular momentum) as the primary rotational variable, half-kicked
  by the torque each step (`L += ½ dt τ`). Angular velocity `ω` is *derived*
  from `angmom` on demand, never integrated directly. This keeps the update
  stable for the asymmetric inertia tensors of non-spherical clumps, where
  integrating `ω` through Euler's equation directly is fragile.
- **Richardson-extrapolated quaternion update.** The orientation quaternion is
  advanced by taking one full step and two half-steps (recomputing `ω` at the
  half-step orientation) and extrapolating `q ← 2 q_half − q_full` for
  second-order accuracy, then renormalizing.

### Two quaternions

A body carries **two** orientation quaternions, and the distinction matters:

- `quaternion` — **body → space**: the live orientation of the body frame in the
  world, updated every step by the Richardson integrator.
- `principal_axes` — **body → principal**: a *fixed* rotation (set at creation)
  from the body frame to the frame in which the inertia tensor is diagonal
  (`principal_moments`).

Inertia work is only simple in the principal frame, so wherever a
principal-frame mapping is needed the two are composed (`q_principal→space =
quaternion * principal_axes`, a Hamilton product). That composed quaternion is
what the angular-velocity helper expects: it rotates `angmom` into the principal
frame, divides by `principal_moments`, and rotates the resulting `ω` back to
space.

## Inertia computation

At body creation the inertia tensor is built one of two ways, auto-selected by
whether the spheres overlap:

- **Non-overlapping spheres** → exact, parallel-axis theorem.
- **Overlapping spheres** → Monte Carlo over the bounding box with a hardcoded
  **100 000 samples**. This double-counts overlap volume correctly but
  introduces ~5 % stochastic noise in the resulting moments at this sample
  count — the price of handling arbitrary overlapping geometry.

Both paths are then diagonalized into `principal_moments` + `principal_axes`.

> **Overlapping-sphere inertia is non-reproducible.** The Monte Carlo path draws
> its 100 000 samples from a *freshly-seeded, unseeded* RNG at each clump's
> creation. The ~5 % noise it leaves in the principal moments is therefore
> different on every run — two runs of the *same* config produce slightly
> different rigid-body inertias. For any overlapping geometry (which includes
> most realistic clumps, e.g. the `sphere7` benchmark) treat the per-body inertia
> as a stochastic quantity, not a fixed number.

> A scalar helper `compute_clump_inertia` is **legacy** — it returns only a
> single averaged moment (the trace ÷ 3) and is kept for backward compatibility.
> New code should use the full-tensor functions.

## Configuration

Clump definitions and insertion live under the `[clump]` TOML section (separate
from `[dem]`). A definition names the body and lists its sub-spheres by offset
from the center of mass; an insert block stamps copies of a named definition
into a region.

```toml
[[clump.definitions]]
name = "dimer"
spheres = [
    { offset = [-0.0003, 0.0, 0.0], radius = 0.001 },
    { offset = [ 0.0003, 0.0, 0.0], radius = 0.001 },
]

[[clump.insert]]
definition = "dimer"
count = 100
density = 2500.0
material = "glass"
velocity = 0.5
region = { type = "block", min = [0.001, 0.001, 0.001], max = [0.019, 0.019, 0.019] }
```

Each definition's `spheres` give body-frame offsets from the center of mass and a
radius. An insert block stamps `count` copies of a named definition into `region`
(default: the domain inset by the **effective clump radius**,
`max over spheres of |offset| + radius`); overlap checks use a 5 % margin. If
`velocity` is set, each component is drawn uniform in `[−v, +v]`.

> **Insertion is rank-0-only.** The clump insert system materializes all bodies
> and sub-spheres on rank 0, then lets the normal MPI exchange distribute them at
> step 0. If you insert clumps programmatically (calling `insert_clump` directly
> from your own setup code) you must follow the same rank-0-only discipline, or
> every rank will create duplicate atoms.

**Neighbor bin size.** A clump's sub-spheres can straddle several neighbor bins,
so set `neighbor.bin_size` well above the clump bounding radius (the `sphere7`
benchmark uses `bin_size = 0.004` for `eff_radius ≈ 0.0011 m`, ~4× the sub-sphere
diameter) to avoid missing contacts at bin boundaries.

The `clump_dimer_drop` example is a minimal end-to-end run; the
`bench_clump_haff_cooling` benchmark validates the rigid-body rotational
dynamics (see the [Validation chapter](../reference/validation.md)).
