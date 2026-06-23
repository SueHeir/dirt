# GPU Granular Plugin — production integration design

Status: **blueprint** (core GPU rewire is complete & validated; this is the
remaining production-integration phase). Nothing here is built yet.

## What's already done (the substrate this builds on)

- `soil_gpu`: generic GPU substrate — `GpuContext`, `CellList`, `DualBuffer`
  (Kokkos-DualView lazy host/device coherence), `NeighborSlots` (ping-pong
  per-neighbor state), `Boundary`/`Plane`, and `GpuState` — the resident
  velocity-Verlet loop with a `GpuForce` hook API + auxiliary-DOF integration
  (generic rotation: state=ω, rate=τ, inv_coeff=1/inertia).
- `dirt_gpu`: DEM constitutive physics as hooks — `GranularForce`
  (Hertz-normal + Mindlin-tangential + spring history + torque) and `WallForce`
  (sphere-plane HM). Validated vs **real dirt** at ~1e-6 (compare_cpu Tier 1),
  faster than the old monolith (compare_cpu Tier 3, ~250× CPU at 216k).

The remaining gap: nothing wires `GpuState`+hooks into the `grass_app`
plugin/schedule system, so real sims can't use it yet.

## The schedule (from soil_core `ParticleSimScheduleSet`)

Per timestep, in order, the CPU granular path runs:

| Phase | CPU systems |
|---|---|
| PreInitialIntegration | `wall_move` |
| InitialIntegration | `soil_verlet::initial_integration`, `dirt_granular::rotational::initial_rotation` |
| Exchange / Neighbor | (MPI migrate, neighbor rebuild) |
| PreForce | `wall_zero_force_accumulators` |
| **Force** | `hertz_mindlin_contact_force`, `wall_contact_force` |
| FinalIntegration | `soil_verlet::final_integration`, `dirt_granular::rotational::final_rotation` |

CPU plugin group: `dirt_granular::GranularDefaultPlugins` (+ `WallPlugin`).

## Decision: force-offload first (LAMMPS GPU-package model), not full KOKKOS residency

Two models were considered:

- **(A) Force-offload** — GPU computes only the **Force** phase; CPU keeps
  integration, quaternion rotation, comm, fixes. Mirrors the existing
  `soil_gpu::VelocityVerletGpuPlugin` (host-authoritative, per-phase sync) and
  the LAMMPS *GPU package*. CPU fixes "just work" (they read/write host arrays).
- **(B) Full residency + coherence** — whole step resident on device, sync to
  host only when a CPU system reads (LAMMPS *KOKKOS*). Fastest, but the schedule
  has CPU systems interleaved between phases (wall_move, exchange, fixes) that
  would each force a sync; correct interleaving is intricate.

**Choose (A) first.** It is correct, composes with every CPU fix, reuses
`GpuState::eval_force_once()`, and offloads the expensive phase. (B) is a later
optimization layered on `DualBuffer` once (A) is proven in real sims.

### Why force-offload is clean here

`GpuState::eval_force_once()` already does exactly the Force phase: build cell
list at current positions, seed (gravity off — integration owns gravity in this
model), run the hooks (contact + wall), leaving device `force` and aux `rate`
(=torque). The contact history (`NeighborSlots`, on device, atom-index-keyed)
persists across calls — no host round-trip for springs.

## Plugin shape (`dirt_granular`, new module e.g. `gpu.rs`)

```
pub struct GpuGranularForcePlugin;   // drop-in for the Force-phase DEM systems

struct GpuGranular {                 // App resource
    ctx: Option<GpuContext>,
    state: Option<GpuState>,         // lazily built when atoms exist
    omega_aux: usize,
    grid: Grid,
    n: usize,
}

// system @ Force phase:
fn gpu_granular_force(
    mut atoms: ResMut<Atom>,
    mut registry: ResMut<AtomDataRegistry>,   // DemAtom: omega(read), torque(write), inv_inertia, radius
    walls: Option<Res<Walls>>,
    mt: Res<MaterialTable>,
    mut g: ResMut<GpuGranular>,
) {
    // 1. lazy (re)build GpuState+GranularForce(+WallForce) on first call / count change
    // 2. upload host pos/vel (Atom) + omega (DemAtom) to GpuState
    //    (gravity OFF: set_params dt, [0,0,0] — CPU integration applies gravity,
    //     OR seed gravity here and skip CPU gravity; pick one, document it)
    // 3. gs.eval_force_once()
    // 4. download force -> atoms.force, aux_rate -> dem.torque
}
```

`build()` registers `gpu_granular_force` at `ParticleSimScheduleSet::Force`
INSTEAD of `hertz_mindlin_contact_force` + `wall_contact_force`. Keep the CPU
integration + rotation plugins unchanged. Falls back to the CPU force systems if
`GpuContext::new()` is `None`.

## The subtleties to handle (and how)

1. **Physics subset.** GPU hooks do plain Hertz-Mindlin only. Assert / warn if
   the `MaterialTable` has nonzero rolling/twisting friction, cohesion, or
   surface energy (the GPU path silently ignores them). Document: GPU plugin is
   for plain-HM material configs.
2. **Grid / domain.** Positions move each step. For the common fixed-box DEM sim
   (walls form a box) the grid is stable — build it once from the box/initial
   AABB with a margin. If a particle leaves the grid, binning clamps (contacts
   near the moved-out region would be missed). v1: build grid from the `Domain`
   bounds (a `Res<Domain>`); rebuild GpuState only if cell count changes —
   **but** a rebuild drops the device contact history (NeighborSlots), so warn /
   minimize. Prefer a grid sized to the domain so it never needs rebuilding.
   Remember: `soil_gpu::Grid::from_positions` takes the **literal** cutoff
   (= 2·r_max for monodisperse); pass `2*r`, not `r`.
3. **Contact history across rebuilds.** Keep `GpuState` alive for the whole run
   (don't rebuild per step). Only pos/vel/omega upload each step; springs stay
   resident. Atom count change (insertion) → rebuild + history reset (acceptable;
   matches a neighbor-list rebuild).
4. **Gravity ownership.** Either CPU applies gravity (a `GravityPlugin` system at
   Force, runs before/after GPU) or GPU seeds it. If GPU seeds gravity in
   `eval_force_once`, a CPU `GravityPlugin` must NOT also add it. Simplest: GPU
   path owns gravity (set_params with real g; eval seeds m*g) and the plugin
   group omits the CPU GravityPlugin. Document.
5. **inv_mass / inv_inertia / radius** come from `Atom` (mass) + `DemAtom`
   (inv_inertia, radius). Build the GPU radius/inv_mass/inv_inertia arrays from
   those on (re)build.
6. **Walls** from `Res<Walls>`: convert `WallPlane{point,normal}` → soil
   `Plane::new`. Rebuild `WallForce` if walls move (or upload wall uniform each
   step — cheap; better: a `WallForce::set_walls` to update the uniform without
   rebuild, since `wall_move` changes planes each step).
   → **add `WallForce::set_walls(&Boundary)`** (update uniform + n_walls) so
   moving walls don't force a rebuild.

## Validation

A plugin-level test mirroring `soil_gpu::plugin::gpu_plugin_matches_cpu_schedule`:
assemble a small granular App with `GpuGranularForcePlugin` and another with the
CPU force systems (plain-HM material, a floor wall), run N steps, assert
positions match within f32 tol (~1e-3). Plus keep `compare_cpu` as the
per-evaluation ground truth.

## Later: (B) full residency + coherence

Wrap device state in `DualBuffer`; CPU systems that read `Atom`/`DemAtom` fields
trip `ensure_host` (sync only those fields, only at that boundary); an all-GPU
chain runs zero syncs. This is the "no surface difference between CPU/GPU modes,
all-GPU is faster" goal. Build it on top of (A) once (A) is validated in a real
sim.
