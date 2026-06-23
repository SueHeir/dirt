# Step 1 — `GpuGranularResidentPlugin` (windowed-resident GPU schedule plugin)

The device-resident counterpart to milestone 1's host-authoritative GPU force
plugin. Instead of a host↔device round-trip every step, it keeps pos/vel/force/
omega **resident on the device** and advances the whole velocity-Verlet loop
(integrate + Hertz–Mindlin contact + planar walls + gravity + rotation) on the
GPU for a **window of K steps per schedule tick** via soil's `run_steps` /
`run_steps_continue`, syncing the host `Atom`/`DemAtom` only at window boundaries.

## How it slots into the schedule
- New file `crates/dirt_granular/src/gpu_resident.rs` — milestone 1's
  `gpu.rs` is byte-for-byte untouched (`git diff main` = 0), so the
  host-authoritative path stays the safe default for MPI / mixed-physics configs.
- It is added **instead of** `VelocityVerletPlugin` + host force +
  `RotationalDynamicsPlugin`: one system at the `Force` phase advances the device
  a whole window and writes back. The device is authoritative between syncs —
  host arrays are never re-uploaded (re-upload would re-prime the force and
  corrupt the contact history; the gate fix made windowing bit-exact).
- Walls are passed as a `Boundary` (planar) and gravity as a constant vector, so
  there's no dependency on host `WallPlugin`/`GravityPlugin`.

## Verification (`examples/resident_plugin_validate.rs`, Apple M5 Pro)
Wall+gravity drop, n=512, 4000 steps:

| run | wall-clock |
|---|---|
| resident plugin, window=500 | 423.8 ms |
| resident plugin, window=1 (host sync every step) | 24925.0 ms |

- **max \|plugin − direct `GpuState`\| = 0.000e0 (bit-identical).** The resident
  plugin reproduces a direct `GpuState` run set up identically, proving the
  schedule integration (resource handling, write-back to `Atom`) is faithful.
  That direct path is validated vs CPU (~1e-4) by `validate_trajectory` and is
  physically correct in `pile`, so by transitivity the resident plugin matches
  CPU within the f32 band.
- **End-to-end speedup (window=1 → 500): 58.8×.** This is the real residency win
  realized *through the schedule plugin*, not a kernel cost model.
- Milestone-1 regression `gpu_plugin_matches_cpu` still passes (~1e-8).

## Honest limitations (the step-2b / generalisation boundary)
- **Single-rank only.** The window boundary is exactly where MPI halo exchange
  would go — that's step 2b. No MPI yet.
- **No host fixes mid-window.** Anything that mutates pos/vel on the host between
  syncs (pin, deform, host-side custom fixes) can't run inside a window; it would
  need its own sync point. Single-rank wall+gravity granular has none.
- **Independent in-schedule CPU+walls comparison not built.** Correctness vs CPU
  is established transitively (plugin == direct `GpuState` == CPU). A fresh
  head-to-head against `GranularDefaultPlugins` *with walls* in one App is the
  remaining nice-to-have; the contact force itself is already CPU-validated.
- Walls are planar and material is single-type (the GPU kernel's scope).
