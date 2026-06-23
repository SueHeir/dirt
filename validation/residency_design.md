# Roadmap Step 1 — Single-GPU residency: design + first increment

**Goal:** keep the granular sim resident on the GPU — data stays on device, the
host touches it only for I/O/checkpoints and (later) MPI halo exchange — instead
of the per-step host↔device round-trip the schedule plugin does today.

## The two GPU paths (current state)

**1. Standalone resident stepper — `soil_gpu::GpuState::run_steps` (gpu_state.rs:278).**
One command encoder, one `submit`: it primes `F(x₀)`, then for each of K steps
does `p_init` (velocity half-kick + position drift + reseed `force = m·g`), the
aux half-kick (rotation), an on-device **cell-list rebuild** (`cell_list.record`),
the force hooks (`GranularForce`, `WallForce`), `p_final` (second half-kick), and
the second aux kick — all on-device, for all K steps, with **zero host
involvement** between them. `pile`/`validate_trajectory` use this. Single-rank
residency therefore already exists and is validated:
- `validate_trajectory`: resident trajectories match the CPU baseline to ~1e-4.
- This is the path the benchmark below measures.

**2. Schedule plugin — `dirt_granular::gpu_granular_force` (dirt_granular/src/gpu.rs:119).**
Host-authoritative force *offload* (the LAMMPS GPU-package model, its own docs say
so at gpu.rs:8–10). Per step: re-upload pos/vel/omega + cell-list rebuild
(`set_state`), `eval_force_once` (force only, no integrate), `download_force` +
`download_aux_rate`, then the **host** does Verlet integrate (`VelocityVerletPlugin`),
quaternion rotation (`RotationalDynamicsPlugin`), neighbor rebuild, pbc/exchange,
pin, deform. This is the path MPI composes with (milestone 1 made it ghost-aware),
and it pays a full host↔device round-trip — including two blocking GPU readbacks —
**every step**.

## First increment (this commit): quantify the residency win, prove batching exact

`crates/dirt_gpu/examples/residency_bench.rs` advances the *same* granular scene
K=500 steps two ways and measures wall-clock on the Apple M5 Pro:

| N | resident `run_steps(K)` | per-step-sync (plugin model) | speedup | batched vs stepwise |
|---|---|---|---|---|
| 8,000 | 39.9 ms | 1517 ms | **38×** | **0.0e0 (bit-exact)** |
| 64,000 | 172 ms | 2674 ms | **15.6×** | **0.0e0 (bit-exact)** |

- **per-step-sync** is the exact cost model of path 2: each step does `set_state`
  (upload + cell rebuild) + `eval_force_once` + `download_force` +
  `download_aux_rate` + host integrate. The two blocking readbacks per step stall
  the pipeline, so the gap is largest at small N (latency-bound) and shrinks as
  per-step compute grows (15.6× at 64k).
- **batched-vs-stepwise = 0.0** confirms `run_steps(K)` (one submit) is bit-identical
  to K×`run_steps(1)` — residency batching introduces no error.

Correctness of the resident integrator itself is already established
(validate_trajectory ~1e-4 vs CPU; the host-authoritative plugin matches CPU
hertz to 1.2e-8). So: resident is both **correct** and **15–38× cheaper** than
per-step sync. That is the motivation for moving the schedule path onto it.

## What full schedule residency requires (the remaining work)

The standalone path (1) bypasses the dirt schedule entirely, so it can't compose
with host-side fixes or MPI. Full residency = make the **schedule** path run K
resident steps between host-sync points. The natural residency window is one
**neighbor-rebuild interval**: between rebuilds the neighbor list and ghosts are
fixed (GpuState rebuilds its *own* cell list on-device each step anyway), so the
host need not touch state.

**Stays on host (sync points — window boundaries):**
- soil neighbor rebuild (displacement-triggered) — re-derives ghosts; force/cell
  topology can change → must re-sync host↔device.
- `pbc` / `exchange` (MPI migration) — moves atoms between ranks.
- `pin` / `deform` — fixes that mutate pos/vel mid-stream; any such fix forces a
  window boundary (or must itself move on-device).
- I/O / checkpoint / diagnostics that read host pos/vel.

**Moves on device (within a window):** force (already), velocity-Verlet integrate,
rotation aux-DOF, zero-force seed, gravity (`set_params`), wall force — all already
exist as `run_steps` does them; the work is *driving* them from the schedule.

**Refactor plan:**
1. A `GpuGranularResidentPlugin` (kept behind a switch; the milestone-1
   host-authoritative `GpuGranularForcePlugin` stays the safe default) that owns
   the integrate+force phases: instead of per-step `eval_force_once`, it holds a
   resident `GpuState` and calls `run_steps(K)` for the steps within a rebuild
   window, replacing the host `VelocityVerletPlugin`/`RotationalDynamicsPlugin`
   for those steps. The driver advances K timesteps per resident call.
2. Sync host pos/vel/omega ⇄ device only at window boundaries (rebuild / exchange
   / I-O), not per step. This is where the 15–38× lands in the real schedule.
3. The window boundary is exactly the **MPI halo-exchange point**, so this is the
   prerequisite for step 2b (GPU-resident halos): pack ghosts from the device
   buffers at the boundary, exchange (f64 wire), unpack into device buffers.

**Tension to resolve:** the schedule is per-step (one `app.run()` = one step),
but the resident stepper wants to own the K-step window. The plugin must either
advance K steps per `app.run()` (changing the driver's step accounting) or gate
the host integrate/force systems off while the resident `GpuState` is authoritative
between rebuilds. Either is a real but bounded change; correctness is checked the
same way as here (resident vs CPU baseline within the f32 band, and decomposition-
invariance once MPI is in).

## Honest boundary

This commit delivers: the design above + a faithful, bit-exact-verified
measurement that residency is 15–38× cheaper than per-step sync. It does **not**
yet wire the resident model into the dirt schedule — that is the
`GpuGranularResidentPlugin` refactor (step 1 proper), which then unlocks
GPU-resident halos (step 2b). The standalone resident path remains validated and
usable today via `GpuState` directly (`pile`, `validate_trajectory`,
`residency_bench`).
