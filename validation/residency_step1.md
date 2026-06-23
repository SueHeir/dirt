# Step 1 (residency) — schedule wiring attempt + the real blocker found

**Outcome: I did NOT land a working `GpuGranularResidentPlugin`. I found why one
can't be correct yet:** windowing the resident GPU stepper corrupts the
stateful tangential contact history at every window boundary. This is the
concrete blocker that has to be fixed in `soil_gpu` before step-1 (and therefore
step-2b MPI halos) can be correct. A reproducing example is committed.

## What residency requires, and why it breaks

Single-GPU residency = run the sim on-device and only sync to the host at
**window boundaries** (neighbor rebuild, MPI halo exchange, I/O). The standalone
resident stepper `GpuState::run_steps(K)` does K full velocity-Verlet steps in
one submit and is the path `pile` / `validate_trajectory` use — validated ~1e-4
vs CPU **as a single uninterrupted window**.

But `run_steps` **primes the force at entry** (cell-list rebuild + `GranularForce`
hook) before its step loop, and the `GranularForce` hook advances the **Mindlin
tangential spring history** every time it runs. So splitting a run into windows —
`run_steps(K)` per window, which is exactly the residency model — re-primes (re-
advances history) at each boundary, double-counting the tangential history there.

## Reproduction (`crates/dirt_gpu/examples/resident_validate.rs`)

Same deterministic wall+gravity frictional pile (n=1000, μ=0.5), advanced 4000
steps three ways on the Apple M5 Pro:

| path | how | wall-clock | max ‖pos − A‖ |
|---|---|---|---|
| **A** single window | `run_steps(4000)` once (1 prime) | 284 ms | — (reference) |
| **B** windowed | `run_steps(250)`×16 (16 primes) | 274 ms | **4.9e-2** |
| **C** per-step | `run_steps(1)`×4000 (+ readback) | 15765 ms | 4.4e-2 |

- **A vs B = 4.9e-2** (≈ a full grain radius): windowing alone — *without any host
  round-trip* — changes the answer. The only difference between A and B is the
  number of `run_steps` calls, i.e. the number of force re-primes. This isolates
  the re-prime/history bug.
- The residency *speed* win is real and large (A is **55× faster** than the
  per-step path C), but it's unusable across rebuild/MPI/I-O boundaries until the
  re-prime issue is fixed.

## The fix (the actual step-1 prerequisite, soil_gpu)

`run_steps` must separate "compute the force needed for the next half-kick" from
"advance the tangential history." Options:
1. A history-neutral prime: the entry force eval computes contact force without
   incrementing the Mindlin spring (the spring is advanced only inside the step
   loop, once per step).
2. A `run_steps_continue(K)` that skips the entry prime and trusts the force
   buffer left valid by the previous window's last step — so stitched windows
   produce the same trajectory as one window.

Either makes windowing idempotent (A == B), which is the correctness gate for a
`GpuGranularResidentPlugin` and for GPU-resident MPI halos (step 2b), whose sync
points ARE window boundaries.

## What is NOT affected

- **Milestone 1** (host-authoritative ghost-aware GPU force, on `main`) is
  untouched and still correct: it uses `eval_force_once` (one force eval per
  step) + host integrate, so it advances history exactly once per step. This
  step-1 work added only a new example; it changed neither the milestone-1 plugin
  nor soil.
- The standalone single-window resident path (`pile`, `validate_trajectory`)
  remains valid — it just can't be chopped into host-sync windows yet.

## Honest boundary

Delivered: the precise blocker + a reproducing benchmark + the soil_gpu fix
needed. NOT delivered: a working resident schedule plugin (it would be incorrect
until the re-prime fix lands). Recommend the soil_gpu `run_steps` history-neutral
prime as the next concrete step before retrying the plugin.

---

## RESOLUTION (gate closed — bit-exact)

The earlier diagnosis (above) was **wrong**: re-priming history was a minor
effect, not the cause. The real divergence had three parts, fixed as follows.

1. **Non-deterministic cell list (the dominant cause).** soil_gpu's cell-list
   scatter used `atomicAdd` to place atoms in cells, so within-cell neighbour
   order was a race — non-deterministic across command submits. A chaotic
   frictional pile amplified that f32 summation-order noise to ~5e-2 over 4000
   steps. Fix: a deterministic per-cell sort (`sort_cells` in cell_list.wgsl).
   The sim is now bit-reproducible run-to-run.

2. **Re-prime velocity mismatch.** Each `run_steps` re-evaluates the force at
   entry using the *full* end-of-step velocity, but the integrator's force uses
   the *mid-step* half-kick velocity, so re-priming at a window boundary diverges
   the velocity-dependent contact damping. Fix: `run_steps_continue(K)` continues
   without re-priming, trusting the resident force buffer.

3. **Test artifact.** `resident_validate`'s `reset()` re-uploaded pos/vel but
   never cleared on-device contact history, so paths B/C started from path A's
   leftover springs. Fix: build a fresh GpuState + hooks per path.

**Result (n=1000 frictional pile, 4000 steps, M5 Pro):** single-window vs
windowed `run_steps(250)+run_steps_continue(250)×15` → **max |Δpos| = 0.000e0
(bit-exact)**, 54× faster than per-step. Per-step (re-prime every step) still
diverges 4.6e-2 — confirming residency-continue is the correct model.
Milestone-1 unit test and validate_trajectory still pass.

**Residency windowing is therefore correct and unblocked.** The remaining step-1
work — the `GpuGranularResidentPlugin` schedule wiring — can now be built on
`run_steps`/`run_steps_continue` at the neighbour-rebuild cadence.

Requires the soil_gpu commit (deterministic cell list + run_steps_continue);
soil `main` must be pushed before this dirt branch builds against it.
