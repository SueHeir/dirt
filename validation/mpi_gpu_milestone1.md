# Milestone 1 — MPI domain decomposition + GPU contact force (correctness)

**Result: a 2-rank MPI run with the GPU contact force now reproduces the 1-rank
and CPU results to the f32 noise floor.** The scale-out architecture (MPI
domain-decomp + GPU per rank) is functionally correct for particle–particle
contact.

## The bug that was fixed

`gpu_granular_force` / `GpuGranular::build` (crates/dirt_granular/src/gpu.rs)
uploaded and binned only **local** atoms (`0..nlocal`), ignoring **ghost** atoms
(`nlocal..atoms.len()`). The CPU `hertz_mindlin_contact_force` uses the neighbour
list, which includes ghost pairs. So under MPI a local atom contacting a
neighbour-rank atom (delivered as a ghost by border/forward comm) had that
contact **silently dropped** → wrong forces at every subdomain boundary.

## The fix (host-authoritative, no residency work)

Upload **local + ghost** (`nall = atoms.len()`), build the GPU cell list over all
of them, evaluate the i-centric / no-newton kernel (each atom's force is complete
from its own neighbour loop), and accumulate force/torque back onto **locals
only** (`0..nlocal`); ghost forces are redundant and discarded — no reverse-comm.

**Contact history:** the on-device tangential-spring history is keyed by atom
slot. Ghost *count* is constant between neighbour rebuilds (forward_comm only
moves ghosts), so keying the GPU-state rebuild on `nall` triggers it only on
neighbour rebuilds — the same cadence at which the host re-establishes its lists.
Local-atom history therefore persists between rebuilds; no soil_gpu change needed.

## Verification (dense settling bed, 600→474 grains, walled box, 12000 steps)

Particle–particle contact on CPU (`GranularDefaultPlugins`) vs GPU
(`GranularGpuPlugins`); walls + gravity are CPU in both, isolating the contact
kernel. The bed packs straddling x = 0.01, the n=2 decomposition split, creating
many cross-boundary contacts. Metric = max per-atom position diff, matched by
global tag. GPU is f32 → compared at the f32 trajectory level (~1e-3).

| comparison | atoms | max pos diff | verdict |
|---|---|---|---|
| CPU n1 vs CPU n2 (harness check) | 474 | 1.4e-11 | exact — CPU MPI is correct |
| **FIXED GPU n2 vs CPU n2** | **474** | **1.5e-3** | **PASS — f32 floor** |
| **FIXED GPU n2 vs GPU n1** | **474** | **1.8e-3** | **PASS — decomposition-invariant** |
| FIXED GPU n1 vs CPU n1 (f32 baseline) | 474 | 1.4e-3 | f32-vs-f64 reference |
| UNFIXED GPU n2 vs CPU n2 | 337 (!) | 2.5e-2 | FAIL — drops boundary contacts |
| UNFIXED GPU n1 vs CPU n1 (no ghosts) | 474 | 1.4e-3 | fine — bug is MPI-boundary only |

The fix is load-bearing: without it, 2-rank GPU loses 137 boundary grains (they
escape with no cross-rank support) and the survivors are **16× outside** the f32
band; with it, all atoms are conserved and the result sits at the f32 noise floor
— the GPU decomposition adds no error beyond f32 itself. The unfixed n=1 run
(no MPI ghosts) is fine, confirming the bug was purely the boundary ghost-drop.

## Notes / caveats

- On this Mac both ranks share the one Metal GPU — this is a **correctness** test,
  not a performance test. Per-rank GPU binding (roadmap step 3) needs multi-GPU
  hardware to show speedup.
- The config uses gz = −40 to settle quickly; ~21% of grains tunnel the floor
  wall under that strong gravity. This loss is **identical** across CPU and GPU
  (deterministic), so the comparison is fair; it's a config-tuning artifact, not
  a correctness issue. A gentler g + longer run would conserve all grains.
- Reproduce: `cargo build --release --example mpi_gpu_validate`, then
  `DIRT_FORCE={cpu,gpu} mpiexec -n {1,2} target/release/examples/mpi_gpu_validate
  examples/mpi_gpu_validate/config_n{1,2}.toml`.
