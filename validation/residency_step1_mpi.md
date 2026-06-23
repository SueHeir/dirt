# Step 1 × MPI: 2-rank resident GPU (correctness)

`GpuGranularResidentMpiPlugin` (crates/dirt_granular/src/gpu_resident_mpi.rs) fuses
milestone-1's ghost-awareness with step-1 residency: each schedule tick it uploads
**local+ghost** (so cross-rank boundary contacts aren't dropped), advances `window`
velocity-Verlet steps on-device, and writes **locals** back. Verified on a dense
settling bed straddling the x=0.01 MPI split (Apple M5 Pro, both ranks share the one
GPU → correctness only, no speedup).

| run | atoms | lost | max|Δpos| vs 1-rank |
|---|---|---|---|
| n=2, **window=1** | 474 | 0 | **0.000e0 — bit-exact, decomposition-invariant** |
| n=2, window>1 | — | — | diverges (see below) |

**Result:** at window=1 the resident plugin reproduces the single-rank run
**bit-for-bit** under 2-rank MPI — the ghost-aware contact + window-boundary
download/exchange/re-upload is correct. This is "2 MPI ranks + GPU, resident,
correct."

**The window>1 boundary (step-2b seam).** `forward_comm_borders` refreshes ghost
positions **every step** (soil_core comm.rs:186, `CommState::CommunicateOnly`), but
this plugin refreshes them once per *tick* (= per `window` device steps). So for
window>1 the device's ghosts are frozen for `window-1` steps while their owner rank
advances them → boundary-local contacts see a stale ghost trajectory → divergence
(and, empirically, instability/atom churn that slows the run). Correct residency
under MPI therefore requires **GPU-resident halos (step 2b)**: pack/unpack
`forward_comm` from device buffers so ghosts refresh on-device each step without a
host round-trip. That is multi-GPU-perf-bound and the next roadmap item.

**Unchanged:** milestone-1 (`gpu_plugin_matches_cpu`) and the single-rank resident
plugin remain the safe defaults; this is an additive plugin.
