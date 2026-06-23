# Roadmap step 4 — interior/boundary compute–comm overlap

**Goal.** Compute interior contact forces *while the ghost halo is in flight*, then
finish the boundary forces once it lands — hiding MPI communication latency behind
useful compute. The split must be exact: the overlapped force has to reproduce the
standard single-pass force bit-for-bit.

## The split

A neighbour pair `(i, j)` is **boundary** iff `j` is a ghost (`j >= nlocal`): those
pairs read positions that the halo exchange refreshes, so they must run *after* the
halo lands. Every other pair is **interior** (`j < nlocal`) and needs no ghosts, so
it can run *during* the exchange.

`contact_force_core(pass)` (`dirt_granular/src/contact.rs`) parameterises the pair
loop by `ForcePass { All, Interior, Boundary }`:

- `Interior` resets the active-contact flags and processes only `j < nlocal` pairs,
  but does **not** prune (the history lifecycle is still open).
- `Boundary` does **not** reset, processes only `j >= nlocal` pairs, then prunes.
- `All` is the original single pass (reset → all pairs → prune), bit-unchanged.

So `Interior` then `Boundary` share one reset/prune lifecycle, exactly as `All` does.

## The overlap mechanism

1. `grass_mpi::sendrecv_batch_overlap_f64_into(ops, overlap)` — posts every
   `Isend`/`Irecv`, runs `overlap()` while they fly, then `wait_all`. (MPI backend
   truly overlaps; the trait default runs `overlap()` then a blocking batch.)
2. `soil_core::forward_comm_overlap(.., overlap)` — mirrors `forward_comm` but never
   takes the aggregated single-`sendrecv` shortcut, so every round uses the
   overlap-capable batch. The caller's interior-force closure runs during the first
   round's in-flight window; remaining rounds then land.
3. `dirt_granular::overlapped_contact_force` — runs `contact_force_core(Interior)` as
   that closure, then `contact_force_core(Boundary)` on the landed ghosts. Selectable
   via `DIRT_OVERLAP_FORCE=1` (defaults off; standard force otherwise).

## Verification

**Unit (split exactness).** `interior_boundary_split_matches_single_pass`: 3 atoms
(2 local + 1 ghost), one interior pair and one boundary pair. `Interior + Boundary`
vs `All` agree to `< 1e-15` (bit-for-bit). Part of the 33-test `dirt_granular` suite,
all green.

**Integration (2-rank MPI).** `mpi_loadbalance_validate` (periodic frictional gas,
2000 atoms, 2000 steps, fixed dt), 2 ranks, standard vs overlapped force, compared by
global tag over the union of both ranks:

```bash
DIRT_LB_EVERY=0                      mpiexec -n 2 target/release/examples/mpi_loadbalance_validate cfg_n2.toml   # standard
DIRT_OVERLAP_FORCE=1 DIRT_LB_EVERY=0 mpiexec -n 2 target/release/examples/mpi_loadbalance_validate cfg_n2.toml   # overlapped
```

**Result: bit-identical.** 2000/2000 atoms match exactly, 0 lost. The interior force
computed during the in-flight halo, plus the boundary force on the landed ghosts,
reproduces the standard single-pass force with zero drift under real MPI.

## Provenance

- `grass_mpi`: `sendrecv_batch_overlap_f64_into` (d59dea3).
- `soil_core::comm`: `forward_comm_overlap` (8a2c627).
- `dirt_granular::contact`: `ForcePass`, `contact_force_core`, `overlapped_contact_force`.
