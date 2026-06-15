//! perf_mpi_scaling — DEM throughput and MPI strong/weak scaling.
//!
//! This is a *performance* benchmark, not a validation one: it measures how fast
//! DIRT integrates a fixed physics setup and how that rate scales with MPI ranks.
//! It records nothing about correctness — `sweep.py` turns the raw timings into
//! strong/weak scaling curves and (optionally) overlays LAMMPS.
//!
//! The recorder times a **steady-state window** only. It skips a warm-up prefix
//! (the first `WARMUP_FRACTION` of the run — insertion, settling, transient
//! rebinning), then `barrier()`s all ranks, starts a wall clock, and at the end
//! reports throughput as **particle-steps per second**:
//!
//!   particle_steps_per_s = N_global * steps_measured / wall_seconds
//!
//! a size-independent figure of merit that compares directly across particle
//! counts, rank counts, and codes. The global particle count is obtained with an
//! all-reduce so it is correct under domain decomposition.
//!
//! Build WITH default features so the `mpi_backend` is compiled in, and launch
//! under mpiexec (1 rank is a valid serial smoke test):
//!
//! ```bash
//! cargo build --release --example perf_mpi_scaling
//! mpiexec -n 4 target/release/examples/perf_mpi_scaling examples/perf_mpi_scaling/config.toml
//! ```

use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;
use std::time::Instant;

/// Fraction of the total run spent warming up (skipped before timing starts).
/// The measured window is the remaining `1 - WARMUP_FRACTION`.
const WARMUP_FRACTION: f64 = 0.4;

/// How often (in steps) the running summary is recomputed and the CSV rewritten.
/// Each write barriers all ranks, so keep it coarse; the final write is the one
/// that covers the full steady-state window.
const MEASURE_EVERY: usize = 2000;

/// Timing state. Single-instance resource, mutated by the recorder.
struct PerfTracker {
    started: bool,
    t0: Option<Instant>,
    step0: usize,
    n_global: f64,
}

impl PerfTracker {
    fn new() -> Self {
        Self { started: false, t0: None, step0: 0, n_global: 0.0 }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin) // bed scenario settles under gravity; gas sets g = 0
        .add_plugins(WallPlugin); // bed scenario stands on a floor wall; gas has none

    app.add_resource(PerfTracker::new());
    app.add_update_system(record_throughput, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// Times the steady-state window and writes the throughput summary.
///
/// Runs every step in `PostFinalIntegration`. Before the warm-up boundary it does
/// nothing; at the boundary it snapshots the global particle count and starts the
/// clock; afterwards, at a coarse cadence (and on the final step), it rewrites a
/// one-row summary CSV on rank 0.
fn record_throughput(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
    mut tracker: ResMut<PerfTracker>,
) {
    let step = run_state.total_cycle;
    // `cycle_remaining` holds the *constant target* step count per stage (= stage
    // `steps`); `cycle_count` is the running per-stage counter. Total run length is
    // the sum of the targets.
    let total: usize = run_state.cycle_remaining.iter().map(|&r| r as usize).sum();
    if total == 0 {
        return;
    }
    let warmup = (total as f64 * WARMUP_FRACTION) as usize;

    // ── Start of the measured window: sync, snapshot N_global, start the clock ──
    if !tracker.started {
        if step < warmup {
            return;
        }
        let n_global = comm.all_reduce_sum_f64(atoms.nlocal as f64);
        comm.barrier();
        tracker.t0 = Some(Instant::now());
        tracker.step0 = step;
        tracker.n_global = n_global;
        tracker.started = true;
        return;
    }

    // ── Periodic / final summary write ──
    // Steps left = Σ(target − counter) across stages; <= 1 means the final step.
    let steps_left: usize = run_state
        .cycle_remaining
        .iter()
        .zip(run_state.cycle_count.iter())
        .map(|(&target, &done)| (target as usize).saturating_sub(done as usize))
        .sum();
    let is_end = steps_left <= 1;
    if step % MEASURE_EVERY != 0 && !is_end {
        return;
    }

    comm.barrier();
    let elapsed = tracker.t0.expect("clock started").elapsed().as_secs_f64();
    let steps_done = step - tracker.step0;

    // ── Conservation diagnostics (collective: all ranks must call) ──
    // Current global particle count — if it differs from the start count, the
    // domain-decomposition migration/halo path lost (or duplicated) atoms.
    // Total kinetic energy guards against a blow-up (NaN/inf) under MPI.
    let n_now = comm.all_reduce_sum_f64(atoms.nlocal as f64);
    let mut ke_local = 0.0_f64;
    for i in 0..atoms.nlocal as usize {
        let m = atoms.mass[i];
        let v = atoms.vel[i];
        ke_local += 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    let ke_total = comm.all_reduce_sum_f64(ke_local);

    if comm.rank() != 0 || steps_done == 0 || elapsed <= 0.0 {
        return;
    }

    let decomp = comm.processor_decomposition();
    let ranks = comm.size();
    let particle_steps = tracker.n_global * steps_done as f64;
    let pstep_per_s = particle_steps / elapsed;
    let us_per_pstep = elapsed * 1.0e6 / particle_steps;
    let us_per_step = elapsed * 1.0e6 / steps_done as f64;

    let output_dir = input
        .output_dir
        .as_deref()
        .unwrap_or("examples/perf_mpi_scaling/data");
    let data_dir = format!("{}/data", output_dir);
    fs::create_dir_all(&data_dir).ok();
    let filepath = format!("{}/perf_results.csv", data_dir);

    // Truncate-and-rewrite: a single header + row, refreshed each cadence so the
    // last write reflects the full measured window.
    let mut f = fs::File::create(&filepath).expect("cannot create perf_results.csv");
    writeln!(
        f,
        "ranks,px,py,pz,n_global_start,n_global_end,ke_total,steps_measured,wall_s,particle_steps_per_s,us_per_particle_step,us_per_step"
    )
    .unwrap();
    writeln!(
        f,
        "{},{},{},{},{},{},{:.6e},{},{:.6e},{:.6e},{:.6e},{:.6e}",
        ranks,
        decomp[0],
        decomp[1],
        decomp[2],
        tracker.n_global as u64,
        n_now as u64,
        ke_total,
        steps_done,
        elapsed,
        pstep_per_s,
        us_per_pstep,
        us_per_step
    )
    .unwrap();
}
