//! Hopper fill/discharge benchmark for the **region_coherence** DEM quiescence
//! optimization: a coarse cell grid classifies cells whose particles move
//! coherently as either *plug* (internal contacts force-freeze to a cached
//! constant force, thawing when the pair drifts) or *sleeping* (integration and
//! internal pair computation skipped entirely).
//!
//! Everything lives in this example — the DIRT crates are unmodified. The
//! stock Hertz-Mindlin contact plugin is replaced by a local copy with the
//! optimization hooks ([`contact`]) plus the region state machine
//! ([`quiescence`]).
//!
//! Stage 1 ("filling") rains particles into the hopper through rate insertion
//! and lets them settle on a blocker wall. Stage 2 ("flowing") removes the
//! blocker on its first step (a fixed, deterministic step count so that all
//! benchmark variants see identical staging) and the bed discharges onto the
//! floor below, where it settles again.
//!
//! ```bash
//! cargo run --release --no-default-features --example hopper_quiescence -- \
//!     examples/hopper_quiescence/config_baseline.toml
//! ```
//!
//! Per-step statistics (kinetic energy, discharged count, asleep/plug/frozen
//! fractions) are appended to `<config>_stats.csv` next to the config file.

mod contact;
mod quiescence;

use std::fs::File;
use std::io::Write as IoWrite;
use std::time::Instant;

use contact::{QcStore, QuiescentContactPlugin, MODE_ACTIVE, MODE_ADJACENT, MODE_ASLEEP};
use dirt_core::prelude::*;
use quiescence::QuiescencePlugin;

/// Height of the funnel outlet plane; particles below it count as discharged.
const GATE_Z: f64 = 0.058;
/// Stats CSV output cadence in steps.
const STATS_EVERY: usize = 1000;

struct Bench {
    t0: Instant,
    csv: File,
    gate_opened: bool,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args
        .get(1)
        .map(|p| p.replace(".toml", "_stats.csv"))
        .unwrap_or_else(|| "hopper_quiescence_stats.csv".into());
    let mut csv = File::create(&csv_path).expect("create stats csv");
    writeln!(
        csv,
        "step,elapsed_s,ke_total,n_atoms,n_discharged,n_active,n_adjacent,n_asleep,pairs,skipped_pairs,top_z"
    )
    .unwrap();

    let mut app = App::new();
    // GranularDefaultPlugins minus HertzMindlinContactPlugin, which is replaced
    // by this example's QuiescentContactPlugin.
    let use_stock = std::env::var("USE_STOCK_CONTACT").is_ok();
    app.add_plugins(CorePlugins)
        .add_plugins(DemAtomPlugin)
        .add_plugins(DemAtomInsertPlugin)
        .add_plugins(VelocityVerletPlugin::new());
    if use_stock {
        app.add_plugins(HertzMindlinContactPlugin);
    } else {
        app.add_plugins(QuiescentContactPlugin);
    }
    app.add_plugins(RotationalDynamicsPlugin)
        .add_plugins(QuiescencePlugin)
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin);

    app.add_resource(Bench {
        t0: Instant::now(),
        csv,
        gate_opened: false,
    });
    app.add_update_system(
        open_gate.run_if(in_stage("flowing")),
        ParticleSimScheduleSet::PostFinalIntegration,
    );
    app.add_update_system(write_stats, ParticleSimScheduleSet::PostFinalIntegration);

    let t0 = Instant::now();
    app.start();
    let total = t0.elapsed().as_secs_f64();

    if let Some(atoms) = app.get_resource_ref::<Atom>() {
        let nlocal = atoms.nlocal as usize;
        let discharged = (0..nlocal).filter(|&i| (atoms.pos[i][2] as f64) < GATE_Z).count();
        println!(
            "FINAL: {} particles, {} discharged below gate",
            nlocal, discharged
        );
        // Diagnostic: z distribution and speed percentiles.
        let mut zs: Vec<f64> = (0..nlocal).map(|i| atoms.pos[i][2] as f64).collect();
        let mut vs: Vec<f64> = (0..nlocal)
            .map(|i| {
                let v = atoms.vel[i];
                (v[0] as f64 * v[0] as f64 + v[1] as f64 * v[1] as f64 + v[2] as f64 * v[2] as f64)
                    .sqrt()
            })
            .collect();
        zs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        vs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if nlocal > 0 {
            println!(
                "FINAL z: min={:.4} p25={:.4} p50={:.4} p75={:.4} max={:.4}",
                zs[0], zs[nlocal / 4], zs[nlocal / 2], zs[3 * nlocal / 4], zs[nlocal - 1]
            );
            println!(
                "FINAL |v|: p50={:.4} p90={:.4} p99={:.4} max={:.4}",
                vs[nlocal / 2], vs[nlocal * 9 / 10], vs[nlocal * 99 / 100], vs[nlocal - 1]
            );
        }
    }
    println!("TOTAL WALL TIME: {:.2} s (stats: {})", total, csv_path);
}

/// Remove the blocker wall on the first step of the "flowing" run stage, and
/// wake every frozen particle. Removing a support is a disturbance the local
/// force-deviation wake can't see (funnel-wall friction seamlessly takes up the
/// load with no net-force change), so a frozen arch would otherwise jam the
/// outlet permanently — we break it by waking the whole bed once.
fn open_gate(
    mut bench: ResMut<Bench>,
    mut walls: ResMut<Walls>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
) {
    if bench.gate_opened {
        return;
    }
    walls.deactivate_by_name("blocker");
    bench.gate_opened = true;

    if let Some(mut store) = registry.get_mut::<QcStore>() {
        for i in 0..store.mode.len() {
            store.mode[i] = MODE_ACTIVE;
            store.still_count[i] = 0;
            store.has_base[i] = false;
            store.accum_fnet[i] = [0.0; 3];
            store.accum_scale[i] = 0.0;
        }
    }

    if comm.rank() == 0 {
        println!(
            "Step {}: gate opened (blocker removed + all particles woken), wall time {:.2} s",
            run_state.total_cycle,
            bench.t0.elapsed().as_secs_f64()
        );
    }
}

/// Append per-step statistics to the CSV every [`STATS_EVERY`] steps.
fn write_stats(
    mut bench: ResMut<Bench>,
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
) {
    let step = run_state.total_cycle;
    if step % STATS_EVERY != 0 {
        return;
    }
    let nlocal = atoms.nlocal as usize;

    let mut ke = 0.0;
    let mut discharged = 0usize;
    // Pile crest: tallest particle, a settled-bed height proxy for the fill /
    // discharge fidelity check (read at the end-of-filling step in bench.sh).
    let mut top_z = f64::MIN;
    for i in 0..nlocal {
        let v = atoms.vel[i];
        ke += 0.5
            * atoms.mass[i] as f64
            * (v[0] as f64 * v[0] as f64
                + v[1] as f64 * v[1] as f64
                + v[2] as f64 * v[2] as f64);
        if (atoms.pos[i][2] as f64) < GATE_Z {
            discharged += 1;
        }
        if atoms.pos[i][2] as f64 > top_z {
            top_z = atoms.pos[i][2] as f64;
        }
    }
    let top_z = if nlocal > 0 { top_z } else { 0.0 };

    let (mut active, mut adjacent, mut asleep, mut pairs, mut skipped) = (0usize, 0, 0, 0, 0);
    if let Some(store) = registry.get::<QcStore>() {
        for i in 0..nlocal.min(store.mode.len()) {
            match store.mode[i] {
                MODE_ASLEEP => asleep += 1,
                MODE_ADJACENT => adjacent += 1,
                _ => active += 1,
            }
        }
        pairs = store.n_pairs + store.n_skipped;
        skipped = store.n_skipped;
    }

    let elapsed = bench.t0.elapsed().as_secs_f64();
    writeln!(
        bench.csv,
        "{},{:.3},{:.6e},{},{},{},{},{},{},{},{:.5}",
        step, elapsed, ke, nlocal, discharged, active, adjacent, asleep, pairs, skipped, top_z
    )
    .unwrap();
}
