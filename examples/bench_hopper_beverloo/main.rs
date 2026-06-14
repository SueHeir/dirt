//! bench_hopper_beverloo — validates the Beverloo discharge law for a hopper.
//!
//! Particles fill a quasi-2D **slot** hopper (periodic in y, depth a few grain
//! diameters), settle on a blocker wall, then the blocker is removed and the bed
//! discharges through a bottom slot orifice of opening width `D` under gravity.
//! The steady mass-flow rate `W` (per unit slot depth) follows Beverloo's law for
//! a 2D slot:  `W = C · ρ_b · √g · (D − k·d)^(3/2)`.
//!
//! This recorder is thin: it removes the blocker once the bed settles, then logs
//! the cumulative count and mass of particles that have fallen below the orifice
//! plane versus time to `<output_dir>/data/hopper_beverloo_results.csv`. The
//! steady-state slope (W) and the Beverloo exponent are fit in `sweep.py`.
//!
//! ```bash
//! cargo run --release --example bench_hopper_beverloo --no-default-features -- examples/bench_hopper_beverloo/config.toml
//! ```

use dirt_core::prelude::*;
use serde::Deserialize;
use std::collections::HashSet;
use std::fs;
use std::io::Write as IoWrite;

/// `[hopper_beverloo]` recorder parameters from the config.
#[derive(Deserialize, Clone)]
#[serde(default)]
struct HopperBeverlooConfig {
    /// z-height of the orifice plane; a particle counts as discharged once its
    /// center drops below this.
    orifice_z: f64,
    /// Sampling interval (steps) for the cumulative-discharge CSV rows.
    sample_interval: usize,
}

impl Default for HopperBeverlooConfig {
    fn default() -> Self {
        Self { orifice_z: 0.0, sample_interval: 2000 }
    }
}

#[derive(Clone, Debug, PartialEq, Default, StageEnum)]
enum Phase {
    #[default]
    #[stage("filling")]
    Filling,
    #[stage("flowing")]
    Flowing,
}

/// Records cumulative discharge (count + mass) of particles crossing the orifice
/// plane `z = orifice_z` downward, sampled every `sample_interval` steps.
struct DischargeTracker {
    orifice_z: f64,
    sample_interval: usize,
    flow_start_step: usize,
    /// Tags already counted as discharged (crossed below the orifice once).
    discharged: HashSet<u32>,
    discharged_mass: f64,
    /// Buffered rows: (time_since_open, cumulative_count, cumulative_mass).
    rows: Vec<(f64, usize, f64)>,
}

impl DischargeTracker {
    fn new(orifice_z: f64, sample_interval: usize) -> Self {
        Self {
            orifice_z,
            sample_interval,
            flow_start_step: 0,
            discharged: HashSet::new(),
            discharged_mass: 0.0,
            rows: Vec::new(),
        }
    }
}

fn main() {
    // The orifice plane and sample interval are read from the config's
    // [hopper_beverloo] table so a single recorder serves the whole sweep.
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin)
        .add_plugins(StatesPlugin::new(Phase::Filling, ParticleSimScheduleSet::PostFinalIntegration))
        .add_plugins(StageAdvancePlugin::<Phase>::new(ParticleSimScheduleSet::PostFinalIntegration));

    let (orifice_z, sample_interval) = read_tracker_params(&mut app);
    app.add_resource(DischargeTracker::new(orifice_z, sample_interval));

    // The two `[[run]]` stages are named "filling" and "flowing"; the stack auto-
    // advances Phase at each stage boundary. On entering "flowing" we pull the
    // blocker wall and arm the discharge recorder. During "flowing" we record.
    app.add_update_system(
        open_orifice.run_if(on_enter_stage("flowing")),
        ParticleSimScheduleSet::PostFinalIntegration,
    );
    app.add_update_system(
        record_discharge.run_if(in_stage("flowing")),
        ParticleSimScheduleSet::PostFinalIntegration,
    );

    app.start();
}

/// Parse the `[hopper_beverloo]` section (with defaults) before the app starts.
fn read_tracker_params(app: &mut App) -> (f64, usize) {
    let config = app
        .get_resource_ref::<Config>()
        .expect("Config resource must exist");
    let cfg = config.section::<HopperBeverlooConfig>("hopper_beverloo");
    (cfg.orifice_z, cfg.sample_interval)
}

/// On entering the flowing stage: remove the blocker wall and mark the time
/// origin for the discharge curve.
fn open_orifice(
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    mut walls: ResMut<Walls>,
    mut tracker: ResMut<DischargeTracker>,
) {
    walls.deactivate_by_name("blocker");
    tracker.flow_start_step = run_state.total_cycle;
    if comm.rank() == 0 {
        println!(
            "Step {}: blocker removed — orifice open, discharge begins",
            run_state.total_cycle
        );
    }
}

/// During flowing: count particles that have fallen below the orifice plane and
/// buffer cumulative (time, count, mass) samples. Writes the CSV once the run ends
/// (when no particles remain above the orifice, or on the final sampled step).
fn record_discharge(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
    mut tracker: ResMut<DischargeTracker>,
) {
    let step = run_state.total_cycle;
    let nlocal = atoms.nlocal as usize;

    // Tally any particle that is now below the orifice and not yet counted.
    let oz = tracker.orifice_z;
    let mut new_mass = 0.0;
    let mut newly: Vec<u32> = Vec::new();
    for i in 0..nlocal {
        if atoms.pos[i][2] < oz {
            let tag = atoms.tag[i];
            if !tracker.discharged.contains(&tag) {
                newly.push(tag);
                new_mass += atoms.mass[i];
            }
        }
    }
    for tag in newly {
        tracker.discharged.insert(tag);
    }
    tracker.discharged_mass += new_mass;

    if step % tracker.sample_interval != 0 {
        return;
    }

    // Globalize the cumulative discharge across ranks for this sample.
    let count = comm.all_reduce_sum_f64(tracker.discharged.len() as f64) as usize;
    let mass = comm.all_reduce_sum_f64(tracker.discharged_mass);
    let t = (step - tracker.flow_start_step) as f64 * atoms.dt;
    tracker.rows.push((t, count, mass));

    if comm.rank() == 0 {
        println!(
            "  [discharge] t={:.4e}s  count={}  mass={:.4e}kg",
            t, count, mass
        );
    }

    // Rewrite the CSV on every sample so the latest cumulative curve is always
    // on disk, whether the run ends by draining or by exhausting its step budget.
    write_csv(&tracker, &input, comm.rank());
}

fn write_csv(tracker: &DischargeTracker, input: &Input, rank: i32) {
    if rank != 0 {
        return;
    }
    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| "examples/bench_hopper_beverloo".to_string());
    let data_dir = format!("{}/data", out_dir);
    fs::create_dir_all(&data_dir).ok();
    let path = format!("{}/hopper_beverloo_results.csv", data_dir);
    let mut f = fs::File::create(&path)
        .unwrap_or_else(|e| panic!("Cannot create {}: {}", path, e));
    writeln!(f, "time,count,mass").unwrap();
    for (t, c, m) in &tracker.rows {
        writeln!(f, "{:.10e},{},{:.10e}", t, c, m).unwrap();
    }
    println!("  discharge results -> {}", path);
}
