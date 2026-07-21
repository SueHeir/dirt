//! bench_column_collapse — validates granular column-collapse runout scaling
//! against the planar experimental aspect-ratio law of Lajeunesse et al.
//! (2004), doi:10.1063/1.1736611. The high-aspect axisymmetric law of Lube et
//! al. is deliberately not used: this is a quasi-2D planar-gate benchmark.
//!
//! A quasi-2D rectangular column of grains (initial width L0, height H) is held
//! against a vertical gate wall on a flat floor. Stage 1 ("settle") lets the
//! loosely-inserted column pack down under gravity into a static column. Stage 2
//! ("collapse") removes the gate on its first step; the column collapses and
//! spreads along +x until it comes to rest. The dimensionless runout
//! (L_f - L0)/L0 is expected to follow:
//!   (L_f - L0)/L0 ~ 1.2 a        (a <~ 2-3, linear regime)
//!   (L_f - L0)/L0 ~ 1.6 a^(2/3)  (a >~ 3,   power-law regime)
//! with a = H/L0.
//!
//! This recorder is analysis-free: it dumps the final (x, y, z, radius) of every
//! particle to `<output_dir>/data/column_collapse_results.csv`. All runout
//! extraction, regime fitting and PASS/FAIL live in `sweep.py`.
//!
//! ```bash
//! cargo run --release --example bench_column_collapse --no-default-features -- examples/bench_column_collapse/config.toml
//! ```

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

/// Name of the removable vertical gate wall (matches `name = "gate"` in config).
const GATE_NAME: &str = "gate";
/// Number of collapse integration steps between quiescence witnesses.
///
/// The benchmark is single-rank, so this records the complete population used
/// by the runout analysis rather than a rank-local proxy.
// The released stage retains the 1 us resolved timestep, so 100,000 steps
// retain the 0.1 s physical witness interval used by the acceptance protocol.
const ARREST_SAMPLE_INTERVAL: u64 = 100_000;
const PREPARATION_SAMPLE_INTERVAL: u64 = 100_000;

/// Tracks gate release so it happens exactly once.
struct CollapseTracker {
    gate_opened: bool,
    collapse_steps: u64,
    settle_steps: u64,
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(FixesPlugin)
        .add_plugins(WallPlugin);

    app.add_resource(CollapseTracker {
        gate_opened: false,
        collapse_steps: 0,
        settle_steps: 0,
    });

    // Start the released dynamics at the first collapse integration boundary.
    // The gate must not remain active for one unassisted force evaluation after
    // the preparation aids have been removed.
    app.add_update_system(
        begin_collapse.run_if(in_stage("collapse")),
        ParticleSimScheduleSet::PreInitialIntegration,
    );
    app.add_update_system(
        record_collapse_quiescence.run_if(in_stage("collapse")),
        ParticleSimScheduleSet::PostFinalIntegration,
    );
    app.add_update_system(
        record_preparation_quiescence.run_if(in_stage("settle")),
        ParticleSimScheduleSet::PostFinalIntegration,
    );

    app.start();

    // Dump the final deposit once the run has finished and the bed is at rest.
    dump_deposit(&app);
}

/// Preserve a sustained, still-gated rest witness.  The fixed 0.8 s protocol
/// is unchanged; analysis rejects a source that only happens to be slow in
/// its final frame.
fn record_preparation_quiescence(
    mut tracker: ResMut<CollapseTracker>,
    atoms: Res<Atom>,
    input: Res<Input>,
    comm: Res<CommResource>,
) {
    tracker.settle_steps += 1;
    if tracker.settle_steps % PREPARATION_SAMPLE_INTERVAL != 0 {
        return;
    }
    let mut vmax = 0.0f64;
    for velocity in atoms.vel.iter().take(atoms.nlocal as usize) {
        vmax = vmax.max(
            ((velocity[0] as f64).powi(2)
                + (velocity[1] as f64).powi(2)
                + (velocity[2] as f64).powi(2))
            .sqrt(),
        );
    }
    if comm.rank() == 0 {
        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_column_collapse".to_string());
        let data_dir = format!("{out_dir}/data");
        // This is the first persistent witness written by a fresh case.  The
        // release callback also creates this directory, but it runs only after
        // the final preparation sample; do not make a valid still-gated
        // witness depend on a later lifecycle phase.
        fs::create_dir_all(&data_dir).unwrap_or_else(|e| {
            panic!("Cannot create preparation-witness directory {data_dir}: {e}")
        });
        let path = format!("{data_dir}/column_collapse_preparation.csv");
        let mut f = if tracker.settle_steps == PREPARATION_SAMPLE_INTERVAL {
            let mut created =
                fs::File::create(&path).unwrap_or_else(|e| panic!("Cannot create {path}: {e}"));
            writeln!(created, "settle_step,particle_count,max_speed_m_s").unwrap();
            created
        } else {
            fs::OpenOptions::new()
                .append(true)
                .open(&path)
                .unwrap_or_else(|e| panic!("Cannot append {path}: {e}"))
        };
        // The last sampled settle frame is the same still-gated state that the
        // following pre-integration callback releases.  Keep population and
        // speed in this one time-series witness rather than creating a second,
        // redundant release-state file at a subtly different callback point.
        writeln!(f, "{},{},{:.10e}", tracker.settle_steps, atoms.nlocal, vmax).unwrap();
    }
}

/// Begin the released stage before its first force/integration evaluation.
///
/// Damping and the displacement limiter are strictly source-preparation aids,
/// while the gate is a retained support.  Removing all three at this common
/// pre-integration boundary prevents a one-step hybrid state (unlimited,
/// undamped, but still gated). The preceding final settle sample is the
/// still-gated release-rest witness for the actual collapse dynamics.
fn begin_collapse(
    mut tracker: ResMut<CollapseTracker>,
    mut fixes: ResMut<FixesRegistry>,
    mut walls: ResMut<Walls>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    input: Res<Input>,
) {
    if tracker.gate_opened {
        return;
    }
    fixes.cundall.clear();
    fixes.nve_limit.clear();
    // The released geometry is evidence, not an inferred nominal particle
    // count.  Record it before support removal so analysis can fit the actual
    // settled aspect ratio of each realization.
    let dem = registry.expect::<DemAtom>("record_release_state");
    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| "examples/bench_column_collapse".to_string());
    let data_dir = format!("{out_dir}/data");
    fs::create_dir_all(&data_dir).expect("create release-state directory");
    let path = format!("{data_dir}/column_collapse_release.csv");
    let mut file = fs::File::create(&path).unwrap_or_else(|e| panic!("Cannot create {path}: {e}"));
    writeln!(file, "x,y,z,radius").unwrap();
    for i in 0..atoms.nlocal as usize {
        writeln!(
            file,
            "{:.10e},{:.10e},{:.10e},{:.10e}",
            atoms.pos[i][0], atoms.pos[i][1], atoms.pos[i][2], dem.radius[i]
        )
        .unwrap();
    }
    walls.deactivate_by_name(GATE_NAME);
    // A new collapse must never inherit a previous run's arrest witnesses.
    let arrest_path = format!("{data_dir}/column_collapse_arrest.csv");
    let mut arrest = fs::File::create(&arrest_path)
        .unwrap_or_else(|e| panic!("Cannot create {arrest_path}: {e}"));
    writeln!(arrest, "collapse_step,max_speed_m_s").unwrap();
    tracker.gate_opened = true;
    if comm.rank() == 0 {
        println!(
            "Step {}: gate removed — column released.",
            run_state.total_cycle
        );
    }
}

/// Record a sparse, sustained-quiescence witness through the collapse stage.
///
/// A terminal speed alone can catch a transient low-velocity phase.  The
/// analysis requires the final four equally spaced samples to meet its Froude
/// threshold, which establishes a bounded interval of rest without changing a
/// runout estimator or an empirical acceptance band.
fn record_collapse_quiescence(
    mut tracker: ResMut<CollapseTracker>,
    atoms: Res<Atom>,
    input: Res<Input>,
    comm: Res<CommResource>,
) {
    if !tracker.gate_opened {
        return;
    }
    tracker.collapse_steps += 1;
    if tracker.collapse_steps % ARREST_SAMPLE_INTERVAL != 0 {
        return;
    }
    let mut vmax = 0.0f64;
    for velocity in atoms.vel.iter().take(atoms.nlocal as usize) {
        let speed = ((velocity[0] as f64) * (velocity[0] as f64)
            + (velocity[1] as f64) * (velocity[1] as f64)
            + (velocity[2] as f64) * (velocity[2] as f64))
            .sqrt();
        vmax = vmax.max(speed);
    }
    if comm.rank() == 0 {
        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_column_collapse".to_string());
        let path = format!("{out_dir}/data/column_collapse_arrest.csv");
        let mut f = fs::OpenOptions::new()
            .append(true)
            .open(&path)
            .unwrap_or_else(|e| panic!("Cannot append {path}: {e}"));
        writeln!(f, "{},{:.10e}", tracker.collapse_steps, vmax).unwrap();
    }
}

/// Write the deposit profile (per-particle x, y, z, radius) so `sweep.py` can
/// extract the final runout L_f. Called after `start()`, so positions are the
/// settled rest state.
fn dump_deposit(app: &App) {
    let atoms = match app.get_resource_ref::<Atom>() {
        Some(a) => a,
        None => return,
    };
    let registry = app
        .get_resource_ref::<AtomDataRegistry>()
        .expect("AtomDataRegistry must exist");
    let dem = registry.expect::<DemAtom>("dump_deposit");
    let nlocal = atoms.nlocal as usize;

    let out_dir = app
        .get_resource_ref::<Input>()
        .and_then(|i| i.output_dir.clone())
        .unwrap_or_else(|| "examples/bench_column_collapse".to_string());
    let data_dir = format!("{}/data", out_dir);
    fs::create_dir_all(&data_dir).ok();
    let results_file = format!("{}/column_collapse_results.csv", data_dir);

    let mut f = fs::File::create(&results_file)
        .unwrap_or_else(|e| panic!("Cannot create {}: {}", results_file, e));
    writeln!(f, "x,y,z,radius").unwrap();
    let mut vmax = 0.0f64;
    for i in 0..nlocal {
        writeln!(
            f,
            "{:.10e},{:.10e},{:.10e},{:.10e}",
            atoms.pos[i][0], atoms.pos[i][1], atoms.pos[i][2], dem.radius[i]
        )
        .unwrap();
        let v = atoms.vel[i];
        let s = ((v[0] as f64) * (v[0] as f64)
            + (v[1] as f64) * (v[1] as f64)
            + (v[2] as f64) * (v[2] as f64))
            .sqrt();
        if s > vmax {
            vmax = s;
        }
    }

    println!(
        "FINAL: {} particles dumped -> {} (max |v| = {:.3e} m/s)",
        nlocal, results_file, vmax
    );
    let state_file = format!("{data_dir}/column_collapse_final_state.csv");
    let mut state =
        fs::File::create(&state_file).unwrap_or_else(|e| panic!("Cannot create {state_file}: {e}"));
    writeln!(state, "particle_count,max_speed_m_s").unwrap();
    writeln!(state, "{nlocal},{vmax:.10e}").unwrap();
}
