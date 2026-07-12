//! Real DIRT bonded-sphere timestep probe.
//!
//! The sweep driver launches this example repeatedly with a fixed-free chain
//! of bonded spheres.  This is deliberately a propagation problem rather than
//! a one-bond oscillator: the short-wavelength chain mode is what makes the
//! Guo et al. axial-wave timestep bound relevant.

use dirt_core::dirt_fixes::FixesPlugin;
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

struct Probe {
    initialized: bool,
    written: bool,
    mode: String,
    drive: f64,
    anchor_tag: Option<u32>,
    mover_tag: Option<u32>,
    x0: f64,
    rows: Vec<(usize, f64, f64, f64)>,
}

impl Probe {
    fn new() -> Self {
        Self {
            initialized: false,
            written: false,
            mode: std::env::var("FIBER_TIMESTEP_MODE").unwrap_or_else(|_| "axial".to_string()),
            drive: std::env::var("FIBER_TIMESTEP_DRIVE")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(1.0e-3),
            anchor_tag: None,
            mover_tag: None,
            x0: 0.0,
            rows: Vec::new(),
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(FixesPlugin)
        .add_plugins(DemBondPlugin);

    app.add_resource(Probe::new());
    app.add_update_system(init_probe, ParticleSimScheduleSet::PreForce);
    app.add_update_system(record_probe, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

fn init_probe(mut atoms: ResMut<Atom>, mut probe: ResMut<Probe>) {
    if probe.initialized || atoms.nlocal < 3 {
        return;
    }
    let anchor = (0..atoms.nlocal as usize)
        .min_by(|&a, &b| atoms.pos[a][0].partial_cmp(&atoms.pos[b][0]).unwrap())
        .unwrap();
    let mover = (0..atoms.nlocal as usize)
        .max_by(|&a, &b| atoms.pos[a][0].partial_cmp(&atoms.pos[b][0]).unwrap())
        .unwrap();

    probe.anchor_tag = Some(atoms.tag[anchor]);
    probe.mover_tag = Some(atoms.tag[mover]);
    probe.x0 = atoms.pos[mover][0] as f64;

    if probe.mode == "bending" {
        // A transverse tip impulse excites the bonded chain's bending/shear
        // dynamics; no analytic single-coordinate reduction is assumed.
        atoms.vel[mover][2] = probe.drive as dirt_core::soil_core::Real;
    } else {
        atoms.vel[mover][0] = probe.drive as dirt_core::soil_core::Real;
    }
    probe.initialized = true;
}

fn record_probe(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    run_config: Res<RunConfig>,
    input: Res<Input>,
    mut probe: ResMut<Probe>,
) {
    if probe.written || !probe.initialized {
        return;
    }
    let (tag_a, tag_b) = match (probe.anchor_tag, probe.mover_tag) {
        (Some(a), Some(b)) => (a, b),
        _ => return,
    };
    if !(0..atoms.nlocal as usize).any(|i| atoms.tag[i] == tag_a) {
        return;
    }
    let ib = match (0..atoms.nlocal as usize).find(|&i| atoms.tag[i] == tag_b) {
        Some(i) => i,
        None => return,
    };

    let q = if probe.mode == "bending" {
        atoms.pos[ib][2] as f64
    } else {
        atoms.pos[ib][0] as f64 - probe.x0
    };
    let finite = if q.is_finite()
        && (0..atoms.nlocal as usize).all(|i| {
            atoms.pos[i].iter().all(|x| (*x as f64).is_finite())
                && atoms.vel[i].iter().all(|x| (*x as f64).is_finite())
        }) {
        1.0
    } else {
        0.0
    };
    let step = run_state.total_cycle;
    probe.rows.push((step, step as f64 * atoms.dt, q, finite));

    let last_step = run_config.current_stage(0).steps as usize;
    if step + 1 >= last_step {
        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_fiber_timestep".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();
        let path = format!("{}/timestep_probe.csv", data_dir);
        let mut f =
            fs::File::create(&path).unwrap_or_else(|e| panic!("Cannot create {}: {}", path, e));
        writeln!(f, "step,t,q,finite").unwrap();
        for (step, t, q, finite) in &probe.rows {
            writeln!(f, "{},{:.10e},{:.10e},{:.0}", step, t, q, finite).unwrap();
        }
        println!("=== Fiber Timestep Probe ===");
        println!("  mode       : {}", probe.mode);
        println!("  rows       : {}", probe.rows.len());
        println!("  output     : {}", path);
        probe.written = true;
    }
}
