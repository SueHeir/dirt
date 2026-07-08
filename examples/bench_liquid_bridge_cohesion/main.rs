//! Pendular liquid-bridge cohesion recorder.
//!
//! A frozen sphere and a slowly separating free sphere exercise the opt-in
//! Willett et al. (2000) liquid-bridge force. The recorder writes the normal
//! force versus surface separation; validation and plotting live in `sweep.py`.

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

struct BridgeTracker {
    moving_tag: Option<u32>,
    frozen_tag: Option<u32>,
    trace: Vec<(usize, f64, f64)>,
    finished: bool,
}

impl BridgeTracker {
    fn new() -> Self {
        Self {
            moving_tag: None,
            frozen_tag: None,
            trace: Vec::new(),
            finished: false,
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(FixesPlugin);

    app.add_resource(BridgeTracker::new());
    app.add_update_system(record_bridge, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

fn index_of_tag(atoms: &Atom, tag: u32) -> Option<usize> {
    (0..atoms.nlocal as usize).find(|&i| atoms.tag[i] == tag)
}

fn record_bridge(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run_state: Res<RunState>,
    input: Res<Input>,
    mut tracker: ResMut<BridgeTracker>,
) {
    if tracker.finished || atoms.nlocal < 2 {
        return;
    }
    if tracker.moving_tag.is_none() {
        for i in 0..atoms.nlocal as usize {
            let speed2 = atoms.vel[i][0] as f64 * atoms.vel[i][0] as f64
                + atoms.vel[i][1] as f64 * atoms.vel[i][1] as f64
                + atoms.vel[i][2] as f64 * atoms.vel[i][2] as f64;
            if speed2 > 0.0 {
                tracker.moving_tag = Some(atoms.tag[i]);
            } else {
                tracker.frozen_tag = Some(atoms.tag[i]);
            }
        }
    }

    let Some(m) = tracker.moving_tag.and_then(|tag| index_of_tag(&atoms, tag)) else {
        return;
    };
    let Some(f) = tracker.frozen_tag.and_then(|tag| index_of_tag(&atoms, tag)) else {
        return;
    };
    let dem = registry.expect::<DemAtom>("record_bridge");
    let d = [
        atoms.pos[m][0] as f64 - atoms.pos[f][0] as f64,
        atoms.pos[m][1] as f64 - atoms.pos[f][1] as f64,
        atoms.pos[m][2] as f64 - atoms.pos[f][2] as f64,
    ];
    let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if dist == 0.0 {
        return;
    }
    let n = [d[0] / dist, d[1] / dist, d[2] / dist];
    let separation = dist - dem.radius[m] - dem.radius[f];
    let fvec = [
        atoms.force[m][0] as f64,
        atoms.force[m][1] as f64,
        atoms.force[m][2] as f64,
    ];
    let normal_force = fvec[0] * n[0] + fvec[1] * n[1] + fvec[2] * n[2];
    tracker
        .trace
        .push((run_state.total_cycle, separation, normal_force));

    let last_step = match (
        run_state.cycle_count.first(),
        run_state.cycle_remaining.first(),
    ) {
        (Some(&done), Some(&total)) => total > 0 && done + 1 >= total,
        _ => false,
    };
    if last_step {
        tracker.finished = true;
        let out_dir = input
            .output_dir
            .clone()
            .unwrap_or_else(|| "examples/bench_liquid_bridge_cohesion".to_string());
        let data_dir = format!("{}/data", out_dir);
        fs::create_dir_all(&data_dir).ok();
        let path = format!("{}/bridge_trace.csv", data_dir);
        let mut file =
            fs::File::create(&path).unwrap_or_else(|e| panic!("Cannot create {}: {}", path, e));
        writeln!(file, "step,separation,normal_force").unwrap();
        for (step, separation, normal_force) in &tracker.trace {
            writeln!(file, "{},{:.12e},{:.12e}", step, separation, normal_force).unwrap();
        }
        println!(
            "wrote {} liquid-bridge samples -> {}",
            tracker.trace.len(),
            path
        );
    }
}
