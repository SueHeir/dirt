//! Slot hopper discharge: particles fill a funnel, then flow through an exit.
//!
//! Demonstrates multi-stage simulation with runtime wall control (removing a
//! blocker wall when particles settle) and KE-based stage transitions.
//!
//! ```bash
//! cargo run --example hopper --no-default-features -- examples/hopper/config.toml
//! ```

use dirt_core::prelude::*;

const FILLING_STAGE: &str = "filling";

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GranularTempPlugin) // hopper/validate.py reads data/GranularTemp.txt
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin);

    app.add_update_system(
        check_settled.run_if(in_stage(FILLING_STAGE)),
        ParticleSimScheduleSet::PostFinalIntegration,
    );

    app.start();
}

/// Check if particles have settled (KE near zero) and remove the blocker wall.
fn check_settled(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    mut walls: ResMut<Walls>,
    mut run_control: ResMut<SchedulerManager>,
) {
    let step = run_state.total_cycle;
    // Wait at least 1000 steps for particles to start moving, then check every 100 steps
    if step < 1000 || step % 100 != 0 {
        return;
    }

    let nlocal = atoms.nlocal as usize;
    let local_ke: f64 = (0..nlocal)
        .map(|i| {
            let vx = atoms.vel[i][0] as f64;
            let vy = atoms.vel[i][1] as f64;
            let vz = atoms.vel[i][2] as f64;
            0.5 * atoms.mass[i] as f64 * (vx * vx + vy * vy + vz * vz)
        })
        .sum();
    let global_ke = comm.all_reduce_sum_f64(local_ke);

    if global_ke < 1e-5 {
        walls.deactivate_by_name("blocker");
        run_control.advance_requested = true;
        if comm.rank() == 0 {
            println!(
                "Step {}: KE = {:.3e} J — particles settled, removing blocker wall",
                step, global_ke
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn referenced_stage_exists_in_all_configs() {
        for source in [
            include_str!("config.toml"),
            include_str!("validate_config.toml"),
            include_str!("validate_long_config.toml"),
        ] {
            let config = Config::from_str(source);
            let stages = RunConfig::from_config(&config);
            assert!(stages
                .stages
                .iter()
                .any(|stage| stage.name.as_deref() == Some(FILLING_STAGE)));
        }
    }
}
