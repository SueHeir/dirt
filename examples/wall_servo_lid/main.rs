//! wall_servo_lid — a servo-controlled lid compresses a particle bed, then is released.
//!
//! A box of glass spheres settles under gravity, then a **servo** plane wall
//! (the "lid") drives its velocity to reach a target downward contact force,
//! compressing the bed to a prescribed pressure. After the bed is compacted the
//! lid is removed at runtime with [`Walls::deactivate_by_name`], letting the bed
//! relax.
//!
//! Demonstrates:
//! - a servo wall (`servo = { target_force, max_velocity, gain }`),
//! - runtime wall control via `deactivate_by_name`,
//! - a multi-stage run keyed on settling.
//!
//! ```bash
//! cargo run --release --example wall_servo_lid --no-default-features -- \
//!     examples/wall_servo_lid/config.toml
//! ```

use dirt_core::prelude::*;

const COMPRESS_STAGE: &str = "compress";

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin);

    // After enough steps under the servo lid, release the lid and relax.
    app.add_update_system(
        release_lid.run_if(in_stage(COMPRESS_STAGE)),
        ParticleSimScheduleSet::PostFinalIntegration,
    );

    app.start();
}

/// Once the compression stage has run long enough, deactivate the servo lid by
/// name and advance to the relax stage.
fn release_lid(
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    mut walls: ResMut<Walls>,
    mut run_control: ResMut<SchedulerManager>,
) {
    let step = run_state.total_cycle;
    if step < 60_000 {
        return;
    }
    walls.deactivate_by_name("lid");
    run_control.advance_requested = true;
    if comm.rank() == 0 {
        println!("Step {step}: servo lid released, bed relaxing");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn referenced_stage_exists_in_default_config() {
        let config = Config::from_str(include_str!("config.toml"));
        let stages = RunConfig::from_config(&config);
        assert!(stages
            .stages
            .iter()
            .any(|stage| stage.name.as_deref() == Some(COMPRESS_STAGE)));
    }
}
