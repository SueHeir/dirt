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

#[derive(Clone, Debug, PartialEq, Default, StageEnum)]
enum Phase {
    #[default]
    #[stage("compress")]
    Compress,
    #[stage("relax")]
    Relax,
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin)
        .add_plugins(StatesPlugin::new(
            Phase::Compress,
            ParticleSimScheduleSet::PostFinalIntegration,
        ))
        .add_plugins(StageAdvancePlugin::<Phase>::new(
            ParticleSimScheduleSet::PostFinalIntegration,
        ));

    // After enough steps under the servo lid, release the lid and relax.
    app.add_update_system(
        release_lid.run_if(in_state(Phase::Compress)),
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
    mut next_state: ResMut<NextState<Phase>>,
) {
    let step = run_state.total_cycle;
    if step < 60_000 {
        return;
    }
    walls.deactivate_by_name("lid");
    next_state.set(Phase::Relax);
    if comm.rank() == 0 {
        println!("Step {step}: servo lid released, bed relaxing");
    }
}
