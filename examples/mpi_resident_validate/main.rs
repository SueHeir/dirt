//! mpi_resident_validate — MPI + the windowed-resident GPU stepper.
//!
//! Drives `GpuGranularResidentMpiPlugin` (ghost-aware, window K from `DIRT_WINDOW`)
//! through the real MPI schedule. Each schedule tick = one host forward_comm + K
//! device steps. Dumps final local state by global tag (mergeable across ranks).
//!
//! `forward_comm` refreshes ghosts EVERY step but this plugin refreshes them once
//! per tick (per K device steps), so window=1 is correct (ghosts fresh each step)
//! and window>1 freezes ghosts for K-1 steps → boundary divergence growing with K.
//!
//! ```bash
//! DIRT_WINDOW=1 mpiexec -n 2 target/release/examples/mpi_resident_validate examples/mpi_gpu_validate/config_n2.toml
//! ```
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

struct DumpState {
    window: usize,
}

fn main() {
    let window: usize = std::env::var("DIRT_WINDOW").ok().and_then(|s| s.parse().ok()).unwrap_or(1);

    // Box walls matching config domain [0,0.02]^2 x z>=0 (floor + 4 sides), on-device.
    let mut boundary = Boundary::new();
    boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]));
    boundary.push(Plane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]));
    boundary.push(Plane::new([0.02, 0.0, 0.0], [-1.0, 0.0, 0.0]));
    boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]));
    boundary.push(Plane::new([0.0, 0.02, 0.0], [0.0, -1.0, 0.0]));

    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(DemAtomPlugin)
        .add_plugins(DemAtomInsertPlugin)
        .add_plugins(GpuGranularResidentMpiPlugin { boundary, gravity: [0.0, 0.0, -40.0], window });
    app.add_resource(DumpState { window });
    app.add_update_system(dump_final_state, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

fn dump_final_state(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    dump: Res<DumpState>,
) {
    let steps_left: usize = run_state
        .cycle_remaining
        .iter()
        .zip(run_state.cycle_count.iter())
        .map(|(&target, &done)| (target as usize).saturating_sub(done as usize))
        .sum();
    if steps_left > 1 {
        return;
    }
    let ranks = comm.size();
    let rank = comm.rank();
    let dir = "examples/mpi_resident_validate/data";
    fs::create_dir_all(dir).ok();
    let path = format!("{dir}/state_w{}_n{ranks}_rank{rank}.csv", dump.window);
    let mut f = fs::File::create(&path).expect("cannot create state csv");
    writeln!(f, "tag,x,y,z").unwrap();
    for i in 0..atoms.nlocal as usize {
        writeln!(
            f,
            "{},{:.9e},{:.9e},{:.9e}",
            atoms.tag[i], atoms.pos[i][0], atoms.pos[i][1], atoms.pos[i][2]
        )
        .unwrap();
    }
}
