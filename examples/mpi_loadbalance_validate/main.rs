//! mpi_loadbalance_validate — verify dynamic load balancing (roadmap step 5) is
//! decomposition-invariant: a 2-rank run with the periodic `LoadBalancePlugin`
//! must reproduce a 2-rank static run (and a 1-rank run) bit-for-bit by global
//! tag. Periodic frictional gas (no walls/gravity) so there's no atom loss and
//! the physics is purely decomposition-dependent.
//!
//! ```bash
//! # static (no rebalance) vs rebalanced — compare by global tag
//! DIRT_LB_EVERY=0   mpiexec -n 2 target/release/examples/mpi_loadbalance_validate <cfg_n2>
//! DIRT_LB_EVERY=50  mpiexec -n 2 target/release/examples/mpi_loadbalance_validate <cfg_n2>
//! ```
use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

struct DumpInfo {
    every: usize,
}

fn main() {
    let every: usize = std::env::var("DIRT_LB_EVERY")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);

    let mut app = App::new();
    app.add_plugins(CorePlugins).add_plugins(GranularDefaultPlugins);
    if every > 0 {
        // Rebalance the x-decomposition by particle count every `every` steps.
        app.add_plugins(LoadBalancePlugin { every, nbins: 64 });
    }
    app.add_resource(DumpInfo { every });
    app.add_update_system(dump_final_state, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

fn dump_final_state(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    info: Res<DumpInfo>,
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
    let dir = "examples/mpi_loadbalance_validate/data";
    fs::create_dir_all(dir).ok();
    let path = format!("{dir}/state_e{}_n{ranks}_rank{rank}.csv", info.every);
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
