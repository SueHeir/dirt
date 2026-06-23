//! mpi_gpu_validate — correctness of MPI domain decomposition + GPU contact force.
//!
//! Runs a dense settling bed under gravity with the particle–particle contact
//! force on either the CPU (`GranularDefaultPlugins`) or the GPU
//! (`GranularGpuPlugins`), selected by the `DIRT_FORCE` env var (`cpu`|`gpu`,
//! default `gpu`). Walls and gravity are CPU in BOTH, so the only thing that
//! differs is the particle–particle contact — isolating the GPU kernel.
//!
//! At the final step each rank writes its LOCAL atoms (global tag + position) to
//! `data/state_<force>_n<ranks>_rank<rank>.csv`. A merge of all ranks' files,
//! sorted by global tag, is directly comparable across rank counts — so a 2-rank
//! decomposition can be diffed against a 1-rank run (decomposition invariance)
//! and against the CPU reference.
//!
//! ```bash
//! cargo build --release --example mpi_gpu_validate
//! DIRT_FORCE=cpu mpiexec -n 1 target/release/examples/mpi_gpu_validate examples/mpi_gpu_validate/config_n1.toml
//! DIRT_FORCE=gpu mpiexec -n 2 target/release/examples/mpi_gpu_validate examples/mpi_gpu_validate/config_n2.toml
//! ```

use dirt_core::prelude::*;
use std::fs;
use std::io::Write as IoWrite;

struct DumpState {
    force: String,
}

fn main() {
    let force = std::env::var("DIRT_FORCE").unwrap_or_else(|_| "gpu".into());

    let mut app = App::new();
    app.add_plugins(CorePlugins);
    if force == "cpu" {
        app.add_plugins(GranularDefaultPlugins);
    } else {
        app.add_plugins(GranularGpuPlugins);
    }
    app.add_plugins(GravityPlugin).add_plugins(WallPlugin);

    app.add_resource(DumpState { force });
    app.add_update_system(dump_final_state, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// On the final step, each rank writes its local atoms (global tag + position).
fn dump_final_state(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
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
    let output_dir = input
        .output_dir
        .as_deref()
        .unwrap_or("examples/mpi_gpu_validate");
    let data_dir = format!("{}/data", output_dir);
    fs::create_dir_all(&data_dir).ok();
    let path = format!("{}/state_{}_n{}_rank{}.csv", data_dir, dump.force, ranks, rank);
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
