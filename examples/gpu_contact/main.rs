//! One DEM run on either contact force, writing the final state for comparison.
//!
//! ```sh
//! cargo build --release --features gpu --example gpu_contact
//! DIRT_WRITE_EVERY=1000 DIRT_CONTACT=cpu target/release/examples/gpu_contact examples/gpu_contact/config.toml
//! DIRT_WRITE_EVERY=1000 DIRT_CONTACT=gpu target/release/examples/gpu_contact examples/gpu_contact/config.toml
//! ```
//!
//! `validate.py` checks results against the CPU, and `speed.py` measures steps per second.
//! ```sh
//! python3 examples/gpu_contact/validate.py
//! ```

use std::fs::File;
use std::io::Write as IoWrite;

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use dirt_gpu::{GranularGpuPlugins, GranularResidentPlugins};

/// `cpu`, `gpu` or `resident`, from the environment since the command line belongs to the config.
fn mode() -> String {
    std::env::var("DIRT_CONTACT").unwrap_or_else(|_| "cpu".to_string())
}

/// Writes every local atom's state to `<mode>.rank<r>.csv` on the last step.
fn write_final(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    run: Res<RunState>,
    comm: Res<CommResource>,
) {
    // The last write of the run is what stays on disk. Off unless asked for, because at
    // large sizes it would swamp a speed measurement.
    let every = std::env::var("DIRT_WRITE_EVERY").map_or(0, |v| v.parse().unwrap_or(0));
    if every == 0 || (run.total_cycle + 1) % every != 0 {
        return;
    }
    let dem = registry.expect::<DemAtom>("write_final");
    let name = format!("{}.rank{}.csv", mode(), comm.rank());
    let mut out = File::create(&name).expect("create output");
    writeln!(out, "tag,x,y,z,vx,vy,vz,wx,wy,wz").unwrap();
    for i in 0..atoms.nlocal() as usize {
        let (p, v, w) = (atoms.pos[i], atoms.vel[i], dem.omega[i]);
        writeln!(
            out,
            "{},{},{},{},{},{},{},{},{},{}",
            atoms.tag[i], p[0], p[1], p[2], v[0], v[1], v[2], w[0], w[1], w[2]
        )
        .unwrap();
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins);
    match mode().as_str() {
        "gpu" => {
            app.add_plugins(GranularGpuPlugins);
        }
        "resident" => {
            let window = std::env::var("DIRT_WINDOW").map_or(100, |v| v.parse().expect("DIRT_WINDOW"));
            app.add_plugins(GranularResidentPlugins { window });
        }
        _ => {
            app.add_plugins(GranularDefaultPlugins);
        }
    }
    app.add_plugins(GravityPlugin).add_plugins(WallPlugin);
    app.add_update_system(write_final, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}
