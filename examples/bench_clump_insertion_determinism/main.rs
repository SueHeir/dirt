//! bench_clump_insertion_determinism — config-level clump insertion fingerprint.
//!
//! Runs the normal DIRT plugin stack through `ClumpPlugin`, including
//! `clump_insert_atoms`, then writes the inserted atom/body state to a CSV file.
//! The companion `sweep.py` runs the same config twice and a changed-seed config
//! once, proving that `[[clump.insert]].seed` makes the production setup path
//! byte-stable while still changing the insertion stream when the seed changes.

use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};

use dirt_core::prelude::*;

fn main() -> std::io::Result<()> {
    let output = env::args()
        .nth(2)
        .unwrap_or_else(|| "examples/bench_clump_insertion_determinism/data/state.csv".to_string());

    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(ClumpPlugin);
    app.start();

    write_state(&app, &output)
}

fn write_state(app: &App, output: &str) -> std::io::Result<()> {
    let atoms = app
        .get_resource_ref::<Atom>()
        .expect("Atom resource should exist after setup");
    let bodies = app
        .get_resource_ref::<MultisphereBodyStore>()
        .expect("MultisphereBodyStore should exist after setup");

    if let Some(parent) = std::path::Path::new(output).parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut out = BufWriter::new(File::create(output)?);
    writeln!(out, "kind,index,field,x,y,z,w")?;
    for i in 0..atoms.nlocal as usize {
        writeln!(
            out,
            "atom,{i},pos,{},{},{},",
            atoms.pos[i][0].to_bits(),
            atoms.pos[i][1].to_bits(),
            atoms.pos[i][2].to_bits()
        )?;
        writeln!(
            out,
            "atom,{i},vel,{},{},{},",
            atoms.vel[i][0].to_bits(),
            atoms.vel[i][1].to_bits(),
            atoms.vel[i][2].to_bits()
        )?;
    }
    for (i, body) in bodies.bodies.iter().enumerate() {
        writeln!(
            out,
            "body,{i},com_pos,{},{},{},",
            body.com_pos[0].to_bits(),
            body.com_pos[1].to_bits(),
            body.com_pos[2].to_bits()
        )?;
        writeln!(
            out,
            "body,{i},com_vel,{},{},{},",
            body.com_vel[0].to_bits(),
            body.com_vel[1].to_bits(),
            body.com_vel[2].to_bits()
        )?;
        writeln!(
            out,
            "body,{i},quat,{},{},{},{}",
            body.quaternion[0].to_bits(),
            body.quaternion[1].to_bits(),
            body.quaternion[2].to_bits(),
            body.quaternion[3].to_bits()
        )?;
    }
    Ok(())
}
