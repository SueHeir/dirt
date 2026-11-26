//! OVITO-native output for rigid sphero-cylinder particles.
//!
//! The generic SOIL dump pipeline handles MPI gathering and output cadence.
//! This module contributes capsule geometry/orientation columns and the
//! `ovito_spherocylinder` format. Its column names match OVITO's LAMMPS importer:
//! `AsphericalShape.X/Y/Z` and `quati/quatj/quatk/quatw`.

use std::fs::{self, File};
use std::io::{BufWriter, Error, ErrorKind, Write};
use std::path::Path;

use dirt_atom::DemAtom;
use grass_app::prelude::*;
use grass_scheduler::prelude::Res;
use soil_core::{CommResource, Input, ScheduleSetupSet};
use soil_print::{DumpConfig, DumpFrame, DumpRegistry};

use crate::CapsuleAtom;

const SHAPE_COLUMN: &str = "__capsule_aspherical_shape";
const QUATERNION_XYZ_COLUMN: &str = "__capsule_quaternion_xyz";
const QUATERNION_W_COLUMN: &str = "__capsule_quaternion_w";
const MOLECULE_COLUMN: &str = "__capsule_molecule";
const INTERNAL_PREFIX: &str = "__capsule_";

/// Register OVITO sphero-cylinder columns and writer with SOIL's dump pipeline.
///
/// `CorePlugins`/`PrintPlugin` must be installed before `CapsulePlugin`, as in
/// the standard DIRT application setup. If no dump registry exists, capsule
/// mechanics still work and only this optional output integration is skipped.
pub(crate) fn register_capsule_ovito_output(app: &mut App) {
    let Some(dump) = app.get_resource_ref::<DumpRegistry>() else {
        return;
    };

    dump.register_vector(SHAPE_COLUMN, |atoms, registry| {
        let dem = registry.expect::<DemAtom>("capsule OVITO shape");
        let capsules = registry.expect::<CapsuleAtom>("capsule OVITO shape");
        let nlocal = atoms.nlocal() as usize;
        (0..nlocal)
            .map(|i| {
                let radius = dem.radius[i];
                // Match TANGLE's native sphero-cylinder convention: X/Y are
                // radius and Z is the full centerline/cylindrical length.
                [radius, radius, 2.0 * capsules.half_length[i]]
            })
            .collect()
    });
    dump.register_vector(QUATERNION_XYZ_COLUMN, |atoms, registry| {
        let dem = registry.expect::<DemAtom>("capsule OVITO orientation");
        let nlocal = atoms.nlocal() as usize;
        (0..nlocal)
            .map(|i| {
                let q = dem.quaternion[i]; // DIRT: [w, x, y, z]
                [q[1], q[2], q[3]]
            })
            .collect()
    });
    dump.register_scalar(QUATERNION_W_COLUMN, |atoms, registry| {
        let dem = registry.expect::<DemAtom>("capsule OVITO orientation");
        (0..atoms.nlocal() as usize)
            .map(|i| dem.quaternion[i][0])
            .collect()
    });
    dump.register_scalar(MOLECULE_COLUMN, |atoms, registry| {
        let dem = registry.expect::<DemAtom>("capsule OVITO molecule");
        dem.molecule_id[..atoms.nlocal() as usize].to_vec()
    });
    dump.register_format("ovito_spherocylinder", write_ovito_spherocylinder);
    drop(dump);

    app.add_setup_system(write_view_script, ScheduleSetupSet::PostSetup);
}

fn required_vector<'a>(frame: &'a DumpFrame, name: &str) -> std::io::Result<&'a [[f64; 3]]> {
    frame
        .vectors
        .iter()
        .find(|(candidate, _)| candidate == name)
        .map(|(_, values)| values.as_slice())
        .ok_or_else(|| {
            Error::new(
                ErrorKind::InvalidData,
                format!("missing dump column {name}"),
            )
        })
}

fn required_scalar<'a>(frame: &'a DumpFrame, name: &str) -> std::io::Result<&'a [f64]> {
    frame
        .scalars
        .iter()
        .find(|(candidate, _)| candidate == name)
        .map(|(_, values)| values.as_slice())
        .ok_or_else(|| {
            Error::new(
                ErrorKind::InvalidData,
                format!("missing dump column {name}"),
            )
        })
}

/// Write one MPI-assembled OVITO/LAMMPS sphero-cylinder frame.
///
/// The output uses OVITO-recognized property names. Selecting the
/// `Spherocylinder` particle shape therefore produces continuous rigid capsules
/// with the same length convention as TANGLE's trajectory writer.
pub fn write_ovito_spherocylinder(frame: &DumpFrame) -> std::io::Result<()> {
    let shape = required_vector(frame, SHAPE_COLUMN)?;
    let quaternion_xyz = required_vector(frame, QUATERNION_XYZ_COLUMN)?;
    let quaternion_w = required_scalar(frame, QUATERNION_W_COLUMN)?;
    let molecule = required_scalar(frame, MOLECULE_COLUMN)?;
    let path = format!("{}.lammpstrj", frame.path_stem);
    let mut writer = BufWriter::new(File::create(path)?);
    let boundary = |axis: usize| {
        if frame.box_info.periodic[axis] {
            "pp"
        } else {
            "ff"
        }
    };

    writeln!(writer, "ITEM: TIMESTEP")?;
    writeln!(writer, "{}", frame.step)?;
    writeln!(writer, "ITEM: NUMBER OF ATOMS")?;
    writeln!(writer, "{}", frame.n())?;
    writeln!(
        writer,
        "ITEM: BOX BOUNDS {} {} {}",
        boundary(0),
        boundary(1),
        boundary(2)
    )?;
    for axis in 0..3 {
        writeln!(
            writer,
            "{} {}",
            frame.box_info.low[axis], frame.box_info.high[axis]
        )?;
    }

    write!(
        writer,
        "ITEM: ATOMS id mol type AsphericalShape.X AsphericalShape.Y AsphericalShape.Z quati quatj quatk quatw x y z"
    )?;
    for (name, _) in &frame.scalars {
        if !name.starts_with(INTERNAL_PREFIX) {
            write!(writer, " {name}")?;
        }
    }
    for (name, _) in &frame.vectors {
        if !name.starts_with(INTERNAL_PREFIX) {
            write!(writer, " {name}.X {name}.Y {name}.Z")?;
        }
    }
    if frame.is_ghost.is_some() {
        write!(writer, " IsGhost")?;
    }
    writeln!(writer)?;

    for i in 0..frame.n() {
        let q = quaternion_xyz[i];
        let s = shape[i];
        let p = frame.pos[i];
        write!(
            writer,
            "{} {} {} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}",
            frame.tag[i],
            molecule[i].round() as u64,
            frame.atom_type[i] + 1,
            s[0],
            s[1],
            s[2],
            q[0],
            q[1],
            q[2],
            quaternion_w[i],
            p[0],
            p[1],
            p[2],
        )?;
        for (name, data) in &frame.scalars {
            if !name.starts_with(INTERNAL_PREFIX) {
                write!(writer, " {:.17e}", data[i])?;
            }
        }
        for (name, data) in &frame.vectors {
            if !name.starts_with(INTERNAL_PREFIX) {
                let value = data[i];
                write!(
                    writer,
                    " {:.17e} {:.17e} {:.17e}",
                    value[0], value[1], value[2]
                )?;
            }
        }
        if let Some(ghosts) = &frame.is_ghost {
            write!(writer, " {}", u8::from(ghosts[i]))?;
        }
        writeln!(writer)?;
    }
    writer.flush()
}

fn write_view_script(dump: Res<DumpConfig>, input: Res<Input>, comm: Res<CommResource>) {
    if comm.rank() != 0 || dump.interval == 0 || dump.format != "ovito_spherocylinder" {
        return;
    }
    let output = input.output_dir.as_deref().unwrap_or(".");
    let directory = Path::new(output).join("dump");
    fs::create_dir_all(&directory).expect("create capsule OVITO dump directory");
    let path = directory.join("view_spherocylinders.py");
    let mut writer = BufWriter::new(File::create(&path).expect("create OVITO view script"));
    writeln!(writer, "# Run with: ovitos {}", path.display()).unwrap();
    writeln!(writer, "from pathlib import Path").unwrap();
    writeln!(writer, "import ovito").unwrap();
    writeln!(writer, "from ovito.io import import_file").unwrap();
    writeln!(writer, "from ovito.vis import ParticlesVis").unwrap();
    writeln!(writer, "HERE = Path(__file__).resolve().parent").unwrap();
    writeln!(
        writer,
        "pipeline = import_file(str(HERE / 'dump_*.lammpstrj'), multiple_frames=True, sort_particles=True)"
    )
    .unwrap();
    writeln!(writer, "pipeline.add_to_scene(name='DIRT capsules')").unwrap();
    writeln!(writer, "data = pipeline.compute()").unwrap();
    writeln!(
        writer,
        "data.particles.vis.shape = ParticlesVis.Shape.Spherocylinder"
    )
    .unwrap();
    writeln!(writer, "data.cell.vis.enabled = True").unwrap();
    writeln!(
        writer,
        "ovito.scene.save(str(HERE / 'spherocylinders.ovito'))"
    )
    .unwrap();
    writeln!(
        writer,
        "print(f\"Saved OVITO session: {{HERE / 'spherocylinders.ovito'}}\")"
    )
    .unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    use soil_print::BoxInfo;

    #[test]
    fn writes_tangle_compatible_spherocylinder_columns_and_full_length() {
        let stem = std::env::temp_dir().join(format!("dirt-capsule-ovito-{}", std::process::id()));
        let frame = DumpFrame {
            step: 12,
            box_info: BoxInfo {
                low: [0.0; 3],
                high: [2.0; 3],
                periodic: [true, false, false],
            },
            path_stem: stem.to_string_lossy().into_owned(),
            tag: vec![7],
            atom_type: vec![1],
            pos: vec![[0.5, 1.0, 1.5]],
            scalars: vec![
                (QUATERNION_W_COLUMN.to_string(), vec![1.0]),
                (MOLECULE_COLUMN.to_string(), vec![3.0]),
            ],
            vectors: vec![
                (SHAPE_COLUMN.to_string(), vec![[0.1, 0.1, 1.5]]),
                (QUATERNION_XYZ_COLUMN.to_string(), vec![[0.0; 3]]),
            ],
            ..Default::default()
        };

        write_ovito_spherocylinder(&frame).unwrap();
        let path = stem.with_extension("lammpstrj");
        let text = fs::read_to_string(&path).unwrap();
        fs::remove_file(path).unwrap();
        assert!(text.contains(
            "AsphericalShape.X AsphericalShape.Y AsphericalShape.Z quati quatj quatk quatw"
        ));
        let row = text
            .lines()
            .last()
            .unwrap()
            .split_whitespace()
            .collect::<Vec<_>>();
        assert_eq!(row[0], "7");
        assert_eq!(row[1], "3");
        assert_eq!(row[2], "2");
        assert!((row[5].parse::<f64>().unwrap() - 1.5).abs() < 1.0e-12);
    }
}
