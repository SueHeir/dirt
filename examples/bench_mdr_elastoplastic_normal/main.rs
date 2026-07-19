//! Quasi-static MDR elastic-plastic normal-contact trace.
//!
//! The executable reads a declarative TOML file, drives two identical spheres
//! through a prescribed loading/unloading overlap path, evaluates DIRT's
//! `contact_model = "mdr"` force through the normal fused contact loop, and
//! writes the measured force-displacement trace for `sweep.py` to validate.

use dirt_core::dirt_atom::DemAtom;
use dirt_core::dirt_granular::contact::{contact_force_core, ForcePass};
use dirt_core::dirt_granular::tangential::ContactHistoryStore;
use dirt_core::prelude::*;
use serde::Deserialize;
use std::env;
use std::fs;
use std::io::Write;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BenchConfig {
    radius: f64,
    density: f64,
    youngs_mod: f64,
    poisson_ratio: f64,
    yield_stress: f64,
    surface_energy: f64,
    damping: f64,
    max_overlap: f64,
    points_per_leg: usize,
    output_csv: String,
}

fn push_sphere(atom: &mut Atom, dem: &mut DemAtom, tag: u32, x: f64, radius: f64, density: f64) {
    let mass = density * 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
    // Each raw core row below is paired with every DemAtom row before use.
    unsafe { atom.push_test_atom(tag, [x, 0.0, 0.0], radius, mass) };
    dem.radius.push(radius);
    dem.density.push(density);
    dem.inv_inertia.push(1.0 / (0.4 * mass * radius * radius));
    dem.quaternion.push([1.0, 0.0, 0.0, 0.0]);
    dem.omega.push([0.0; 3]);
    dem.ang_mom.push([0.0; 3]);
    dem.torque.push([0.0; 3]);
    dem.body_id.push(0.0);
}

fn main() {
    let config_path = env::args()
        .nth(1)
        .unwrap_or_else(|| "examples/bench_mdr_elastoplastic_normal/config.toml".to_string());
    let text = fs::read_to_string(&config_path)
        .unwrap_or_else(|e| panic!("cannot read {config_path}: {e}"));
    let cfg: BenchConfig = toml::from_str(&text).unwrap_or_else(|e| panic!("bad config: {e}"));

    let mut atom = Atom::new();
    atom.dt = 1.0;
    atom.nlocal = 2;
    atom.natoms = 2;
    let mut dem = DemAtom::new();
    push_sphere(&mut atom, &mut dem, 1, 0.0, cfg.radius, cfg.density);
    push_sphere(
        &mut atom,
        &mut dem,
        2,
        2.0 * cfg.radius,
        cfg.radius,
        cfg.density,
    );

    let mut history = ContactHistoryStore::new();
    history.contacts.resize_with(2, Vec::new);

    let mut registry = AtomDataRegistry::new();
    registry.try_register(dem, atom.len()).unwrap();
    registry.try_register(history, atom.len()).unwrap();

    let mut neighbor = Neighbor::new();
    neighbor.neighbor_offsets = vec![0, 1, 1];
    neighbor.neighbor_indices = vec![1];
    neighbor.newton = true;

    let mut mt = MaterialTable::new();
    mt.contact_model = "mdr".to_string();
    mt.limit_damping = false;
    mt.add(
        Material::new(
            "powder",
            Elastic::new(cfg.youngs_mod, cfg.poisson_ratio, 1.0),
        )
        .with_friction(Friction {
            sliding: 0.0,
            ..Friction::default()
        })
        .with_adhesion(Adhesion::SurfaceEnergy {
            energy: cfg.surface_energy,
        })
        .with_mdr(Mdr {
            yield_stress: cfg.yield_stress,
            psi_b: 0.5,
            damping: cfg.damping,
        }),
    )
    .unwrap();
    mt.build_pair_tables();

    if let Some(parent) = std::path::Path::new(&cfg.output_csv).parent() {
        fs::create_dir_all(parent).unwrap();
    }
    let mut out = fs::File::create(&cfg.output_csv).unwrap();
    writeln!(out, "phase,delta,force").unwrap();

    let n = cfg.points_per_leg.max(2);
    for leg in 0..2 {
        for k in 0..n {
            let frac = k as f64 / (n - 1) as f64;
            let delta = if leg == 0 {
                frac * cfg.max_overlap
            } else {
                (1.0 - frac) * cfg.max_overlap
            };
            atom.pos[0] = [0.0, 0.0, 0.0];
            atom.pos[1] = [(2.0 * cfg.radius - delta) as _, 0.0, 0.0];
            for force in atom.force.iter_mut() {
                *force = [0.0; 3];
            }
            {
                let mut dem = registry.expect_mut::<DemAtom>("bench_mdr");
                dem.torque.fill([0.0; 3]);
            }
            contact_force_core(&mut atom, &neighbor, &registry, &mt, None, ForcePass::All);
            let force = -(atom.force[0][0] as f64);
            let phase = if leg == 0 { "loading" } else { "unloading" };
            writeln!(out, "{phase},{delta:.12e},{force:.12e}").unwrap();
        }
    }
    println!("MDR trace written to {}", cfg.output_csv);
}
