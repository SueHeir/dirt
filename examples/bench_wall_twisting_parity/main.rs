//! bench_wall_twisting_parity — one-step wall twisting torque parity check.
//!
//! The example places one sphere at the same overlap and local normal against a
//! plane wall, cylinder wall, sphere wall, and spherical region wall. It gives
//! the sphere a pure spin about that normal and records the wall twisting torque.
//! The expected torque is the plane-wall law `tau = mu_tw |F_n| R*`.

use dirt_core::dirt_atom::DemAtom;
use dirt_core::dirt_wall::{wall_contact_force, WallCylinder, WallRegion, WallSphere};
use dirt_core::prelude::*;
use dirt_core::soil_core::region::Region;
use serde::Deserialize;
use std::fs;
use std::io::Write as IoWrite;

#[derive(Deserialize)]
struct ConfigFile {
    material: MaterialCfg,
    scenario: ScenarioCfg,
}

#[derive(Deserialize)]
struct MaterialCfg {
    name: String,
    youngs_mod: f64,
    poisson_ratio: f64,
    restitution: f64,
    twisting_friction: f64,
}

#[derive(Clone, Deserialize)]
struct ScenarioCfg {
    radius: f64,
    wall_radius: f64,
    overlap: f64,
    omega: f64,
    dt: f64,
    out_dir: String,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let config_path = args
        .get(1)
        .map(String::as_str)
        .unwrap_or("examples/bench_wall_twisting_parity/config.toml");
    let cfg: ConfigFile = toml::from_str(&fs::read_to_string(config_path).expect("read config"))
        .expect("parse config");

    let cases = [
        ("plane", make_plane_walls(&cfg.scenario)),
        ("cylinder", make_cylinder_walls(&cfg.scenario)),
        ("sphere", make_sphere_walls(&cfg.scenario)),
        ("region", make_region_walls(&cfg.scenario)),
    ];

    let data_dir = format!("{}/data", cfg.scenario.out_dir);
    fs::create_dir_all(&data_dir).expect("create data dir");
    let path = format!("{data_dir}/wall_twisting_parity.csv");
    let mut f = fs::File::create(&path).expect("create results");
    writeln!(
        f,
        "geometry,torque_x,torque_y,torque_z,expected_tau_x,rel_err"
    )
    .unwrap();

    let mut max_rel_err: f64 = 0.0;
    for (name, walls) in cases {
        let torque = run_case(&cfg.material, &cfg.scenario, walls);
        let expected = expected_tau_x(&cfg.material, &cfg.scenario);
        let rel_err = ((torque[0] - expected) / expected).abs();
        max_rel_err = max_rel_err.max(rel_err);
        writeln!(
            f,
            "{name},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            torque[0], torque[1], torque[2], expected, rel_err
        )
        .unwrap();
    }

    println!("wall_twisting_parity results: {path}");
    println!("max_rel_err = {:.6e}", max_rel_err);
}

fn run_case(mat: &MaterialCfg, sc: &ScenarioCfg, walls: Walls) -> [f64; 3] {
    let mut app = App::new();
    app.add_resource(make_atom(sc));
    app.add_resource(make_registry(sc));
    app.add_resource(make_material_table(mat));
    app.add_resource(walls);
    app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let torque = registry
        .expect::<DemAtom>("bench_wall_twisting_parity")
        .torque[0];
    torque
}

fn make_atom(sc: &ScenarioCfg) -> Atom {
    let density = 2500.0;
    let mass = density * 4.0 / 3.0 * std::f64::consts::PI * sc.radius.powi(3);
    let mut atom = Atom::new();
    atom.dt = sc.dt;
    atom.push_test_atom(
        1,
        [sc.wall_radius - sc.radius + sc.overlap, 0.0, 0.0],
        sc.radius,
        mass,
    );
    atom.nlocal = 1;
    atom.natoms = 1;
    atom
}

fn make_registry(sc: &ScenarioCfg) -> AtomDataRegistry {
    let density = 2500.0;
    let mass = density * 4.0 / 3.0 * std::f64::consts::PI * sc.radius.powi(3);
    let mut dem = DemAtom::new();
    dem.radius.push(sc.radius);
    dem.density.push(density);
    dem.inv_inertia
        .push(1.0 / (0.4 * mass * sc.radius * sc.radius));
    dem.quaternion.push([1.0, 0.0, 0.0, 0.0]);
    dem.omega.push([sc.omega, 0.0, 0.0]);
    dem.ang_mom.push([0.0; 3]);
    dem.torque.push([0.0; 3]);
    dem.body_id.push(0.0);
    let mut registry = AtomDataRegistry::new();
    registry.register(dem);
    registry
}

fn make_material_table(mat: &MaterialCfg) -> MaterialTable {
    let mut mt = MaterialTable::new();
    mt.add_material_extended(
        &mat.name,
        mat.youngs_mod,
        mat.poisson_ratio,
        mat.restitution,
        0.0,
        0.0,
        0.0,
        0.0,
        mat.twisting_friction,
        0.0,
        0.0,
    );
    mt.build_pair_tables();
    mt
}

fn expected_tau_x(mat: &MaterialCfg, sc: &ScenarioCfg) -> f64 {
    let e_eff = mat.youngs_mod / (2.0 * (1.0 - mat.poisson_ratio * mat.poisson_ratio));
    let sdr = (sc.overlap * sc.radius).sqrt();
    let k_n = 4.0 / 3.0 * e_eff * sdr;
    mat.twisting_friction * (k_n * sc.overlap).abs() * sc.radius
}

fn make_empty_walls() -> Walls {
    Walls {
        planes: Vec::new(),
        active: Vec::new(),
        cylinders: Vec::new(),
        cylinder_active: Vec::new(),
        spheres: Vec::new(),
        sphere_active: Vec::new(),
        regions: Vec::new(),
        region_active: Vec::new(),
        time: 0.0,
        tangential_springs: std::collections::HashMap::new(),
        rolling_springs: std::collections::HashMap::new(),
    }
}

fn make_plane_walls(sc: &ScenarioCfg) -> Walls {
    let mut walls = make_empty_walls();
    walls.planes.push(WallPlane {
        point_x: sc.wall_radius,
        point_y: 0.0,
        point_z: 0.0,
        normal_x: -1.0,
        normal_y: 0.0,
        normal_z: 0.0,
        material_index: 0,
        name: None,
        bound_x_low: f64::NEG_INFINITY,
        bound_x_high: f64::INFINITY,
        bound_y_low: f64::NEG_INFINITY,
        bound_y_high: f64::INFINITY,
        bound_z_low: f64::NEG_INFINITY,
        bound_z_high: f64::INFINITY,
        velocity: [0.0; 3],
        motion: WallMotion::Static,
        origin: [sc.wall_radius, 0.0, 0.0],
        force_accumulator: 0.0,
        temperature: None,
    });
    walls.active.push(true);
    walls
}

fn make_cylinder_walls(sc: &ScenarioCfg) -> Walls {
    let mut walls = make_empty_walls();
    walls.cylinders.push(WallCylinder {
        axis: 2,
        center: [0.0, 0.0],
        radius: sc.wall_radius,
        lo: -0.01,
        hi: 0.01,
        inside: true,
        material_index: 0,
        name: None,
        force_accumulator: 0.0,
        temperature: None,
    });
    walls.cylinder_active.push(true);
    walls
}

fn make_sphere_walls(sc: &ScenarioCfg) -> Walls {
    let mut walls = make_empty_walls();
    walls.spheres.push(WallSphere {
        center: [0.0, 0.0, 0.0],
        radius: sc.wall_radius,
        inside: true,
        material_index: 0,
        name: None,
        force_accumulator: 0.0,
        temperature: None,
    });
    walls.sphere_active.push(true);
    walls
}

fn make_region_walls(sc: &ScenarioCfg) -> Walls {
    let mut walls = make_empty_walls();
    walls.regions.push(WallRegion {
        region: Region::Sphere {
            center: [0.0, 0.0, 0.0],
            radius: sc.wall_radius,
        },
        inside: true,
        material_index: 0,
        name: None,
        force_accumulator: 0.0,
        temperature: None,
    });
    walls.region_active.push(true);
    walls
}
