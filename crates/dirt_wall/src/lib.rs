//! Wall contact forces for DEM simulations using Hertz or Hooke normal contact
//! with viscous damping, tangential/rolling/twisting friction, and optional
//! adhesion (JKR, DMT, SJKR cohesion).
//!
//! # Contact mechanics
//!
//! Wall contacts reuse the same per-pair mixing tables as particle–particle
//! contacts: the particle's material and the wall's `material` index a row of
//! [`MaterialTable`] (`e_eff_ij`, `g_eff_ij`, `beta_ij`, `friction_ij`,
//! `rolling_friction_ij`, `twisting_friction_ij`, …). Because a wall is treated
//! as infinitely massive and infinitely flat, the **effective radius is the
//! particle radius** (`R* = particle_radius`) and the reduced mass is the
//! particle mass.
//!
//! Beyond the normal force + damping, walls apply:
//!
//! - **Tangential (Mindlin sliding) friction** — incremental spring-history
//!   model with a per-contact tangential spring, Coulomb-capped at `μ |F_n|`.
//!   Supported by **all** wall types (plane, cylinder, sphere, region).
//! - **Rolling resistance** — `constant` (default) or `sds` (spring–dashpot–
//!   slider) model, mirroring the particle–particle rolling model with the wall
//!   as a zero-spin second body. Supported by **all** wall types.
//! - **Twisting friction** — constant-torque model opposing spin about the
//!   local contact normal. Supported by **all** wall types.
//!
//! Frictionless walls (`friction = 0`) are byte-for-byte unchanged from a
//! pure-normal contact. Tangential and rolling spring histories are stored on
//! the [`Walls`] resource (`tangential_springs`, `rolling_springs`).
//!
//! ## Adhesion-model asymmetry
//!
//! Adhesion support differs by wall geometry:
//!
//! - **Plane walls** support JKR and DMT (`surface_energy`) *and* SJKR cohesion
//!   (`cohesion_energy`), including the JKR extended-range pull-off regime.
//! - **Cylinder, sphere, and region walls** support **SJKR cohesion only**
//!   (`cohesion_energy`); their `surface_energy` is not consulted, so JKR/DMT
//!   pull-off is unavailable on curved/region walls. Use a plane wall if you
//!   need JKR/DMT against a wall.
//!
//! ## Wall temperature
//!
//! Every wall config accepts an optional `temperature` (K). This crate
//! **stores** it on the wall but never reads it — it is a hook for an external
//! heat-transfer system (e.g. a thermal-conduction plugin) to consult a wall's
//! temperature. It has no effect on the contact force computed here.
//!
//! # Wall Types
//!
//! | Type | Description | Config key |
//! |------|-------------|------------|
//! | **Plane** | Infinite flat plane defined by a point and unit normal | `type = "plane"` |
//! | **Cylinder** | Infinite cylinder along X/Y/Z axis with finite axial bounds | `type = "cylinder"` |
//! | **Sphere** | Sphere defined by center and radius | `type = "sphere"` |
//! | **Region** | Any [`Region`] shape used as a wall surface | `type = "region"` |
//!
//! All wall types treat the wall as having infinite mass and infinite radius
//! for contact mechanics, so the effective radius equals the particle radius
//! and the reduced mass equals the particle mass.
//!
//! # Motion Types
//!
//! | Motion | Description |
//! |--------|-------------|
//! | **Static** | Wall does not move (default) |
//! | **Constant velocity** | Wall translates at a fixed velocity each timestep |
//! | **Oscillating** | Sinusoidal displacement along the wall normal |
//! | **Servo** | Proportional controller adjusting velocity to reach a target force |
//!
//! Motion is currently supported only for plane walls.
//!
//! # TOML Configuration
//!
//! Walls are defined as `[[wall]]` array-of-tables entries. Each entry requires
//! a `material` field matching a name in `[[dem.materials]]`.
//!
//! ```toml
//! # Plane wall (floor at z=0, normal pointing up)
//! [[wall]]
//! type = "plane"
//! point_z = 0.0
//! normal_z = 1.0
//! material = "glass"
//! name = "floor"                  # optional, for runtime enable/disable
//!
//! # Cylinder wall (particles confined inside a z-aligned cylinder)
//! [[wall]]
//! type = "cylinder"
//! axis = "z"
//! center = [0.005, 0.005]         # center in the XY plane
//! radius = 0.004
//! lo = 0.0                        # axial lo bound (default: -inf)
//! hi = 0.01                       # axial hi bound (default: +inf)
//! inside = true                   # particles live inside the cylinder
//! material = "glass"
//!
//! # Sphere wall (particles confined inside a sphere)
//! [[wall]]
//! type = "sphere"
//! center = [0.005, 0.005, 0.005]
//! radius = 0.004
//! inside = true
//! material = "glass"
//!
//! # Region wall (any Region shape as a wall surface)
//! [[wall]]
//! type = "region"
//! inside = true
//! material = "glass"
//! region = { type = "cone", center = [0.005, 0.005], axis = "z",
//!            rad_lo = 0.004, rad_hi = 0.002, lo = 0.0, hi = 0.01 }
//!
//! # Moving wall with constant velocity
//! [[wall]]
//! type = "plane"
//! normal_z = 1.0
//! material = "glass"
//! velocity = [0.0, 0.0, -0.01]    # [vx, vy, vz]
//!
//! # Oscillating wall (sinusoidal along normal)
//! [[wall]]
//! type = "plane"
//! point_z = 0.1
//! normal_z = 1.0
//! material = "glass"
//! oscillate = { amplitude = 0.001, frequency = 50.0 }
//!
//! # Servo-controlled wall (adjusts velocity to reach target force)
//! [[wall]]
//! type = "plane"
//! point_z = 0.1
//! normal_z = -1.0
//! material = "glass"
//! servo = { target_force = 100.0, max_velocity = 0.1, gain = 0.001 }
//! ```
//!
//! # Plugin Registration
//!
//! Add [`WallPlugin`] to your app. It depends on `DemAtomPlugin` (for
//! [`MaterialTable`] and [`DemAtom`] data).
//!
//! Wall configuration is preflighted through [`Plugin::try_build`]. A malformed
//! `[[wall]]` entry is returned as an [`AppError`] to the outer runner, which
//! owns formatting and MPI-consistent termination.
//!
//! # Systems
//!
//! | System | Schedule | Purpose |
//! |--------|----------|---------|
//! | [`wall_move`] | `PreInitialIntegration` | Updates wall positions from motion modes |
//! | [`wall_zero_force_accumulators`] | `PreForce` | Zeros per-wall force accumulators |
//! | [`wall_contact_force`] | `Force` | Computes normal contact + damping + adhesion |

// Public API documentation-completeness gate: every public item in this crate
// must carry a doc comment. Enforced on both `cargo build` (rustc) and
// `cargo doc` (rustdoc; e.g. `RUSTDOCFLAGS="-D missing_docs"`). Document real
// API intent here — do not add empty doc comments just to satisfy the lint.
#![deny(missing_docs)]

mod config;
mod contact;
mod geometry;
mod motion;
mod plugin;

pub use config::{OscillateDef, ServoDef, WallDef};
pub use contact::wall_contact_force;
pub use geometry::{WallCylinder, WallPlane, WallRegion, WallSphere, Walls};
pub use motion::{wall_move, wall_zero_force_accumulators, WallMotion};
pub use plugin::WallPlugin;

#[cfg(test)]
use dirt_atom::MaterialTable;
#[cfg(test)]
use grass_app::prelude::*;
#[cfg(test)]
use soil_core::ParticleSimScheduleSet;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::curved_or_region_wall_surface_energy_warning;
    use crate::contact::wall_normal_force;
    use dirt_atom::DemAtom;
    use dirt_test_utils::{make_material_table, push_dem_test_atom, ParticleFixture, ParticleSpec};
    use soil_core::region::Region;
    use soil_core::{Atom, AtomDataRegistry};

    fn make_wall_plane(
        point_x: f64,
        point_y: f64,
        point_z: f64,
        normal_x: f64,
        normal_y: f64,
        normal_z: f64,
    ) -> WallPlane {
        let mag = (normal_x * normal_x + normal_y * normal_y + normal_z * normal_z).sqrt();
        WallPlane {
            point_x,
            point_y,
            point_z,
            normal_x: normal_x / mag,
            normal_y: normal_y / mag,
            normal_z: normal_z / mag,
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
            origin: [point_x, point_y, point_z],
            force_accumulator: 0.0,
            temperature: None,
        }
    }

    fn make_walls(planes: Vec<WallPlane>) -> Walls {
        let n = planes.len();
        Walls {
            planes,
            active: vec![true; n],
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

    fn wall_def(wall_type: &str, material: &str) -> WallDef {
        WallDef {
            wall_type: wall_type.to_string(),
            point_x: 0.0,
            point_y: 0.0,
            point_z: 0.0,
            normal_x: 0.0,
            normal_y: 0.0,
            normal_z: 1.0,
            axis: None,
            center: None,
            radius: None,
            lo: None,
            hi: None,
            inside: None,
            material: material.to_string(),
            name: None,
            bound_x_low: f64::NEG_INFINITY,
            bound_x_high: f64::INFINITY,
            bound_y_low: f64::NEG_INFINITY,
            bound_y_high: f64::INFINITY,
            bound_z_low: f64::NEG_INFINITY,
            bound_z_high: f64::INFINITY,
            velocity: None,
            oscillate: None,
            servo: None,
            region: None,
            temperature: None,
        }
    }

    fn material_table_with_surface_energy() -> MaterialTable {
        let mut mt = MaterialTable::new();
        mt.add_material_full("dry", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0, 0.0);
        mt.add_material_full("sticky", 8.7e9, 0.3, 0.95, 0.4, 0.0, 0.0, 0.25);
        mt.add_material("sjkr", 8.7e9, 0.3, 0.95, 0.4, 0.0, 1.0);
        mt
    }

    fn material_table_with_twisting() -> MaterialTable {
        let mut mt = MaterialTable::new();
        mt.add_material_extended("glass", 8.7e9, 0.3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0);
        mt.build_pair_tables();
        mt
    }

    fn run_wall_twist_case(walls: Walls, pos: [f64; 3], omega: [f64; 3]) -> [f64; 3] {
        let fixture = ParticleFixture::single(ParticleSpec::new(0, pos, 0.001)).build();
        fixture
            .registry
            .expect_mut::<DemAtom>("run_wall_twist_case")
            .omega[0] = omega;
        let (atom, _neighbor, registry, _materials) = fixture.into_parts();
        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(material_table_with_twisting());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
        let torque = registry.expect::<DemAtom>("run_wall_twist_case").torque[0];
        torque
    }

    #[test]
    fn twisting_friction_matches_plane_for_equivalent_curved_wall_normals() {
        let pos = [0.0095, 0.0, 0.0];
        let omega = [-5.0, 0.0, 0.0];

        let plane = run_wall_twist_case(
            make_walls(vec![make_wall_plane(0.01, 0.0, 0.0, -1.0, 0.0, 0.0)]),
            pos,
            omega,
        );
        let cylinder = run_wall_twist_case(
            make_walls_with_cylinder(WallCylinder {
                axis: 2,
                center: [0.0, 0.0],
                radius: 0.01,
                lo: -0.01,
                hi: 0.01,
                inside: true,
                material_index: 0,
                name: None,
                force_accumulator: 0.0,
                temperature: None,
            }),
            pos,
            omega,
        );
        let sphere = run_wall_twist_case(
            make_walls_with_sphere(WallSphere {
                center: [0.0, 0.0, 0.0],
                radius: 0.01,
                inside: true,
                material_index: 0,
                name: None,
                force_accumulator: 0.0,
                temperature: None,
            }),
            pos,
            omega,
        );

        assert!(
            plane[0] > 0.0 && plane[1].abs() < 1e-20 && plane[2].abs() < 1e-20,
            "plane twist torque should oppose spin about -x, got {plane:?}"
        );
        for (name, torque) in [("cylinder", cylinder), ("sphere", sphere)] {
            for axis in 0..3 {
                assert!(
                    (torque[axis] - plane[axis]).abs() < 1e-12,
                    "{name} twisting torque should match plane axis {axis}: plane={plane:?} {name}={torque:?}"
                );
            }
        }
    }

    #[test]
    fn twisting_friction_matches_plane_for_equivalent_region_wall_normal() {
        let pos = [0.0095, 0.0, 0.0];
        let omega = [-5.0, 0.0, 0.0];

        let plane = run_wall_twist_case(
            make_walls(vec![make_wall_plane(0.01, 0.0, 0.0, -1.0, 0.0, 0.0)]),
            pos,
            omega,
        );
        let region = run_wall_twist_case(
            make_walls_with_region(WallRegion {
                region: Region::Sphere {
                    center: [0.0, 0.0, 0.0],
                    radius: 0.01,
                },
                inside: true,
                material_index: 0,
                name: None,
                force_accumulator: 0.0,
                temperature: None,
            }),
            pos,
            omega,
        );

        for axis in 0..3 {
            assert!(
                (region[axis] - plane[axis]).abs() < 1e-12,
                "region twisting torque should match plane axis {axis}: plane={plane:?} region={region:?}"
            );
        }
    }

    #[test]
    fn hooke_wall_normal_uses_kn_and_beta_tables() {
        let mut mt = MaterialTable::new();
        mt.contact_model = "hooke".to_string();
        mt.limit_damping = false;
        mt.add_material_extended(
            "grain", 8.7e9, 0.3, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e5, 0.0,
        );
        mt.build_pair_tables();

        let delta = 2.5e-5;
        let m_r = 0.002;
        let v_n = -0.2;
        let beta = mt.beta_ij[0][0];
        let expected = 1.0e5 * delta - 2.0 * beta * (1.0e5_f64 * m_r).sqrt() * v_n;

        let force = wall_normal_force(&mt, 0, 0, 0.005, delta, v_n, m_r, false);
        assert!(
            (force - expected).abs() / expected.abs() < 1.0e-12,
            "Hooke wall force should use kn_ij and beta_ij: got {force}, expected {expected}"
        );
    }

    #[test]
    fn warns_for_curved_and_region_wall_surface_energy() {
        let mt = material_table_with_surface_energy();

        for geometry in ["cylinder", "sphere", "region"] {
            let mut wall = wall_def(geometry, "sticky");
            wall.name = Some(format!("{geometry}-guard"));
            let msg = curved_or_region_wall_surface_energy_warning(&wall, &mt, 1)
                .expect("curved/region walls with surface_energy must warn");

            assert!(
                msg.contains(&format!("{geometry} wall")),
                "warning must name the wall geometry: {msg}"
            );
            assert!(
                msg.contains(&format!("{geometry}-guard")),
                "warning must name the configured wall when available: {msg}"
            );
            assert!(
                msg.contains("sticky") && msg.contains("surface_energy"),
                "warning must name the material and ignored field: {msg}"
            );
            assert!(
                msg.contains("plane-wall-only") && msg.contains("JKR/DMT"),
                "warning must explain the plane-wall-only JKR/DMT limitation: {msg}"
            );
        }
    }

    #[test]
    fn surface_energy_warning_keeps_valid_wall_configs_silent() {
        let mt = material_table_with_surface_energy();

        assert!(
            curved_or_region_wall_surface_energy_warning(&wall_def("plane", "sticky"), &mt, 1)
                .is_none(),
            "plane walls with surface_energy remain accepted without warning"
        );
        assert!(
            curved_or_region_wall_surface_energy_warning(&wall_def("cylinder", "dry"), &mt, 0)
                .is_none(),
            "curved walls without surface_energy must not warn"
        );
        assert!(
            curved_or_region_wall_surface_energy_warning(&wall_def("region", "sjkr"), &mt, 2)
                .is_none(),
            "region walls using cohesion_energy for SJKR behavior must not warn"
        );
    }

    #[test]
    fn wall_repulsive_for_overlap() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        push_dem_test_atom(&mut atom, &mut dem, 0, [0.01, 0.01, 0.0005], radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls(vec![make_wall_plane(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)]);

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(
            atom.force[0][2] > 0.0,
            "atom should be pushed away from wall, got {}",
            atom.force[0][2]
        );
        assert!((atom.force[0][0]).abs() < 1e-15);
        assert!((atom.force[0][1]).abs() < 1e-15);
    }

    #[test]
    fn wall_zero_for_no_overlap() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        push_dem_test_atom(&mut atom, &mut dem, 0, [0.01, 0.01, 0.002], radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls(vec![make_wall_plane(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)]);

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!((atom.force[0][2]).abs() < 1e-15);
    }

    #[test]
    fn inactive_wall_applies_no_force() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        push_dem_test_atom(&mut atom, &mut dem, 0, [0.01, 0.01, 0.0005], radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let mut plane = make_wall_plane(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        plane.name = Some("blocker".into());
        let mut walls = make_walls(vec![plane]);
        walls.active[0] = false;

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(
            (atom.force[0][2]).abs() < 1e-15,
            "inactive wall should apply no force"
        );
    }

    #[test]
    fn angled_wall_force_direction() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        push_dem_test_atom(&mut atom, &mut dem, 0, [0.0003, 0.0, 0.0003], radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls(vec![make_wall_plane(0.0, 0.0, 0.0, 1.0, 0.0, 1.0)]);

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(atom.force[0][0] > 0.0, "force_x should be positive");
        assert!(atom.force[0][2] > 0.0, "force_z should be positive");
        assert!(
            (atom.force[0][0] - atom.force[0][2]).abs() < 1e-10,
            "force_x and force_z should be equal for 45-degree wall"
        );
        assert!((atom.force[0][1]).abs() < 1e-15);
    }

    #[test]
    fn bounded_wall_ignores_out_of_bounds_atom() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        push_dem_test_atom(&mut atom, &mut dem, 0, [0.05, 0.01, 0.0005], radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let mut wall = make_wall_plane(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        wall.bound_x_low = 0.0;
        wall.bound_x_high = 0.04;

        let walls = make_walls(vec![wall]);

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(
            (atom.force[0][2]).abs() < 1e-15,
            "out-of-bounds atom should get no wall force"
        );
    }

    #[test]
    fn wall_cohesion_attractive_for_small_overlap() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        push_dem_test_atom(&mut atom, &mut dem, 0, [0.01, 0.01, 0.000999], radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls(vec![make_wall_plane(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)]);

        let mut mt = dirt_atom::MaterialTable::new();
        mt.add_material("glass", 8.7e9, 0.3, 0.95, 0.4, 0.0, 1e9);
        mt.build_pair_tables();

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(mt);
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(
            atom.force[0][2] < 0.0,
            "wall cohesion should produce attractive force, got {}",
            atom.force[0][2]
        );
    }

    // ── Moving wall tests ───────────────────────────────────────────────────

    #[test]
    fn constant_velocity_wall_moves() {
        let mut atom = Atom::new();
        atom.dt = 0.001;
        atom.nlocal = 0;
        atom.natoms = 0;

        let mut plane = make_wall_plane(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        plane.velocity = [0.0, 0.0, -0.01];
        plane.motion = WallMotion::ConstantVelocity;

        let walls = make_walls(vec![plane]);

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(walls);
        app.add_update_system(wall_move, ParticleSimScheduleSet::PreInitialIntegration);
        app.organize_systems();
        app.run();

        let walls = app.get_resource_ref::<Walls>().unwrap();
        assert!(
            (walls.planes[0].point_z - (-0.00001)).abs() < 1e-15,
            "wall should move, got {}",
            walls.planes[0].point_z
        );
        assert!((walls.time - 0.001).abs() < 1e-15);
    }

    #[test]
    fn oscillating_wall_follows_sine() {
        let mut atom = Atom::new();
        atom.dt = 0.001;
        atom.nlocal = 0;
        atom.natoms = 0;

        let amplitude = 0.002;
        let frequency = 50.0;
        let mut plane = make_wall_plane(0.0, 0.0, 0.1, 0.0, 0.0, 1.0);
        plane.motion = WallMotion::Oscillate {
            amplitude,
            frequency,
        };

        let walls = make_walls(vec![plane]);

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(walls);
        app.add_update_system(wall_move, ParticleSimScheduleSet::PreInitialIntegration);
        app.organize_systems();
        app.run();

        let walls = app.get_resource_ref::<Walls>().unwrap();
        let expected_phase = 2.0 * std::f64::consts::PI * frequency * 0.001;
        let expected_disp = amplitude * expected_phase.sin();
        assert!(
            (walls.planes[0].point_z - (0.1 + expected_disp)).abs() < 1e-12,
            "oscillating wall z={}, expected {}",
            walls.planes[0].point_z,
            0.1 + expected_disp
        );
    }

    #[test]
    fn servo_wall_adjusts_velocity() {
        let mut atom = Atom::new();
        atom.dt = 0.001;
        atom.nlocal = 0;
        atom.natoms = 0;

        let mut plane = make_wall_plane(0.0, 0.0, 0.1, 0.0, 0.0, -1.0);
        plane.motion = WallMotion::Servo {
            target_force: 100.0,
            max_velocity: 0.1,
            gain: 0.001,
        };
        // Simulate accumulated force = 50 (below target)
        plane.force_accumulator = 50.0;

        let walls = make_walls(vec![plane]);

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(walls);
        app.add_update_system(wall_move, ParticleSimScheduleSet::PreInitialIntegration);
        app.organize_systems();
        app.run();

        let walls = app.get_resource_ref::<Walls>().unwrap();
        // error = 100 - 50 = 50, vel_mag = 0.001 * 50 = 0.05 (within max)
        // velocity along normal (-z): vel = 0.05 * [0, 0, -1] = [0, 0, -0.05]
        assert!(
            (walls.planes[0].velocity[2] - (-0.05)).abs() < 1e-10,
            "servo velocity z={}, expected -0.05",
            walls.planes[0].velocity[2]
        );
        // Position should move
        assert!(
            walls.planes[0].point_z < 0.1,
            "servo wall should have moved"
        );
    }

    #[test]
    fn moving_wall_relative_velocity_affects_damping() {
        // A wall moving toward a stationary atom should produce higher force
        // than a static wall with the same overlap
        let radius = 0.001;

        let run_with_wall_vel = |wall_vel: [f64; 3]| -> f64 {
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            push_dem_test_atom(&mut atom, &mut dem, 0, [0.01, 0.01, 0.0005], radius);
            atom.nlocal = 1;
            atom.natoms = 1;

            let mut registry = AtomDataRegistry::new();
            registry.try_register(dem, atom.len()).unwrap();

            let mut plane = make_wall_plane(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
            plane.velocity = wall_vel;
            plane.motion = WallMotion::ConstantVelocity;
            let walls = make_walls(vec![plane]);

            let mut app = App::new();
            app.add_resource(atom);
            app.add_resource(registry);
            app.add_resource(make_material_table());
            app.add_resource(walls);
            app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();

            let atom = app.get_resource_ref::<Atom>().unwrap();
            atom.force[0][2]
        };

        let f_static = run_with_wall_vel([0.0, 0.0, 0.0]);
        let f_approaching = run_with_wall_vel([0.0, 0.0, 1.0]); // wall moving toward atom

        // Wall approaching means relative velocity is negative (approaching)
        // which increases damping force, so total repulsion should be higher
        assert!(
            f_approaching > f_static,
            "approaching wall should increase repulsive force: static={}, approaching={}",
            f_static,
            f_approaching
        );
    }

    #[test]
    fn wall_history_initializes_without_particle_contact_plugin() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        push_dem_test_atom(&mut atom, &mut dem, 0, [0.0, 0.0, 0.0005], radius);
        atom.vel[0][0] = 0.25;
        dem.omega[0][1] = 40.0;
        atom.dt = 1.0e-5;
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let mut material_table = dirt_atom::MaterialTable::new();
        material_table.rolling_model = "sds".to_string();
        material_table.add_material_with_sds(
            "glass", 8.7e9, 0.3, 0.95, 0.5, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-7, 0.0, 0.0, 0.0,
        );
        material_table.build_pair_tables();

        let walls = make_walls(vec![make_wall_plane(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)]);

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(material_table);
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(
            atom.force[0][0] < 0.0,
            "wall tangential friction should oppose positive x slip, got fx={}",
            atom.force[0][0]
        );

        let walls = app.get_resource_ref::<Walls>().unwrap();
        let key = (0u8, 0usize, 0u32);
        let tangential = walls
            .tangential_springs
            .get(&key)
            .expect("dirt_wall must initialize plane tangential history itself");
        let rolling = walls
            .rolling_springs
            .get(&key)
            .expect("dirt_wall must initialize plane SDS rolling history itself");

        let tangential_mag = (tangential[0] * tangential[0]
            + tangential[1] * tangential[1]
            + tangential[2] * tangential[2])
            .sqrt();
        let rolling_mag =
            (rolling[0] * rolling[0] + rolling[1] * rolling[1] + rolling[2] * rolling[2]).sqrt();
        assert!(
            tangential_mag > 0.0,
            "stored wall tangential spring should advance from zero"
        );
        assert!(
            rolling_mag > 0.0,
            "stored wall rolling spring should advance from zero"
        );
    }

    // ── Cylinder wall tests ────────────────────────────────────────────────

    fn make_walls_with_cylinder(cyl: WallCylinder) -> Walls {
        Walls {
            planes: Vec::new(),
            active: Vec::new(),
            cylinders: vec![cyl],
            cylinder_active: vec![true],
            spheres: Vec::new(),
            sphere_active: Vec::new(),
            regions: Vec::new(),
            region_active: Vec::new(),
            time: 0.0,
            tangential_springs: std::collections::HashMap::new(),
            rolling_springs: std::collections::HashMap::new(),
        }
    }

    fn make_walls_with_sphere(sph: WallSphere) -> Walls {
        Walls {
            planes: Vec::new(),
            active: Vec::new(),
            cylinders: Vec::new(),
            cylinder_active: Vec::new(),
            spheres: vec![sph],
            sphere_active: vec![true],
            regions: Vec::new(),
            region_active: Vec::new(),
            time: 0.0,
            tangential_springs: std::collections::HashMap::new(),
            rolling_springs: std::collections::HashMap::new(),
        }
    }

    fn make_walls_with_region(reg: WallRegion) -> Walls {
        Walls {
            planes: Vec::new(),
            active: Vec::new(),
            cylinders: Vec::new(),
            cylinder_active: Vec::new(),
            spheres: Vec::new(),
            sphere_active: Vec::new(),
            regions: vec![reg],
            region_active: vec![true],
            time: 0.0,
            tangential_springs: std::collections::HashMap::new(),
            rolling_springs: std::collections::HashMap::new(),
        }
    }

    fn make_named_wall_set(name: &str) -> Walls {
        let mut plane = make_wall_plane(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        plane.name = Some(name.to_string());

        Walls {
            planes: vec![plane],
            active: vec![true],
            cylinders: vec![WallCylinder {
                axis: 2,
                center: [0.005, 0.005],
                radius: 0.004,
                lo: 0.0,
                hi: 0.01,
                inside: true,
                material_index: 0,
                name: Some(name.to_string()),
                force_accumulator: 0.0,
                temperature: None,
            }],
            cylinder_active: vec![true],
            spheres: vec![WallSphere {
                center: [0.005, 0.005, 0.005],
                radius: 0.004,
                inside: true,
                material_index: 0,
                name: Some(name.to_string()),
                force_accumulator: 0.0,
                temperature: None,
            }],
            sphere_active: vec![true],
            regions: vec![WallRegion {
                region: Region::Sphere {
                    center: [0.005, 0.005, 0.005],
                    radius: 0.004,
                },
                inside: true,
                material_index: 0,
                name: Some(name.to_string()),
                force_accumulator: 0.0,
                temperature: None,
            }],
            region_active: vec![true],
            time: 0.0,
            tangential_springs: std::collections::HashMap::new(),
            rolling_springs: std::collections::HashMap::new(),
        }
    }

    #[test]
    fn activate_by_name_reactivates_all_wall_geometries() {
        let mut walls = make_named_wall_set("gate");

        walls.deactivate_by_name("gate");
        assert_eq!(walls.active, vec![false]);
        assert_eq!(walls.cylinder_active, vec![false]);
        assert_eq!(walls.sphere_active, vec![false]);
        assert_eq!(walls.region_active, vec![false]);

        walls.activate_by_name("gate");
        assert_eq!(walls.active, vec![true]);
        assert_eq!(walls.cylinder_active, vec![true]);
        assert_eq!(walls.sphere_active, vec![true]);
        assert_eq!(walls.region_active, vec![true]);
    }

    #[test]
    fn activate_by_name_unknown_name_is_noop() {
        let mut walls = make_named_wall_set("gate");

        walls.deactivate_by_name("gate");
        walls.activate_by_name("missing");

        assert_eq!(walls.active, vec![false]);
        assert_eq!(walls.cylinder_active, vec![false]);
        assert_eq!(walls.sphere_active, vec![false]);
        assert_eq!(walls.region_active, vec![false]);
    }

    #[test]
    fn cylinder_inside_repels_toward_center() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        // Particle near the wall of a z-cylinder centered at (0.005, 0.005), radius 0.004
        // Place particle at radial distance 0.0035 from axis (gap = 0.004 - 0.0035 = 0.0005 < radius)
        push_dem_test_atom(
            &mut atom,
            &mut dem,
            0,
            [0.005 + 0.0035, 0.005, 0.005],
            radius,
        );
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls_with_cylinder(WallCylinder {
            axis: 2, // Z
            center: [0.005, 0.005],
            radius: 0.004,
            lo: 0.0,
            hi: 0.01,
            inside: true,
            material_index: 0,
            name: None,
            force_accumulator: 0.0,
            temperature: None,
        });

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        // Force should push particle toward center (negative x direction)
        assert!(
            atom.force[0][0] < 0.0,
            "cylinder should push particle toward center, got fx={}",
            atom.force[0][0]
        );
        assert!((atom.force[0][1]).abs() < 1e-15, "no y force");
        assert!((atom.force[0][2]).abs() < 1e-15, "no z force");
    }

    #[test]
    fn cylinder_no_force_when_not_touching() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        // Particle well inside cylinder (far from wall)
        push_dem_test_atom(&mut atom, &mut dem, 0, [0.005, 0.005, 0.005], radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls_with_cylinder(WallCylinder {
            axis: 2,
            center: [0.005, 0.005],
            radius: 0.004,
            lo: 0.0,
            hi: 0.01,
            inside: true,
            material_index: 0,
            name: None,
            force_accumulator: 0.0,
            temperature: None,
        });

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        let f_mag =
            (atom.force[0][0].powi(2) + atom.force[0][1].powi(2) + atom.force[0][2].powi(2)).sqrt();
        assert!(
            f_mag < 1e-15,
            "no force when not touching cylinder wall, got {}",
            f_mag
        );
    }

    #[test]
    fn sphere_inside_repels_toward_center() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        // Particle near the wall of a sphere centered at (0.005, 0.005, 0.005), radius 0.004
        push_dem_test_atom(
            &mut atom,
            &mut dem,
            0,
            [0.005 + 0.0035, 0.005, 0.005],
            radius,
        );
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls_with_sphere(WallSphere {
            center: [0.005, 0.005, 0.005],
            radius: 0.004,
            inside: true,
            material_index: 0,
            name: None,
            force_accumulator: 0.0,
            temperature: None,
        });

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        // Force should push particle toward center (negative x direction)
        assert!(
            atom.force[0][0] < 0.0,
            "sphere should push particle toward center, got fx={}",
            atom.force[0][0]
        );
        assert!((atom.force[0][1]).abs() < 1e-15, "no y force");
        assert!((atom.force[0][2]).abs() < 1e-15, "no z force");
    }

    #[test]
    fn sphere_no_force_when_not_touching() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        // Particle at center of sphere
        push_dem_test_atom(&mut atom, &mut dem, 0, [0.005, 0.005, 0.005], radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls_with_sphere(WallSphere {
            center: [0.005, 0.005, 0.005],
            radius: 0.004,
            inside: true,
            material_index: 0,
            name: None,
            force_accumulator: 0.0,
            temperature: None,
        });

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        let f_mag =
            (atom.force[0][0].powi(2) + atom.force[0][1].powi(2) + atom.force[0][2].powi(2)).sqrt();
        assert!(
            f_mag < 1e-15,
            "no force when not touching sphere wall, got {}",
            f_mag
        );
    }

    // ══════════════════════════════════════════════════════════════════════
    // VALIDATION: Cylinder wall force direction always points radially
    // For an inside cylinder, the force should always point toward the axis
    // regardless of where the particle is around the circumference.
    // Test at multiple angular positions around a Z-cylinder.
    // ══════════════════════════════════════════════════════════════════════
    #[test]
    fn cylinder_force_always_points_radially_inward() {
        let particle_radius = 0.001;
        let cyl_radius = 0.01;
        let center = [0.005, 0.005];
        // Place particles near the wall at different angles
        let angles: Vec<f64> = vec![0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0];

        for angle in &angles {
            let r = cyl_radius - 0.5 * particle_radius; // overlap = 0.5 * particle_radius
            let px = center[0] + r * angle.cos();
            let py = center[1] + r * angle.sin();
            let pz = 0.005;

            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            push_dem_test_atom(&mut atom, &mut dem, 0, [px, py, pz], particle_radius);
            atom.nlocal = 1;
            atom.natoms = 1;

            let mut registry = AtomDataRegistry::new();
            registry.try_register(dem, atom.len()).unwrap();

            let walls = make_walls_with_cylinder(WallCylinder {
                axis: 2,
                center,
                radius: cyl_radius,
                lo: 0.0,
                hi: 0.01,
                inside: true,
                material_index: 0,
                name: None,
                force_accumulator: 0.0,
                temperature: None,
            });

            let mut app = App::new();
            app.add_resource(atom);
            app.add_resource(registry);
            app.add_resource(make_material_table());
            app.add_resource(walls);
            app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();

            let atom = app.get_resource_ref::<Atom>().unwrap();
            let fx = atom.force[0][0];
            let fy = atom.force[0][1];
            let fz = atom.force[0][2];

            // Force should be purely radial (no z component)
            assert!(
                fz.abs() < 1e-12,
                "angle={:.1}: no z force expected, got {:.6e}",
                angle,
                fz
            );

            // Force direction should point toward axis center
            let dx = px - center[0];
            let dy = py - center[1];
            let r_actual = (dx * dx + dy * dy).sqrt();
            // Radial unit vector (outward): (dx/r, dy/r)
            // Force should oppose this (inward): dot(f, r_hat) < 0
            let f_dot_r = fx * dx / r_actual + fy * dy / r_actual;
            assert!(
                f_dot_r < 0.0,
                "angle={:.1}: force should point inward, f·r_hat={:.6e}",
                angle,
                f_dot_r
            );

            // Force magnitude should be nonzero
            let f_mag = (fx * fx + fy * fy).sqrt();
            assert!(
                f_mag > 0.0,
                "angle={:.1}: force magnitude should be nonzero",
                angle
            );
        }
    }

    // ══════════════════════════════════════════════════════════════════════
    // VALIDATION: Cylinder axial bounds are enforced
    // Particles outside lo/hi should not get any force from the cylinder.
    // ══════════════════════════════════════════════════════════════════════
    #[test]
    fn cylinder_axial_bounds_enforced() {
        let particle_radius = 0.001;
        let cyl_radius = 0.01;
        let center = [0.005, 0.005];

        // Place particle near the wall but outside axial bounds (below lo)
        let r = cyl_radius - 0.5 * particle_radius;
        let px = center[0] + r;
        let py = center[1];
        let pz = -0.001; // below lo=0.0

        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        push_dem_test_atom(&mut atom, &mut dem, 0, [px, py, pz], particle_radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls_with_cylinder(WallCylinder {
            axis: 2,
            center,
            radius: cyl_radius,
            lo: 0.0,
            hi: 0.01,
            inside: true,
            material_index: 0,
            name: None,
            force_accumulator: 0.0,
            temperature: None,
        });

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        let f_mag =
            (atom.force[0][0].powi(2) + atom.force[0][1].powi(2) + atom.force[0][2].powi(2)).sqrt();
        assert!(
            f_mag < 1e-15,
            "No force outside axial bounds: f_mag={:.6e}",
            f_mag
        );
    }

    // ══════════════════════════════════════════════════════════════════════
    // VALIDATION: Sphere wall force direction at multiple positions
    // For an inside sphere, force should always point toward the center.
    // ══════════════════════════════════════════════════════════════════════
    #[test]
    fn sphere_force_points_toward_center_at_multiple_positions() {
        let particle_radius = 0.001;
        let sph_radius = 0.01;
        let sph_center = [0.005, 0.005, 0.005];

        // Test positions along different axes
        let directions: Vec<[f64; 3]> = vec![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];

        for dir in &directions {
            let mag = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
            let nd = [dir[0] / mag, dir[1] / mag, dir[2] / mag];
            let r = sph_radius - 0.5 * particle_radius;
            let px = sph_center[0] + r * nd[0];
            let py = sph_center[1] + r * nd[1];
            let pz = sph_center[2] + r * nd[2];

            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            push_dem_test_atom(&mut atom, &mut dem, 0, [px, py, pz], particle_radius);
            atom.nlocal = 1;
            atom.natoms = 1;

            let mut registry = AtomDataRegistry::new();
            registry.try_register(dem, atom.len()).unwrap();

            let walls = make_walls_with_sphere(WallSphere {
                center: sph_center,
                radius: sph_radius,
                inside: true,
                material_index: 0,
                name: None,
                force_accumulator: 0.0,
                temperature: None,
            });

            let mut app = App::new();
            app.add_resource(atom);
            app.add_resource(registry);
            app.add_resource(make_material_table());
            app.add_resource(walls);
            app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();

            let atom = app.get_resource_ref::<Atom>().unwrap();
            let fx = atom.force[0][0];
            let fy = atom.force[0][1];
            let fz = atom.force[0][2];
            let f_mag = (fx * fx + fy * fy + fz * fz).sqrt();

            assert!(f_mag > 0.0, "dir={:?}: force should be nonzero", dir);

            // Force should point toward center: dot(f, r_hat) < 0
            let dx = px - sph_center[0];
            let dy = py - sph_center[1];
            let dz = pz - sph_center[2];
            let r_actual = (dx * dx + dy * dy + dz * dz).sqrt();
            let f_dot_r = fx * dx / r_actual + fy * dy / r_actual + fz * dz / r_actual;
            assert!(
                f_dot_r < 0.0,
                "dir={:?}: force should point toward center, f·r_hat={:.6e}",
                dir,
                f_dot_r
            );
        }
    }

    // ══════════════════════════════════════════════════════════════════════
    // VALIDATION: Wall force at exact contact (delta=0) should be zero
    // When a particle's surface just touches the wall (no overlap),
    // the elastic force should be zero.
    // ══════════════════════════════════════════════════════════════════════
    #[test]
    fn wall_zero_force_at_exact_contact() {
        let radius = 0.001;

        // Place particle center at exactly radius from wall -> delta = 0
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        push_dem_test_atom(&mut atom, &mut dem, 0, [0.01, 0.01, radius], radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls(vec![make_wall_plane(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)]);

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        let f_mag =
            (atom.force[0][0].powi(2) + atom.force[0][1].powi(2) + atom.force[0][2].powi(2)).sqrt();
        assert!(
            f_mag < 1e-10,
            "Force at exact contact should be ~zero, got {:.6e}",
            f_mag
        );
    }

    // ══════════════════════════════════════════════════════════════════════
    // VALIDATION: Wall Hertz force scales as delta^(3/2) for plane walls
    // Same as particle-particle Hertz, but with R_eff = R_particle.
    // ══════════════════════════════════════════════════════════════════════
    #[test]
    fn wall_hertz_force_scales_as_delta_three_halves() {
        let radius = 0.001;

        let wall_force_at = |delta: f64| -> f64 {
            let distance = radius - delta; // signed distance from wall to center
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            push_dem_test_atom(&mut atom, &mut dem, 0, [0.01, 0.01, distance], radius);
            atom.nlocal = 1;
            atom.natoms = 1;

            let mut registry = AtomDataRegistry::new();
            registry.try_register(dem, atom.len()).unwrap();

            let walls = make_walls(vec![make_wall_plane(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)]);

            let mut app = App::new();
            app.add_resource(atom);
            app.add_resource(registry);
            app.add_resource(make_material_table());
            app.add_resource(walls);
            app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();

            let atom = app.get_resource_ref::<Atom>().unwrap();
            atom.force[0][2].abs()
        };

        let deltas = [1e-5, 2e-5, 4e-5, 6e-5, 8e-5];
        let forces: Vec<f64> = deltas.iter().map(|d| wall_force_at(*d)).collect();

        for i in 1..deltas.len() {
            let expected_ratio = (deltas[i] / deltas[0]).powf(1.5);
            let actual_ratio = forces[i] / forces[0];
            let rel_err = ((actual_ratio - expected_ratio) / expected_ratio).abs();
            assert!(
                rel_err < 0.01,
                "Wall Hertz scaling: delta ratio {:.1}, expected F ratio {:.4}, got {:.4} (rel err {:.4})",
                deltas[i] / deltas[0], expected_ratio, actual_ratio, rel_err
            );
        }
    }

    // ── Region wall tests ─────────────────────────────────────────────────

    #[test]
    fn region_sphere_inside_repels() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        // Particle near sphere wall surface (inside sphere of radius 0.004)
        push_dem_test_atom(
            &mut atom,
            &mut dem,
            0,
            [0.005 + 0.0035, 0.005, 0.005],
            radius,
        );
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls_with_region(WallRegion {
            region: Region::Sphere {
                center: [0.005, 0.005, 0.005],
                radius: 0.004,
            },
            inside: true,
            material_index: 0,
            name: None,
            force_accumulator: 0.0,
            temperature: None,
        });

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        assert!(
            atom.force[0][0] < 0.0,
            "region sphere wall should push particle toward center, got fx={}",
            atom.force[0][0]
        );
        assert!((atom.force[0][1]).abs() < 1e-15, "no y force");
        assert!((atom.force[0][2]).abs() < 1e-15, "no z force");
    }

    #[test]
    fn region_sphere_no_force_when_far() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        // Particle at center of sphere (far from wall)
        push_dem_test_atom(&mut atom, &mut dem, 0, [0.005, 0.005, 0.005], radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls_with_region(WallRegion {
            region: Region::Sphere {
                center: [0.005, 0.005, 0.005],
                radius: 0.004,
            },
            inside: true,
            material_index: 0,
            name: None,
            force_accumulator: 0.0,
            temperature: None,
        });

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        let f_mag =
            (atom.force[0][0].powi(2) + atom.force[0][1].powi(2) + atom.force[0][2].powi(2)).sqrt();
        assert!(
            f_mag < 1e-15,
            "no force when far from region wall, got {}",
            f_mag
        );
    }

    #[test]
    fn region_block_inside_repels() {
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.001;

        // Particle near the +z face of a block (inside, close to top)
        push_dem_test_atom(&mut atom, &mut dem, 0, [0.005, 0.005, 0.0095], radius);
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls_with_region(WallRegion {
            region: Region::Block {
                min: [0.0, 0.0, 0.0],
                max: [0.01, 0.01, 0.01],
            },
            inside: true,
            material_index: 0,
            name: None,
            force_accumulator: 0.0,
            temperature: None,
        });

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        // Should be pushed away from the +z face (downward)
        assert!(
            atom.force[0][2] < 0.0,
            "region block wall should push particle away from +z face, got fz={}",
            atom.force[0][2]
        );
    }

    #[test]
    fn region_cone_inside_repels() {
        use soil_core::region::Axis;
        let mut atom = Atom::new();
        let mut dem = DemAtom::new();
        let radius = 0.0005;

        // Cone: z-axis, rad_lo=0.004 at z=0, rad_hi=0.002 at z=0.01
        // At z=0.005, radius = 0.003
        // Place particle at radial distance 0.0028 from axis (gap = 0.003 - 0.0028 = 0.0002 < radius)
        push_dem_test_atom(
            &mut atom,
            &mut dem,
            0,
            [0.005 + 0.0028, 0.005, 0.005],
            radius,
        );
        atom.nlocal = 1;
        atom.natoms = 1;

        let mut registry = AtomDataRegistry::new();
        registry.try_register(dem, atom.len()).unwrap();

        let walls = make_walls_with_region(WallRegion {
            region: Region::Cone {
                center: [0.005, 0.005],
                axis: Axis::Z,
                rad_lo: 0.004,
                rad_hi: 0.002,
                lo: 0.0,
                hi: 0.01,
            },
            inside: true,
            material_index: 0,
            name: None,
            force_accumulator: 0.0,
            temperature: None,
        });

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(registry);
        app.add_resource(make_material_table());
        app.add_resource(walls);
        app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
        app.organize_systems();
        app.run();

        let atom = app.get_resource_ref::<Atom>().unwrap();
        // Force should push toward center (negative x direction)
        assert!(
            atom.force[0][0] < 0.0,
            "cone wall should push particle toward center, got fx={}",
            atom.force[0][0]
        );
    }

    #[test]
    fn region_wall_force_matches_dedicated_sphere() {
        // Verify that a region sphere wall produces the same force as the dedicated sphere wall
        let radius = 0.001;
        let sphere_center = [0.005, 0.005, 0.005];
        let sphere_radius = 0.004;
        let particle_pos = [0.005 + 0.0035, 0.005, 0.005];

        // Run with dedicated sphere wall
        let f_dedicated = {
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            push_dem_test_atom(&mut atom, &mut dem, 0, particle_pos, radius);
            atom.nlocal = 1;
            atom.natoms = 1;
            let mut registry = AtomDataRegistry::new();
            registry.try_register(dem, atom.len()).unwrap();
            let walls = make_walls_with_sphere(WallSphere {
                center: sphere_center,
                radius: sphere_radius,
                inside: true,
                material_index: 0,
                name: None,
                force_accumulator: 0.0,
                temperature: None,
            });
            let mut app = App::new();
            app.add_resource(atom);
            app.add_resource(registry);
            app.add_resource(make_material_table());
            app.add_resource(walls);
            app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();
            let atom = app.get_resource_ref::<Atom>().unwrap();
            [atom.force[0][0], atom.force[0][1], atom.force[0][2]]
        };

        // Run with region sphere wall
        let f_region = {
            let mut atom = Atom::new();
            let mut dem = DemAtom::new();
            push_dem_test_atom(&mut atom, &mut dem, 0, particle_pos, radius);
            atom.nlocal = 1;
            atom.natoms = 1;
            let mut registry = AtomDataRegistry::new();
            registry.try_register(dem, atom.len()).unwrap();
            let walls = make_walls_with_region(WallRegion {
                region: Region::Sphere {
                    center: sphere_center,
                    radius: sphere_radius,
                },
                inside: true,
                material_index: 0,
                name: None,
                force_accumulator: 0.0,
                temperature: None,
            });
            let mut app = App::new();
            app.add_resource(atom);
            app.add_resource(registry);
            app.add_resource(make_material_table());
            app.add_resource(walls);
            app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
            app.organize_systems();
            app.run();
            let atom = app.get_resource_ref::<Atom>().unwrap();
            [atom.force[0][0], atom.force[0][1], atom.force[0][2]]
        };

        for d in 0..3 {
            assert!(
                (f_dedicated[d] - f_region[d]).abs() < 1e-6 * f_dedicated[d].abs().max(1e-15),
                "force mismatch in dim {}: dedicated={}, region={}",
                d,
                f_dedicated[d],
                f_region[d]
            );
        }
    }

    #[test]
    fn static_wall_unaffected_by_motion_systems() {
        let mut atom = Atom::new();
        atom.dt = 0.001;
        atom.nlocal = 0;
        atom.natoms = 0;

        let plane = make_wall_plane(0.0, 0.0, 0.5, 0.0, 0.0, 1.0);
        let walls = make_walls(vec![plane]);

        let mut app = App::new();
        app.add_resource(atom);
        app.add_resource(walls);
        app.add_update_system(wall_move, ParticleSimScheduleSet::PreInitialIntegration);
        app.organize_systems();
        app.run();

        let walls = app.get_resource_ref::<Walls>().unwrap();
        assert!(
            (walls.planes[0].point_z - 0.5).abs() < 1e-15,
            "static wall should not move"
        );
    }
}
