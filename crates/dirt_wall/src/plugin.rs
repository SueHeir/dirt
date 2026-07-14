//! Wall plugin wiring and configuration validation.

use crate::config::{curved_or_region_wall_surface_energy_warning, WallDef};
use crate::contact::wall_contact_force;
use crate::geometry::{WallCylinder, WallPlane, WallRegion, WallSphere, Walls};
use crate::motion::{wall_move, wall_zero_force_accumulators, WallMotion};
use dirt_atom::MaterialTable;
use dirt_schedule::WALL_CONTACT;
use grass_app::prelude::*;
use grass_scheduler::prelude::*;
use soil_core::{Config, ParticleSimScheduleSet};

// ── Plugin ──────────────────────────────────────────────────────────────────

/// Plugin that registers wall contact force systems from `[[wall]]` TOML config.
///
/// Parses all `[[wall]]` entries, resolves material indices, and creates the
/// [`Walls`] resource. Registers three systems:
///
/// - [`wall_move`] — updates wall positions (oscillation, servo, constant velocity)
/// - [`wall_zero_force_accumulators`] — zeros per-wall force accumulators before force pass
/// - [`wall_contact_force`] — computes Hertz + damping + adhesion forces for all wall types
///
/// # Dependencies
///
/// Requires `DemAtomPlugin` (provides [`MaterialTable`] and [`DemAtom`]).
pub struct WallPlugin;

impl Plugin for WallPlugin {
    fn dependencies(&self) -> Vec<std::any::TypeId> {
        grass_app::type_ids![dirt_atom::DemAtomPlugin]
    }

    fn default_config(&self) -> Option<&str> {
        Some(
            r#"# Wall definitions (uncomment to add walls)
# [[wall]]
# point_x = 0.0
# point_y = 0.0
# point_z = 0.0
# normal_x = 0.0
# normal_y = 0.0
# normal_z = 1.0
# material = "glass"        # must match a [[dem.materials]] name
# name = "floor"            # optional name for runtime enable/disable
# velocity = [0.0, 0.0, -0.01]  # constant velocity (optional)
# oscillate = { amplitude = 0.001, frequency = 50.0 }  # sinusoidal (optional)
# servo = { target_force = 100.0, max_velocity = 0.1, gain = 0.001 }  # servo (optional)"#,
        )
    }

    fn build(&self, app: &mut App) {
        validate_wall_plugin_config(app)
            .unwrap_or_else(|error| panic!("WallPlugin failed to build: {error}"));
        let walls = {
            let config = app
                .get_resource_ref::<Config>()
                .expect("Config resource must exist");
            let wall_defs: Vec<WallDef> = if let Some(val) = config.table.get("wall") {
                match val {
                    toml::Value::Array(arr) => arr
                        .iter()
                        .enumerate()
                        .map(|(idx, v)| match v.clone().try_into::<WallDef>() {
                            Ok(w) => w,
                            Err(e) => {
                                eprintln!("ERROR: failed to parse [[wall]] entry {}: {}", idx, e);
                                panic!(
                                    "WallPlugin preflight should reject malformed [[wall]] entries"
                                );
                            }
                        })
                        .collect(),
                    toml::Value::Table(t) => {
                        match toml::Value::Table(t.clone()).try_into::<WallDef>() {
                            Ok(w) => vec![w],
                            Err(e) => {
                                eprintln!("ERROR: failed to parse [wall] entry: {}", e);
                                panic!(
                                    "WallPlugin preflight should reject malformed [wall] entries"
                                );
                            }
                        }
                    }
                    _ => {
                        eprintln!("ERROR: [wall] must be a table or array of tables");
                        panic!("WallPlugin preflight should reject a non-table wall value");
                    }
                }
            } else {
                Vec::new()
            };
            drop(config);

            let material_table = app
                .get_resource_ref::<MaterialTable>()
                .expect("MaterialTable must exist before WallPlugin — add DemAtomPlugin first");

            let mut planes = Vec::new();
            let mut cylinders = Vec::new();
            let mut spheres = Vec::new();
            let mut regions: Vec<WallRegion> = Vec::new();

            for w in &wall_defs {
                let mat_idx = match material_table.find_material(&w.material) {
                    Some(idx) => idx as usize,
                    None => {
                        eprintln!(
                            "ERROR: wall material '{}' not found in [[dem.materials]]. Available: {:?}",
                            w.material, material_table.names
                        );
                        panic!("WallPlugin preflight should reject unknown wall materials");
                    }
                };

                if let Some(msg) =
                    curved_or_region_wall_surface_energy_warning(w, &material_table, mat_idx)
                {
                    eprintln!("{msg}");
                }

                match w.wall_type.as_str() {
                    "cylinder" => {
                        let axis_str = w.axis.as_deref().unwrap_or("z");
                        let axis = match axis_str {
                            "x" | "X" => 0,
                            "y" | "Y" => 1,
                            "z" | "Z" => 2,
                            _ => {
                                eprintln!(
                                    "ERROR: cylinder wall axis must be x, y, or z, got '{}'",
                                    axis_str
                                );
                                panic!("WallPlugin preflight should reject unknown cylinder axes");
                            }
                        };
                        let center_vec = w.center.as_ref().unwrap_or_else(|| {
                            eprintln!("ERROR: cylinder wall requires 'center' [c0, c1]");
                            panic!("WallPlugin preflight should require cylinder center");
                        });
                        if center_vec.len() != 2 {
                            eprintln!("ERROR: cylinder wall 'center' must have 2 elements");
                            panic!("WallPlugin preflight should validate cylinder center length");
                        }
                        let center = [center_vec[0], center_vec[1]];
                        let radius = w.radius.unwrap_or_else(|| {
                            eprintln!("ERROR: cylinder wall requires 'radius'");
                            panic!("WallPlugin preflight should require cylinder radius");
                        });
                        let lo = w.lo.unwrap_or(f64::NEG_INFINITY);
                        let hi = w.hi.unwrap_or(f64::INFINITY);
                        let inside = w.inside.unwrap_or(false);
                        cylinders.push(WallCylinder {
                            axis,
                            center,
                            radius,
                            lo,
                            hi,
                            inside,
                            material_index: mat_idx,
                            name: w.name.clone(),
                            force_accumulator: 0.0,
                            temperature: w.temperature,
                        });
                    }
                    "sphere" => {
                        let center_vec = w.center.as_ref().unwrap_or_else(|| {
                            eprintln!("ERROR: sphere wall requires 'center' [x, y, z]");
                            panic!("WallPlugin preflight should require sphere center");
                        });
                        if center_vec.len() != 3 {
                            eprintln!("ERROR: sphere wall 'center' must have 3 elements");
                            panic!("WallPlugin preflight should validate sphere center length");
                        }
                        let center = [center_vec[0], center_vec[1], center_vec[2]];
                        let radius = w.radius.unwrap_or_else(|| {
                            eprintln!("ERROR: sphere wall requires 'radius'");
                            panic!("WallPlugin preflight should require sphere radius");
                        });
                        let inside = w.inside.unwrap_or(false);
                        spheres.push(WallSphere {
                            center,
                            radius,
                            inside,
                            material_index: mat_idx,
                            name: w.name.clone(),
                            force_accumulator: 0.0,
                            temperature: w.temperature,
                        });
                    }
                    "region" => {
                        let region = w.region.clone().unwrap_or_else(|| {
                            eprintln!("ERROR: region wall requires a 'region' field");
                            panic!("WallPlugin preflight should require region wall region");
                        });
                        let inside = w.inside.unwrap_or(false);
                        regions.push(WallRegion {
                            region,
                            inside,
                            material_index: mat_idx,
                            name: w.name.clone(),
                            force_accumulator: 0.0,
                            temperature: w.temperature,
                        });
                    }
                    "plane" => {
                        let mag = (w.normal_x * w.normal_x
                            + w.normal_y * w.normal_y
                            + w.normal_z * w.normal_z)
                            .sqrt();
                        if mag <= 1e-15 {
                            eprintln!(
                                "ERROR: wall normal vector must be non-zero (wall material '{}')",
                                w.material
                            );
                            panic!("WallPlugin preflight should reject zero wall normals");
                        }
                        let nx = w.normal_x / mag;
                        let ny = w.normal_y / mag;
                        let nz = w.normal_z / mag;

                        let (motion, velocity) = if let Some(ref osc) = w.oscillate {
                            (
                                WallMotion::Oscillate {
                                    amplitude: osc.amplitude,
                                    frequency: osc.frequency,
                                },
                                [0.0; 3],
                            )
                        } else if let Some(ref srv) = w.servo {
                            (
                                WallMotion::Servo {
                                    target_force: srv.target_force,
                                    max_velocity: srv.max_velocity,
                                    gain: srv.gain,
                                },
                                [0.0; 3],
                            )
                        } else if let Some(vel) = w.velocity {
                            (WallMotion::ConstantVelocity, vel)
                        } else {
                            (WallMotion::Static, [0.0; 3])
                        };

                        planes.push(WallPlane {
                            point_x: w.point_x,
                            point_y: w.point_y,
                            point_z: w.point_z,
                            normal_x: nx,
                            normal_y: ny,
                            normal_z: nz,
                            material_index: mat_idx,
                            name: w.name.clone(),
                            bound_x_low: w.bound_x_low,
                            bound_x_high: w.bound_x_high,
                            bound_y_low: w.bound_y_low,
                            bound_y_high: w.bound_y_high,
                            bound_z_low: w.bound_z_low,
                            bound_z_high: w.bound_z_high,
                            velocity,
                            motion,
                            origin: [w.point_x, w.point_y, w.point_z],
                            force_accumulator: 0.0,
                            temperature: w.temperature,
                        });
                    }
                    other => {
                        eprintln!(
                            "ERROR: unknown wall type in [[wall]]: '{}'. Expected 'plane', 'cylinder', 'sphere', or 'region'",
                            other
                        );
                        panic!("WallPlugin preflight should reject unknown wall types");
                    }
                }
            }
            drop(material_table);

            let np = planes.len();
            let nc = cylinders.len();
            let ns = spheres.len();
            let nr = regions.len();
            Walls {
                planes,
                active: vec![true; np],
                cylinders,
                cylinder_active: vec![true; nc],
                spheres,
                sphere_active: vec![true; ns],
                regions,
                region_active: vec![true; nr],
                time: 0.0,
                tangential_springs: std::collections::HashMap::new(),
                rolling_springs: std::collections::HashMap::new(),
            }
        };

        app.add_resource(walls);
        app.add_update_system(wall_move, ParticleSimScheduleSet::PreInitialIntegration);
        app.add_update_system(
            wall_zero_force_accumulators,
            ParticleSimScheduleSet::PreForce,
        );
        app.add_update_system(
            wall_contact_force.label(WALL_CONTACT),
            ParticleSimScheduleSet::Force,
        );
    }

    fn try_build(&self, app: &mut App) -> Result<(), AppError> {
        validate_wall_plugin_config(app)?;
        self.build(app);
        Ok(())
    }
}

/// Validate wall input before the plugin mutates the app.  This is deliberately
/// shared with the legacy `build` entry point so callers using
/// `App::try_add_plugins` receive one typed, rank-consistent setup failure.
fn validate_wall_plugin_config(app: &mut App) -> Result<(), AppError> {
    let defs = {
        let config = app
            .get_resource_ref::<Config>()
            .ok_or_else(|| AppError::message("WallPlugin requires Config"))?;
        match config.table.get("wall") {
            Some(toml::Value::Array(entries)) => entries
                .iter()
                .enumerate()
                .map(|(index, entry)| {
                    entry.clone().try_into::<WallDef>().map_err(|error| {
                        AppError::message(format!(
                            "failed to parse [[wall]] entry {index}: {error}"
                        ))
                    })
                })
                .collect::<Result<Vec<_>, _>>()?,
            Some(toml::Value::Table(entry)) => vec![toml::Value::Table(entry.clone())
                .try_into::<WallDef>()
                .map_err(|error| {
                    AppError::message(format!("failed to parse [wall] entry: {error}"))
                })?],
            Some(_) => {
                return Err(AppError::message(
                    "[wall] must be a table or array of tables",
                ))
            }
            None => Vec::new(),
        }
    };
    let materials = app
        .get_resource_ref::<MaterialTable>()
        .ok_or_else(|| AppError::message("WallPlugin requires DemAtomPlugin"))?;

    for wall in defs {
        if materials.find_material(&wall.material).is_none() {
            return Err(AppError::message(format!(
                "wall material '{}' not found in [[dem.materials]]. Available: {:?}",
                wall.material, materials.names
            )));
        }
        match wall.wall_type.as_str() {
            "plane" => {
                let magnitude = (wall.normal_x * wall.normal_x
                    + wall.normal_y * wall.normal_y
                    + wall.normal_z * wall.normal_z)
                    .sqrt();
                if magnitude <= 1e-15 {
                    return Err(AppError::message(format!(
                        "wall normal vector must be non-zero (wall material '{}')",
                        wall.material
                    )));
                }
            }
            "cylinder" => {
                let axis = wall.axis.as_deref().unwrap_or("z");
                if !matches!(axis, "x" | "X" | "y" | "Y" | "z" | "Z") {
                    return Err(AppError::message(format!(
                        "cylinder wall axis must be x, y, or z, got '{axis}'"
                    )));
                }
                let center = wall.center.as_ref().ok_or_else(|| {
                    AppError::message("cylinder wall requires 'center' [c0, c1]")
                })?;
                if center.len() != 2 {
                    return Err(AppError::message("cylinder wall 'center' must have 2 elements"));
                }
                if wall.radius.is_none() {
                    return Err(AppError::message("cylinder wall requires 'radius'"));
                }
            }
            "sphere" => {
                let center = wall.center.as_ref().ok_or_else(|| {
                    AppError::message("sphere wall requires 'center' [x, y, z]")
                })?;
                if center.len() != 3 {
                    return Err(AppError::message("sphere wall 'center' must have 3 elements"));
                }
                if wall.radius.is_none() {
                    return Err(AppError::message("sphere wall requires 'radius'"));
                }
            }
            "region" => {
                if wall.region.is_none() {
                    return Err(AppError::message("region wall requires a 'region' field"));
                }
            }
            other => return Err(AppError::message(format!(
                "unknown wall type in [[wall]]: '{other}'. Expected 'plane', 'cylinder', 'sphere', or 'region'"
            ))),
        }
    }
    Ok(())
}
