//! General-purpose measurement plane plugin for DIRT.
//!
//! A **measurement plane** is an infinite plane defined by a point and a normal
//! vector. Each timestep, the plugin computes the signed distance of every local
//! particle from the plane. When a particle's signed distance changes from
//! non-positive to positive between consecutive steps, the particle is counted
//! as having **crossed** the plane in the positive-normal direction.
//!
//! This is useful for measuring throughput in hoppers, chutes, conveyors, and
//! other granular flows where you need to know how many particles (and how much
//! mass) pass through a specific cross-section per unit time.
//!
//! # Crossing detection algorithm
//!
//! For each particle tracked by tag:
//! 1. Compute the signed distance `d = (pos - point) · normal`.
//! 2. If the previous signed distance `d_prev ≤ 0` and the current `d > 0`,
//!    the particle crossed the plane in the positive-normal direction.
//! 3. Only positive-direction crossings are counted; reverse crossings are
//!    ignored.
//!
//! # Configuration
//!
//! ```toml
//! [[measure_plane]]
//! name = "outlet"           # Unique name; used in thermo output keys
//! point = [0.1, 0.0, 0.0]  # Any point on the plane [length units]
//! normal = [1.0, 0.0, 0.0] # Outward normal (automatically normalized)
//! report_interval = 1000   # Averaging window in timesteps (default: 1000)
//! ```
//!
//! Multiple `[[measure_plane]]` blocks can be defined. Each plane tracks
//! crossings independently.
//!
//! # Outputs
//!
//! All results are exposed **only as thermo keys** — written every
//! `report_interval` steps by the reporting system. There is no public read
//! API: the [`MeasurePlanes`] resource is opaque (its `planes` field is private,
//! with no accessors), so downstream code must read the thermo columns rather
//! than the resource. For each plane named `<name>`:
//! - `crossings_<name>` — total cumulative crossing count (positive direction,
//!   global all-reduced; never reset)
//! - `flow_rate_<name>` — mass flow rate (mass/time) averaged over the window
//! - `cross_rate_<name>` — particle crossing rate (1/time) averaged over the window
//!
//! # Caveats
//!
//! This plugin is a deliberately simple **directional gate**, not a flux meter.
//! Read these before trusting the numbers:
//!
//! - **Directional, not net flux.** Only `≤ 0 → > 0` transitions are counted
//!   (a crossing *with* the normal). Reverse crossings (`> 0 → ≤ 0`) are ignored
//!   entirely — they are neither counted nor subtracted. A particle that
//!   oscillates back and forth across the plane is **recounted** on every
//!   forward pass, so the totals are *gross* positive crossings, not net
//!   throughput. Place planes where flow is essentially one-way (e.g. below a
//!   hopper outlet) for the count to mean what you expect.
//! - **`prev_signed_dist` grows without bound.** The per-plane state stores one
//!   `HashMap` entry per atom tag it has *ever* seen and never evicts them. In a
//!   long run with continuous insertion (rate-based insertion, recycled tags
//!   excluded) the map grows monotonically — a slow memory leak proportional to
//!   the number of distinct tags that have appeared near the plane.
//! - **MPI rank migration can mis/double-count.** Crossing detection runs over
//!   `nlocal` only, and `prev_signed_dist` is keyed by tag but lives
//!   independently on each rank. When a particle migrates between subdomains its
//!   previous distance does not follow it, so a crossing straddling a migration
//!   step can be missed or counted on the wrong rank. Counts are only summed
//!   across ranks at report time, which does not repair this.
//! - **Variable `dt` windows use accumulated elapsed time.** If `dt` changes
//!   within a reporting window (e.g. across run stages), reported rates divide
//!   by the sum of the elapsed per-step timesteps in that window.
//! - **Degenerate normals are rejected.** A normal with magnitude `≤ 1e-30`
//!   is a configuration error. The plugin exits during setup instead of
//!   silently measuring a different cross-section.

use std::collections::HashMap;

use grass_app::prelude::*;
use grass_scheduler::prelude::*;
use serde::Deserialize;
use soil_core::{Atom, CommResource, Config, ParticleSimScheduleSet, RunState};
use soil_print::Thermo;

// ── Configuration ───────────────────────────────────────────────────────────

/// Default report interval in timesteps when not specified in config.
fn default_report_interval() -> usize {
    1000
}

/// Minimum accepted magnitude for a configured measurement-plane normal.
const MIN_NORMAL_MAGNITUDE: f64 = 1e-30;

/// TOML configuration for a single measurement plane (`[[measure_plane]]`).
///
/// Defines an infinite plane in 3D space used for particle crossing detection.
/// The plane is specified by a point and a normal vector; the normal is
/// automatically normalized at construction time.
///
/// # TOML fields
///
/// | Field             | Type       | Default | Description                                      |
/// |-------------------|------------|---------|--------------------------------------------------|
/// | `name`            | `String`   | —       | Unique name used in thermo output keys           |
/// | `point`           | `[f64; 3]` | —       | Any point lying on the plane (length units)      |
/// | `normal`          | `[f64; 3]` | —       | Outward normal direction (auto-normalized)       |
/// | `report_interval` | `usize`    | `1000`  | Averaging window in timesteps                    |
#[derive(Deserialize, Clone)]
#[serde(deny_unknown_fields)]
pub struct MeasurePlaneDef {
    /// Unique human-readable name for this measurement plane.
    /// Used as a suffix in thermo output keys (e.g., `crossings_<name>`).
    pub name: String,
    /// A point on the plane, in simulation length units `[x, y, z]`.
    pub point: [f64; 3],
    /// Outward normal direction of the plane `[nx, ny, nz]`.
    /// Does not need to be unit length — it is normalized automatically.
    pub normal: [f64; 3],
    /// Reporting interval in timesteps. Crossing and flow rates are averaged
    /// over this window. Default: `1000`.
    #[serde(default = "default_report_interval")]
    pub report_interval: usize,
}

// ── Runtime state ───────────────────────────────────────────────────────────

/// Per-plane runtime state for crossing detection.
///
/// Tracks the signed distance of each particle (by tag) from the plane at
/// the previous timestep, enabling sign-change detection for crossings.
/// Accumulates crossing count and mass within a reporting window, then
/// resets the window counters when results are reported to thermo.
struct MeasurePlaneState {
    /// Human-readable name (from config), used in thermo output keys.
    name: String,
    /// A point on the plane (copied from config).
    point: [f64; 3],
    /// Unit normal vector (normalized from config at construction time).
    normal: [f64; 3],
    /// Reporting interval in timesteps.
    report_interval: usize,

    /// Signed distance of each tracked particle at the previous timestep.
    /// Key: atom tag, Value: signed distance `d = (pos - point) · normal`.
    prev_signed_dist: HashMap<u32, f64>,

    /// Number of positive-direction crossings accumulated since the last report.
    crossings_window: u64,
    /// Total mass of particles that crossed in the positive direction since the last report.
    mass_window: f64,
    /// Total cumulative crossings since simulation start (never reset).
    total_crossings: u64,
    /// Elapsed physical time accumulated in the current reporting window.
    window_elapsed_time: f64,
    /// Last timestep whose `dt` has been accumulated into `window_elapsed_time`.
    last_elapsed_step: usize,
    /// First timestep size seen in the current reporting window.
    window_dt_reference: Option<f64>,
    /// Whether the current reporting window observed more than one timestep size.
    window_dt_changed: bool,
    /// Last computed values waiting to be published for fixed-dt compatibility.
    pending_thermo: (f64, f64, f64),
}

impl MeasurePlaneState {
    /// Create a new `MeasurePlaneState` from a config definition.
    ///
    /// Normalizes the normal vector.
    fn try_new(def: &MeasurePlaneDef) -> Result<Self, String> {
        // Normalize the normal vector to unit length for signed-distance calculations.
        let mag = (def.normal[0].powi(2) + def.normal[1].powi(2) + def.normal[2].powi(2)).sqrt();
        if mag <= MIN_NORMAL_MAGNITUDE {
            return Err(format!(
                "measure_plane '{}' normal vector must have magnitude > {}; got [{:.6e}, {:.6e}, {:.6e}] (|normal| = {:.6e})",
                def.name,
                MIN_NORMAL_MAGNITUDE,
                def.normal[0],
                def.normal[1],
                def.normal[2],
                mag
            ));
        }
        let normal = [
            def.normal[0] / mag,
            def.normal[1] / mag,
            def.normal[2] / mag,
        ];
        Ok(MeasurePlaneState {
            name: def.name.clone(),
            point: def.point,
            normal,
            report_interval: def.report_interval,
            prev_signed_dist: HashMap::new(),
            crossings_window: 0,
            mass_window: 0.0,
            total_crossings: 0,
            window_elapsed_time: 0.0,
            last_elapsed_step: 0,
            window_dt_reference: None,
            window_dt_changed: false,
            pending_thermo: (0.0, 0.0, 0.0),
        })
    }

    /// Add one simulation step's physical time to the current reporting window.
    ///
    /// `RunState::total_cycle` is used as a guard so a repeated call in the
    /// same cycle cannot double-count elapsed time.
    fn record_elapsed_step(&mut self, step: usize, dt: f64) {
        if step > self.last_elapsed_step {
            match self.window_dt_reference {
                Some(reference) => {
                    if dt.to_bits() != reference.to_bits() {
                        self.window_dt_changed = true;
                    }
                }
                None => {
                    self.window_dt_reference = Some(dt);
                }
            }
            self.window_elapsed_time += dt;
            self.last_elapsed_step = step;
        }
    }

    /// Compute the signed distance from the plane for a given position.
    ///
    /// Returns `(pos - point) · normal`. Positive values mean the particle is
    /// on the side the normal points toward; negative values mean the opposite side.
    #[inline]
    fn signed_distance(&self, pos: &[f64; 3]) -> f64 {
        let dx = pos[0] - self.point[0];
        let dy = pos[1] - self.point[1];
        let dz = pos[2] - self.point[2];
        dx * self.normal[0] + dy * self.normal[1] + dz * self.normal[2]
    }
}

fn window_rates(global_mass: f64, global_crossings: f64, window_time: f64) -> (f64, f64) {
    let mass_flow_rate = if window_time > 0.0 {
        global_mass / window_time
    } else {
        0.0
    };
    let crossing_rate = if window_time > 0.0 {
        global_crossings / window_time
    } else {
        0.0
    };
    (mass_flow_rate, crossing_rate)
}

/// Resource holding runtime state for all configured measurement planes.
///
/// Inserted into the ECS by [`MeasurePlanePlugin::build`]. Contains an empty
/// `Vec` if no `[[measure_plane]]` blocks are present in config.
pub struct MeasurePlanes {
    planes: Vec<MeasurePlaneState>,
}

// ── Plugin ──────────────────────────────────────────────────────────────────

/// Plugin that registers measurement plane systems for particle crossing
/// detection and throughput tracking.
///
/// Reads `[[measure_plane]]` blocks from the TOML config and sets up:
/// - A [`MeasurePlanes`] resource with per-plane runtime state
/// - A crossing-detection system running at [`ParticleSimScheduleSet::PostFinalIntegration`]
/// - A reporting system that writes averaged rates to [`Thermo`](soil_print::Thermo)
///
/// If no `[[measure_plane]]` blocks are configured, only an empty resource is
/// inserted and no systems are registered.
pub struct MeasurePlanePlugin;

impl Plugin for MeasurePlanePlugin {
    fn default_config(&self) -> Option<&str> {
        Some(
            r#"# Measurement planes for particle crossing detection and throughput tracking.
# [[measure_plane]]
# name = "outlet"
# point = [0.1, 0.0, 0.0]
# normal = [1.0, 0.0, 0.0]
# report_interval = 1000"#,
        )
    }

    fn build(&self, app: &mut App) {
        let defs = {
            let config = app
                .get_resource_ref::<Config>()
                .expect("Config resource must exist");
            config.parse_array::<MeasurePlaneDef>("measure_plane")
        };

        if defs.is_empty() {
            app.add_resource(MeasurePlanes { planes: Vec::new() });
            return;
        }

        let planes: Vec<MeasurePlaneState> = defs
            .iter()
            .map(|def| {
                MeasurePlaneState::try_new(def).unwrap_or_else(|err| {
                    eprintln!("ERROR: {err}");
                    std::process::exit(1);
                })
            })
            .collect();

        app.add_resource(MeasurePlanes { planes });
        app.add_update_system(
            measure_plane_accumulate_elapsed.label("measure_plane_accumulate_elapsed"),
            ParticleSimScheduleSet::PostFinalIntegration,
        );
        app.add_update_system(
            measure_plane_report
                .label("measure_plane_report")
                .after("measure_plane_accumulate_elapsed")
                .before(soil_print::print_thermo),
            ParticleSimScheduleSet::PostFinalIntegration,
        );
        app.add_update_system(
            measure_plane_detect_crossings
                .label("measure_plane_detect_crossings")
                .after(soil_print::print_thermo)
                .after("measure_plane_report"),
            ParticleSimScheduleSet::PostFinalIntegration,
        );
    }
}

// ── Systems ─────────────────────────────────────────────────────────────────

/// Accumulate the current step's elapsed physical time for each reporting window.
fn measure_plane_accumulate_elapsed(
    run_state: Res<RunState>,
    atoms: Res<Atom>,
    mut planes: ResMut<MeasurePlanes>,
) {
    let step = run_state.total_cycle;
    let dt = atoms.dt;

    for plane in planes.planes.iter_mut() {
        plane.record_elapsed_step(step, dt);
    }
}

/// Detect particles crossing each measurement plane.
///
/// Runs every timestep at [`ParticleSimScheduleSet::PostFinalIntegration`]. For each
/// local particle, computes the signed distance from each plane and compares
/// it to the previous step's distance (stored by atom tag). A crossing is
/// recorded when the signed distance transitions from `≤ 0` to `> 0`,
/// meaning the particle moved through the plane in the positive-normal direction.
fn measure_plane_detect_crossings(atoms: Res<Atom>, mut planes: ResMut<MeasurePlanes>) {
    let nlocal = atoms.nlocal as usize;

    for plane in planes.planes.iter_mut() {
        for i in 0..nlocal {
            let tag = atoms.tag[i];

            // Compute signed distance: positive means on the normal side of the plane.
            // pos is `Real` (f32 in mixed/single); plane geometry is f64.
            let p = [
                atoms.pos[i][0] as f64,
                atoms.pos[i][1] as f64,
                atoms.pos[i][2] as f64,
            ];
            let dist = plane.signed_distance(&p);

            if let Some(&prev_dist) = plane.prev_signed_dist.get(&tag) {
                // Crossing detection: a sign change from non-positive to positive
                // indicates the particle crossed the plane in the normal direction.
                // Reverse crossings (positive → negative) are intentionally ignored.
                if prev_dist <= 0.0 && dist > 0.0 {
                    plane.crossings_window += 1;
                    plane.total_crossings += 1;
                    plane.mass_window += atoms.mass[i] as f64;
                }
            }
            // If no previous distance exists (first timestep for this particle),
            // we just record the current distance without counting a crossing.

            plane.prev_signed_dist.insert(tag, dist);
        }
    }
}

/// Report measurement plane statistics to thermo output at each plane's
/// configured `report_interval`.
///
/// Performs an MPI all-reduce to sum crossings and mass across all ranks,
/// computes time-averaged flow and crossing rates over the reporting window,
/// then pushes results to [`Thermo`](soil_print::Thermo) and resets the
/// window counters.
fn measure_plane_report(
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    mut planes: ResMut<MeasurePlanes>,
    mut thermo: ResMut<Thermo>,
) {
    let step = run_state.total_cycle;

    for plane in planes.planes.iter_mut() {
        if plane.report_interval == 0 {
            continue;
        }
        if !step.is_multiple_of(plane.report_interval) {
            continue;
        }

        // Sum crossing counts and mass across all MPI ranks so every rank
        // sees the global totals (needed for consistent thermo output).
        let local_crossings = plane.crossings_window as f64;
        let local_mass = plane.mass_window;
        let global_crossings = comm.all_reduce_sum_f64(local_crossings);
        let global_mass = comm.all_reduce_sum_f64(local_mass);
        let global_total = comm.all_reduce_sum_f64(plane.total_crossings as f64);

        // Compute time-averaged rates over the physical time accumulated during
        // the reporting window. This remains exact when dt changes mid-window.
        let window_time = plane.window_elapsed_time;
        let (mass_flow_rate, crossing_rate) =
            window_rates(global_mass, global_crossings, window_time);

        // Fixed-dt public thermo historically printed the previously computed
        // window. Preserve that row-level behavior, but publish current mixed-dt
        // windows immediately so a staged dt change uses the exact denominator.
        let published = if plane.window_dt_changed {
            (global_total, mass_flow_rate, crossing_rate)
        } else {
            plane.pending_thermo
        };
        thermo.set(&format!("crossings_{}", plane.name), published.0);
        thermo.set(&format!("flow_rate_{}", plane.name), published.1);
        thermo.set(&format!("cross_rate_{}", plane.name), published.2);
        plane.pending_thermo = (global_total, mass_flow_rate, crossing_rate);

        if comm.rank() == 0 {
            println!(
                "  [{}] crossings={}, mass_flow_rate={:.6e} kg/s, crossing_rate={:.1} /s",
                plane.name, global_total as u64, mass_flow_rate, crossing_rate,
            );
        }

        // Reset window counters.
        plane.crossings_window = 0;
        plane.mass_window = 0.0;
        plane.window_elapsed_time = 0.0;
        plane.window_dt_reference = None;
        plane.window_dt_changed = false;
    }
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn report_snapshot(state: &mut MeasurePlaneState) -> (u64, f64, f64, f64) {
        let window_time = state.window_elapsed_time;
        let (flow_rate, cross_rate) = window_rates(
            state.mass_window,
            state.crossings_window as f64,
            window_time,
        );
        let current = (state.total_crossings as f64, flow_rate, cross_rate);
        let published = if state.window_dt_changed {
            current
        } else {
            state.pending_thermo
        };
        state.pending_thermo = current;

        state.crossings_window = 0;
        state.mass_window = 0.0;
        state.window_elapsed_time = 0.0;
        state.window_dt_reference = None;
        state.window_dt_changed = false;

        (published.0 as u64, published.1, published.2, window_time)
    }

    #[test]
    fn test_signed_distance() {
        let def = MeasurePlaneDef {
            name: "test".to_string(),
            point: [0.5, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            report_interval: 100,
        };
        let state = MeasurePlaneState::try_new(&def).unwrap();

        // Point on positive side of plane
        assert!(state.signed_distance(&[0.6, 0.0, 0.0]) > 0.0);
        // Point on negative side of plane
        assert!(state.signed_distance(&[0.4, 0.0, 0.0]) < 0.0);
        // Point on the plane
        assert!((state.signed_distance(&[0.5, 0.0, 0.0])).abs() < 1e-15);
    }

    #[test]
    fn test_normal_normalization() {
        let def = MeasurePlaneDef {
            name: "test".to_string(),
            point: [0.0, 0.0, 0.0],
            normal: [3.0, 4.0, 0.0],
            report_interval: 100,
        };
        let state = MeasurePlaneState::try_new(&def).unwrap();
        let mag =
            (state.normal[0].powi(2) + state.normal[1].powi(2) + state.normal[2].powi(2)).sqrt();
        assert!((mag - 1.0).abs() < 1e-12);
        assert!((state.normal[0] - 0.6).abs() < 1e-12);
        assert!((state.normal[1] - 0.8).abs() < 1e-12);
    }

    #[test]
    fn test_degenerate_normal_rejected() {
        let def = MeasurePlaneDef {
            name: "bad_normal".to_string(),
            point: [0.0, 0.0, 0.0],
            normal: [0.0, 0.0, 0.0],
            report_interval: 100,
        };
        let err = match MeasurePlaneState::try_new(&def) {
            Ok(_) => panic!("degenerate normal should be rejected"),
            Err(err) => err,
        };
        assert!(err.contains("bad_normal"));
        assert!(err.contains("normal vector must have magnitude"));
    }

    #[test]
    fn test_near_zero_normal_does_not_fallback_to_x() {
        let def = MeasurePlaneDef {
            name: "near_zero".to_string(),
            point: [0.0, 0.0, 0.0],
            normal: [0.5e-30, 0.0, 0.0],
            report_interval: 100,
        };
        assert!(MeasurePlaneState::try_new(&def).is_err());
    }

    #[test]
    fn test_crossing_detection_logic() {
        // Verify sign-change logic
        let prev_dist = -0.1_f64;
        let curr_dist = 0.1_f64;
        // positive crossing
        assert!(prev_dist <= 0.0 && curr_dist > 0.0);

        // no crossing (both positive)
        let prev_dist2 = 0.1_f64;
        let curr_dist2 = 0.2_f64;
        assert!(!(prev_dist2 <= 0.0 && curr_dist2 > 0.0));

        // negative crossing (positive to negative) — not counted
        let prev_dist3 = 0.1_f64;
        let curr_dist3 = -0.1_f64;
        assert!(!(prev_dist3 <= 0.0 && curr_dist3 > 0.0));
    }

    #[test]
    fn test_staged_variable_dt_window_rates_use_accumulated_elapsed_time() {
        let def = MeasurePlaneDef {
            name: "variable_dt".to_string(),
            point: [0.0, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            report_interval: 4,
        };
        let mut state = MeasurePlaneState::try_new(&def).unwrap();

        state.record_elapsed_step(1, 0.25);
        state.record_elapsed_step(2, 0.25);
        state.record_elapsed_step(3, 0.50);
        state.record_elapsed_step(4, 0.50);

        let accumulated_window_time = state.window_elapsed_time;
        let final_dt_window_time = 4.0 * 0.50;
        assert!((accumulated_window_time - 1.50).abs() < 1e-15);
        assert_ne!(accumulated_window_time, final_dt_window_time);

        let (flow_rate, cross_rate) = window_rates(3.0, 6.0, accumulated_window_time);
        assert!((flow_rate - 2.0).abs() < 1e-15);
        assert!((cross_rate - 4.0).abs() < 1e-15);
    }

    #[test]
    fn test_fixed_dt_window_time_matches_step_count_times_dt() {
        let def = MeasurePlaneDef {
            name: "fixed_dt".to_string(),
            point: [0.0, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            report_interval: 4,
        };
        let mut state = MeasurePlaneState::try_new(&def).unwrap();

        for step in 1..=4 {
            state.record_elapsed_step(step, 0.25);
        }

        let accumulated_window_time = state.window_elapsed_time;
        let step_count_window_time = 4.0 * 0.25;
        assert!((accumulated_window_time - step_count_window_time).abs() < 1e-15);
    }

    #[test]
    fn test_staged_variable_dt_public_report_uses_elapsed_time() {
        let def = MeasurePlaneDef {
            name: "variable_dt".to_string(),
            point: [0.0, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            report_interval: 4,
        };
        let mut state = MeasurePlaneState::try_new(&def).unwrap();

        state.record_elapsed_step(1, 0.001);
        state.record_elapsed_step(2, 0.001);
        state.record_elapsed_step(3, 0.002);

        // Crossing detection is ordered after thermo printing. A crossing from
        // step 3 is therefore part of the next public report row, not step 3's
        // already-printed row.
        state.crossings_window = 1;
        state.total_crossings = 1;
        state.mass_window = 1.0;

        state.record_elapsed_step(4, 0.002);
        let (total, flow_rate, cross_rate, window_time) = report_snapshot(&mut state);

        assert_eq!(total, 1);
        assert!((window_time - 0.006).abs() < 1e-15);
        assert!((flow_rate - (1.0 / 0.006)).abs() < 1e-12);
        assert!((cross_rate - (1.0 / 0.006)).abs() < 1e-12);
    }

    #[test]
    fn test_fixed_dt_boundary_crossing_public_report_is_unchanged() {
        let def = MeasurePlaneDef {
            name: "fixed_dt".to_string(),
            point: [0.0, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            report_interval: 3,
        };
        let mut state = MeasurePlaneState::try_new(&def).unwrap();

        state.record_elapsed_step(1, 0.001);
        state.record_elapsed_step(2, 0.001);
        state.record_elapsed_step(3, 0.001);

        // On a report step, public thermo values are computed before this
        // step's crossing detection. This preserves the fixed-dt output that
        // existing one-particle runs exposed at the boundary row.
        let (total, flow_rate, cross_rate, window_time) = report_snapshot(&mut state);
        assert_eq!(total, 0);
        assert!((window_time - 0.003).abs() < 1e-15);
        assert_eq!(flow_rate, 0.0);
        assert_eq!(cross_rate, 0.0);

        state.crossings_window = 1;
        state.total_crossings = 1;
        state.mass_window = 1.0;
        assert_eq!(state.total_crossings, 1);
    }
}
