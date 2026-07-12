//! Normal-stress-fixed, wall-driven flexible-fiber shear cell.
//!
//! The lower wall is held during preparation and then translated in `x`; the
//! upper wall retains its force-servo in `y`.  Measurements are made from the
//! live wall reactions and current particle positions, never from a closure.

use std::collections::VecDeque;
use std::fs;
use std::io::Write as IoWrite;
use std::sync::Once;

use dirt_core::dirt_bond::DemBondPlugin;
use dirt_core::prelude::*;

// The normal force controller needs a finite physical interval to compact this
// population.  This interval is intentionally separate from the shearing
// stage: the lower wall must not move while the prescribed normal load is
// still being approached.
/// Let the initially loose packing fall onto the lower wall before the lid is
/// allowed to move. Starting a downward force servo while the packing is in
/// free fall makes it chase the bed and can report a misleading zero reaction.
const SETTLE_STEPS: usize = 100_000;
const SAMPLE_EVERY: usize = 1_000;
const QUALIFICATION_SAMPLES: usize = 40;
const NORMAL_REL_TOL: f64 = 0.15;
const LOWER_WALL: &str = "lower_drive";
const LID: &str = "normal_load_lid";
const BEADS_PER_FIBRE: f64 = 17.0;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin)
        .add_plugins(DemBondPlugin);
    app.add_resource(LoadQualification::default());
    app.add_resource(ShearStart::default());
    app.add_resource(LidStage::default());
    app.add_update_system(stage_lid, ParticleSimScheduleSet::Setup);
    app.add_update_system(enable_shear, ParticleSimScheduleSet::PreInitialIntegration);
    app.add_update_system(record_cell, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

/// Captured servo settings while the lid is deliberately held during gravity
/// settling. This makes the two-stage protocol executable rather than relying
/// on a lucky initial lid height.
#[derive(Default)]
struct LidStage {
    servo: Option<(f64, f64, f64)>,
}

fn stage_lid(run_state: Res<RunState>, mut walls: ResMut<Walls>, mut stage: ResMut<LidStage>) {
    let lid = walls
        .planes
        .iter_mut()
        .find(|wall| wall.name.as_deref() == Some(LID))
        .expect("Guo cell requires normal_load_lid");
    if run_state.total_cycle < SETTLE_STEPS {
        if stage.servo.is_none() {
            stage.servo = match lid.motion {
                WallMotion::Servo {
                    target_force,
                    max_velocity,
                    gain,
                } => Some((target_force, max_velocity, gain)),
                WallMotion::Static => None,
                _ => panic!("Guo cell lid must be a force-servo wall"),
            };
        }
        lid.motion = WallMotion::Static;
        lid.velocity = [0.0; 3];
    } else if let Some((target_force, max_velocity, gain)) = stage.servo.take() {
        lid.motion = WallMotion::Servo {
            target_force,
            max_velocity,
            gain,
        };
    }
}

/// Starts the tangential lower-wall drive only after gravitational settling and
/// force-servo compaction.  The lid's normal controller remains active.
fn enable_shear(
    run_state: Res<RunState>,
    qualification: Res<LoadQualification>,
    mut walls: ResMut<Walls>,
    mut shear_start: ResMut<ShearStart>,
) {
    // Do not turn a nominal number of settling steps into an implicit claim
    // that the prescribed normal load was attained.  Shearing starts only
    // after the same measured-load band that the post-processing gate uses.
    if run_state.total_cycle < SETTLE_STEPS || !qualification.is_qualified() {
        return;
    }
    let wall = walls
        .planes
        .iter_mut()
        .find(|w| w.name.as_deref() == Some(LOWER_WALL))
        .expect("Guo cell requires a lower_drive wall");
    if wall.velocity[0] == 0.0 {
        wall.velocity[0] = 0.020; // paper's 20 mm/s case; y remains immobile.
        shear_start.step = Some(run_state.total_cycle);
    }
}

/// Actual cycle on which the qualified tangential drive starts.
#[derive(Default)]
struct ShearStart {
    step: Option<usize>,
}

/// Bounded recorder history used to qualify the actual lid load before shear.
/// This is deliberately benchmark-local: it is a protocol decision, not a
/// general DEM data model.
#[derive(Default)]
struct LoadQualification {
    stresses: VecDeque<f64>,
    target: Option<f64>,
}

impl LoadQualification {
    fn observe(&mut self, stress: f64, target: f64) {
        self.target = Some(target);
        self.stresses.push_back(stress);
        if self.stresses.len() > QUALIFICATION_SAMPLES {
            self.stresses.pop_front();
        }
    }

    fn is_qualified(&self) -> bool {
        if self.stresses.len() != QUALIFICATION_SAMPLES {
            return false;
        }
        let target = self.target.expect("qualification target is recorded");
        let mean = self.stresses.iter().sum::<f64>() / QUALIFICATION_SAMPLES as f64;
        (mean - target).abs() / target <= NORMAL_REL_TOL
    }
}

/// Recorder for the two paper observables: lid shear reaction / area and
/// solid fraction.  The normal reaction is retained as a servo-convergence
/// diagnostic, not substituted by its requested set point.
fn record_cell(
    atoms: Res<Atom>,
    domain: Res<Domain>,
    walls: Res<Walls>,
    run_state: Res<RunState>,
    input: Res<Input>,
    comm: Res<CommResource>,
    mut qualification: ResMut<LoadQualification>,
    shear_start: Res<ShearStart>,
) {
    let step = run_state.total_cycle;
    if step % SAMPLE_EVERY != 0 {
        return;
    }
    let lid = walls
        .planes
        .iter()
        .find(|w| w.name.as_deref() == Some(LID))
        .expect("Guo cell requires normal_load_lid");
    let lower = walls
        .planes
        .iter()
        .find(|w| w.name.as_deref() == Some(LOWER_WALL))
        .expect("Guo cell requires lower_drive");
    // Plane walls are replicated on every rank whereas contacts are local.
    // Use the complete-wall reactions for both the normal-load protocol and
    // the reported stress; otherwise rank zero would qualify and report a
    // different physical cell from the one the other ranks advance.
    let lid_force = comm.all_reduce_sum_f64(lid.force_accumulator);
    let lower_x_reaction = comm.all_reduce_sum_f64(lower.force_vector[0]);
    // Table 2 represents a 21.6-mm, 2.4-mm-diameter cord by 17 spheres at
    // 1.2-mm spacing. Summing sphere volumes counts their overlaps twice.
    // Record the physical spherocylinder volume so this is comparable with
    // the paper's rubber-cord solid fraction, while retaining the actual
    // DEM contact geometry for the solver itself.
    let local_fibres = atoms.nlocal as f64 / BEADS_PER_FIBRE;
    let radius = 0.0012_f64;
    let fibre_length = 0.0216_f64;
    let physical_fibre_volume =
        std::f64::consts::PI * radius.powi(2) * (fibre_length - 2.0 * radius)
            + 4.0 * std::f64::consts::PI * radius.powi(3) / 3.0;
    let solid = comm.all_reduce_sum_f64(local_fibres * physical_fibre_volume);
    // `nlocal` is per-rank in a decomposed cell.  The campaign contract is
    // about the physical fibre population, so record the global count.
    let global_atoms = comm.all_reduce_sum_f64(atoms.nlocal as f64).round() as usize;
    // The published numerical comparator is periodic in x/z: use its actual
    // rectangular planform, never the obsolete circular-cup area.
    let area = (domain.boundaries_high[0] - domain.boundaries_low[0])
        * (domain.boundaries_high[2] - domain.boundaries_low[2]);
    // During the explicit gravity-settle stage the lid is intentionally static;
    // only servo-stage measurements may qualify the normal load.
    let target_force = match &lid.motion {
        WallMotion::Servo { target_force, .. } => Some(*target_force),
        WallMotion::Static if step < SETTLE_STEPS => None,
        _ => panic!("Guo cell lid must be static only during settle or force-servo controlled"),
    };
    let height = lid.point_y - lower.point_y;
    if let Some(target_force) = target_force {
        // Every rank must acquire the same qualification state before its
        // local lower-wall copy may begin translating.
        qualification.observe(lid_force / area, target_force / area);
    }
    if comm.rank() != 0 {
        return;
    }
    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| "examples/bench_guo2018_fiber_shear_cell/generated".to_owned());
    fs::create_dir_all(&out_dir).expect("cannot create output directory");
    let path = format!("{out_dir}/cell_history.csv");
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        let mut f = fs::File::create(&path).expect("cannot create cell history");
        writeln!(f, "step,time,stage,shear_strain,lid_force_n,lid_x_reaction_n,normal_stress_pa,shear_stress_pa,solid_fraction,height_m,n_atoms")
            .unwrap();
    });
    let mut f = fs::OpenOptions::new()
        .append(true)
        .open(&path)
        .expect("cannot append cell history");
    writeln!(
        f,
        "{},{},{},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{}",
        step,
        step as f64 * atoms.dt,
        if step < SETTLE_STEPS {
            "settle"
        } else if shear_start.step.is_some() {
            "shear"
        } else {
            "normal_load"
        },
        shear_start.step.map_or(0.0, |start| {
            lower.velocity[0] * step.saturating_sub(start) as f64 * atoms.dt / height
        }),
        lid_force,
        lower_x_reaction,
        lid_force / area,
        lower_x_reaction.abs() / area,
        solid / (area * height),
        height,
        global_atoms
    )
    .unwrap();
}
