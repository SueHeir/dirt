//! Through-thickness tensile probe for TANGLE's polished fake felt.
//!
//! The bottom and top 15% of capsule centers are captured once as persistent,
//! migration-safe grip labels. The grips then move symmetrically in `z`, while
//! DIRT reports reaction force, through-thickness stress, and engineering strain.
//!
//! ```bash
//! cargo build --release --example fake_felted_tension
//! mpirun -np 8 target/release/examples/fake_felted_tension \
//!   examples/fake_felted_tension/config_smoke.toml
//! ```

use std::fs::{self, File};
use std::io::{BufWriter, Write as IoWrite};

use dirt_capsule::{CapsuleContactPlugin, CapsulePlugin};
use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use dirt_core::soil_core::{register_atom_data, AtomData, BondStore, REVERSE_SEND_FORCE};
use soil_derive::AtomData;

const GRIP_FRACTION: f64 = 0.15;
const DEFAULT_SETTLE_STEPS: u64 = 5;
const DEFAULT_RAMP_STEPS: u64 = 500;
const DEFAULT_GRIP_SPEED: f64 = 0.5;
const DEFAULT_RECORD_EVERY: u64 = 1;

/// Persistent grip membership: -1 bottom, 0 free, +1 top.
///
/// This is registered atom data rather than a local tag set, so membership is
/// permuted with neighbor sorting and migrates with its owning capsule in MPI.
#[derive(AtomData, Default)]
struct GripLabel {
    side: Vec<f64>,
}

#[derive(Default)]
struct TensionState {
    initialized: bool,
    initial_grip_center_gap: f64,
    initial_gauge_length: f64,
    area: f64,
    grip_speed: f64,
    writer: Option<BufWriter<File>>,
}

fn env_u64(name: &str, fallback: u64) -> u64 {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(fallback)
}

fn env_f64(name: &str, fallback: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(fallback)
}

fn settle_steps() -> u64 {
    env_u64("TENSION_SETTLE_STEPS", DEFAULT_SETTLE_STEPS)
}

fn grip_speed() -> f64 {
    env_f64("TENSION_GRIP_SPEED", DEFAULT_GRIP_SPEED)
}

fn ramp_steps() -> u64 {
    env_u64("TENSION_RAMP_STEPS", DEFAULT_RAMP_STEPS)
}

fn commanded_grip_speed(step: u64, maximum_speed: f64) -> f64 {
    let pull_step = step.saturating_sub(settle_steps());
    let ramp = ramp_steps();
    let fraction = if step < settle_steps() {
        0.0
    } else if ramp == 0 {
        1.0
    } else {
        ((pull_step + 1) as f64 / ramp as f64).min(1.0)
    };
    fraction * maximum_speed
}

fn requested_final_strain_percent() -> Option<f64> {
    std::env::var("TENSION_FINAL_STRAIN_PERCENT")
        .ok()
        .map(|value| {
            value
                .parse::<f64>()
                .expect("TENSION_FINAL_STRAIN_PERCENT must be a number")
        })
}

fn speed_for_final_strain(
    final_strain_percent: f64,
    initial_gap: f64,
    dt: f64,
    total_steps: u64,
) -> f64 {
    assert!(
        final_strain_percent >= 0.0,
        "TENSION_FINAL_STRAIN_PERCENT must be non-negative"
    );
    assert!(dt > 0.0, "timestep must be positive");

    let pull_steps = total_steps.saturating_sub(settle_steps());
    assert!(pull_steps > 0, "run must include steps after grip settling");
    let ramp = ramp_steps();
    let ramped_steps = pull_steps.min(ramp);
    let ramp_weight = if ramp == 0 {
        0.0
    } else {
        ramped_steps as f64 * (ramped_steps + 1) as f64 / (2.0 * ramp as f64)
    };
    let full_speed_weight = (pull_steps - ramped_steps) as f64;
    let effective_steps = ramp_weight + full_speed_weight;

    (final_strain_percent / 100.0) * initial_gap / (2.0 * dt * effective_steps)
}

fn record_every() -> u64 {
    env_u64("TENSION_RECORD_EVERY", DEFAULT_RECORD_EVERY).max(1)
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(DemAtomPlugin)
        .add_plugins(DemAtomInsertPlugin)
        .add_plugins(dirt_core::soil_verlet::VelocityVerletPlugin::new())
        .add_plugins(CapsuleContactPlugin)
        .add_plugins(CapsulePlugin)
        .add_plugins(DemBondPlugin);

    register_atom_data!(&mut app, GripLabel::default());
    {
        let dump = app
            .get_resource_ref::<DumpRegistry>()
            .expect("CorePlugins should install DumpRegistry");
        dump.register_scalar("grip_side", |atoms, registry| {
            let grips = registry.expect::<GripLabel>("tension OVITO grip labels");
            grips.side[..atoms.nlocal() as usize].to_vec()
        });
    }
    app.add_resource(TensionState::default());
    app.add_update_system(
        initialize_and_drive_grips,
        ParticleSimScheduleSet::PreInitialIntegration,
    );
    app.add_update_system(
        measure_and_clamp_grips.after(REVERSE_SEND_FORCE),
        ParticleSimScheduleSet::PostForce,
    );
    app.start();
}

fn initialize_and_drive_grips(
    mut atoms: ResMut<Atom>,
    particles: ParticlesWith<'_, (Write<GripLabel>, Write<DemAtom>)>,
    comm: Res<CommResource>,
    domain: Res<Domain>,
    input: Res<Input>,
    run_config: Res<RunConfig>,
    run_state: Res<RunState>,
    mut state: ResMut<TensionState>,
) {
    let nlocal = atoms.nlocal() as usize;
    let mut local_min = f64::INFINITY;
    let mut local_max = f64::NEG_INFINITY;
    for i in 0..nlocal {
        let z = atoms.pos[i][2] as f64;
        local_min = local_min.min(z);
        local_max = local_max.max(z);
    }
    let global_min = comm.all_reduce_min_f64(local_min);
    let global_max = -comm.all_reduce_min_f64(-local_max);

    particles.with(|(mut grips, mut dem)| {
        if !state.initialized {
            let depth = GRIP_FRACTION * (global_max - global_min);
            let bottom_limit = global_min + depth;
            let top_limit = global_max - depth;
            let mut bottom_sum = 0.0;
            let mut top_sum = 0.0;
            let mut bottom_count = 0.0;
            let mut top_count = 0.0;

            for i in 0..nlocal {
                let z = atoms.pos[i][2] as f64;
                let side = if z <= bottom_limit {
                    bottom_sum += z;
                    bottom_count += 1.0;
                    -1.0
                } else if z >= top_limit {
                    top_sum += z;
                    top_count += 1.0;
                    1.0
                } else {
                    0.0
                };
                grips.side[i] = side;
            }

            bottom_sum = comm.all_reduce_sum_f64(bottom_sum);
            top_sum = comm.all_reduce_sum_f64(top_sum);
            bottom_count = comm.all_reduce_sum_f64(bottom_count);
            top_count = comm.all_reduce_sum_f64(top_count);
            assert!(bottom_count > 0.0 && top_count > 0.0, "both grips need particles");
            state.initial_grip_center_gap = top_sum / top_count - bottom_sum / bottom_count;
            // Engineering strain is referenced to the material between the
            // inner grip boundaries, not to the distance between grip centers.
            state.initial_gauge_length = top_limit - bottom_limit;
            state.area = domain.size[0] * domain.size[1];
            let total_steps = run_config
                .stages
                .iter()
                .map(|stage| stage.steps as u64)
                .sum();
            state.grip_speed = requested_final_strain_percent()
                .map(|strain| {
                    speed_for_final_strain(
                        strain,
                        state.initial_gauge_length,
                        atoms.dt as f64,
                        total_steps,
                    )
                })
                .unwrap_or_else(grip_speed);

            if comm.rank() == 0 {
                let output = input
                    .output_dir
                    .as_deref()
                    .unwrap_or("examples/fake_felted_tension/output");
                fs::create_dir_all(output).expect("create tension output directory");
                let path = format!("{output}/tension.csv");
                let mut writer = BufWriter::new(File::create(&path).expect("create tension.csv"));
                writeln!(
                    writer,
                    "step,time,phase,gap,displacement,gauge_length,strain,bottom_reaction_z,top_reaction_z,tensile_force,engineering_stress,bond_count,bonds_broken"
                )
                .expect("write tension.csv header");
                state.writer = Some(writer);
                println!("=== TANGLE fake-felted through-thickness tension ===");
                println!(
                    "  specimen z extent {:.3} um; grip depth {:.3} um ({:.0}% each)",
                    1.0e6 * (global_max - global_min),
                    1.0e6 * depth,
                    100.0 * GRIP_FRACTION
                );
                println!(
                    "  bottom grip {} capsules; top grip {} capsules; initial mean-center gap {:.3} um",
                    bottom_count as usize,
                    top_count as usize,
                    1.0e6 * state.initial_grip_center_gap
                );
                println!(
                    "  ungripped gauge length {:.3} um (engineering-strain reference)",
                    1.0e6 * state.initial_gauge_length
                );
                println!(
                    "  symmetric grip speed +/-{:.3} m/s after {} settling steps and a {}-step ramp",
                    state.grip_speed,
                    settle_steps(),
                    ramp_steps()
                );
                if let Some(strain) = requested_final_strain_percent() {
                    println!("  prescribed final engineering strain: {strain:.4}%");
                }
                println!("  force-displacement output: {path}");
            }
            state.initialized = true;
        }

        let speed = commanded_grip_speed(run_state.total_cycle as u64, state.grip_speed);
        for i in 0..nlocal {
            let side = grips.side[i];
            if side == 0.0 {
                continue;
            }
            atoms.vel[i] = [0.0, 0.0, (side * speed) as Real];
            dem.omega[i] = [0.0; 3];
            dem.ang_mom[i] = [0.0; 3];
        }
    });
}

fn measure_and_clamp_grips(
    mut atoms: ResMut<Atom>,
    particles: ParticlesWith<'_, (Read<GripLabel>, Write<DemAtom>, Read<BondStore>)>,
    comm: Res<CommResource>,
    run_state: Res<RunState>,
    metrics: Option<Res<dirt_core::dirt_bond::BondMetrics>>,
    mut state: ResMut<TensionState>,
) {
    if !state.initialized {
        return;
    }
    let nlocal = atoms.nlocal() as usize;
    let should_record = (run_state.total_cycle as u64) % record_every() == 0;
    let mut bottom_z = 0.0;
    let mut top_z = 0.0;
    let mut bottom_count = 0.0;
    let mut top_count = 0.0;
    let mut bottom_reaction = 0.0;
    let mut top_reaction = 0.0;
    let mut local_bonds = 0usize;

    particles.with(|(grips, mut dem, bonds)| {
        for i in 0..nlocal {
            let side = grips.side[i];
            if should_record {
                if side < 0.0 {
                    bottom_z += atoms.pos[i][2] as f64;
                    bottom_count += 1.0;
                    bottom_reaction += atoms.force[i][2] as f64;
                } else if side > 0.0 {
                    top_z += atoms.pos[i][2] as f64;
                    top_count += 1.0;
                    top_reaction += atoms.force[i][2] as f64;
                }
            }
            if side != 0.0 {
                atoms.force[i] = [0.0; 3];
                dem.torque[i] = [0.0; 3];
            }
        }
        if should_record {
            // Each bond is stored at both endpoints; count only the lower-tag owner.
            for i in 0..nlocal.min(bonds.bonds.len()) {
                local_bonds += bonds.bonds[i]
                    .iter()
                    .filter(|bond| atoms.tag[i] < bond.partner_tag)
                    .count();
            }
        }
    });

    if !should_record {
        return;
    }

    bottom_z = comm.all_reduce_sum_f64(bottom_z);
    top_z = comm.all_reduce_sum_f64(top_z);
    bottom_count = comm.all_reduce_sum_f64(bottom_count);
    top_count = comm.all_reduce_sum_f64(top_count);
    bottom_reaction = comm.all_reduce_sum_f64(bottom_reaction);
    top_reaction = comm.all_reduce_sum_f64(top_reaction);
    let global_bonds = comm.all_reduce_sum_f64(local_bonds as f64) as usize;
    let local_broken = metrics
        .as_deref()
        .map(|value| value.total_bonds_broken)
        .unwrap_or(0);
    let global_broken = comm.all_reduce_sum_f64(local_broken as f64) as usize;

    if comm.rank() != 0 {
        return;
    }
    let gap = top_z / top_count - bottom_z / bottom_count;
    let displacement = gap - state.initial_grip_center_gap;
    let strain = displacement / state.initial_gauge_length;
    let tensile_force = 0.5 * (bottom_reaction.abs() + top_reaction.abs());
    let stress = tensile_force / state.area;
    let phase = if (run_state.total_cycle as u64) < settle_steps() {
        "settle"
    } else {
        "pull"
    };
    let time = run_state.total_cycle as f64 * atoms.dt;
    let gauge_length = state.initial_gauge_length;

    if let Some(writer) = state.writer.as_mut() {
        writeln!(
            writer,
            "{},{:.9e},{},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{},{}",
            run_state.total_cycle,
            time,
            phase,
            gap,
            displacement,
            gauge_length,
            strain,
            bottom_reaction,
            top_reaction,
            tensile_force,
            stress,
            global_bonds,
            global_broken
        )
        .expect("write tension row");
        writer.flush().expect("flush tension.csv");
    }
    println!(
        "[tension] step={} phase={} strain={:.5}% force={:.4e} N stress={:.4e} Pa bonds={} broken={}",
        run_state.total_cycle,
        phase,
        100.0 * strain,
        tensile_force,
        stress,
        global_bonds,
        global_broken
    );
}
