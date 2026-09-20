//! Clump energy-conservation diagnostic (localizes the multisphere comm energy bug).
//!
//! Two identical sphere-7 clumps collide head-on, centrally, with a nearly
//! elastic frictionless material. The SAME collision is run two ways:
//!
//! - `interior` — the clumps meet in the middle of the box.
//! - `boundary` — the clumps meet exactly on the periodic seam (x = 0 ≡ x_high),
//!   so at contact one clump's leading sub-sphere interacts with the *periodic
//!   ghost image* of the other clump's sub-sphere.
//!
//! Both runs are physically identical, so they must dissipate an identical amount
//! of kinetic energy. If the `boundary` run injects (or loses) energy relative to
//! `interior`, the cross-seam clump communication path is the culprit.
//!
//! ```bash
//! for m in interior boundary; do CLUMP_ENERGY_MODE=$m cargo run --release \
//!   --example bench_clump_energy_conservation --no-default-features \
//!   --features precision-double -- examples/bench_clump_energy_conservation/config.toml; done
//! ```

use std::fs;
use std::io::Write;
use std::sync::OnceLock;

use dirt_core::prelude::*;

fn mode() -> String {
    std::env::var("CLUMP_ENERGY_MODE").unwrap_or_else(|_| "interior".to_string())
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(DemAtomPlugin)
        .add_plugins(DemAtomInsertPlugin)
        .add_plugins(HertzMindlinContactPlugin)
        .add_plugins(RotationalDynamicsPlugin)
        .add_plugins(ClumpPlugin);

    // Place the two clumps deterministically (runs at Setup, after ClumpPlugin
    // registered the "sphere7" definition; there is no [[clump.insert]]).
    app.add_setup_system(place_two_clumps, ScheduleSetupSet::Setup);

    app.add_update_system(measure_energy, ParticleSimScheduleSet::PostFinalIntegration);

    app.start();
}

/// Deterministically place two sphere-7 clumps for a head-on central collision.
fn place_two_clumps(
    comm: Res<CommResource>,
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    mut body_store: ResMut<MultisphereBodyStore>,
    clump_registry: Res<ClumpRegistry>,
    material_table: Res<dirt_atom::MaterialTable>,
) {
    if comm.rank() != 0 {
        return;
    }
    // Guard against setup re-running per run-stage (there is only one stage here).
    if !body_store.bodies.is_empty() {
        return;
    }

    let def = clump_registry
        .find("sphere7")
        .expect("clump definition 'sphere7' must exist");
    let mat_idx = material_table
        .names
        .iter()
        .position(|n| n == "glass")
        .expect("material 'glass' must exist") as u32;
    let density = 2500.0;

    let yc = 0.01;
    let zc = 0.01;
    let v = 0.5; // each clump's speed → closing speed 1.0 m/s
                 // COM offset from the collision point (COM separation = 2*half_gap).
    let half_gap: f64 = std::env::var("CLUMP_ENERGY_GAP")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.002);
    // Impact parameter (transverse offset of clump B): >0 makes the collision
    // glancing so, with friction, it induces spin — stressing the torque and
    // tangential cross-seam aggregation. Default 0.0008 m (~1.6 sub-radii).
    let b: f64 = std::env::var("CLUMP_ENERGY_B")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.0008);

    // interior: A left of center moving +x, B right of center moving -x → collide at x=0.01.
    // boundary: the SAME configuration rigidly translated by +half the box in x, so the
    // collision lands exactly on the periodic seam (x=0 ≡ 0.02). By periodic translational
    // invariance the two runs MUST be bit-identical; any divergence is the seam bug.
    let (mut xa, va, mut xb, vb) = (0.01 - half_gap, v, 0.01 + half_gap, -v);
    if mode() == "boundary" {
        let box_x = 0.02;
        xa = (xa + 0.5 * box_x).rem_euclid(box_x);
        xb = (xb + 0.5 * box_x).rem_euclid(box_x);
    }

    dirt_clump::insert_clump(
        &mut atoms,
        &registry,
        &mut body_store,
        def,
        [xa, yc, zc],
        [va, 0.0, 0.0],
        density,
        mat_idx,
        1,
    );
    dirt_clump::insert_clump(
        &mut atoms,
        &registry,
        &mut body_store,
        def,
        [xb, yc + b, zc],
        [vb, 0.0, 0.0],
        density,
        mat_idx,
        2,
    );
    body_store.generate_map();

    println!(
        "[energy_test] mode={} placed 2 clumps: A x={:.4} vx={:+.2}, B x={:.4} vx={:+.2}",
        mode(),
        xa,
        va,
        xb,
        vb
    );
}

/// Total body kinetic energy (translational + rotational), logged every 100 steps.
fn measure_energy(
    atoms: Res<Atom>,
    bodies: Res<MultisphereBodyStore>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
) {
    let step = run_state.total_cycle;
    if step % 100 != 0 {
        return;
    }
    if bodies.bodies.is_empty() {
        return;
    }

    let mut ke_trans = 0.0_f64;
    let mut ke_rot = 0.0_f64;
    for body in &bodies.bodies {
        let m = body.total_mass;
        let v2 = body.com_vel[0] * body.com_vel[0]
            + body.com_vel[1] * body.com_vel[1]
            + body.com_vel[2] * body.com_vel[2];
        ke_trans += 0.5 * m * v2;

        // rotational KE = 1/2 ω·I·ω, evaluated in the principal frame.
        let q = body.quaternion;
        let qc = [q[0], -q[1], -q[2], -q[3]];
        let omega_body = dirt_clump::quat_rotate(qc, body.omega);
        let qp = body.principal_axes;
        let qpc = [qp[0], -qp[1], -qp[2], -qp[3]];
        let wp = dirt_clump::quat_rotate(qpc, omega_body);
        ke_rot += 0.5
            * (body.principal_moments[0] * wp[0] * wp[0]
                + body.principal_moments[1] * wp[1] * wp[1]
                + body.principal_moments[2] * wp[2] * wp[2]);
    }
    let ke_trans = comm.all_reduce_sum_f64(ke_trans);
    let ke_rot = comm.all_reduce_sum_f64(ke_rot);
    let ke_total = ke_trans + ke_rot;

    if comm.rank() != 0 {
        return;
    }

    let time = step as f64 * atoms.dt;
    let output_dir = match input.output_dir.as_deref() {
        Some(dir) => dir.to_string(),
        None => "examples/bench_clump_energy_conservation/data".to_string(),
    };
    let filepath = format!("{}/energy_{}.csv", output_dir, mode());

    static INIT: OnceLock<()> = OnceLock::new();
    let _ = INIT.get_or_init(|| {
        fs::create_dir_all(&output_dir).ok();
        let mut f = fs::File::create(&filepath).expect("cannot create energy csv");
        writeln!(f, "step,time,ke_trans,ke_rot,ke_total").unwrap();
        ()
    });

    let mut f = fs::OpenOptions::new()
        .append(true)
        .open(&filepath)
        .expect("cannot open energy csv");
    writeln!(
        f,
        "{},{:.8e},{:.10e},{:.10e},{:.10e}",
        step, time, ke_trans, ke_rot, ke_total
    )
    .unwrap();
}
