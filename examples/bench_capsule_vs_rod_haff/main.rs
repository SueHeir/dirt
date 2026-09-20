//! Haff cooling: single-body sphero-cylinder (capsule) gas vs. rod-of-spheres
//! (multisphere clump) gas.
//!
//! Both are elongated non-spherical particles cooling in a periodic box with no
//! gravity; a correct capsule implementation should follow Haff's law just like
//! the validated clump gas, and the two cooling curves should broadly agree.
//! Mode is selected by `CAPSULE_VS_ROD_MODE` (`capsule` | `rod`).
//!
//! ```bash
//! CAPSULE_VS_ROD_MODE=capsule cargo run --release --example bench_capsule_vs_rod_haff \
//!   --no-default-features --features precision-double -- \
//!   examples/bench_capsule_vs_rod_haff/config_capsule.toml
//! CAPSULE_VS_ROD_MODE=rod cargo run --release --example bench_capsule_vs_rod_haff \
//!   --no-default-features --features precision-double -- \
//!   examples/bench_capsule_vs_rod_haff/config_rod.toml
//! ```

use std::fs;
use std::io::Write;
use std::sync::OnceLock;

use dirt_capsule::*;
use dirt_core::prelude::*;

fn mode() -> String {
    std::env::var("CAPSULE_VS_ROD_MODE").unwrap_or_else(|_| "capsule".to_string())
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(DemAtomPlugin)
        .add_plugins(DemAtomInsertPlugin);

    if mode() == "rod" {
        app.add_plugins(HertzMindlinContactPlugin)
            .add_plugins(RotationalDynamicsPlugin)
            .add_plugins(ClumpPlugin);
        app.add_update_system(measure_rod, ParticleSimScheduleSet::PostFinalIntegration);
    } else {
        // Capsules are free atoms: velocity-Verlet handles translation, CapsulePlugin
        // handles the full-tensor orientation, CapsuleContactPlugin the contact.
        app.add_plugins(dirt_core::soil_verlet::VelocityVerletPlugin::new())
            .add_plugins(CapsuleContactPlugin)
            .add_plugins(CapsulePlugin);
        app.add_update_system(
            measure_capsule,
            ParticleSimScheduleSet::PostFinalIntegration,
        );
    }

    app.start();
}

fn write_row(step: usize, time: f64, t_trans: f64, t_rot: f64) {
    let dir = "examples/bench_capsule_vs_rod_haff/data".to_string();
    let path = format!("{}/cooling_{}.csv", dir, mode());
    static INIT: OnceLock<()> = OnceLock::new();
    let _ = INIT.get_or_init(|| {
        fs::create_dir_all(&dir).ok();
        let mut f = fs::File::create(&path).expect("create cooling csv");
        writeln!(f, "step,time,T_trans,T_rot,T_total").unwrap();
    });
    let mut f = fs::OpenOptions::new()
        .append(true)
        .open(&path)
        .expect("open cooling csv");
    writeln!(
        f,
        "{},{:.8e},{:.8e},{:.8e},{:.8e}",
        step,
        time,
        t_trans,
        t_rot,
        t_trans + t_rot
    )
    .unwrap();
}

/// Translational + rotational granular temperature of the capsule gas.
fn measure_capsule(
    atoms: Res<Atom>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    particles: ParticlesWith<'_, (Read<dirt_core::dirt_atom::DemAtom>, Read<CapsuleAtom>)>,
) {
    let step = run_state.total_cycle;
    if step % 20000 != 0 {
        return;
    }
    let nlocal = atoms.nlocal() as usize;

    // Mean (mass-weighted) velocity first.
    let (mut m_tot, mut mv) = (0.0_f64, [0.0_f64; 3]);
    particles.with(|(_dem, caps)| {
        for i in 0..nlocal.min(caps.is_capsule.len()) {
            if caps.is_capsule[i] == 0.0 {
                continue;
            }
            let m = atoms.mass[i] as f64;
            m_tot += m;
            for d in 0..3 {
                mv[d] += m * atoms.vel[i][d] as f64;
            }
        }
    });
    let m_tot = comm.all_reduce_sum_f64(m_tot);
    let mv = [
        comm.all_reduce_sum_f64(mv[0]),
        comm.all_reduce_sum_f64(mv[1]),
        comm.all_reduce_sum_f64(mv[2]),
    ];
    let v_mean = if m_tot > 0.0 {
        [mv[0] / m_tot, mv[1] / m_tot, mv[2] / m_tot]
    } else {
        [0.0; 3]
    };

    let (mut ke_t, mut ke_r) = (0.0_f64, 0.0_f64);
    particles.with(|(dem, caps)| {
        for i in 0..nlocal.min(caps.is_capsule.len()) {
            if caps.is_capsule[i] == 0.0 {
                continue;
            }
            let m = atoms.mass[i] as f64;
            let mut dv2 = 0.0;
            for d in 0..3 {
                let dv = atoms.vel[i][d] as f64 - v_mean[d];
                dv2 += dv * dv;
            }
            ke_t += m * dv2;
            // Rotational: ω·I·ω in the principal frame (principal_axes = identity).
            let q = dem.quaternion[i];
            let qc = [q[0], -q[1], -q[2], -q[3]];
            let wb = dirt_capsule::quat_rotate(qc, dem.omega[i]);
            let im = caps.principal_moments[i];
            ke_r += im[0] * wb[0] * wb[0] + im[1] * wb[1] * wb[1] + im[2] * wb[2] * wb[2];
        }
    });
    let ke_t = comm.all_reduce_sum_f64(ke_t);
    let ke_r = comm.all_reduce_sum_f64(ke_r);
    if comm.rank() != 0 || m_tot <= 0.0 {
        return;
    }
    write_row(
        step,
        step as f64 * atoms.dt,
        ke_t / (3.0 * m_tot),
        ke_r / (3.0 * m_tot),
    );
}

/// Translational + rotational granular temperature of the clump gas.
fn measure_rod(
    atoms: Res<Atom>,
    bodies: Res<MultisphereBodyStore>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
) {
    let step = run_state.total_cycle;
    if step % 20000 != 0 || bodies.bodies.is_empty() {
        return;
    }
    let (mut m_tot, mut mv) = (0.0_f64, [0.0_f64; 3]);
    for b in &bodies.bodies {
        m_tot += b.total_mass;
        for d in 0..3 {
            mv[d] += b.total_mass * b.com_vel[d];
        }
    }
    let m_tot = comm.all_reduce_sum_f64(m_tot);
    let mv = [
        comm.all_reduce_sum_f64(mv[0]),
        comm.all_reduce_sum_f64(mv[1]),
        comm.all_reduce_sum_f64(mv[2]),
    ];
    let v_mean = if m_tot > 0.0 {
        [mv[0] / m_tot, mv[1] / m_tot, mv[2] / m_tot]
    } else {
        [0.0; 3]
    };

    let (mut ke_t, mut ke_r) = (0.0_f64, 0.0_f64);
    for b in &bodies.bodies {
        let mut dv2 = 0.0;
        for d in 0..3 {
            let dv = b.com_vel[d] - v_mean[d];
            dv2 += dv * dv;
        }
        ke_t += b.total_mass * dv2;
        let q = b.quaternion;
        let qc = [q[0], -q[1], -q[2], -q[3]];
        let ob = dirt_core::dirt_clump::quat_rotate(qc, b.omega);
        let qp = b.principal_axes;
        let qpc = [qp[0], -qp[1], -qp[2], -qp[3]];
        let op = dirt_core::dirt_clump::quat_rotate(qpc, ob);
        let im = b.principal_moments;
        ke_r += im[0] * op[0] * op[0] + im[1] * op[1] * op[1] + im[2] * op[2] * op[2];
    }
    let ke_t = comm.all_reduce_sum_f64(ke_t);
    let ke_r = comm.all_reduce_sum_f64(ke_r);
    if comm.rank() != 0 || m_tot <= 0.0 {
        return;
    }
    write_row(
        step,
        step as f64 * atoms.dt,
        ke_t / (3.0 * m_tot),
        ke_r / (3.0 * m_tot),
    );
}
