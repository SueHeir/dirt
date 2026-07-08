//! Rod-like glued-sphere Lees-Edwards shear-flow aspect-ratio benchmark.
//!
//! This is a small DIRT counterpart to Guo et al. (JFM 713, 2012): frictionless
//! glued-sphere rods in homogeneous simple shear. The recorder writes
//! body-level stress, solid fraction, apparent friction, and alignment metrics to
//! `rod_shear_results.csv`; `sweep.py` compares their aspect-ratio trends with
//! the published DEM trends.

use std::fs;
use std::io::Write as IoWrite;
use std::sync::Once;

use dirt_core::prelude::*;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(DemAtomPlugin)
        .add_plugins(DemAtomInsertPlugin)
        .add_plugins(HertzMindlinContactPlugin)
        .add_plugins(RotationalDynamicsPlugin)
        .add_plugins(ClumpPlugin)
        .add_plugins(DeformPlugin);

    app.add_update_system(
        record_rod_shear,
        ParticleSimScheduleSet::PostFinalIntegration,
    );
    app.start();
}

fn record_rod_shear(
    atoms: Res<Atom>,
    bodies: Res<MultisphereBodyStore>,
    domain: Res<Domain>,
    virial: Option<Res<VirialStress>>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
) {
    let step = run_state.total_cycle;
    if step % 1000 != 0 {
        return;
    }

    let vol = domain.volume;
    if vol <= 0.0 || bodies.bodies.is_empty() {
        return;
    }

    let ly = domain.size[1];
    let gdot = if ly > 0.0 {
        domain.boundary_vel[0] / ly
    } else {
        0.0
    };
    let yc = 0.5 * (domain.boundaries_low[1] + domain.boundaries_high[1]);

    let mut kin = [0.0f64; 6];
    let mut ke_fluct = 0.0f64;
    let mut m_total = 0.0f64;
    let mut vol_solid = 0.0f64;
    let mut orient_x2 = 0.0f64;
    let mut align_abs_x = 0.0f64;
    let mut nbodies = 0.0f64;

    for body in &bodies.bodies {
        let m = body.total_mass;
        m_total += m;
        let vx = body.com_vel[0] - gdot * (body.com_pos[1] - yc);
        let vy = body.com_vel[1];
        let vz = body.com_vel[2];
        kin[0] += m * vx * vx;
        kin[1] += m * vy * vy;
        kin[2] += m * vz * vz;
        kin[3] += m * vx * vy;
        kin[4] += m * vx * vz;
        kin[5] += m * vy * vz;
        ke_fluct += m * (vx * vx + vy * vy + vz * vz);

        for r in &body.sub_sphere_radii {
            vol_solid += (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
        }

        let axis_body = if body.body_offsets.len() >= 2 {
            let first = body.body_offsets[0];
            let last = body.body_offsets[body.body_offsets.len() - 1];
            let v = [last[0] - first[0], last[1] - first[1], last[2] - first[2]];
            let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if n > 0.0 {
                [v[0] / n, v[1] / n, v[2] / n]
            } else {
                [1.0, 0.0, 0.0]
            }
        } else {
            [1.0, 0.0, 0.0]
        };
        let axis = dirt_clump::quat_rotate(body.quaternion, axis_body);
        orient_x2 += axis[0] * axis[0];
        align_abs_x += axis[0].abs();
        nbodies += 1.0;
    }

    let vir = match virial.as_ref() {
        Some(v) => [v.xx, v.yy, v.zz, v.xy, v.xz, v.yz],
        None => [0.0; 6],
    };

    let mut acc = [
        kin[0],
        kin[1],
        kin[2],
        kin[3],
        kin[4],
        kin[5],
        vir[0],
        vir[1],
        vir[2],
        vir[3],
        vir[4],
        vir[5],
        ke_fluct,
        m_total,
        vol_solid,
        orient_x2,
        align_abs_x,
        nbodies,
    ];
    for a in acc.iter_mut() {
        *a = comm.all_reduce_sum_f64(*a);
    }
    let kin = [acc[0], acc[1], acc[2], acc[3], acc[4], acc[5]];
    let vir = [acc[6], acc[7], acc[8], acc[9], acc[10], acc[11]];
    let ke_fluct = acc[12];
    let m_total = acc[13];
    let vol_solid = acc[14];
    let orient_x2 = acc[15];
    let align_abs_x = acc[16];
    let nbodies = acc[17];

    let mut sig = [0.0f64; 6];
    for k in 0..6 {
        sig[k] = (kin[k] - vir[k]) / vol;
    }
    let p = (sig[0] + sig[1] + sig[2]) / 3.0;
    let mu = if p.abs() > 0.0 {
        sig[3].abs() / p.abs()
    } else {
        0.0
    };
    let t_gran = if m_total > 0.0 {
        ke_fluct / (3.0 * m_total)
    } else {
        0.0
    };
    let phi = vol_solid / vol;
    let order_x = if nbodies > 0.0 {
        0.5 * (3.0 * orient_x2 / nbodies - 1.0)
    } else {
        0.0
    };
    let align_x = if nbodies > 0.0 {
        align_abs_x / nbodies
    } else {
        0.0
    };

    if comm.rank() != 0 {
        return;
    }

    let time = step as f64 * atoms.dt;
    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| "examples/bench_rod_shear_aspect_ratio/data".to_string());
    let path = format!("{}/rod_shear_results.csv", out_dir);

    static INIT: Once = Once::new();
    INIT.call_once(|| {
        fs::create_dir_all(&out_dir).ok();
        let mut f = fs::File::create(&path).expect("cannot create rod_shear_results.csv");
        writeln!(
            f,
            "step,time,gdot,sxx,syy,szz,sxy,sxz,syz,p,mu,T,phi,order_x,align_x"
        )
        .unwrap();
    });

    let mut f = fs::OpenOptions::new()
        .append(true)
        .open(&path)
        .expect("cannot open rod_shear_results.csv");
    writeln!(
        f,
        "{},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e}",
        step,
        time,
        gdot,
        sig[0],
        sig[1],
        sig[2],
        sig[3],
        sig[4],
        sig[5],
        p,
        mu,
        t_gran,
        phi,
        order_x,
        align_x
    )
    .unwrap();
}
