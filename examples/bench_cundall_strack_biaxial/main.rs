//! bench_cundall_strack_biaxial — compact dense-assembly biaxial compression check.
//!
//! The example runs a closed, walled dense granular assembly under a moving
//! loading platen. It records wall-reaction resultants, internal stress,
//! density, coordination, and contact-fabric anisotropy.

use dirt_core::dirt_atom::DemAtom;
use dirt_core::prelude::*;
use dirt_core::soil_core::Neighbor;
use std::fs;
use std::io::Write as IoWrite;

const RECORD_INTERVAL: usize = 100;
// The initial random insertion is intentionally loose.  Record only after
// the prescribed walls have compacted it into the dense late-loading regime;
// otherwise a handful of transient contacts could be mislabeled "biaxial".
const BIAXIAL_START_STEP: usize = 13_000;

#[derive(Clone, Copy)]
struct InitialBox {
    lx: f64,
    ly: f64,
    volume: f64,
    loading_top_z: f64,
}

struct BiaxialRecorder {
    initial: Option<InitialBox>,
    written_header: bool,
    reaction_sum: [f64; 3],
    reaction_samples: usize,
}

impl BiaxialRecorder {
    fn new() -> Self {
        Self {
            initial: None,
            written_header: false,
            reaction_sum: [0.0; 3],
            reaction_samples: 0,
        }
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin)
        .add_plugins(FixesPlugin)
        .add_plugins(DeformPlugin);

    app.add_resource(BiaxialRecorder::new());
    app.add_update_system(record_biaxial, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

#[allow(clippy::too_many_arguments)]
fn record_biaxial(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    domain: Res<Domain>,
    neighbor: Res<Neighbor>,
    virial: Option<Res<VirialStress>>,
    walls: Res<Walls>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
    mut recorder: ResMut<BiaxialRecorder>,
) {
    let step = run_state.total_cycle;
    if step < BIAXIAL_START_STEP || domain.volume <= 0.0 {
        return;
    }

    // `Domain` owns the fixed neighbour/decomposition bounds.  It is not the
    // specimen cell once plane walls move, so using it here silently reported
    // zero lateral strain and an incorrect volume.  The measured cell is
    // bounded by the live, named specimen walls.
    let wall_position_x = |name: &str| {
        walls
            .planes
            .iter()
            .find(|wall| wall.name.as_deref() == Some(name))
            .map(|wall| wall.point_x)
            .expect("biaxial benchmark requires named x walls")
    };
    let x_low = wall_position_x("x_low");
    let x_high = wall_position_x("x_high");
    let lx = x_high - x_low;
    if lx <= 0.0 {
        return;
    }
    let ly = (domain.boundaries_high[1] - domain.boundaries_low[1]) as f64;
    let loading_top_z = walls
        .planes
        .iter()
        .find(|wall| wall.name.as_deref() == Some("loading_top"))
        .map(|wall| wall.point_z)
        .expect("biaxial benchmark requires loading_top wall");
    let wall_force = |name: &str| {
        walls
            .planes
            .iter()
            .find(|wall| wall.name.as_deref() == Some(name))
            .map(|wall| wall.force_accumulator)
            .expect("biaxial benchmark requires named walls")
    };
    // Cundall--Strack's source observable is the in-plane horizontal resultant
    // of a two-dimensional specimen.  This specimen is one grain thick in y;
    // the y walls only constrain that suppressed direction and must *not* be
    // folded into F_H.  Doing so was the earlier 3-D proxy error.
    recorder.reaction_sum[0] += 0.5 * (wall_force("x_low") + wall_force("x_high"));
    recorder.reaction_sum[2] += wall_force("loading_top");
    recorder.reaction_samples += 1;
    if step % RECORD_INTERVAL != 0 {
        return;
    }
    let reaction_samples = recorder.reaction_samples as f64;
    let f_h = recorder.reaction_sum[0] / reaction_samples;
    let f_v = recorder.reaction_sum[2] / reaction_samples;
    // Both resultants are integrated force reactions on the source-equivalent
    // in-plane walls.  The common one-particle out-of-plane depth cancels.
    let wall_force_ratio = if f_v > 0.0 { f_h / f_v } else { 0.0 };
    recorder.reaction_sum = [0.0; 3];
    recorder.reaction_samples = 0;
    let cell_volume = lx * ly * loading_top_z;
    let initial = *recorder.initial.get_or_insert(InitialBox {
        lx,
        ly,
        volume: cell_volume,
        loading_top_z,
    });

    let nlocal = atoms.nlocal as usize;
    let dem = registry.expect::<DemAtom>("record_biaxial");
    let mut kin = [0.0f64; 6];
    let mut solid_volume = 0.0f64;
    for i in 0..nlocal {
        let m = atoms.mass[i] as f64;
        let vx = atoms.vel[i][0] as f64;
        let vy = atoms.vel[i][1] as f64;
        let vz = atoms.vel[i][2] as f64;
        kin[0] += m * vx * vx;
        kin[1] += m * vy * vy;
        kin[2] += m * vz * vz;
        kin[3] += m * vx * vy;
        kin[4] += m * vx * vz;
        kin[5] += m * vy * vz;

        let r = dem.radius[i] as f64;
        solid_volume += (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
    }

    let vir = match virial.as_ref() {
        Some(v) => [v.xx, v.yy, v.zz, v.xy, v.xz, v.yz],
        None => [0.0; 6],
    };

    let (contacts, coord_sum, fabric) = contact_metrics(&atoms, &neighbor, &dem);

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
        solid_volume,
        contacts,
        coord_sum,
        fabric[0],
        fabric[1],
        fabric[2],
        fabric[3],
        fabric[4],
        fabric[5],
    ];
    for a in acc.iter_mut() {
        *a = comm.all_reduce_sum_f64(*a);
    }

    if comm.rank() != 0 {
        return;
    }

    let mut sig = [0.0f64; 6];
    for k in 0..6 {
        sig[k] = (acc[k] - acc[6 + k]) / cell_volume;
    }
    let p = (sig[0] + sig[1] + sig[2]) / 3.0;
    let q = 0.5 * (sig[0] + sig[1]) - sig[2];
    let stress_ratio = if p.abs() > 0.0 { q / p } else { 0.0 };
    let lateral_axial_stress_ratio = if sig[2].abs() > 0.0 {
        0.5 * (sig[0] + sig[1]) / sig[2]
    } else {
        0.0
    };
    let axial_strain = if initial.loading_top_z > 0.0 {
        (initial.loading_top_z - loading_top_z) / initial.loading_top_z
    } else {
        0.0
    };
    let lateral_strain =
        0.5 * (((lx - initial.lx) / initial.lx) + ((ly - initial.ly) / initial.ly));
    let volumetric_strain = if initial.volume > 0.0 {
        (initial.volume - cell_volume) / initial.volume
    } else {
        0.0
    };
    let contacts = acc[13];
    let coordination = if atoms.natoms > 0 {
        acc[14] / atoms.natoms as f64
    } else {
        0.0
    };
    let fxx = if contacts > 0.0 {
        acc[15] / contacts
    } else {
        0.0
    };
    let fyy = if contacts > 0.0 {
        acc[16] / contacts
    } else {
        0.0
    };
    let fzz = if contacts > 0.0 {
        acc[17] / contacts
    } else {
        0.0
    };
    let fxy = if contacts > 0.0 {
        acc[18] / contacts
    } else {
        0.0
    };
    let fxz = if contacts > 0.0 {
        acc[19] / contacts
    } else {
        0.0
    };
    let fyz = if contacts > 0.0 {
        acc[20] / contacts
    } else {
        0.0
    };
    let fabric_anisotropy = fzz - 0.5 * (fxx + fyy);
    let phi = acc[12] / cell_volume;
    let time = step as f64 * atoms.dt;

    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| "examples/bench_cundall_strack_biaxial".to_string());
    let data_dir = format!("{}/data", out_dir);
    fs::create_dir_all(&data_dir).ok();
    let path = format!("{}/biaxial_results.csv", data_dir);

    if !recorder.written_header {
        let mut f = fs::File::create(&path).expect("cannot create biaxial_results.csv");
        writeln!(
            f,
            "step,time,axial_strain,lateral_strain,volumetric_strain,phi,sxx,syy,szz,sxy,sxz,syz,p,q,stress_ratio,lateral_axial_stress_ratio,f_h_mean,f_v_mean,wall_force_ratio,contacts,coordination,fxx,fyy,fzz,fxy,fxz,fyz,fabric_anisotropy"
        )
        .unwrap();
        recorder.written_header = true;
    }

    let mut f = fs::OpenOptions::new()
        .append(true)
        .open(&path)
        .expect("cannot open biaxial_results.csv");
    writeln!(
        f,
        "{},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e}",
        step, time, axial_strain, lateral_strain, volumetric_strain, phi, sig[0], sig[1], sig[2],
        sig[3], sig[4], sig[5], p, q, stress_ratio, lateral_axial_stress_ratio, f_h, f_v, wall_force_ratio, contacts,
        coordination, fxx, fyy, fzz, fxy, fxz, fyz, fabric_anisotropy
    )
    .unwrap();
}

fn contact_metrics(atoms: &Atom, neighbor: &Neighbor, dem: &DemAtom) -> (f64, f64, [f64; 6]) {
    let nlocal = atoms.nlocal as usize;
    let mut contacts = 0.0;
    let mut coord_sum = 0.0;
    let mut fabric = [0.0f64; 6];
    let weight = if neighbor.newton { 1.0 } else { 0.5 };

    for (i, j) in neighbor.pairs(nlocal) {
        let dx = atoms.pos[j][0] as f64 - atoms.pos[i][0] as f64;
        let dy = atoms.pos[j][1] as f64 - atoms.pos[i][1] as f64;
        let dz = atoms.pos[j][2] as f64 - atoms.pos[i][2] as f64;
        let dist_sq = dx * dx + dy * dy + dz * dz;
        let sum_r = dem.radius[i] + dem.radius[j];
        if dist_sq <= 0.0 || dist_sq >= sum_r * sum_r {
            continue;
        }
        let inv_dist = 1.0 / dist_sq.sqrt();
        let nx = dx * inv_dist;
        let ny = dy * inv_dist;
        let nz = dz * inv_dist;
        contacts += weight;
        coord_sum += 2.0 * weight;
        fabric[0] += nx * nx * weight;
        fabric[1] += ny * ny * weight;
        fabric[2] += nz * nz * weight;
        fabric[3] += nx * ny * weight;
        fabric[4] += nx * nz * weight;
        fabric[5] += ny * nz * weight;
    }

    (contacts, coord_sum, fabric)
}
