//! bench_granular_conductivity — granular-temperature conductivity from a
//! gas-free vibro-fluidized bed.
//!
//! A bed of glass beads sits under gravity on an **oscillating base wall**
//! (periodic in x,z). The vibrating base injects fluctuation (granular-
//! temperature) energy at the bottom; it conducts upward and is removed by
//! inelastic dissipation, producing steady **Φ(y)** and **T(y)** profiles — the
//! canonical inhomogeneous kinetic-theory benchmark, and a gas-free fluidized
//! bed. From the steady profiles `sweep.py` extracts the KT **conductivity κ(Φ)**
//! (the one transport coefficient the homogeneous LEBC rig can't see) via the
//! conduction–dissipation balance. Turning the base static (`amplitude=0`) makes
//! the same rig a **de-fluidization** test: T decays, Φ rises, the bed consolidates.
//!
//! The recorder streams horizontal-slab profiles each output interval:
//! per y-bin the solid fraction Φ(y), the granular temperature
//! `T(y) = ⅓⟨|v − v̄_bin|²⟩` (per-bin mean removed, so the coherent vibration is
//! not counted), and the kinetic fluctuation-energy flux
//! `q_y(y) = (1/V_bin) Σ ½m|v'|² v'_y` to `data/conductivity_profiles.csv`.
//!
//! ```bash
//! cargo run --release --example bench_granular_conductivity --no-default-features -- examples/bench_granular_conductivity/config.toml
//! ```

use std::fs;
use std::io::Write as IoWrite;
use std::sync::Once;

use dirt_core::prelude::*;
use dirt_core::dirt_atom::DemAtom;

const NBINS: usize = 40;

fn main() {
    let mut app = App::new();
    app.add_plugins(CorePlugins)
        .add_plugins(GranularDefaultPlugins)
        .add_plugins(GravityPlugin)
        .add_plugins(WallPlugin);

    app.add_update_system(record_profile, ParticleSimScheduleSet::PostFinalIntegration);
    app.start();
}

/// Stream per-y-bin Φ(y), T(y) (per-bin streaming velocity removed), and the
/// kinetic fluctuation-energy flux q_y(y).
fn record_profile(
    atoms: Res<Atom>,
    registry: Res<AtomDataRegistry>,
    domain: Res<Domain>,
    run_state: Res<RunState>,
    comm: Res<CommResource>,
    input: Res<Input>,
) {
    let step = run_state.total_cycle;
    if step % 5000 != 0 {
        return;
    }
    let nlocal = atoms.nlocal as usize;
    let ylo = domain.boundaries_low[1];
    let yhi = domain.boundaries_high[1];
    let dy = (yhi - ylo) / NBINS as f64;
    if dy <= 0.0 {
        return;
    }
    let bin_vol = domain.size[0] * domain.size[2] * dy;
    let dem = registry.get::<DemAtom>();

    let bin_of = |y: f64| -> usize {
        (((y - ylo) / dy).floor() as i64).clamp(0, NBINS as i64 - 1) as usize
    };

    // Pass 1: per-bin mass, momentum, solid volume.
    let mut m_sum = [0.0f64; NBINS];
    let mut mv = [[0.0f64; 3]; NBINS];
    let mut vol_solid = [0.0f64; NBINS];
    for i in 0..nlocal {
        let b = bin_of(atoms.pos[i][1]);
        let m = atoms.mass[i];
        m_sum[b] += m;
        for d in 0..3 {
            mv[b][d] += m * atoms.vel[i][d];
        }
        if let Some(ref dem) = dem {
            let r = dem.radius[i];
            vol_solid[b] += (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
        }
    }

    // Pass 2: per-bin fluctuation KE and kinetic heat flux, about the per-bin mean.
    let mut ke = [0.0f64; NBINS];
    let mut qy = [0.0f64; NBINS];
    for i in 0..nlocal {
        let b = bin_of(atoms.pos[i][1]);
        if m_sum[b] <= 0.0 {
            continue;
        }
        let vbar = [mv[b][0] / m_sum[b], mv[b][1] / m_sum[b], mv[b][2] / m_sum[b]];
        let m = atoms.mass[i];
        let vp = [
            atoms.vel[i][0] - vbar[0],
            atoms.vel[i][1] - vbar[1],
            atoms.vel[i][2] - vbar[2],
        ];
        let v2 = vp[0] * vp[0] + vp[1] * vp[1] + vp[2] * vp[2];
        ke[b] += m * v2;
        qy[b] += 0.5 * m * v2 * vp[1];
    }

    // Reduce across ranks (no-op on a single process).
    let reduce = |arr: &mut [f64; NBINS]| {
        for v in arr.iter_mut() {
            *v = comm.all_reduce_sum_f64(*v);
        }
    };
    reduce(&mut m_sum);
    reduce(&mut vol_solid);
    reduce(&mut ke);
    reduce(&mut qy);

    if comm.rank() != 0 {
        return;
    }

    let time = step as f64 * atoms.dt;
    let out_dir = input
        .output_dir
        .clone()
        .unwrap_or_else(|| "examples/bench_granular_conductivity/data".to_string());
    let path = format!("{}/conductivity_profiles.csv", out_dir);
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        fs::create_dir_all(&out_dir).ok();
        let mut f = fs::File::create(&path).expect("cannot create conductivity_profiles.csv");
        writeln!(f, "step,time,y,phi,T,qy").unwrap();
    });
    let mut f = fs::OpenOptions::new().append(true).open(&path).expect("cannot open profiles csv");
    for b in 0..NBINS {
        let yc = ylo + (b as f64 + 0.5) * dy;
        let phi = vol_solid[b] / bin_vol;
        let t = if m_sum[b] > 0.0 { ke[b] / (3.0 * m_sum[b]) } else { 0.0 };
        let qy_density = qy[b] / bin_vol;
        writeln!(f, "{},{:.8e},{:.8e},{:.8e},{:.8e},{:.8e}", step, time, yc, phi, t, qy_density).unwrap();
    }
}
