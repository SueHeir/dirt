//! coherence_validate — Phase 4 of coherence_plan.md. Proves the scheduler-mediated
//! host↔device coherence layer: a CPU system added to a GPU-resident config now
//! syncs transparently instead of silently dropping, and the lazy device→host pull
//! reproduces the eager per-window download bit-for-bit.
//!
//! Four resident-GPU runs of the same drop (window W, T ticks):
//!   1. eager_base   — coherence OFF, no CPU systems
//!   2. coh_base     — coherence ON, a CPU reader each tick (no writer)
//!   3. coh_writer   — coherence ON, a CPU velocity-damp writer + reader
//!   4. eager_writer — coherence OFF, the same writer
//!
//! Assertions:
//!   A. coh_base ≈ eager_base       (lazy pull == eager download, bit-faithful)
//!   B. coh_base forced T syncs      (each host read attributed + counted)
//!   C. coh_writer ≠ coh_base        (the CPU writer is RESPECTED under coherence)
//!   D. eager_writer == eager_base   (the CPU writer is SILENTLY DROPPED without it)
//! C vs D is the headline: the identical CPU writer changes the trajectory under
//! coherence but is dropped by the eager resident path.
//!
//! Run: cargo run -p dirt_granular --example coherence_validate \
//!        --no-default-features --features precision-double,gpu_coherence --release

use std::any::TypeId;

use grass_app::prelude::*;
use grass_scheduler::prelude::*;

use dirt_atom::{DemAtom, MaterialTable};
use dirt_granular::{gpu_granular_resident_step, resident_coherence_registry, ResidentGpu};
use dirt_gpu::{Boundary, GpuContext, Plane};
use soil_core::{Atom, AtomDataRegistry, ParticleSimScheduleSet, Real};

const SIDE: usize = 6; // 6^3 = 216 grains
const R: f32 = 0.05;
const DENSITY: f64 = 2500.0;
const WINDOW: usize = 50;
const TICKS: usize = 20; // 1000 device steps total
const GRAVITY: [f32; 3] = [0.0, 0.0, -9.81];
const DAMP: f64 = 0.9; // CPU writer: scale velocity each tick

struct Scene {
    pos: Vec<[f32; 3]>,
    n: usize,
    mass: f64,
    inv_inertia: f32,
    boundary: Boundary,
    dt: f32,
}

fn make_mt() -> MaterialTable {
    let mut mt = MaterialTable::new();
    mt.add_material("soft", 1.0e6, 0.3, 0.9, 0.4, 0.0, 0.0);
    mt.build_pair_tables();
    mt
}

fn build_scene() -> Scene {
    let spacing = 2.05 * R;
    let mut pos = Vec::new();
    for ix in 0..SIDE {
        for iy in 0..SIDE {
            for iz in 0..SIDE {
                let f = (ix + iy * SIDE + iz * SIDE * SIDE) as f64;
                pos.push([
                    1.5 * R + ix as f32 * spacing + (0.13 * f).sin() as f32 * 0.03 * R,
                    1.5 * R + iy as f32 * spacing + (0.27 * f).cos() as f32 * 0.03 * R,
                    1.5 * R + iz as f32 * spacing,
                ]);
            }
        }
    }
    let n = pos.len();
    let mass = DENSITY * 4.0 / 3.0 * std::f64::consts::PI * (R as f64).powi(3);
    let inv_inertia = (1.0 / (0.4 * mass * (R as f64).powi(2))) as f32;

    let mt = make_mt();
    let e_eff = mt.e_eff_ij[0][0] as f32;
    let delta = 0.05 * R;
    let k_n = (4.0 / 3.0) * e_eff * (delta * R).sqrt();
    let tc = 2.0 * std::f32::consts::PI * (mass as f32 / k_n).sqrt();
    let dt = tc / 40.0;

    let box_w = SIDE as f32 * spacing + R;
    let mut boundary = Boundary::new();
    boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]));
    boundary.push(Plane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]));
    boundary.push(Plane::new([box_w, 0.0, 0.0], [-1.0, 0.0, 0.0]));
    boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]));
    boundary.push(Plane::new([0.0, box_w, 0.0], [0.0, -1.0, 0.0]));

    Scene { pos, n, mass, inv_inertia, boundary, dt }
}

fn make_app_state(s: &Scene) -> (Atom, AtomDataRegistry) {
    let mut atom = Atom::new();
    atom.dt = s.dt as f64;
    for (i, p) in s.pos.iter().enumerate() {
        atom.push_test_atom(i as u32, [p[0] as f64, p[1] as f64, p[2] as f64], R as f64, s.mass);
    }
    atom.nlocal = s.n as u32;
    atom.natoms = s.n as u64;

    let mut dem = DemAtom::new();
    for _ in 0..s.n {
        dem.radius.push(R as f64);
        dem.density.push(DENSITY);
        dem.inv_inertia.push(s.inv_inertia as f64);
        dem.quaternion.push([1.0, 0.0, 0.0, 0.0]);
        dem.omega.push([0.0; 3]);
        dem.ang_mom.push([0.0; 3]);
        dem.torque.push([0.0; 3]);
        dem.body_id.push(0.0);
    }
    let mut registry = AtomDataRegistry::new();
    registry.register(dem);
    (atom, registry)
}

// ── CPU systems added to the resident config ─────────────────────────────────

/// Reader (PostFinalIntegration): reads `Atom` after the resident step → forces a
/// device→host pull under coherence, and keeps the host mirror current for readout.
struct Probe {
    pos0_z: f64,
}
fn reader_system(atoms: Res<Atom>, mut probe: ResMut<Probe>) {
    probe.pos0_z = atoms.pos[0][2] as f64;
}

/// Writer (InitialIntegration): damps every particle's velocity. Honestly declared
/// `ResMut<Atom>`, so the scheduler treats it as a write to the mirror trigger.
struct Damp(f64);
fn writer_system(mut atoms: ResMut<Atom>, damp: Res<Damp>) {
    let n = atoms.nlocal as usize;
    let f = damp.0 as Real;
    for i in 0..n {
        atoms.vel[i][0] *= f;
        atoms.vel[i][1] *= f;
        atoms.vel[i][2] *= f;
    }
}

struct RunOpts {
    coherence: bool,
    writer: bool,
    reader: bool,
}

/// Run the scene through the resident step `TICKS` times. Returns (final positions,
/// forced device→host sync count).
fn run(s: &Scene, ctx: GpuContext, opts: RunOpts) -> (Vec<[f64; 3]>, u64) {
    let (atom, registry) = make_app_state(s);
    let mut app = App::new();
    app.add_resource(atom);
    app.add_resource(registry);
    app.add_resource(make_mt());
    app.add_resource(ResidentGpu::new(Some(ctx), WINDOW, GRAVITY, s.boundary.clone()));
    if opts.coherence {
        app.add_resource(resident_coherence_registry());
    }
    app.add_update_system(gpu_granular_resident_step, ParticleSimScheduleSet::Force);
    if opts.writer {
        app.add_resource(Damp(DAMP));
        app.add_update_system(writer_system, ParticleSimScheduleSet::InitialIntegration);
    }
    if opts.reader {
        app.add_resource(Probe { pos0_z: 0.0 });
        app.add_update_system(reader_system, ParticleSimScheduleSet::PostFinalIntegration);
    }
    app.organize_systems();

    for _ in 0..TICKS {
        app.run();
    }

    let a = app.get_resource_ref::<Atom>().unwrap();
    let pos: Vec<[f64; 3]> =
        (0..s.n).map(|i| [a.pos[i][0] as f64, a.pos[i][1] as f64, a.pos[i][2] as f64]).collect();
    drop(a);
    let syncs = if opts.coherence {
        app.get_resource_ref::<CoherenceRegistry>().unwrap().syncs(TypeId::of::<Atom>())
    } else {
        0
    };
    (pos, syncs)
}

fn max_diff(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    a.iter().zip(b).flat_map(|(x, y)| (0..3).map(move |d| (x[d] - y[d]).abs())).fold(0.0, f64::max)
}

fn main() {
    let Some(ctx) = GpuContext::new() else {
        eprintln!("no GPU adapter; skipping coherence_validate");
        return;
    };
    let s = build_scene();
    println!(
        "coherence_validate: n={} window={} ticks={} dt={:.3e} adapter={}",
        s.n, WINDOW, TICKS, s.dt, ctx.adapter_info
    );

    let (eager_base, _) = run(&s, ctx.clone(), RunOpts { coherence: false, writer: false, reader: false });
    let (coh_base, coh_syncs) = run(&s, ctx.clone(), RunOpts { coherence: true, writer: false, reader: true });
    let (coh_writer, _) = run(&s, ctx.clone(), RunOpts { coherence: true, writer: true, reader: true });
    let (eager_writer, _) = run(&s, ctx.clone(), RunOpts { coherence: false, writer: true, reader: false });

    let a = max_diff(&coh_base, &eager_base); // lazy pull vs eager download
    let c = max_diff(&coh_writer, &coh_base); // writer respected under coherence
    let d = max_diff(&eager_writer, &eager_base); // writer dropped without coherence

    println!("\n  A  max|coh_base - eager_base|     : {a:.3e}   (lazy pull == eager; want < 1e-6)");
    println!("  B  forced device->host syncs      : {coh_syncs}        (one per reader tick; want {TICKS})");
    println!("  C  max|coh_writer - coh_base|      : {c:.3e}   (CPU writer RESPECTED; want > 1e-4)");
    println!("  D  max|eager_writer - eager_base|  : {d:.3e}   (CPU writer DROPPED;  want < 1e-9)");

    let pass_a = a < 1e-6;
    let pass_b = coh_syncs == TICKS as u64;
    let pass_c = c > 1e-4;
    let pass_d = d < 1e-9;
    let pass = pass_a && pass_b && pass_c && pass_d;

    println!(
        "\n  => {} A={} B={} C={} D={}\n     The same CPU velocity-damp writer changes the trajectory under coherence\n     (C) but is silently dropped by the eager resident path (D) — coherence makes\n     host systems safe to add, and the lazy device->host pull is bit-faithful (A).",
        if pass { "PASS:" } else { "FAIL:" },
        pass_a, pass_b, pass_c, pass_d
    );
    assert!(pass, "coherence_validate failed: A={a:.3e} B={coh_syncs} C={c:.3e} D={d:.3e}");
}
