//! Rigorous comparison of the GPU DEM (`dirt_gpu`) against dirt's REAL CPU force
//! code (`dirt_granular::hertz_mindlin_contact_force` + `dirt_wall`), plus a fair
//! full-step benchmark.
//!
//! Build/run in precision-double so the host `Atom` is f64, exactly matching
//! dirt's CPU. The GPU kernels are always f32 (data cast on upload), so f32-vs-f64
//! differences of ~1e-3 are expected and asserted.
//!
//! ```text
//! cargo run --release --example compare_cpu \
//!     -p dirt_gpu --no-default-features --features precision-double
//! ```
//!
//! Tier 1: per-evaluation force/torque correctness vs real dirt (configs a–d).
//! Tier 3: fair full-step CPU-vs-GPU benchmark for several N.

use grass_app::prelude::*;
use grass_scheduler::prelude::{Res, ResMut};

use dirt_atom::{DemAtom, MaterialTable};
use dirt_gpu::GpuContext;
use dirt_granular::contact::hertz_mindlin_contact_force;
use dirt_granular::rotational::{final_rotation, initial_rotation};
use dirt_granular::tangential::ContactHistoryStore;
use dirt_test_utils::make_material_table;
use dirt_wall::{wall_contact_force, WallPlane, WallMotion, Walls};
use soil_core::{Atom, AtomDataRegistry, Neighbor, ParticleSimScheduleSet};
use soil_verlet::{final_integration, initial_integration};
use std::time::Instant;

const DENSITY: f64 = 2500.0;

fn sphere_mass(r: f64) -> f64 {
    DENSITY * 4.0 / 3.0 * std::f64::consts::PI * r.powi(3)
}
fn sphere_inv_inertia(r: f64) -> f64 {
    let m = sphere_mass(r);
    1.0 / (0.4 * m * r * r)
}

/// A single material with rolling/twisting/cohesion/surface-energy forced to 0,
/// so dirt takes the plain Hertz + Mindlin path that the GPU implements.
fn matched_material_table() -> MaterialTable {
    let mut mt = make_material_table();
    // make_material_table already sets rolling=0 and cohesion=0, but be defensive:
    // zero every field the GPU omits, on the pre-built pair tables.
    for row in mt.rolling_friction_ij.iter_mut() { row.iter_mut().for_each(|v| *v = 0.0); }
    for row in mt.twisting_friction_ij.iter_mut() { row.iter_mut().for_each(|v| *v = 0.0); }
    for row in mt.cohesion_energy_ij.iter_mut() { row.iter_mut().for_each(|v| *v = 0.0); }
    for row in mt.surface_energy_ij.iter_mut() { row.iter_mut().for_each(|v| *v = 0.0); }
    mt
}

/// Extract the scalar GPU params from the (single-material) table.
fn gpu_scalars(mt: &MaterialTable) -> (f32, f32, f32, f32) {
    (
        mt.e_eff_ij[0][0] as f32,
        mt.beta_ij[0][0] as f32,
        mt.g_eff_ij[0][0] as f32,
        mt.friction_ij[0][0] as f32,
    )
}

/// Push one DEM atom (atom + dem in lock-step) with full state.
#[allow(clippy::too_many_arguments)]
fn push_atom(
    atom: &mut Atom,
    dem: &mut DemAtom,
    tag: u32,
    pos: [f64; 3],
    vel: [f64; 3],
    omega: [f64; 3],
    radius: f64,
) {
    let mass = sphere_mass(radius);
    atom.push_test_atom(tag, pos, radius, mass);
    let i = atom.pos.len() - 1;
    atom.vel[i] = [vel[0] as _, vel[1] as _, vel[2] as _];
    dem.radius.push(radius);
    dem.density.push(DENSITY);
    dem.inv_inertia.push(1.0 / (0.4 * mass * radius * radius));
    dem.quaternion.push([1.0, 0.0, 0.0, 0.0]);
    dem.omega.push(omega);
    dem.ang_mom.push([0.0; 3]);
    dem.torque.push([0.0; 3]);
    dem.body_id.push(0.0);
}

// ── Tier 1: per-evaluation force/torque correctness ─────────────────────────

struct EvalResult {
    force: Vec<[f64; 3]>,
    torque: Vec<[f64; 3]>,
}

/// Run dirt's REAL particle contact force ONCE (no integration) on a full
/// neighbor list. Returns per-atom force and torque.
fn cpu_particle_force_once(
    pos: &[[f64; 3]],
    vel: &[[f64; 3]],
    omega: &[[f64; 3]],
    radius: f64,
    mt: MaterialTable,
) -> EvalResult {
    let n = pos.len();
    let mut app = App::new();
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    atom.dt = 1e-6;
    for i in 0..n {
        push_atom(&mut atom, &mut dem, i as u32, pos[i], vel[i], omega[i], radius);
        hist.contacts.push(Vec::new());
    }
    atom.nlocal = n as u32;
    atom.natoms = n as u64;

    // Full neighbor list (newton=false: every ordered pair, each i sees all j).
    let mut neighbor = Neighbor::new();
    neighbor.newton = false;
    let mut offsets = vec![0u32];
    let mut indices = Vec::new();
    for i in 0..n {
        for j in 0..n {
            if i != j {
                indices.push(j as u32);
            }
        }
        offsets.push(indices.len() as u32);
    }
    neighbor.neighbor_offsets = offsets;
    neighbor.neighbor_indices = indices;

    let mut registry = AtomDataRegistry::new();
    registry.register(dem);
    registry.register(hist);

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(mt);
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let dem = registry.expect::<DemAtom>("read");
    EvalResult {
        force: (0..n).map(|i| [atom.force[i][0] as f64, atom.force[i][1] as f64, atom.force[i][2] as f64]).collect(),
        torque: (0..n).map(|i| dem.torque[i]).collect(),
    }
}

/// Run dirt's REAL wall contact force ONCE on a set of plane walls.
fn cpu_wall_force_once(
    pos: &[[f64; 3]],
    vel: &[[f64; 3]],
    omega: &[[f64; 3]],
    radius: f64,
    walls_geo: &[([f64; 3], [f64; 3])],
    mt: MaterialTable,
) -> EvalResult {
    let n = pos.len();
    let mut app = App::new();
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    atom.dt = 1e-6;
    for i in 0..n {
        push_atom(&mut atom, &mut dem, i as u32, pos[i], vel[i], omega[i], radius);
    }
    atom.nlocal = n as u32;
    atom.natoms = n as u64;

    let mut planes = Vec::new();
    let mut active = Vec::new();
    for (point, normal) in walls_geo {
        let mag = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        planes.push(WallPlane {
            point_x: point[0], point_y: point[1], point_z: point[2],
            normal_x: normal[0] / mag, normal_y: normal[1] / mag, normal_z: normal[2] / mag,
            material_index: 0,
            name: None,
            bound_x_low: f64::NEG_INFINITY, bound_x_high: f64::INFINITY,
            bound_y_low: f64::NEG_INFINITY, bound_y_high: f64::INFINITY,
            bound_z_low: f64::NEG_INFINITY, bound_z_high: f64::INFINITY,
            velocity: [0.0; 3],
            motion: WallMotion::Static,
            origin: *point,
            force_accumulator: 0.0,
            temperature: None,
        });
        active.push(true);
    }
    let walls = Walls {
        planes, active,
        cylinders: Vec::new(), cylinder_active: Vec::new(),
        spheres: Vec::new(), sphere_active: Vec::new(),
        regions: Vec::new(), region_active: Vec::new(),
        time: 0.0,
        tangential_springs: std::collections::HashMap::new(),
        rolling_springs: std::collections::HashMap::new(),
    };

    let mut registry = AtomDataRegistry::new();
    registry.register(dem);

    app.add_resource(atom);
    app.add_resource(walls);
    app.add_resource(registry);
    app.add_resource(mt);
    app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
    app.organize_systems();
    app.run();

    let atom = app.get_resource_ref::<Atom>().unwrap();
    let registry = app.get_resource_ref::<AtomDataRegistry>().unwrap();
    let dem = registry.expect::<DemAtom>("read");
    EvalResult {
        force: (0..n).map(|i| [atom.force[i][0] as f64, atom.force[i][1] as f64, atom.force[i][2] as f64]).collect(),
        torque: (0..n).map(|i| dem.torque[i]).collect(),
    }
}

/// GPU force/torque ONCE via the NEW Force-hook stack (soil GpuState + dirt
/// GranularForce + WallForce), i.e. the path that will replace the monolith.
fn gpu_force_once(
    ctx: GpuContext,
    pos: &[[f64; 3]],
    vel: &[[f64; 3]],
    omega: &[[f64; 3]],
    radius: f64,
    e_eff: f32, beta: f32, g_eff: f32, mu: f32,
    walls_geo: &[([f64; 3], [f64; 3])],
) -> EvalResult {
    use dirt_gpu::{GranularConfig, GranularForce, WallForce};
    use soil_gpu::{Boundary, GpuState, Grid as SoilGrid, Plane};

    let n = pos.len();
    let posf: Vec<[f32; 3]> = pos.iter().map(|p| [p[0] as f32, p[1] as f32, p[2] as f32]).collect();
    let velf: Vec<[f32; 3]> = vel.iter().map(|p| [p[0] as f32, p[1] as f32, p[2] as f32]).collect();
    let omf: Vec<[f32; 3]> = omega.iter().map(|p| [p[0] as f32, p[1] as f32, p[2] as f32]).collect();
    let radiusf = vec![radius as f32; n];
    let inv_mass = vec![(1.0 / sphere_mass(radius)) as f32; n];
    let inv_inertia = vec![sphere_inv_inertia(radius) as f32; n];

    let mut posf_grid = posf.clone();
    for (point, _) in walls_geo {
        posf_grid.push([point[0] as f32, point[1] as f32, point[2] as f32]);
    }
    // soil's Grid takes the literal cutoff (= sum of radii = 2*r), unlike dirt's
    // Grid which doubles r_max internally. Pass 2*r so binsize matches (bin=2r,
    // ±1 stencil = exact) and the cell count isn't 8x inflated.
    let grid = SoilGrid::from_positions(&posf_grid, (2.0 * radius) as f32);

    let mut gs = GpuState::new(ctx, n, grid.total_cells);
    gs.set_params(1e-6, [0.0, 0.0, 0.0]); // dt drives spring integration; gravity off
    gs.set_state(&posf, &velf, &inv_mass, grid);
    let om = gs.add_aux_dof();
    gs.set_aux_inv_coeff(om, &inv_inertia);
    gs.set_aux_state(om, &omf);

    let cfg = GranularConfig { e_eff, beta, g_eff, mu, dt: 1e-6 };
    gs.add_force_hook(Box::new(GranularForce::new(&gs, &grid, om, &radiusf, cfg)));
    if !walls_geo.is_empty() {
        let mut b = Boundary::new();
        for (point, normal) in walls_geo {
            b.push(Plane::new(
                [point[0] as f32, point[1] as f32, point[2] as f32],
                [normal[0] as f32, normal[1] as f32, normal[2] as f32],
            ));
        }
        gs.add_force_hook(Box::new(WallForce::new(&gs, om, &radiusf, &b, e_eff, beta, g_eff, mu, 1e-6)));
    }

    gs.eval_force_once();
    let gf = gs.download_force();
    let gt = gs.download_aux_rate(om);
    EvalResult {
        force: (0..n).map(|i| [gf[i][0] as f64, gf[i][1] as f64, gf[i][2] as f64]).collect(),
        torque: (0..n).map(|i| [gt[i][0] as f64, gt[i][1] as f64, gt[i][2] as f64]).collect(),
    }
}

/// Component-wise max relative difference, with an absolute floor on the denom.
fn max_rel(a: &[[f64; 3]], b: &[[f64; 3]], floor: f64) -> (f64, f64) {
    let mut max = 0.0f64;
    let mut max_mag = 0.0f64;
    for i in 0..a.len() {
        for c in 0..3 {
            let denom = b[i][c].abs().max(floor);
            max = max.max((a[i][c] - b[i][c]).abs() / denom);
            max_mag = max_mag.max(b[i][c].abs());
        }
    }
    (max, max_mag)
}

// ── Tier 1: per-evaluation force/torque vs REAL dirt, via the Force-hook stack
// (soil GpuState + dirt GranularForce/WallForce) ────────────────────────────
fn tier1(ctx: &GpuContext) {
    println!("\n=== Tier 1: per-evaluation force/torque vs REAL dirt (Force-hook stack) ===");
    let mt = matched_material_table();
    let (e_eff, beta, g_eff, mu) = gpu_scalars(&mt);

    let r = 0.001_f64;
    let overlap = 0.0001_f64;
    let d = 2.0 * r - overlap;
    let f_floor = 1.0;
    let t_floor = 1e-6;

    // (a) head-on normal.
    {
        let pos = vec![[0.0, 0.0, 0.0], [d, 0.0, 0.0]];
        let vel = vec![[0.0, 0.0, 0.0], [-0.1, 0.0, 0.0]];
        let omega = vec![[0.0; 3]; 2];
        let cpu = cpu_particle_force_once(&pos, &vel, &omega, r, matched_material_table());
        let gpu = gpu_force_once(ctx.clone(), &pos, &vel, &omega, r, e_eff, beta, g_eff, mu, &[]);
        let (fr, fm) = max_rel(&gpu.force, &cpu.force, f_floor);
        let (tr, _) = max_rel(&gpu.torque, &cpu.torque, t_floor);
        println!("(a) head-on overlap (normal): max rel force={fr:.2e} torque={tr:.2e}  |F|~{fm:.3e}");
        assert!(fr < 5e-3, "(a) force rel diff too large: {fr:.2e}");
    }
    // (b) overlap + sliding (tangential).
    {
        let pos = vec![[0.0, 0.0, 0.0], [d, 0.0, 0.0]];
        let vel = vec![[0.0, 0.05, 0.0], [0.0, -0.05, 0.0]];
        let omega = vec![[0.0; 3]; 2];
        let cpu = cpu_particle_force_once(&pos, &vel, &omega, r, matched_material_table());
        let gpu = gpu_force_once(ctx.clone(), &pos, &vel, &omega, r, e_eff, beta, g_eff, mu, &[]);
        let (fr, fm) = max_rel(&gpu.force, &cpu.force, f_floor);
        let (tr, tm) = max_rel(&gpu.torque, &cpu.torque, t_floor);
        println!("(b) overlap + sliding (tangential): max rel force={fr:.2e} torque={tr:.2e}  |F|~{fm:.3e} |T|~{tm:.3e}");
        assert!(fr < 5e-3, "(b) force rel diff too large: {fr:.2e}");
        assert!(tr < 5e-3, "(b) torque rel diff too large: {tr:.2e}");
    }
    // (c) overlap + spin.
    {
        let pos = vec![[0.0, 0.0, 0.0], [d, 0.0, 0.0]];
        let vel = vec![[0.0; 3]; 2];
        let omega = vec![[0.0, 0.0, 50.0], [0.0, 0.0, -30.0]];
        let cpu = cpu_particle_force_once(&pos, &vel, &omega, r, matched_material_table());
        let gpu = gpu_force_once(ctx.clone(), &pos, &vel, &omega, r, e_eff, beta, g_eff, mu, &[]);
        let (fr, fm) = max_rel(&gpu.force, &cpu.force, f_floor);
        let (tr, tm) = max_rel(&gpu.torque, &cpu.torque, t_floor);
        println!("(c) overlap + spin (omega x r + torque): max rel force={fr:.2e} torque={tr:.2e}  |F|~{fm:.3e} |T|~{tm:.3e}");
        assert!(fr < 5e-3, "(c) force rel diff too large: {fr:.2e}");
        assert!(tr < 5e-3, "(c) torque rel diff too large: {tr:.2e}");
    }
    // (d) floor wall + slide + spin.
    {
        let pos = vec![[0.0, 0.0, r - overlap]];
        let vel = vec![[0.2, -0.1, -0.3]];
        let omega = vec![[1.5, -0.8, 0.6]];
        let walls_geo = vec![([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])];
        let cpu = cpu_wall_force_once(&pos, &vel, &omega, r, &walls_geo, matched_material_table());
        let gpu = gpu_force_once(ctx.clone(), &pos, &vel, &omega, r, e_eff, beta, g_eff, mu, &walls_geo);
        let (fr, fm) = max_rel(&gpu.force, &cpu.force, f_floor);
        let (tr, tm) = max_rel(&gpu.torque, &cpu.torque, t_floor);
        println!("(d) floor wall + slide + spin: max rel force={fr:.2e} torque={tr:.2e}  |F|~{fm:.3e} |T|~{tm:.3e}");
        assert!(fr < 5e-3, "(d) wall force rel diff too large: {fr:.2e}");
        assert!(tr < 5e-3, "(d) wall torque rel diff too large: {tr:.2e}");
    }
    println!("Tier 1: Force-hook stack matches real dirt within 5e-3 (f32 vs f64).");
}

// ── Tier 3: fair full-step benchmark ────────────────────────────────────────

/// Build a packing of N particles (lightly overlapping grains) above a floor.
fn make_packing(n: usize, r: f64) -> Vec<[f64; 3]> {
    let side = (n as f64).cbrt().ceil() as usize;
    let spacing = 1.9 * r; // slight overlap so many contacts are live
    let mut pos = Vec::with_capacity(n);
    for k in 0..n {
        let (ix, iy, iz) = (k % side, (k / side) % side, k / (side * side));
        let f = k as f64;
        pos.push([
            r + ix as f64 * spacing + (0.13 * f).sin() * 0.02 * r,
            r + iy as f64 * spacing + (0.27 * f).cos() * 0.02 * r,
            r + iz as f64 * spacing, // sitting just above floor z=0
        ]);
    }
    pos
}

/// Full neighbor list (newton=false), built ONCE and reused across the timed
/// steps (the particles move only slightly, mirroring the GPU which also keeps
/// the same grid). Uses a uniform cell list so the build is O(N), not O(N^2).
fn build_full_neighbor_list(pos: &[[f64; 3]], r: f64) -> Neighbor {
    let n = pos.len();
    let cutoff = 2.0 * r;
    let cut_sq = cutoff * cutoff;
    let bin = cutoff;
    let mut lo = [f64::MAX; 3];
    let mut hi = [f64::MIN; 3];
    for p in pos {
        for d in 0..3 { lo[d] = lo[d].min(p[d]); hi[d] = hi[d].max(p[d]); }
    }
    let mut nb = [1i64; 3];
    for d in 0..3 {
        nb[d] = (((hi[d] - lo[d]) / bin).floor() as i64 + 1).max(1);
    }
    let cell_of = |p: &[f64; 3]| -> [i64; 3] {
        [
            (((p[0] - lo[0]) / bin) as i64).clamp(0, nb[0] - 1),
            (((p[1] - lo[1]) / bin) as i64).clamp(0, nb[1] - 1),
            (((p[2] - lo[2]) / bin) as i64).clamp(0, nb[2] - 1),
        ]
    };
    let idx = |c: [i64; 3]| -> usize { ((c[0] * nb[1] + c[1]) * nb[2] + c[2]) as usize };
    let total = (nb[0] * nb[1] * nb[2]) as usize;
    let mut cells: Vec<Vec<usize>> = vec![Vec::new(); total];
    for (i, p) in pos.iter().enumerate() {
        cells[idx(cell_of(p))].push(i);
    }
    let mut offsets = vec![0u32];
    let mut indices = Vec::new();
    for i in 0..n {
        let c = cell_of(&pos[i]);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let cc = [c[0] + dx, c[1] + dy, c[2] + dz];
                    if cc[0] < 0 || cc[0] >= nb[0] || cc[1] < 0 || cc[1] >= nb[1] || cc[2] < 0 || cc[2] >= nb[2] {
                        continue;
                    }
                    for &j in &cells[idx(cc)] {
                        if j == i { continue; }
                        let dxp = pos[j][0] - pos[i][0];
                        let dyp = pos[j][1] - pos[i][1];
                        let dzp = pos[j][2] - pos[i][2];
                        if dxp * dxp + dyp * dyp + dzp * dzp < cut_sq {
                            indices.push(j as u32);
                        }
                    }
                }
            }
        }
        offsets.push(indices.len() as u32);
    }
    let mut neighbor = Neighbor::new();
    neighbor.newton = false;
    neighbor.neighbor_offsets = offsets;
    neighbor.neighbor_indices = indices;
    neighbor
}

/// Zero force and torque each step (PreForce), mirroring soil's force-zeroing.
fn zero_force_torque(mut atoms: ResMut<Atom>, registry: Res<AtomDataRegistry>) {
    let n = atoms.len();
    for i in 0..n { atoms.force[i] = [0.0; 3]; }
    let mut dem = registry.expect_mut::<DemAtom>("zero");
    for i in 0..n { dem.torque[i] = [0.0; 3]; }
}

fn bench_cpu_full_step(n: usize, r: f64, steps: usize, mt: MaterialTable) -> f64 {
    let pos = make_packing(n, r);
    let neighbor = build_full_neighbor_list(&pos, r);

    let mut app = App::new();
    let mut atom = Atom::new();
    let mut dem = DemAtom::new();
    let mut hist = ContactHistoryStore::new();
    // Stable dt ~ 1/40 of contact period.
    let mass = sphere_mass(r);
    let (e_eff, _, _, _) = gpu_scalars(&mt);
    let k_n = 4.0 / 3.0 * e_eff as f64 * (0.05 * r * r).sqrt();
    let dt = (2.0 * std::f64::consts::PI * (mass / k_n).sqrt()) / 40.0;
    atom.dt = dt;
    for i in 0..n {
        push_atom(&mut atom, &mut dem, i as u32, pos[i], [0.0; 3], [0.0; 3], r);
        hist.contacts.push(Vec::new());
    }
    atom.nlocal = n as u32;
    atom.natoms = n as u64;

    // Floor + four side walls.
    let walls = box_walls(&pos, r);

    let mut registry = AtomDataRegistry::new();
    registry.register(dem);
    registry.register(hist);

    app.add_resource(atom);
    app.add_resource(neighbor);
    app.add_resource(registry);
    app.add_resource(mt);
    app.add_resource(walls);

    // One full velocity-Verlet DEM step, all REAL dirt systems, single-threaded:
    app.add_update_system(initial_integration, ParticleSimScheduleSet::InitialIntegration);
    app.add_update_system(initial_rotation, ParticleSimScheduleSet::InitialIntegration);
    app.add_update_system(zero_force_torque, ParticleSimScheduleSet::PreForce);
    app.add_update_system(hertz_mindlin_contact_force, ParticleSimScheduleSet::Force);
    app.add_update_system(wall_contact_force, ParticleSimScheduleSet::Force);
    app.add_update_system(final_integration, ParticleSimScheduleSet::FinalIntegration);
    app.add_update_system(final_rotation, ParticleSimScheduleSet::FinalIntegration);
    app.organize_systems();

    // Warm-up (build caches, JIT-free but stabilize allocs).
    app.run();
    let t0 = Instant::now();
    for _ in 0..steps { app.run(); }
    let elapsed = t0.elapsed().as_secs_f64();
    elapsed * 1000.0 / steps as f64
}

fn box_walls(pos: &[[f64; 3]], r: f64) -> Walls {
    let mut hi = [f64::MIN; 3];
    for p in pos { for d in 0..3 { hi[d] = hi[d].max(p[d]); } }
    let geo = [
        ([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        ([hi[0] + r, 0.0, 0.0], [-1.0, 0.0, 0.0]),
        ([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        ([0.0, hi[1] + r, 0.0], [0.0, -1.0, 0.0]),
    ];
    let mut planes = Vec::new();
    let mut active = Vec::new();
    for (point, normal) in geo {
        planes.push(WallPlane {
            point_x: point[0], point_y: point[1], point_z: point[2],
            normal_x: normal[0], normal_y: normal[1], normal_z: normal[2],
            material_index: 0, name: None,
            bound_x_low: f64::NEG_INFINITY, bound_x_high: f64::INFINITY,
            bound_y_low: f64::NEG_INFINITY, bound_y_high: f64::INFINITY,
            bound_z_low: f64::NEG_INFINITY, bound_z_high: f64::INFINITY,
            velocity: [0.0; 3], motion: WallMotion::Static, origin: point,
            force_accumulator: 0.0, temperature: None,
        });
        active.push(true);
    }
    Walls {
        planes, active,
        cylinders: Vec::new(), cylinder_active: Vec::new(),
        spheres: Vec::new(), sphere_active: Vec::new(),
        regions: Vec::new(), region_active: Vec::new(),
        time: 0.0,
        tangential_springs: std::collections::HashMap::new(),
        rolling_springs: std::collections::HashMap::new(),
    }
}

/// Same full-step benchmark, but via the NEW Force-hook stack (soil GpuState +
/// dirt GranularForce + WallForce). Proves the resident loop runs at scale and
/// at monolith-equivalent speed through the hook path.
fn bench_gpu_full_step(ctx: GpuContext, n: usize, r: f64, steps: usize, mt: &MaterialTable) -> f64 {
    use dirt_gpu::{GranularConfig, GranularForce, WallForce};
    use soil_gpu::{Boundary, GpuState, Grid as SoilGrid, Plane};

    let pos = make_packing(n, r);
    let posf: Vec<[f32; 3]> = pos.iter().map(|p| [p[0] as f32, p[1] as f32, p[2] as f32]).collect();
    let velf = vec![[0.0f32; 3]; n];
    let omf = vec![[0.0f32; 3]; n];
    let radiusf = vec![r as f32; n];
    let inv_mass = vec![(1.0 / sphere_mass(r)) as f32; n];
    let inv_inertia = vec![sphere_inv_inertia(r) as f32; n];
    let (e_eff, beta, g_eff, mu) = gpu_scalars(mt);

    let mass = sphere_mass(r);
    let k_n = 4.0 / 3.0 * e_eff as f64 * (0.05 * r * r).sqrt();
    let dt = ((2.0 * std::f64::consts::PI * (mass / k_n).sqrt()) / 40.0) as f32;

    let grid = SoilGrid::from_positions(&posf, (2.0 * r) as f32);
    let mut hi = [f32::MIN; 3];
    for p in &posf { for d in 0..3 { hi[d] = hi[d].max(p[d]); } }

    let mut gs = GpuState::new(ctx, n, grid.total_cells);
    gs.set_params(dt, [0.0, 0.0, -9.81]);
    gs.set_state(&posf, &velf, &inv_mass, grid);
    let om = gs.add_aux_dof();
    gs.set_aux_inv_coeff(om, &inv_inertia);
    gs.set_aux_state(om, &omf);

    let cfg = GranularConfig { e_eff, beta, g_eff, mu, dt };
    gs.add_force_hook(Box::new(GranularForce::new(&gs, &grid, om, &radiusf, cfg)));
    let mut b = Boundary::new();
    b.push(Plane::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]));
    b.push(Plane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]));
    b.push(Plane::new([hi[0] + r as f32, 0.0, 0.0], [-1.0, 0.0, 0.0]));
    b.push(Plane::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]));
    b.push(Plane::new([0.0, hi[1] + r as f32, 0.0], [0.0, -1.0, 0.0]));
    gs.add_force_hook(Box::new(WallForce::new(&gs, om, &radiusf, &b, e_eff, beta, g_eff, mu, dt)));

    gs.run_steps(1); // warm-up
    gs.wait();
    let t0 = Instant::now();
    gs.run_steps(steps);
    gs.wait();
    let elapsed = t0.elapsed().as_secs_f64();
    elapsed * 1000.0 / steps as f64
}

fn tier3(ctx: &GpuContext) {
    println!("\n=== Tier 3: fair full-step benchmark (CPU real dirt vs GPU resident) ===");
    println!("CPU step = soil_verlet initial_integration + rotational initial + dirt hertz_mindlin");
    println!("           + dirt_wall + final_integration + rotational final (single-threaded).");
    println!("GPU step = run_steps via soil GpuState + dirt GranularForce/WallForce hooks, resident");
    println!("           (integrate + aux-rotation + cell-list + Hertz/Mindlin + wall, one submit).");
    let r = 0.001_f64;
    let mt = matched_material_table();

    let configs = [(8_000usize, 200usize), (64_000, 60), (216_000, 30)];
    println!("\n{:>10} | {:>14} | {:>14} | {:>10}", "N", "CPU ms/step", "GPU ms/step", "speedup");
    println!("{}", "-".repeat(60));
    for (n, steps) in configs {
        let cpu_ms = bench_cpu_full_step(n, r, steps, matched_material_table());
        let gpu_ms = bench_gpu_full_step(ctx.clone(), n, r, steps, &mt);
        println!("{:>10} | {:>14.4} | {:>14.4} | {:>9.1}x", n, cpu_ms, gpu_ms, cpu_ms / gpu_ms);
    }
}

fn main() {
    let Some(ctx) = GpuContext::new() else {
        eprintln!("No GPU adapter available; cannot run comparison.");
        std::process::exit(1);
    };
    println!("GPU adapter: {}", ctx.adapter_info);
    tier1(&ctx);
    tier3(&ctx);
    println!("\nDone.");
}
