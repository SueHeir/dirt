//! Full beam bonded-pair force as a `GpuForce` hook (Phase B of gpu_bonds_plan.md).
//!
//! Ports the CPU `dirt_bond` beam: normal + history shear + twist + bending
//! moments, per-channel critical damping, and beam-stress breakage. Built on the
//! Phase-A persistent-bond layout ([`crate::BondTopology`]) extended with per-bond
//! shear (`delta_t`) and rotation (`delta_theta`) history, updated in place.
//!
//! Atomic-free i-centric, like the contact kernel: a bond `(i,j)` is double-stored,
//! and each endpoint computes its own half as the exact mirror image of the other
//! (n̂ flips, v_rel flips, ω_rel flips, the histories stay negatives of each other),
//! applying force/torque only to itself. This reproduces the CPU's single-owner
//! `+f on i / −f on j`, `+M on i / −M on j`, shear-torque-same-on-both result —
//! verified term by term against `dirt_bond::bond_force`.
//!
//! Scope (this file): elastic beam (no plasticity) with the constant beam-stress
//! breakage criterion (σ = Fₙ/A + 2|M_bend|r_b/J, τ = |F_t|/A + |M_tor|r_b/J). Per-
//! bond Weibull thresholds and plasticity are follow-ons. Material mode only
//! (stiffness from E, G + per-bond geometry).
//!
//! Torque: this hook OWNS the aux-rate (torque) buffer for the atoms it touches
//! (writes, not accumulates) — correct for a bonds-only resident loop. Composing
//! with the granular hook needs a torque-seed kernel in soil_gpu (a later step).

use bytemuck::{Pod, Zeroable};
use soil_gpu::{GpuForce, GpuState};

use crate::BondTopology;

/// Elastic beam bond parameters (material mode + constant breakage).
#[derive(Clone, Copy, Debug)]
pub struct BeamBondConfig {
    /// Bond radius = `bond_radius_ratio · min(Rᵢ, Rⱼ)`.
    pub bond_radius_ratio: f32,
    pub youngs_modulus: f32,
    pub shear_modulus: f32,
    pub beta_normal: f32,
    pub beta_shear: f32,
    pub beta_twist: f32,
    pub beta_bending: f32,
    /// Tensile beam-stress limit σ_max; set huge for unbreakable.
    pub sigma_max: f32,
    /// Shear beam-stress limit τ_max; set huge for unbreakable.
    pub tau_max: f32,
    pub dt: f32,
    /// Periodic box lengths (0 on an axis = non-periodic). When set, the bond
    /// vector uses the triclinic minimum image — so a BPM aggregate spanning a
    /// periodic boundary stays bonded. Bonds don't use the cell list, so this needs
    /// no cell-list change.
    pub lx: f32,
    pub ly: f32,
    pub lz: f32,
    /// Lees–Edwards xy tilt (x shift per y-image). 0 = orthogonal box.
    pub tilt_xy: f32,
    /// If true, ACCUMULATE torque (`+=`) onto the aux-rate buffer instead of owning
    /// it (`=`). Set true when composing with a contact hook that seeds torque first
    /// (the resident BPM loop); false for a bonds-only run where this hook owns it.
    pub accumulate_torque: bool,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct BeamParams {
    n: u32,
    dt: f32,
    ratio: f32,
    e_mod: f32,
    g_mod: f32,
    beta_n: f32,
    beta_t: f32,
    beta_tor: f32,
    beta_bend: f32,
    sigma_max: f32,
    tau_max: f32,
    lx: f32,
    ly: f32,
    lz: f32,
    tilt_xy: f32,
    accum: u32,
}

/// Elastic beam bond force hook. Group 0 binds the resident pos/vel/force/omega/
/// torque + radius/inv_mass/inv_inertia + params; group 1 the persistent bond CSR
/// (offsets/partner/r0/broken) plus per-bond `delta_t` and `delta_theta` history.
pub struct BeamBondForce {
    n: u32,
    pipeline: wgpu::ComputePipeline,
    g0: wgpu::BindGroup,
    g1: wgpu::BindGroup,
    broken_buf: wgpu::Buffer,
    staging: wgpu::Buffer,
    num_bonds: usize,
}

impl BeamBondForce {
    pub fn new(
        gs: &GpuState,
        omega_aux: usize,
        radius: &[f32],
        topo: &BondTopology,
        cfg: BeamBondConfig,
    ) -> Self {
        let ctx = gs.context();
        let device = &ctx.device;
        let n = gs.n();
        assert_eq!(radius.len(), n);
        assert_eq!(topo.offsets.len(), n + 1);
        let num_bonds = topo.num_bonds();

        let storage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST;
        let mk_u32 = |label: &str, data: &[u32], src: bool| {
            let bytes = (data.len().max(1) * 4) as u64;
            let usage = if src { storage | wgpu::BufferUsages::COPY_SRC } else { storage };
            let buf = device.create_buffer(&wgpu::BufferDescriptor { label: Some(label), size: bytes, usage, mapped_at_creation: false });
            if !data.is_empty() { ctx.queue.write_buffer(&buf, 0, bytemuck::cast_slice(data)); }
            buf
        };
        let mk_f32 = |label: &str, data: &[f32]| {
            let bytes = (data.len().max(1) * 4) as u64;
            let buf = device.create_buffer(&wgpu::BufferDescriptor { label: Some(label), size: bytes, usage: storage, mapped_at_creation: false });
            if !data.is_empty() { ctx.queue.write_buffer(&buf, 0, bytemuck::cast_slice(data)); }
            buf
        };

        let radius_buf = mk_f32("beam_radius", radius);
        let offsets_buf = mk_u32("beam_offsets", &topo.offsets, false);
        let partner_buf = mk_u32("beam_partner", &topo.partner, false);
        let r0_buf = mk_f32("beam_r0", &topo.r0);
        let broken_buf = mk_u32("beam_broken", &vec![0u32; num_bonds], true);
        let delta_t_buf = mk_f32("beam_delta_t", &vec![0.0f32; num_bonds * 3]);
        let delta_theta_buf = mk_f32("beam_delta_theta", &vec![0.0f32; num_bonds * 3]);

        let staging = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("beam_broken_staging"), size: (num_bonds.max(1) * 4) as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
        });

        let params = BeamParams {
            n: n as u32, dt: cfg.dt, ratio: cfg.bond_radius_ratio,
            e_mod: cfg.youngs_modulus, g_mod: cfg.shear_modulus,
            beta_n: cfg.beta_normal, beta_t: cfg.beta_shear, beta_tor: cfg.beta_twist, beta_bend: cfg.beta_bending,
            sigma_max: cfg.sigma_max, tau_max: cfg.tau_max,
            lx: cfg.lx, ly: cfg.ly, lz: cfg.lz, tilt_xy: cfg.tilt_xy,
            accum: u32::from(cfg.accumulate_torque),
        };
        let params_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("beam_params"), size: std::mem::size_of::<BeamParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
        });
        ctx.queue.write_buffer(&params_buf, 0, bytemuck::bytes_of(&params));

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("beam_bond_force"), source: wgpu::ShaderSource::Wgsl(BEAM_WGSL.into()),
        });
        let st = |binding, read_only| wgpu::BindGroupLayoutEntry {
            binding, visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only }, has_dynamic_offset: false, min_binding_size: None },
            count: None,
        };
        let uni = |binding| wgpu::BindGroupLayoutEntry {
            binding, visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Uniform, has_dynamic_offset: false, min_binding_size: None },
            count: None,
        };

        let g0l = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("beam g0 bgl"),
            entries: &[st(0, true), st(1, true), st(2, false), st(3, true), st(4, false), st(5, true), st(6, true), st(7, true), uni(8)],
        });
        let g0 = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("beam g0 bg"), layout: &g0l,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: gs.pos_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: gs.vel_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: gs.force_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: gs.aux_state_buffer(omega_aux).as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: gs.aux_rate_buffer(omega_aux).as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: radius_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: gs.inv_mass_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 7, resource: gs.aux_inv_coeff_buffer(omega_aux).as_entire_binding() },
                wgpu::BindGroupEntry { binding: 8, resource: params_buf.as_entire_binding() },
            ],
        });

        let g1l = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("beam g1 bgl"),
            entries: &[st(0, true), st(1, true), st(2, true), st(3, false), st(4, false), st(5, false)],
        });
        let g1 = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("beam g1 bg"), layout: &g1l,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: offsets_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: partner_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: r0_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: broken_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: delta_t_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: delta_theta_buf.as_entire_binding() },
            ],
        });

        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("beam pl"), bind_group_layouts: &[Some(&g0l), Some(&g1l)], immediate_size: 0,
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("beam_bond_force"), layout: Some(&layout), module: &shader, entry_point: Some("beam_bond_force"),
            compilation_options: wgpu::PipelineCompilationOptions::default(), cache: None,
        });

        let _ = (radius_buf, offsets_buf, partner_buf, r0_buf, delta_t_buf, delta_theta_buf, params_buf);
        BeamBondForce { n: n as u32, pipeline, g0, g1, broken_buf, staging, num_bonds }
    }

    pub fn download_broken(&self, gs: &GpuState) -> Vec<u32> {
        if self.num_bonds == 0 { return Vec::new(); }
        let ctx = gs.context();
        let bytes = (self.num_bonds * 4) as u64;
        let mut enc = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("beam dl") });
        enc.copy_buffer_to_buffer(&self.broken_buf, 0, &self.staging, 0, bytes);
        ctx.queue.submit(Some(enc.finish()));
        let slice = self.staging.slice(0..bytes);
        slice.map_async(wgpu::MapMode::Read, |_| {});
        ctx.device.poll(wgpu::PollType::wait_indefinitely()).expect("poll");
        let data = slice.get_mapped_range();
        let v: Vec<u32> = bytemuck::cast_slice(&data).to_vec();
        drop(data);
        self.staging.unmap();
        v
    }
}

impl GpuForce for BeamBondForce {
    fn record(&self, pass: &mut wgpu::ComputePass) {
        pass.set_bind_group(0, &self.g0, &[]);
        pass.set_bind_group(1, &self.g1, &[]);
        pass.set_pipeline(&self.pipeline);
        pass.dispatch_workgroups(self.n.div_ceil(64).max(1), 1, 1);
    }
}

const BEAM_WGSL: &str = r#"
const PI: f32 = 3.14159265358979323846;

struct BeamParams {
    n: u32, dt: f32, ratio: f32, e_mod: f32, g_mod: f32,
    beta_n: f32, beta_t: f32, beta_tor: f32, beta_bend: f32,
    sigma_max: f32, tau_max: f32,
    lx: f32, ly: f32, lz: f32, tilt_xy: f32, accum: u32,
};

@group(0) @binding(0) var<storage, read>       pos: array<f32>;
@group(0) @binding(1) var<storage, read>       vel: array<f32>;
@group(0) @binding(2) var<storage, read_write> force_out: array<f32>;
@group(0) @binding(3) var<storage, read>       omega: array<f32>;
@group(0) @binding(4) var<storage, read_write> torque: array<f32>;
@group(0) @binding(5) var<storage, read>       radius: array<f32>;
@group(0) @binding(6) var<storage, read>       inv_mass: array<f32>;
@group(0) @binding(7) var<storage, read>       inv_inertia: array<f32>;
@group(0) @binding(8) var<uniform>             params: BeamParams;

@group(1) @binding(0) var<storage, read>       bond_offsets: array<u32>;
@group(1) @binding(1) var<storage, read>       bond_partner: array<u32>;
@group(1) @binding(2) var<storage, read>       bond_r0: array<f32>;
@group(1) @binding(3) var<storage, read_write> bond_broken: array<u32>;
@group(1) @binding(4) var<storage, read_write> delta_t: array<f32>;       // 3 per bond
@group(1) @binding(5) var<storage, read_write> delta_theta: array<f32>;   // 3 per bond

@compute @workgroup_size(64)
fn beam_bond_force(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n) { return; }
    let bi = 3u * i;
    let xi = vec3<f32>(pos[bi], pos[bi + 1u], pos[bi + 2u]);
    let vi = vec3<f32>(vel[bi], vel[bi + 1u], vel[bi + 2u]);
    let wi = vec3<f32>(omega[bi], omega[bi + 1u], omega[bi + 2u]);
    let imi = inv_mass[i];
    let iii = inv_inertia[i];
    let ri = radius[i];

    var f = vec3<f32>(force_out[bi], force_out[bi + 1u], force_out[bi + 2u]); // accumulate
    var t = vec3<f32>(0.0, 0.0, 0.0);                                          // own torque

    let start = bond_offsets[i];
    let end = bond_offsets[i + 1u];
    for (var k = start; k < end; k = k + 1u) {
        if (bond_broken[k] != 0u) { continue; }
        let j = bond_partner[k];
        let bj = 3u * j;
        let xj = vec3<f32>(pos[bj], pos[bj + 1u], pos[bj + 2u]);
        // Triclinic minimum image (Lees–Edwards): wrap z, then y (a y-image shifts
        // x by the tilt), then x. Lets a BPM aggregate span a periodic boundary.
        // (NOTE: the LE velocity offset Δv for the damping term on a y-wrapped bond
        // is not yet applied — elastic position term only; orthogonal-periodic is
        // fully correct.)
        var d = xj - xi;
        if (params.lz > 0.0) { d.z = d.z - params.lz * round(d.z / params.lz); }
        if (params.ly > 0.0) {
            let ny = round(d.y / params.ly);
            d.y = d.y - params.ly * ny;
            d.x = d.x - params.tilt_xy * ny;
        }
        if (params.lx > 0.0) { d.x = d.x - params.lx * round(d.x / params.lx); }
        let dist = length(d);
        if (dist < 1.0e-20) { continue; }
        let nhat = d / dist;
        let r0 = bond_r0[k];
        let delta = dist - r0;

        // Beam geometry (cylindrical).
        let r_b = params.ratio * min(ri, radius[j]);
        let area = PI * r_b * r_b;
        let jpol = 0.5 * PI * r_b * r_b * r_b * r_b;
        let iben = 0.5 * jpol;
        let len = r0;
        let k_n = params.e_mod * area / len;
        let k_t = params.g_mod * area / len;
        let k_tor = params.g_mod * jpol / len;
        let k_bend = params.e_mod * iben / len;

        // Reduced mass / MOI for critical damping.
        let m_red = 1.0 / (imi + inv_mass[j]);
        let moi_red = 1.0 / (iii + inv_inertia[j]);
        let gamma_n = 2.0 * params.beta_n * sqrt(m_red * max(k_n, 0.0));
        let gamma_t = 2.0 * params.beta_t * sqrt(m_red * max(k_t, 0.0));
        let gamma_tor = 2.0 * params.beta_tor * sqrt(moi_red * max(k_tor, 0.0));
        let gamma_bend = 2.0 * params.beta_bend * sqrt(moi_red * max(k_bend, 0.0));

        // Kinematics at bond mid-point (lever r1 = (L/2) n̂ from i).
        let half_l = 0.5 * len;
        let r1 = half_l * nhat;
        let wj = vec3<f32>(omega[bj], omega[bj + 1u], omega[bj + 2u]);
        let vj = vec3<f32>(vel[bj], vel[bj + 1u], vel[bj + 2u]);
        let v_i_c = vi + cross(wi, r1);
        let v_j_c = vj - cross(wj, r1);   // r2 = -r1
        let v_rel = v_j_c - v_i_c;
        let v_n_s = dot(v_rel, nhat);
        let v_n = v_n_s * nhat;
        let v_t = v_rel - v_n;

        // Normal force (elastic + damping).
        let f_n_mag = k_n * delta + gamma_n * v_n_s;
        let f_n = f_n_mag * nhat;

        // Shear history: reproject ⊥ to new n̂, integrate.
        let kt3 = 3u * k;
        var ds = vec3<f32>(delta_t[kt3], delta_t[kt3 + 1u], delta_t[kt3 + 2u]);
        ds = ds - dot(ds, nhat) * nhat + v_t * params.dt;
        let f_t = k_t * ds + gamma_t * v_t;

        // Rotation kinematics + Δθ split into twist (∥n̂) and bend (⊥n̂).
        let w_rel = wj - wi;
        let w_rel_n_s = dot(w_rel, nhat);
        let w_n = w_rel_n_s * nhat;
        let w_t = w_rel - w_n;
        let dt3 = 3u * k;
        var dth = vec3<f32>(delta_theta[dt3], delta_theta[dt3 + 1u], delta_theta[dt3 + 2u]) + w_rel * params.dt;
        let dth_n_s = dot(dth, nhat);
        let dth_twist = dth_n_s * nhat;
        let dth_bend = dth - dth_twist;

        let m_tor = k_tor * dth_twist + gamma_tor * w_n;
        let m_bend = k_bend * dth_bend + gamma_bend * w_t;

        // Breakage: constant beam-stress criterion at the extreme fibre.
        let f_t_mag = length(f_t);
        let m_bend_mag = length(m_bend);
        let m_tor_mag = length(m_tor);
        let sigma = f_n_mag / area + 2.0 * m_bend_mag * r_b / jpol;
        let tau = f_t_mag / area + m_tor_mag * r_b / jpol;
        if (sigma > params.sigma_max || tau > params.tau_max) {
            bond_broken[k] = 1u;
            continue;
        }

        // Persist the advanced history (in place; the partner stores the mirror).
        delta_t[kt3] = ds.x; delta_t[kt3 + 1u] = ds.y; delta_t[kt3 + 2u] = ds.z;
        delta_theta[dt3] = dth.x; delta_theta[dt3 + 1u] = dth.y; delta_theta[dt3 + 2u] = dth.z;

        // Apply to self (mirror image gives the partner its equal/opposite half).
        f = f + f_n + f_t;
        let tau_shear = cross(r1, f_t);
        t = t + tau_shear + m_tor + m_bend;
    }

    force_out[bi] = f.x; force_out[bi + 1u] = f.y; force_out[bi + 2u] = f.z;
    // Own (=) when standalone, or accumulate (+=) onto a contact hook's seeded
    // torque when composing in the resident BPM loop.
    if (params.accum != 0u) {
        torque[bi] = torque[bi] + t.x;
        torque[bi + 1u] = torque[bi + 1u] + t.y;
        torque[bi + 2u] = torque[bi + 2u] + t.z;
    } else {
        torque[bi] = t.x; torque[bi + 1u] = t.y; torque[bi + 2u] = t.z;
    }
}
"#;

#[cfg(test)]
mod tests {
    use super::*;
    use soil_gpu::{GpuContext, Grid};

    // A single beam bond between two unit spheres. Returns final pos, vel, omega.
    fn run(p0: [[f32; 3]; 2], v0: [[f32; 3]; 2], w0: [[f32; 3]; 2], r0: f32, steps: usize)
        -> (Vec<[f32; 3]>, Vec<[f32; 3]>, Vec<[f32; 3]>) {
        let ctx = GpuContext::new().expect("gpu");
        let rad = 0.5f32;
        let mass = 1.0f32;
        let inertia = 0.4 * mass * rad * rad;
        let inv_mass = vec![1.0 / mass, 1.0 / mass];
        let inv_inertia = vec![1.0 / inertia, 1.0 / inertia];
        let grid = Grid::from_positions(&p0, 4.0 * rad);
        let dt = 1.0e-4f32;
        let mut gs = GpuState::new(ctx, 2, grid.total_cells);
        gs.set_params(dt, [0.0, 0.0, 0.0]);
        gs.set_state(&p0, &v0, &inv_mass, grid);
        let omega = gs.add_aux_dof();
        gs.set_aux_inv_coeff(omega, &inv_inertia);
        gs.set_aux_state(omega, &w0);

        let topo = BondTopology { offsets: vec![0, 1, 2], partner: vec![1, 0], r0: vec![r0, r0] };
        let cfg = BeamBondConfig {
            bond_radius_ratio: 1.0, youngs_modulus: 1.0e7, shear_modulus: 4.0e6,
            beta_normal: 0.1, beta_shear: 0.1, beta_twist: 0.1, beta_bending: 0.1,
            sigma_max: 1.0e30, tau_max: 1.0e30, dt,
            lx: 0.0, ly: 0.0, lz: 0.0, tilt_xy: 0.0, accumulate_torque: false,
        };
        gs.add_force_hook(Box::new(BeamBondForce::new(&gs, omega, &[rad, rad], &topo, cfg)));
        gs.run_steps(steps);
        (gs.download_pos(), gs.download_vel(), gs.download_aux_state(omega))
    }

    #[test]
    fn beam_normal_restores_symmetric_momentum_conserved() {
        if GpuContext::new().is_none() { eprintln!("no GPU; skipping"); return; }
        // Stretched along x: r0=0.9, placed at ±0.5 (L=1.0). Should contract.
        let (p, v, _w) = run([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]], [[0.0; 3]; 2], [[0.0; 3]; 2], 0.9, 300);
        assert!(p[0][0] > -0.5 && p[1][0] < 0.5, "did not contract: {p:?}");
        assert!((p[0][0] + p[1][0]).abs() < 1e-4, "asymmetric: {p:?}");
        assert!((v[0][0] + v[1][0]).abs() < 1e-4, "linear momentum not conserved: {v:?}");
    }

    /// f64 reference for one bond evaluation (history starts at zero, so this is
    /// the first-step force/torque on atom 0). Mirrors the CPU `dirt_bond::bond_force`
    /// formula term for term — the parity oracle.
    #[allow(clippy::too_many_arguments)]
    fn beam_ref(p: [[f64; 3]; 2], v: [[f64; 3]; 2], w: [[f64; 3]; 2], r0: f64, rad: f64,
                mass: f64, inertia: f64, e: f64, g: f64, beta: f64, dt: f64)
        -> ([f64; 3], [f64; 3]) {
        use std::f64::consts::PI;
        let sub = |a: [f64;3], b: [f64;3]| [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
        let add = |a: [f64;3], b: [f64;3]| [a[0]+b[0], a[1]+b[1], a[2]+b[2]];
        let scale = |a: [f64;3], s: f64| [a[0]*s, a[1]*s, a[2]*s];
        let dot = |a: [f64;3], b: [f64;3]| a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
        let cross = |a: [f64;3], b: [f64;3]| [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
        let d = sub(p[1], p[0]);
        let dist = dot(d, d).sqrt();
        let nhat = scale(d, 1.0/dist);
        let delta = dist - r0;
        let r_b = rad.min(rad);
        let area = PI*r_b*r_b;
        let jpol = 0.5*PI*r_b.powi(4);
        let iben = 0.5*jpol;
        let len = r0;
        let k_n = e*area/len; let k_t = g*area/len; let k_tor = g*jpol/len; let k_bend = e*iben/len;
        let m_red = mass*mass/(mass+mass);
        let moi_red = inertia*inertia/(inertia+inertia);
        let gamma_n = 2.0*beta*(m_red*k_n).sqrt();
        let gamma_t = 2.0*beta*(m_red*k_t).sqrt();
        let gamma_tor = 2.0*beta*(moi_red*k_tor).sqrt();
        let gamma_bend = 2.0*beta*(moi_red*k_bend).sqrt();
        let r1 = scale(nhat, 0.5*len);
        let v_i_c = add(v[0], cross(w[0], r1));
        let v_j_c = sub(v[1], cross(w[1], r1));
        let v_rel = sub(v_j_c, v_i_c);
        let v_n_s = dot(v_rel, nhat);
        let v_t = sub(v_rel, scale(nhat, v_n_s));
        let f_n_mag = k_n*delta + gamma_n*v_n_s;
        let f_n = scale(nhat, f_n_mag);
        let ds = scale(v_t, dt); // history starts 0
        let f_t = add(scale(ds, k_t), scale(v_t, gamma_t));
        let w_rel = sub(w[1], w[0]);
        let w_rel_n_s = dot(w_rel, nhat);
        let w_n = scale(nhat, w_rel_n_s);
        let w_t = sub(w_rel, w_n);
        let dth = scale(w_rel, dt);
        let dth_n_s = dot(dth, nhat);
        let dth_twist = scale(nhat, dth_n_s);
        let dth_bend = sub(dth, dth_twist);
        let m_tor = add(scale(dth_twist, k_tor), scale(w_n, gamma_tor));
        let m_bend = add(scale(dth_bend, k_bend), scale(w_t, gamma_bend));
        let force = add(f_n, f_t);
        let tau_shear = cross(r1, f_t);
        let torque = add(tau_shear, add(m_tor, m_bend));
        (force, torque)
    }

    #[test]
    fn beam_force_torque_matches_cpu_formula() {
        if GpuContext::new().is_none() { eprintln!("no GPU; skipping"); return; }
        let ctx = GpuContext::new().unwrap();
        let (rad, mass) = (0.5f64, 1.0f64);
        let inertia = 0.4*mass*rad*rad;
        let (e, g, beta, dt, r0) = (1.0e7f64, 4.0e6f64, 0.1f64, 1.0e-4f64, 0.9f64);
        let p = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let v = [[0.0, 0.2, 0.0], [0.0, 0.0, 0.0]];
        let w = [[0.0, 0.0, 0.3], [0.0, 0.0, 0.0]];

        let (ef, et) = beam_ref(p, v, w, r0, rad, mass, inertia, e, g, beta, dt);

        // GPU: one force evaluation at this config (gravity 0).
        let p0: [[f32;3];2] = [[0.0;3], [1.0,0.0,0.0]];
        let grid = Grid::from_positions(&p0, 4.0*rad as f32);
        let mut gs = GpuState::new(ctx, 2, grid.total_cells);
        gs.set_params(dt as f32, [0.0; 3]);
        gs.set_state(&p0, &[[0.0,0.2,0.0],[0.0;3]], &[1.0/mass as f32, 1.0/mass as f32], grid);
        let omega = gs.add_aux_dof();
        gs.set_aux_inv_coeff(omega, &[1.0/inertia as f32, 1.0/inertia as f32]);
        gs.set_aux_state(omega, &[[0.0,0.0,0.3],[0.0;3]]);
        let topo = BondTopology { offsets: vec![0,1,2], partner: vec![1,0], r0: vec![r0 as f32, r0 as f32] };
        let cfg = BeamBondConfig {
            bond_radius_ratio: 1.0, youngs_modulus: e as f32, shear_modulus: g as f32,
            beta_normal: beta as f32, beta_shear: beta as f32, beta_twist: beta as f32, beta_bending: beta as f32,
            sigma_max: 1.0e30, tau_max: 1.0e30, dt: dt as f32,
            lx: 0.0, ly: 0.0, lz: 0.0, tilt_xy: 0.0, accumulate_torque: false,
        };
        gs.add_force_hook(Box::new(BeamBondForce::new(&gs, omega, &[rad as f32, rad as f32], &topo, cfg)));
        gs.eval_force_once();
        let gf = gs.download_force();
        let gt = gs.download_aux_rate(omega);

        let rel = |a: f32, b: f64| (a as f64 - b).abs() / b.abs().max(1.0);
        for c in 0..3 {
            assert!(rel(gf[0][c], ef[c]) < 2e-3, "force[{c}] gpu={} cpu={} ", gf[0][c], ef[c]);
            assert!(rel(gt[0][c], et[c]) < 2e-3, "torque[{c}] gpu={} cpu={}", gt[0][c], et[c]);
        }
        eprintln!("beam parity: F gpu={:?} cpu={:?}\n             T gpu={:?} cpu={:?}", gf[0], ef, gt[0], et);
    }

    /// End-to-end resident stack: a bonded pair in a periodic box, advanced by the
    /// resident loop with BOTH the contact hook (periodic, owns torque) and the beam
    /// bond hook (periodic, accumulate_torque) composing. Validates the residency
    /// steps 2–4 together: periodic contact + on-device PBC remap + resident bond +
    /// torque composition. The bond must hold the pair together under a transverse
    /// kick, with momentum conserved.
    #[test]
    fn resident_contact_plus_bond_periodic_compose() {
        use crate::{GranularConfig, GranularForce};
        let Some(ctx) = GpuContext::new() else { eprintln!("no GPU; skipping"); return; };
        let (r, l) = (0.1f32, 1.0f32);
        let mass = (4.0 / 3.0) * std::f32::consts::PI * r * r * r;
        let inertia = 0.4 * mass * r * r;
        let nc = 6i32;
        let grid = Grid { n: [nc, nc, nc], origin: [0.0; 3], bin_size: l / nc as f32, total_cells: (nc * nc * nc) as usize };
        let p0 = [[0.4f32, 0.5, 0.5], [0.6, 0.5, 0.5]]; // dist 0.2 = r0 (bond at rest)
        let dt = 1.0e-5f32;
        let mut gs = GpuState::new(ctx, 2, grid.total_cells);
        gs.set_params(dt, [0.0; 3]);
        gs.set_box([l, l, l], [0.0; 3], 0.0, 0.0); // periodic box, on-device remap
        gs.set_state(&p0, &[[0.0, 0.5, 0.0], [0.0; 3]], &[1.0 / mass, 1.0 / mass], grid.clone());
        let omega = gs.add_aux_dof();
        gs.set_aux_inv_coeff(omega, &[1.0 / inertia, 1.0 / inertia]);
        gs.set_aux_state(omega, &[[0.0; 3]; 2]);

        // Contact hook FIRST (owns torque, periodic), then bond hook (accumulates).
        let mut cc = GranularConfig::new(1.0e6, 0.1, 4.0e5, 0.5, dt);
        cc.lx = l; cc.ly = l; cc.lz = l;
        gs.add_force_hook(Box::new(GranularForce::new(&gs, &grid, omega, &[r, r], cc)));
        let bcfg = BeamBondConfig {
            bond_radius_ratio: 1.0, youngs_modulus: 1.0e7, shear_modulus: 4.0e6,
            beta_normal: 0.1, beta_shear: 0.1, beta_twist: 0.1, beta_bending: 0.1,
            sigma_max: 1.0e30, tau_max: 1.0e30, dt,
            lx: l, ly: l, lz: l, tilt_xy: 0.0, accumulate_torque: true,
        };
        let topo = BondTopology { offsets: vec![0, 1, 2], partner: vec![1, 0], r0: vec![0.2, 0.2] };
        gs.add_force_hook(Box::new(BeamBondForce::new(&gs, omega, &[r, r], &topo, bcfg)));

        gs.run_steps(500);
        let p = gs.download_pos();
        let v = gs.download_vel();
        // Bond holds the pair: separation stays near r0 (didn't fly apart or collapse).
        let d = [(p[1][0] - p[0][0]), (p[1][1] - p[0][1]), (p[1][2] - p[0][2])];
        let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        assert!((dist - 0.2).abs() < 0.05, "bond did not hold the pair: dist={dist}");
        // Linear momentum conserved at its initial value (atom 0 started vy=0.5,
        // equal masses → total vy stays 0.5; internal contact+bond forces only).
        let mom = [v[0][0] + v[1][0], v[0][1] + v[1][1], v[0][2] + v[1][2]];
        assert!((mom[1] - 0.5).abs() < 1e-3, "y-momentum not conserved: {mom:?}");
        assert!(mom[0].abs() < 1e-3 && mom[2].abs() < 1e-3, "spurious x/z momentum: {mom:?}");
        eprintln!("resident contact+bond periodic: dist={dist} mom={mom:?}");
    }

    #[test]
    fn beam_bond_periodic_minimum_image() {
        if GpuContext::new().is_none() { eprintln!("no GPU; skipping"); return; }
        // One force eval for a bonded pair; returns force on atom 0.
        let eval = |p0: [[f32; 3]; 2], lx: f32| -> [f32; 3] {
            let ctx = GpuContext::new().unwrap();
            let (rad, mass) = (0.5f32, 1.0f32);
            let inertia = 0.4 * mass * rad * rad;
            let grid = Grid::from_positions(&p0, 2.0 * rad);
            let mut gs = GpuState::new(ctx, 2, grid.total_cells);
            gs.set_params(1.0e-4, [0.0; 3]);
            gs.set_state(&p0, &[[0.0; 3]; 2], &[1.0 / mass, 1.0 / mass], grid);
            let omega = gs.add_aux_dof();
            gs.set_aux_inv_coeff(omega, &[1.0 / inertia, 1.0 / inertia]);
            gs.set_aux_state(omega, &[[0.0; 3]; 2]);
            let topo = BondTopology { offsets: vec![0, 1, 2], partner: vec![1, 0], r0: vec![0.2, 0.2] };
            let cfg = BeamBondConfig {
                bond_radius_ratio: 1.0, youngs_modulus: 1.0e7, shear_modulus: 4.0e6,
                beta_normal: 0.0, beta_shear: 0.0, beta_twist: 0.0, beta_bending: 0.0,
                sigma_max: 1.0e30, tau_max: 1.0e30, dt: 1.0e-4,
                lx, ly: 0.0, lz: 0.0, tilt_xy: 0.0, accumulate_torque: false,
            };
            gs.add_force_hook(Box::new(BeamBondForce::new(&gs, omega, &[rad, rad], &topo, cfg)));
            gs.eval_force_once();
            gs.download_force()[0]
        };
        // Periodic in x (L=1): atoms straddle the boundary → min-image separation 0.1.
        let fp = eval([[0.05, 0.0, 0.0], [0.95, 0.0, 0.0]], 1.0);
        // Non-periodic equivalent: the same 0.1 separation, no wrap.
        let fd = eval([[0.1, 0.0, 0.0], [0.0, 0.0, 0.0]], 0.0);
        let rel = |a: f32, b: f32| (a - b).abs() / b.abs().max(1.0);
        for c in 0..3 {
            assert!(rel(fp[c], fd[c]) < 1e-3, "periodic[{c}]={} != direct={}", fp[c], fd[c]);
        }
        assert!(fp[0].abs() > 1.0, "no force across the boundary: {fp:?}");
        eprintln!("periodic bond min-image: F_periodic={fp:?} F_direct={fd:?}");
    }

    #[test]
    fn beam_bending_generates_restoring_torque() {
        if GpuContext::new().is_none() { eprintln!("no GPU; skipping"); return; }
        // At rest length, give atom 0 a spin about z; bending should torque the
        // pair toward alignment and induce motion. Check angular response is
        // nonzero and linear momentum stays conserved.
        let (_p, v, w) = run([[-0.45, 0.0, 0.0], [0.45, 0.0, 0.0]], [[0.0; 3]; 2], [[0.0, 0.0, 5.0], [0.0; 3]], 0.9, 200);
        assert!((v[0][0] + v[1][0]).abs() < 1e-4, "momentum not conserved: {v:?}");
        // Atom 1 should have been spun up by the bond moment (was 0).
        assert!(w[1][2].abs() > 1e-6, "bending did not transmit moment: {w:?}");
    }
}
