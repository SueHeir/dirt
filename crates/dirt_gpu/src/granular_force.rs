//! DEM granular contact force as a `GpuForce` hook (DEM-specific; lives in dirt).
//!
//! This is the Hertz-normal + Mindlin-tangential contact law with persistent
//! per-contact spring history and rotational torque, expressed as a Force hook
//! that plugs into soil's resident loop (`soil_gpu::GpuState`). It binds the
//! resident buffers soil exposes — pos/vel/force, the cell-list outputs
//! (atom_cell/cell_start/sorted_atoms), and an auxiliary DOF for rotation
//! (state = angular velocity ω, rate = torque τ) — in group 0, plus the generic
//! `NeighborSlots` ping-pong (soil) in group 1 for the tangential spring history.
//!
//! Physics is identical to the standalone CPU `hertz_mindlin` / GPU monolith
//! (validated ~1e-6): i-centric, no atomics; each atom owns its half of every
//! spring and the two halves stay exact mirror images (n flips, v_rel flips), so
//! no canonical frame is needed. Soil never sees any of this — it only knows
//! "a force hook" and "an aux DOF". Walls are a separate hook (added next).

use bytemuck::{Pod, Zeroable};
use soil_gpu::{GpuForce, GpuState, Grid, NeighborSlots, SLOTS_WGSL};

/// Material / contact parameters (effective moduli already reduced for the pair).
#[derive(Clone, Copy, Debug)]
pub struct GranularConfig {
    pub e_eff: f32,
    pub beta: f32,
    pub g_eff: f32,
    pub mu: f32,
    pub dt: f32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct GParams {
    n: u32,
    nx: i32,
    ny: i32,
    nz: i32,
    e_eff: f32,
    beta: f32,
    g_eff: f32,
    mu: f32,
    dt: f32,
    _p0: f32,
    _p1: f32,
    _p2: f32,
}

/// Granular contact-force hook. Owns the contact kernel, its group-0 bind group
/// over the resident buffers, the `radius` buffer, and the `NeighborSlots`
/// history (group 1, ping-ponged each step inside `record`).
pub struct GranularForce {
    n: u32,
    pipeline: wgpu::ComputePipeline,
    bind_group: wgpu::BindGroup,
    slots: NeighborSlots,
}

impl GranularForce {
    /// Build the hook over `gs`'s resident buffers. `omega_aux` is the aux-DOF
    /// index registered on `gs` for rotation (state = ω, rate = τ). `radius` is
    /// per-atom. `grid` is the (fixed) cell-list grid `gs` was set up with.
    pub fn new(gs: &GpuState, grid: &Grid, omega_aux: usize, radius: &[f32], cfg: GranularConfig) -> Self {
        let ctx = gs.context();
        let device = &ctx.device;
        let n = gs.n();
        assert_eq!(radius.len(), n);

        let radius_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("gf_radius"), size: (n.max(1) * 4) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
        });
        ctx.queue.write_buffer(&radius_buf, 0, bytemuck::cast_slice(radius));

        let params = GParams {
            n: n as u32, nx: grid.n[0], ny: grid.n[1], nz: grid.n[2],
            e_eff: cfg.e_eff, beta: cfg.beta, g_eff: cfg.g_eff, mu: cfg.mu, dt: cfg.dt,
            _p0: 0.0, _p1: 0.0, _p2: 0.0,
        };
        let params_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("gf_params"), size: std::mem::size_of::<GParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
        });
        ctx.queue.write_buffer(&params_buf, 0, bytemuck::bytes_of(&params));

        let slots = NeighborSlots::new(ctx.clone(), n);
        slots.clear();

        // group 0: resident buffers (from gs) + radius + params; group 1: slots.
        let src = format!("{SLOTS_WGSL}\n{CONTACT_WGSL}");
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("granular_force"), source: wgpu::ShaderSource::Wgsl(src.into()),
        });
        let st = |binding, read_only| wgpu::BindGroupLayoutEntry {
            binding, visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only }, has_dynamic_offset: false, min_binding_size: None },
            count: None,
        };
        let g0 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("gf g0 bgl"),
            entries: &[
                st(0, true), st(1, true), st(2, false), st(3, true), st(4, false),
                st(5, true), st(6, true), st(7, true), st(8, true), st(9, true),
                wgpu::BindGroupLayoutEntry {
                    binding: 10, visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Uniform, has_dynamic_offset: false, min_binding_size: None },
                    count: None,
                },
            ],
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("gf g0 bg"), layout: &g0,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: gs.pos_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: gs.vel_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: gs.force_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: gs.aux_state_buffer(omega_aux).as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: gs.aux_rate_buffer(omega_aux).as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: radius_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: gs.inv_mass_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 7, resource: gs.atom_cell_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 8, resource: gs.cell_start_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 9, resource: gs.sorted_atoms_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 10, resource: params_buf.as_entire_binding() },
            ],
        });
        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("gf pl"), bind_group_layouts: &[Some(&g0), Some(slots.layout())], immediate_size: 0,
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("contact_force"), layout: Some(&layout), module: &shader, entry_point: Some("contact_force"),
            compilation_options: wgpu::PipelineCompilationOptions::default(), cache: None,
        });

        // radius_buf must outlive the bind group; the bind group holds a ref
        // (wgpu refcounts the buffer), so dropping the local handle is fine.
        let _ = radius_buf;
        let _ = params_buf;
        GranularForce { n: n as u32, pipeline, bind_group, slots }
    }
}

impl GpuForce for GranularForce {
    fn record(&self, pass: &mut wgpu::ComputePass) {
        pass.set_bind_group(0, &self.bind_group, &[]);
        pass.set_bind_group(1, self.slots.current_bind_group(), &[]);
        pass.set_pipeline(&self.pipeline);
        pass.dispatch_workgroups(self.n.div_ceil(64).max(1), 1, 1);
        // OLD/NEW history swap for the next step (interior-mutable, &self).
        self.slots.swap();
    }
}

const CONTACT_WGSL: &str = r#"
const SQRT_5_6: f32 = 0.9128709291752769;
const TANGENTIAL_EPSILON: f32 = 1.0e-10;

struct GParams {
    n: u32,
    nx: i32,
    ny: i32,
    nz: i32,
    e_eff: f32,
    beta: f32,
    g_eff: f32,
    mu: f32,
    dt: f32,
    _p0: f32,
    _p1: f32,
    _p2: f32,
};

@group(0) @binding(0)  var<storage, read>       pos: array<f32>;
@group(0) @binding(1)  var<storage, read>       vel: array<f32>;
@group(0) @binding(2)  var<storage, read_write> force_out: array<f32>;
@group(0) @binding(3)  var<storage, read>       omega: array<f32>;
@group(0) @binding(4)  var<storage, read_write> torque: array<f32>;
@group(0) @binding(5)  var<storage, read>       radius: array<f32>;
@group(0) @binding(6)  var<storage, read>       inv_mass: array<f32>;
@group(0) @binding(7)  var<storage, read>       atom_cell: array<u32>;
@group(0) @binding(8)  var<storage, read>       cell_start: array<u32>;
@group(0) @binding(9)  var<storage, read>       sorted_atoms: array<u32>;
@group(0) @binding(10) var<uniform>             params: GParams;

@compute @workgroup_size(64)
fn contact_force(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n) { return; }
    let bi = 3u * i;
    let xi = pos[bi];
    let yi = pos[bi + 1u];
    let zi = pos[bi + 2u];
    let vix = vel[bi];
    let viy = vel[bi + 1u];
    let viz = vel[bi + 2u];
    let wix = omega[bi];
    let wiy = omega[bi + 1u];
    let wiz = omega[bi + 2u];
    let ri = radius[i];
    let imi = inv_mass[i];

    let c = atom_cell[i];
    let cz0 = i32(c) % params.nz;
    let cy0 = (i32(c) / params.nz) % params.ny;
    let cx0 = i32(c) / (params.ny * params.nz);

    // Start from the gravity seed already written into force_out by seed_gravity;
    // torque is fully owned by this hook (no body torque), so start it at zero.
    var fx = force_out[bi];
    var fy = force_out[bi + 1u];
    var fz = force_out[bi + 2u];
    var tx = 0.0;
    var ty = 0.0;
    var tz = 0.0;

    var nslot: u32 = 0u;

    for (var dx = -1; dx <= 1; dx = dx + 1) {
        let cx = cx0 + dx;
        if (cx < 0 || cx >= params.nx) { continue; }
        for (var dy = -1; dy <= 1; dy = dy + 1) {
            let cy = cy0 + dy;
            if (cy < 0 || cy >= params.ny) { continue; }
            for (var dz = -1; dz <= 1; dz = dz + 1) {
                let cz = cz0 + dz;
                if (cz < 0 || cz >= params.nz) { continue; }
                let cc = u32((cx * params.ny + cy) * params.nz + cz);
                let start = cell_start[cc];
                let end = cell_start[cc + 1u];
                for (var m = start; m < end; m = m + 1u) {
                    let j = sorted_atoms[m];
                    if (j == i) { continue; }
                    let bj = 3u * j;
                    let ddx = pos[bj]      - xi;
                    let ddy = pos[bj + 1u] - yi;
                    let ddz = pos[bj + 2u] - zi;
                    let dist_sq = ddx * ddx + ddy * ddy + ddz * ddz;
                    let rj = radius[j];
                    let sum_r = ri + rj;
                    if (dist_sq >= sum_r * sum_r) { continue; }
                    let distance = sqrt(dist_sq);
                    if (distance == 0.0) { continue; }
                    let r_min = min(ri, rj);
                    let delta = min(sum_r - distance, 0.5 * r_min);
                    if (delta <= 0.0) { continue; }
                    let inv_dist = 1.0 / distance;
                    let nx = ddx * inv_dist;
                    let ny = ddy * inv_dist;
                    let nz = ddz * inv_dist;

                    let r_eff = (ri * rj) / sum_r;
                    let sdr = sqrt(delta * r_eff);
                    let s_n = 2.0 * params.e_eff * sdr;
                    let k_n = (4.0 / 3.0) * params.e_eff * sdr;
                    let k_t = 8.0 * params.g_eff * sdr;
                    let m_r = 1.0 / (imi + inv_mass[j]);

                    let r1n_x = ri * nx;
                    let r1n_y = ri * ny;
                    let r1n_z = ri * nz;
                    let vc_ix = vix + (wiy * r1n_z - wiz * r1n_y);
                    let vc_iy = viy + (wiz * r1n_x - wix * r1n_z);
                    let vc_iz = viz + (wix * r1n_y - wiy * r1n_x);
                    let wjx = omega[bj];
                    let wjy = omega[bj + 1u];
                    let wjz = omega[bj + 2u];
                    let r2n_x = rj * nx;
                    let r2n_y = rj * ny;
                    let r2n_z = rj * nz;
                    let vc_jx = vel[bj]      + (-wjy * r2n_z + wjz * r2n_y);
                    let vc_jy = vel[bj + 1u] + (-wjz * r2n_x + wjx * r2n_z);
                    let vc_jz = vel[bj + 2u] + (-wjx * r2n_y + wjy * r2n_x);

                    let vr_x = vc_jx - vc_ix;
                    let vr_y = vc_jy - vc_iy;
                    let vr_z = vc_jz - vc_iz;
                    let v_n = vr_x * nx + vr_y * ny + vr_z * nz;

                    let f_diss_n = 2.0 * params.beta * SQRT_5_6 * sqrt(s_n * m_r) * v_n;
                    let f_n_mag = max(k_n * delta - f_diss_n, 0.0);
                    fx = fx - f_n_mag * nx;
                    fy = fy - f_n_mag * ny;
                    fz = fz - f_n_mag * nz;

                    let vt_x = vr_x - v_n * nx;
                    let vt_y = vr_y - v_n * ny;
                    let vt_z = vr_z - v_n * nz;

                    var s = slot_lookup(i, j);
                    let s_dot_n = s.x * nx + s.y * ny + s.z * nz;
                    s.x = s.x - s_dot_n * nx + vt_x * params.dt;
                    s.y = s.y - s_dot_n * ny + vt_y * params.dt;
                    s.z = s.z - s_dot_n * nz + vt_z * params.dt;

                    let f_t_max = params.mu * abs(f_n_mag);

                    let s_mag = sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
                    let f_t_spring_mag = k_t * s_mag;
                    if (f_t_spring_mag > f_t_max && f_t_spring_mag > TANGENTIAL_EPSILON) {
                        let scale = f_t_max / f_t_spring_mag;
                        s.x = s.x * scale;
                        s.y = s.y * scale;
                        s.z = s.z * scale;
                    }

                    let gamma_t = 2.0 * SQRT_5_6 * params.beta * sqrt(k_t * m_r);
                    var ft_x = k_t * s.x + gamma_t * vt_x;
                    var ft_y = k_t * s.y + gamma_t * vt_y;
                    var ft_z = k_t * s.z + gamma_t * vt_z;
                    let f_t_mag = sqrt(ft_x * ft_x + ft_y * ft_y + ft_z * ft_z);
                    if (f_t_mag > f_t_max && f_t_mag > TANGENTIAL_EPSILON) {
                        let scale = f_t_max / f_t_mag;
                        ft_x = ft_x * scale;
                        ft_y = ft_y * scale;
                        ft_z = ft_z * scale;
                    }

                    fx = fx + ft_x;
                    fy = fy + ft_y;
                    fz = fz + ft_z;
                    tx = tx + (r1n_y * ft_z - r1n_z * ft_y);
                    ty = ty + (r1n_z * ft_x - r1n_x * ft_z);
                    tz = tz + (r1n_x * ft_y - r1n_y * ft_x);

                    if (nslot < SLOT_MAX) {
                        slot_write(i, nslot, j, s);
                        nslot = nslot + 1u;
                    }
                }
            }
        }
    }

    slot_clear_from(i, nslot);

    force_out[bi]      = fx;
    force_out[bi + 1u] = fy;
    force_out[bi + 2u] = fz;
    torque[bi]      = tx;
    torque[bi + 1u] = ty;
    torque[bi + 2u] = tz;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;
    use soil_gpu::GpuContext;

    #[test]
    fn granular_hook_two_spheres_repel_symmetric() {
        let Some(ctx) = GpuContext::new() else {
            eprintln!("no GPU adapter; skipping");
            return;
        };
        // Two equal unit-density spheres (r = 0.5) overlapping head-on by 0.1.
        let r = 0.5f32;
        let mass = (4.0 / 3.0) * std::f32::consts::PI * r * r * r;
        let inertia = 0.4 * mass * r * r;
        let inv_mass = vec![1.0 / mass, 1.0 / mass];
        let inv_inertia = vec![1.0 / inertia, 1.0 / inertia];
        let radius = vec![r, r];

        let p0 = [[-0.45f32, 0.0, 0.0], [0.45, 0.0, 0.0]]; // separation 0.9, overlap 0.1
        let grid = Grid::from_positions(&p0, 2.0 * r); // cutoff = sum of radii
        let mut gs = GpuState::new(ctx, 2, grid.total_cells);
        let dt = 1.0e-5f32;
        gs.set_params(dt, [0.0, 0.0, 0.0]); // no gravity
        gs.set_state(&p0, &[[0.0; 3]; 2], &inv_mass, grid);

        let omega = gs.add_aux_dof();
        gs.set_aux_inv_coeff(omega, &inv_inertia);
        gs.set_aux_state(omega, &[[0.0; 3]; 2]);

        let cfg = GranularConfig { e_eff: 1.0e6, beta: 0.1, g_eff: 4.0e5, mu: 0.5, dt };
        gs.add_force_hook(Box::new(GranularForce::new(&gs, &grid, omega, &radius, cfg)));

        gs.run_steps(3000);
        let p = gs.download_pos();
        let v = gs.download_vel();

        // Normal repulsion pushes them apart, symmetric about the origin, momentum
        // conserved (COM at rest). Head-on, so motion stays on the x-axis.
        assert!(p[0][0] < -0.45 && p[1][0] > 0.45, "did not separate: {p:?}");
        assert!((p[0][0] + p[1][0]).abs() < 1e-4, "asymmetric: {p:?}");
        assert!((v[0][0] + v[1][0]).abs() < 1e-4, "momentum not conserved: {v:?}");
        for a in 0..2 {
            assert!(p[a][1].abs() < 1e-5 && p[a][2].abs() < 1e-5, "off-axis: {p:?}");
        }
        eprintln!("granular hook: x=[{}, {}] v=[{}, {}] (Hertz repulsion, symmetric)",
            p[0][0], p[1][0], v[0][0], v[1][0]);
    }
}
