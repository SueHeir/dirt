//! DEM planar-wall contact response as a `GpuForce` hook (DEM-specific; dirt).
//!
//! Sphere-plane Hertz normal + Mindlin tangential (Coulomb-capped, friction
//! opposing motion) against a set of planar walls. Registered AFTER
//! [`GranularForce`](crate::GranularForce), it read-modify-writes the resident
//! `force`/`torque` (i-centric, each thread owns atom i — no atomics, no race),
//! adding the wall contribution on top of inter-particle contacts.
//!
//! Geometry (the plane set, signed distance) is generic and comes from soil
//! ([`soil_gpu::Boundary`] / `BOUNDARY_WGSL`); the *response* (Hertz/Mindlin
//! against the plane) is DEM physics and lives here. Per-(particle, wall)
//! tangential spring history is owned in place by thread i — no ping-pong needed.

use bytemuck::{Pod, Zeroable};
use soil_gpu::{Boundary, GpuContext, GpuForce, GpuState, BOUNDARY_WGSL};

/// Maximum planar walls (matches the fixed WGSL uniform array size).
pub const MAX_WALLS: usize = 8;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct WallParams {
    n: u32,
    n_walls: u32,
    _p0: u32,
    _p1: u32,
    e_eff: f32,
    beta: f32,
    g_eff: f32,
    mu: f32,
    dt: f32,
    advance_hist: f32,
    _p3: f32,
    _p4: f32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Default)]
struct GpuPlane {
    point: [f32; 4],
    normal: [f32; 4],
}

/// DEM wall-contact hook. Owns its kernel, group-0 bind group over the resident
/// buffers + radius + the wall uniform + the per-(particle,wall) spring history.
pub struct WallForce {
    ctx: GpuContext,
    n: u32,
    n_walls: usize,
    pipeline: wgpu::ComputePipeline,
    bind_group: wgpu::BindGroup,
    bind_group_prime: wgpu::BindGroup,
    walls_buf: wgpu::Buffer,
}

impl WallForce {
    /// Build over `gs`'s resident buffers. `omega_aux` is the rotation aux-DOF
    /// (state = ω, rate = τ) registered on `gs` — the same index passed to
    /// [`GranularForce`](crate::GranularForce). `boundary` is the plane set.
    pub fn new(
        gs: &GpuState,
        omega_aux: usize,
        radius: &[f32],
        boundary: &Boundary,
        e_eff: f32,
        beta: f32,
        g_eff: f32,
        mu: f32,
        dt: f32,
    ) -> Self {
        let ctx = gs.context();
        let device = &ctx.device;
        let n = gs.n();
        assert_eq!(radius.len(), n);
        assert!(boundary.len() <= MAX_WALLS, "at most {MAX_WALLS} walls supported");

        let radius_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("wf_radius"), size: (n.max(1) * 4) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
        });
        ctx.queue.write_buffer(&radius_buf, 0, bytemuck::cast_slice(radius));

        // Pack the boundary into a fixed-size uniform array (unused slots zeroed).
        let planes = pack_planes(boundary);
        let walls_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("wf_walls"), size: std::mem::size_of::<[GpuPlane; MAX_WALLS]>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
        });
        ctx.queue.write_buffer(&walls_buf, 0, bytemuck::cast_slice(&planes));

        let params = WallParams {
            n: n as u32, n_walls: boundary.len() as u32, _p0: 0, _p1: 0,
            e_eff, beta, g_eff, mu, dt, advance_hist: 1.0, _p3: 0.0, _p4: 0.0,
        };
        let params_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("wf_params"), size: std::mem::size_of::<WallParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
        });
        ctx.queue.write_buffer(&params_buf, 0, bytemuck::bytes_of(&params));

        // History-neutral params for the entry prime (advance_hist = 0).
        let params_prime = WallParams { advance_hist: 0.0, ..params };
        let params_buf_prime = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("wf_params_prime"), size: std::mem::size_of::<WallParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
        });
        ctx.queue.write_buffer(&params_buf_prime, 0, bytemuck::bytes_of(&params_prime));

        // Per-(particle, wall) tangential spring history, owned in place; zeroed.
        let wsprings = vec![0.0f32; n.max(1) * MAX_WALLS * 3];
        let wall_spring = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("wf_wall_spring"), size: (wsprings.len() * 4) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC, mapped_at_creation: false,
        });
        ctx.queue.write_buffer(&wall_spring, 0, bytemuck::cast_slice(&wsprings));

        let src = format!("{BOUNDARY_WGSL}\n{WALL_WGSL}");
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("wall_force"), source: wgpu::ShaderSource::Wgsl(src.into()),
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
        let g0 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("wf g0 bgl"),
            entries: &[
                st(0, true), st(1, true), st(2, false), st(3, true), st(4, false),
                st(5, true), st(6, true), uni(7), uni(8), st(9, false),
            ],
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("wf g0 bg"), layout: &g0,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: gs.pos_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: gs.vel_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: gs.force_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: gs.aux_state_buffer(omega_aux).as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: gs.aux_rate_buffer(omega_aux).as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: radius_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: gs.inv_mass_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 7, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 8, resource: walls_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 9, resource: wall_spring.as_entire_binding() },
            ],
        });
        // Same bindings, params → advance_hist=0 buffer, for the prime.
        let bind_group_prime = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("wf g0 bg prime"), layout: &g0,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: gs.pos_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: gs.vel_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: gs.force_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: gs.aux_state_buffer(omega_aux).as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: gs.aux_rate_buffer(omega_aux).as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: radius_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: gs.inv_mass_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 7, resource: params_buf_prime.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 8, resource: walls_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 9, resource: wall_spring.as_entire_binding() },
            ],
        });
        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("wf pl"), bind_group_layouts: &[Some(&g0)], immediate_size: 0,
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("wall_force"), layout: Some(&layout), module: &shader, entry_point: Some("wall_force"),
            compilation_options: wgpu::PipelineCompilationOptions::default(), cache: None,
        });

        WallForce {
            ctx: ctx.clone(), n: n as u32, n_walls: boundary.len(),
            pipeline, bind_group, bind_group_prime, walls_buf,
        }
    }

    /// Update the wall plane geometry in place (e.g. after a `wall_move` step),
    /// without rebuilding the pipeline/bind group. The wall *count* must be
    /// unchanged (only positions/normals move); rebuild via `new` to change count.
    pub fn set_walls(&self, boundary: &Boundary) {
        assert_eq!(
            boundary.len(), self.n_walls,
            "set_walls: wall count changed ({} -> {}); rebuild WallForce instead",
            self.n_walls, boundary.len(),
        );
        let planes = pack_planes(boundary);
        self.ctx.queue.write_buffer(&self.walls_buf, 0, bytemuck::cast_slice(&planes));
    }
}

/// Pack a boundary's planes into the fixed-size uniform array (unused slots 0).
fn pack_planes(boundary: &Boundary) -> [GpuPlane; MAX_WALLS] {
    let mut planes = [GpuPlane::default(); MAX_WALLS];
    for (slot, p) in planes.iter_mut().zip(boundary.to_gpu().iter()) {
        slot.point = p.point;
        slot.normal = p.normal;
    }
    planes
}

impl GpuForce for WallForce {
    fn record(&self, pass: &mut wgpu::ComputePass) {
        pass.set_bind_group(0, &self.bind_group, &[]);
        pass.set_pipeline(&self.pipeline);
        pass.dispatch_workgroups(self.n.div_ceil(64).max(1), 1, 1);
    }

    /// History-neutral prime: evaluate the wall contact force from the current
    /// in-place spring without integrating or writing it back (advance_hist=0),
    /// so re-priming at window boundaries doesn't mutate wall history.
    fn record_prime(&self, pass: &mut wgpu::ComputePass) {
        pass.set_bind_group(0, &self.bind_group_prime, &[]);
        pass.set_pipeline(&self.pipeline);
        pass.dispatch_workgroups(self.n.div_ceil(64).max(1), 1, 1);
    }
}

const WALL_WGSL: &str = r#"
const SQRT_5_6: f32 = 0.9128709291752769;
const TANGENTIAL_EPSILON: f32 = 1.0e-10;
const MAX_WALLS: u32 = 8u;

struct WallParams {
    n: u32,
    n_walls: u32,
    _p0: u32,
    _p1: u32,
    e_eff: f32,
    beta: f32,
    g_eff: f32,
    mu: f32,
    dt: f32,
    advance_hist: f32,
    _p3: f32,
    _p4: f32,
};

struct WallArray { walls: array<BoundaryPlane, 8> };

@group(0) @binding(0) var<storage, read>       pos: array<f32>;
@group(0) @binding(1) var<storage, read>       vel: array<f32>;
@group(0) @binding(2) var<storage, read_write> force_out: array<f32>;
@group(0) @binding(3) var<storage, read>       omega: array<f32>;
@group(0) @binding(4) var<storage, read_write> torque: array<f32>;
@group(0) @binding(5) var<storage, read>       radius: array<f32>;
@group(0) @binding(6) var<storage, read>       inv_mass: array<f32>;
@group(0) @binding(7) var<uniform>             params: WallParams;
@group(0) @binding(8) var<uniform>             walls: WallArray;
@group(0) @binding(9) var<storage, read_write> wall_spring: array<f32>;

@compute @workgroup_size(64)
fn wall_force(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n) { return; }
    let bi = 3u * i;
    let px = pos[bi];
    let py = pos[bi + 1u];
    let pz = pos[bi + 2u];
    let vix = vel[bi];
    let viy = vel[bi + 1u];
    let viz = vel[bi + 2u];
    let wix = omega[bi];
    let wiy = omega[bi + 1u];
    let wiz = omega[bi + 2u];
    let ri = radius[i];
    let imi = inv_mass[i];
    var m_i = 0.0;
    if (imi > 0.0) { m_i = 1.0 / imi; }

    var fx = force_out[bi];
    var fy = force_out[bi + 1u];
    var fz = force_out[bi + 2u];
    var tx = torque[bi];
    var ty = torque[bi + 1u];
    var tz = torque[bi + 2u];

    for (var w: u32 = 0u; w < params.n_walls && w < MAX_WALLS; w = w + 1u) {
        let pt = walls.walls[w].point;
        let nrm = walls.walls[w].normal;
        let nx = nrm.x;
        let ny = nrm.y;
        let nz = nrm.z;
        let signed_dist = (px - pt.x) * nx + (py - pt.y) * ny + (pz - pt.z) * nz;
        let overlap = ri - signed_dist;

        let sb = 3u * (i * MAX_WALLS + w);
        if (overlap <= 0.0) {
            wall_spring[sb] = 0.0;
            wall_spring[sb + 1u] = 0.0;
            wall_spring[sb + 2u] = 0.0;
            continue;
        }
        let delta = min(overlap, 0.5 * ri);

        let r_eff = ri;
        let sdr = sqrt(delta * r_eff);
        let s_n = 2.0 * params.e_eff * sdr;
        let k_n = (4.0 / 3.0) * params.e_eff * sdr;
        let k_t = 8.0 * params.g_eff * sdr;
        let m_r = m_i;

        let oxn_x = wiy * nz - wiz * ny;
        let oxn_y = wiz * nx - wix * nz;
        let oxn_z = wix * ny - wiy * nx;
        let vs_x = vix - ri * oxn_x;
        let vs_y = viy - ri * oxn_y;
        let vs_z = viz - ri * oxn_z;
        let v_n = vs_x * nx + vs_y * ny + vs_z * nz;

        let f_diss_n = 2.0 * params.beta * SQRT_5_6 * sqrt(s_n * m_r) * v_n;
        let f_n_mag = max(k_n * delta - f_diss_n, 0.0);
        fx = fx + f_n_mag * nx;
        fy = fy + f_n_mag * ny;
        fz = fz + f_n_mag * nz;

        let vt_x = vs_x - v_n * nx;
        let vt_y = vs_y - v_n * ny;
        let vt_z = vs_z - v_n * nz;

        var s = vec3<f32>(wall_spring[sb], wall_spring[sb + 1u], wall_spring[sb + 2u]);
        let s_dot_n = s.x * nx + s.y * ny + s.z * nz;
        // advance_hist = 1 in the step loop, 0 in the entry prime: the prime
        // reprojects but does NOT integrate the wall spring (no +vt·dt).
        s.x = s.x - s_dot_n * nx + vt_x * params.dt * params.advance_hist;
        s.y = s.y - s_dot_n * ny + vt_y * params.dt * params.advance_hist;
        s.z = s.z - s_dot_n * nz + vt_z * params.dt * params.advance_hist;

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
        var ft_x = -(k_t * s.x + gamma_t * vt_x);
        var ft_y = -(k_t * s.y + gamma_t * vt_y);
        var ft_z = -(k_t * s.z + gamma_t * vt_z);
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
        let nxf_x = ny * ft_z - nz * ft_y;
        let nxf_y = nz * ft_x - nx * ft_z;
        let nxf_z = nx * ft_y - ny * ft_x;
        tx = tx - ri * nxf_x;
        ty = ty - ri * nxf_y;
        tz = tz - ri * nxf_z;

        // The wall spring is stored in-place (no ping-pong), so the prime must
        // NOT write it back — otherwise the entry prime mutates history.
        if (params.advance_hist != 0.0) {
            wall_spring[sb] = s.x;
            wall_spring[sb + 1u] = s.y;
            wall_spring[sb + 2u] = s.z;
        }
    }

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
    use soil_gpu::{Grid, GpuContext, Plane};
    use crate::{GranularConfig, GranularForce};

    #[test]
    fn wall_hook_particle_settles_on_floor() {
        let Some(ctx) = GpuContext::new() else {
            eprintln!("no GPU adapter; skipping");
            return;
        };
        // One sphere dropped onto a floor at z=0 (normal +z); should settle near
        // resting contact (z ~ r) under gravity + wall repulsion, not fall through.
        let r = 0.5f32;
        let mass = (4.0 / 3.0) * std::f32::consts::PI * r * r * r;
        let inertia = 0.4 * mass * r * r;
        let p0 = [[0.0f32, 0.0, 0.7]]; // starts above the floor
        let grid = Grid::from_positions(&[[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]], 2.0 * r);
        let mut gs = GpuState::new(ctx, 1, grid.total_cells);
        let dt = 1.0e-5f32;
        gs.set_params(dt, [0.0, 0.0, -9.81]);
        gs.set_state(&p0, &[[0.0; 3]], &[1.0 / mass], grid);

        let omega = gs.add_aux_dof();
        gs.set_aux_inv_coeff(omega, &[1.0 / inertia]);
        gs.set_aux_state(omega, &[[0.0; 3]]);

        let radius = vec![r];
        let cfg = GranularConfig { e_eff: 1.0e6, beta: 0.3, g_eff: 4.0e5, mu: 0.5, dt };
        gs.add_force_hook(Box::new(GranularForce::new(&gs, &grid, omega, &radius, cfg)));

        let mut boundary = Boundary::new();
        boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])); // floor
        gs.add_force_hook(Box::new(WallForce::new(&gs, omega, &radius, &boundary, 1.0e6, 0.3, 4.0e5, 0.5, dt)));

        gs.run_steps(20000);
        let p = gs.download_pos()[0];

        // Did not tunnel through the floor, and rests near z = r (small overlap).
        assert!(p[2] > 0.3 && p[2] < 0.6, "did not settle on floor: z={}", p[2]);
        assert!(p[0].abs() < 1e-3 && p[1].abs() < 1e-3, "drifted laterally: {p:?}");
        eprintln!("wall hook: settled at z={} (floor contact)", p[2]);
    }
}
