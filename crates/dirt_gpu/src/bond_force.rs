//! Persistent bonded-pair force as a `GpuForce` hook (DEM/peridynamics; dirt).
//!
//! Unlike the contact force (transient neighbours re-found from the cell list each
//! step), bonds are a **persistent topology**: a fixed per-point family built once,
//! changed only by breakage. So the per-bond state lives in plain buffers and is
//! updated **in place** — no OLD/NEW ping-pong. Layout is CSR (`bond_offsets[i]
//! .. bond_offsets[i+1]` indexes point `i`'s bonds in the flat `bond_partner` /
//! `bond_r0` / `bond_broken` arrays), mirroring the CPU `BondStore`'s per-local-atom
//! `Vec<Vec<BondEntry>>`.
//!
//! Accumulation is i-centric and atomic-free, exactly like the contact kernel: a
//! bond `(i,j)` is double-stored (point `i` owns a copy, point `j` owns a copy),
//! and each thread computes its own endpoint's force as the exact mirror image of
//! the other — so the two halves stay consistent with no canonical frame and no
//! atomics. Breakage is identical on both copies (same geometry, same per-bond
//! threshold → same f32 decision), so they stay in sync.
//!
//! Phase A (this file): the central-force / bond-based-peridynamics law —
//! `f = k_n·(L − r0)` along the bond, break at a critical |strain|. The full beam
//! law (shear/twist/bend moments, beam-stress breakage, plasticity) is Phase B.

use bytemuck::{Pod, Zeroable};
use soil_gpu::{GpuForce, GpuState};

/// Central-bond parameters (Phase A). `break_strain` is the critical |L−r0|/r0; set
/// it large (e.g. `f32::INFINITY` is not allowed in a uniform — use 1e30) for an
/// effectively unbreakable bond.
#[derive(Clone, Copy, Debug)]
pub struct BondConfig {
    pub k_n: f32,
    pub break_strain: f32,
    pub dt: f32,
}

/// Host-side persistent bond topology (CSR). `offsets.len() == n+1`; `partner` and
/// `r0` are flat, indexed by `offsets[i]..offsets[i+1]` for point `i`. `partner`
/// holds the partner's **local atom index** (resolved from global tag at build;
/// re-resolution on migration is a later phase).
#[derive(Clone, Debug, Default)]
pub struct BondTopology {
    pub offsets: Vec<u32>,
    pub partner: Vec<u32>,
    pub r0: Vec<f32>,
}

impl BondTopology {
    pub fn num_bonds(&self) -> usize {
        self.partner.len()
    }

    /// Build a persistent GPU bond topology from the CPU `BondStore`, resolving each
    /// bond's `partner_tag` to its current local/ghost atom index. `offsets` spans
    /// all `atoms.len()` GPU atoms (ghosts own no bonds). Bonds whose partner is not
    /// present (cut by an MPI ghost boundary) are dropped — same `missing_partner`
    /// semantics as the CPU. The flat order matches `BondStore.bonds[i]` per point so
    /// per-bond GPU state lines up with a host download for migration.
    pub fn from_bond_store(bonds: &soil_core::BondStore, atoms: &soil_core::Atom) -> Self {
        let nall = atoms.len();
        let nlocal = atoms.nlocal as usize;
        let mut tag_to_index = std::collections::HashMap::with_capacity(nall);
        for idx in 0..nall {
            tag_to_index.insert(atoms.tag[idx], idx);
        }
        let mut offsets = Vec::with_capacity(nall + 1);
        let mut partner = Vec::new();
        let mut r0 = Vec::new();
        offsets.push(0u32);
        for i in 0..nall {
            if i < nlocal && i < bonds.bonds.len() {
                for e in &bonds.bonds[i] {
                    if let Some(&j) = tag_to_index.get(&e.partner_tag) {
                        partner.push(j as u32);
                        r0.push(e.r0 as f32);
                    }
                }
            }
            offsets.push(partner.len() as u32);
        }
        BondTopology { offsets, partner, r0 }
    }
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct BParams {
    n: u32,
    k_n: f32,
    break_strain: f32,
    dt: f32,
}

/// Persistent bonded-pair force hook. Owns the bond CSR buffers (group 1) and a
/// small group-0 bind over the resident pos/force buffers + params.
pub struct BondForce {
    n: u32,
    pipeline: wgpu::ComputePipeline,
    g0: wgpu::BindGroup,
    g1: wgpu::BindGroup,
    broken_buf: wgpu::Buffer,
    staging: wgpu::Buffer,
    num_bonds: usize,
}

impl BondForce {
    /// Build the hook over `gs`'s resident buffers from a persistent `BondTopology`.
    pub fn new(gs: &GpuState, topo: &BondTopology, cfg: BondConfig) -> Self {
        let ctx = gs.context();
        let device = &ctx.device;
        let n = gs.n();
        assert_eq!(topo.offsets.len(), n + 1, "BondTopology.offsets must be n+1");
        let num_bonds = topo.num_bonds();
        assert_eq!(topo.r0.len(), num_bonds);

        let storage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST;
        let mk_u32 = |label: &str, data: &[u32]| {
            let bytes = (data.len().max(1) * 4) as u64;
            let buf = device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label), size: bytes,
                usage: storage | wgpu::BufferUsages::COPY_SRC, mapped_at_creation: false,
            });
            if !data.is_empty() {
                ctx.queue.write_buffer(&buf, 0, bytemuck::cast_slice(data));
            }
            buf
        };
        let mk_f32 = |label: &str, data: &[f32]| {
            let bytes = (data.len().max(1) * 4) as u64;
            let buf = device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label), size: bytes, usage: storage, mapped_at_creation: false,
            });
            if !data.is_empty() {
                ctx.queue.write_buffer(&buf, 0, bytemuck::cast_slice(data));
            }
            buf
        };

        let offsets_buf = mk_u32("bond_offsets", &topo.offsets);
        let partner_buf = mk_u32("bond_partner", &topo.partner);
        let r0_buf = mk_f32("bond_r0", &topo.r0);
        let broken_buf = mk_u32("bond_broken", &vec![0u32; num_bonds]);

        let staging = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("bond_broken_staging"), size: (num_bonds.max(1) * 4) as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
        });

        let params = BParams { n: n as u32, k_n: cfg.k_n, break_strain: cfg.break_strain, dt: cfg.dt };
        let params_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("bond_params"), size: std::mem::size_of::<BParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false,
        });
        ctx.queue.write_buffer(&params_buf, 0, bytemuck::bytes_of(&params));

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("bond_force"), source: wgpu::ShaderSource::Wgsl(BOND_WGSL.into()),
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
            label: Some("bond g0 bgl"), entries: &[st(0, true), st(1, false), uni(2)],
        });
        let g0 = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("bond g0 bg"), layout: &g0l,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: gs.pos_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: gs.force_buffer().as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: params_buf.as_entire_binding() },
            ],
        });

        let g1l = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("bond g1 bgl"), entries: &[st(0, true), st(1, true), st(2, true), st(3, false)],
        });
        let g1 = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("bond g1 bg"), layout: &g1l,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: offsets_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: partner_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: r0_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: broken_buf.as_entire_binding() },
            ],
        });

        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("bond pl"), bind_group_layouts: &[Some(&g0l), Some(&g1l)], immediate_size: 0,
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("bond_force"), layout: Some(&layout), module: &shader, entry_point: Some("bond_force"),
            compilation_options: wgpu::PipelineCompilationOptions::default(), cache: None,
        });

        // buffers are refcounted by the bind groups; drop local handles.
        let _ = (offsets_buf, partner_buf, r0_buf, params_buf);
        BondForce { n: n as u32, pipeline, g0, g1, broken_buf, staging, num_bonds }
    }

    /// Download the per-bond broken flags (1 = broken). For validation.
    pub fn download_broken(&self, gs: &GpuState) -> Vec<u32> {
        if self.num_bonds == 0 {
            return Vec::new();
        }
        let ctx = gs.context();
        let bytes = (self.num_bonds * 4) as u64;
        let mut enc = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("bond dl") });
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

impl GpuForce for BondForce {
    fn record(&self, pass: &mut wgpu::ComputePass) {
        pass.set_bind_group(0, &self.g0, &[]);
        pass.set_bind_group(1, &self.g1, &[]);
        pass.set_pipeline(&self.pipeline);
        pass.dispatch_workgroups(self.n.div_ceil(64).max(1), 1, 1);
    }
    // record_prime: bonds carry no transient history that double-advances across a
    // window boundary (central force is a pure function of geometry), so the
    // default (== record) is correct. Breakage is monotonic and idempotent here.
}

const BOND_WGSL: &str = r#"
struct BParams { n: u32, k_n: f32, break_strain: f32, dt: f32 };

@group(0) @binding(0) var<storage, read>       pos: array<f32>;
@group(0) @binding(1) var<storage, read_write> force_out: array<f32>;
@group(0) @binding(2) var<uniform>             params: BParams;

@group(1) @binding(0) var<storage, read>       bond_offsets: array<u32>;
@group(1) @binding(1) var<storage, read>       bond_partner: array<u32>;
@group(1) @binding(2) var<storage, read>       bond_r0: array<f32>;
@group(1) @binding(3) var<storage, read_write> bond_broken: array<u32>;

@compute @workgroup_size(64)
fn bond_force(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.n) { return; }
    let bi = 3u * i;
    let xi = pos[bi];
    let yi = pos[bi + 1u];
    let zi = pos[bi + 2u];

    // Accumulate onto the gravity seed (read-modify-write, like the contact hook).
    var fx = force_out[bi];
    var fy = force_out[bi + 1u];
    var fz = force_out[bi + 2u];

    let start = bond_offsets[i];
    let end = bond_offsets[i + 1u];
    for (var k = start; k < end; k = k + 1u) {
        if (bond_broken[k] != 0u) { continue; }
        let j = bond_partner[k];
        let bj = 3u * j;
        let dx = pos[bj]      - xi;
        let dy = pos[bj + 1u] - yi;
        let dz = pos[bj + 2u] - zi;
        let len = sqrt(dx * dx + dy * dy + dz * dz);
        if (len == 0.0) { continue; }
        let r0 = bond_r0[k];
        let stretch = len - r0;

        // Breakage: critical |strain|. Both endpoints see the same len/r0/threshold
        // → the same f32 decision, so the two copies break together.
        if (abs(stretch) > params.break_strain * r0) {
            bond_broken[k] = 1u;
            continue;
        }

        let inv = 1.0 / len;
        let nx = dx * inv;
        let ny = dy * inv;
        let nz = dz * inv;
        // Central force on i: +k·(L−r0)·n̂ (toward j when stretched). The partner's
        // own thread computes the exact mirror (n flips) → equal and opposite.
        let f = params.k_n * stretch;
        fx = fx + f * nx;
        fy = fy + f * ny;
        fz = fz + f * nz;
    }

    force_out[bi]      = fx;
    force_out[bi + 1u] = fy;
    force_out[bi + 2u] = fz;
}
"#;

#[cfg(test)]
mod tests {
    use super::*;
    use soil_gpu::{GpuContext, Grid};

    /// Two equal masses joined by one central bond, initially stretched. With a
    /// finite break_strain above the initial strain the bond restores (the pair
    /// moves together, symmetric, momentum conserved); below it the bond breaks on
    /// step 1 and the pair stays put.
    fn run_pair(break_strain: f32) -> (Vec<[f32; 3]>, Vec<[f32; 3]>, Vec<u32>) {
        let ctx = GpuContext::new().expect("gpu");
        let r0 = 0.8f32;
        let half = 0.45f32; // positions ±0.45 → L=0.9, stretch 0.1, strain 0.125
        let mass = 1.0f32;
        let inv_mass = vec![1.0 / mass, 1.0 / mass];
        let p0 = [[-half, 0.0, 0.0], [half, 0.0, 0.0]];
        let grid = Grid::from_positions(&p0, 2.0 * r0);
        let dt = 1.0e-4f32;
        let mut gs = GpuState::new(ctx, 2, grid.total_cells);
        gs.set_params(dt, [0.0, 0.0, 0.0]);
        gs.set_state(&p0, &[[0.0; 3]; 2], &inv_mass, grid);

        // bond stored on both endpoints
        let topo = BondTopology { offsets: vec![0, 1, 2], partner: vec![1, 0], r0: vec![r0, r0] };
        let cfg = BondConfig { k_n: 1.0e4, break_strain, dt };
        let bond = BondForce::new(&gs, &topo, cfg);
        let handle = std::rc::Rc::new(bond);
        gs.add_force_hook(Box::new(HookRef(std::rc::Rc::clone(&handle))));
        gs.run_steps(500);
        let p = gs.download_pos();
        let v = gs.download_vel();
        let broken = handle.download_broken(&gs);
        (p, v, broken)
    }

    // Local newtype so we can register the hook while keeping a handle for
    // `download_broken` (orphan rule forbids `impl GpuForce for Rc<BondForce>`).
    struct HookRef(std::rc::Rc<BondForce>);
    impl GpuForce for HookRef {
        fn record(&self, pass: &mut wgpu::ComputePass) { self.0.record(pass); }
    }

    #[test]
    fn bond_restores_when_intact() {
        if GpuContext::new().is_none() { eprintln!("no GPU; skipping"); return; }
        let (p, v, broken) = run_pair(0.5); // 0.5 > 0.125 → intact
        assert_eq!(broken, vec![0, 0], "bond should be intact");
        // Stretched bond pulls them together, symmetric, momentum conserved.
        assert!(p[0][0] > -0.45 && p[1][0] < 0.45, "did not contract: {p:?}");
        assert!((p[0][0] + p[1][0]).abs() < 1e-4, "asymmetric: {p:?}");
        assert!((v[0][0] + v[1][0]).abs() < 1e-4, "momentum not conserved: {v:?}");
        eprintln!("bond intact: x=[{}, {}]", p[0][0], p[1][0]);
    }

    #[test]
    fn bond_breaks_past_threshold() {
        if GpuContext::new().is_none() { eprintln!("no GPU; skipping"); return; }
        let (p, _v, broken) = run_pair(0.05); // 0.05 < 0.125 → breaks step 1
        assert_eq!(broken, vec![1, 1], "bond should be broken on both copies");
        // Broken → no force → no gravity → pair stays at ±0.45.
        assert!((p[0][0] + 0.45).abs() < 1e-4 && (p[1][0] - 0.45).abs() < 1e-4, "moved after break: {p:?}");
        eprintln!("bond broken: x=[{}, {}]", p[0][0], p[1][0]);
    }
}
