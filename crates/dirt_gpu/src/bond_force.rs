//! Beam bonds as a [`GpuForce`] hook: normal, history shear, twist and bending, with damping and
//! the `combined_stress` breakage criterion.
//!
//! A bond is stored from both ends, and each end computes its own half as the exact mirror of
//! the other: the axis, relative velocity, relative spin and histories all flip. So each end
//! applies force and torque only to itself, with no atomics, and the result is the CPU's
//! `+f on i, -f on j`, `+M on i, -M on j`, with the same shear torque on both. Both ends reach the
//! same breakage decision, from mirrored arithmetic. Register after
//! [`GranularForce`](crate::GranularForce), which is what starts each torque at zero.

use std::cell::RefCell;
use std::rc::Rc;

use cubecl::prelude::*;
use grass_gpu::cubecl;
use soil_gpu::neighbor_slots::SENTINEL;
use soil_gpu::{DeviceColumn, GpuDevice, GpuForce, GpuReal, GpuState, StepColumns};

/// Fixed per-bond values, in table order.
pub const BOND_FIELDS: usize = 15;
const R0: usize = 0;
const K_N: usize = 1;
const K_T: usize = 2;
const K_TOR: usize = 3;
const K_BEND: usize = 4;
const GAMMA_N: usize = 5;
const GAMMA_T: usize = 6;
const GAMMA_TOR: usize = 7;
const GAMMA_BEND: usize = 8;
const AREA: usize = 9;
const R_B: usize = 10;
const IBEN: usize = 11;
const JPOL: usize = 12;
const LIMIT_TENSILE: usize = 13;
const LIMIT_SHEAR: usize = 14;

/// History per bond end: shear displacement, then relative rotation.
pub const BOND_STATE: usize = 6;

#[cube(launch_unchecked)]
#[allow(clippy::too_many_arguments)]
fn bond_force<F: Float>(
    sizes: &[u32],
    params: &[F],
    start: &[u32],
    table: &[F],
    pos: &[F],
    vel: &[F],
    omega: &[F],
    force: &mut [F],
    torque: &mut [F],
    partner: &mut [u32],
    state: &mut [F],
) {
    let i = ABSOLUTE_POS;
    if (i as u32) >= sizes[0] {
        terminate!();
    }
    let zero = F::new(0.0f32);
    let (dt, advance, lx, ly, lz, tilt_xy) = (
        params[0], params[1], params[2], params[3], params[4], params[5],
    );
    let advancing = advance != zero;
    let (xi, yi, zi) = (pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
    let (vix, viy, viz) = (vel[3 * i], vel[3 * i + 1], vel[3 * i + 2]);
    let (wix, wiy, wiz) = (omega[3 * i], omega[3 * i + 1], omega[3 * i + 2]);
    let mut fx = force[3 * i];
    let mut fy = force[3 * i + 1];
    let mut fz = force[3 * i + 2];
    let mut tx = torque[3 * i];
    let mut ty = torque[3 * i + 1];
    let mut tz = torque[3 * i + 2];

    let end = start[i + 1];
    let mut b = start[i];
    while b < end {
        let slot = b as usize;
        b += 1u32;
        let other = partner[slot];
        if other != comptime!(SENTINEL) {
            let j = other as usize;
            let row = slot * comptime!(BOND_FIELDS);
            let s = slot * comptime!(BOND_STATE);
            let mut dx = pos[3 * j] - xi;
            let mut dy = pos[3 * j + 1] - yi;
            let mut dz = pos[3 * j + 2] - zi;
            if lz > zero {
                dz -= lz * F::round(dz / lz);
            }
            if ly > zero {
                let image = F::round(dy / ly);
                dy -= ly * image;
                dx -= tilt_xy * image;
            }
            if lx > zero {
                dx -= lx * F::round(dx / lx);
            }
            let dist = F::sqrt(dx * dx + dy * dy + dz * dz);
            if dist > zero {
                let (nx, ny, nz) = (dx / dist, dy / dist, dz / dist);
                let r0 = table[row + comptime!(R0)];

                // Kinematics at the mid-point, a lever of L/2 along the axis.
                let half = F::new(0.5f32) * r0;
                let (r1x, r1y, r1z) = (half * nx, half * ny, half * nz);
                let (wjx, wjy, wjz) = (omega[3 * j], omega[3 * j + 1], omega[3 * j + 2]);
                let vr_x = (vel[3 * j] - (wjy * r1z - wjz * r1y)) - (vix + wiy * r1z - wiz * r1y);
                let vr_y =
                    (vel[3 * j + 1] - (wjz * r1x - wjx * r1z)) - (viy + wiz * r1x - wix * r1z);
                let vr_z =
                    (vel[3 * j + 2] - (wjx * r1y - wjy * r1x)) - (viz + wix * r1y - wiy * r1x);
                let v_n = vr_x * nx + vr_y * ny + vr_z * nz;
                let vt_x = vr_x - v_n * nx;
                let vt_y = vr_y - v_n * ny;
                let vt_z = vr_z - v_n * nz;

                let f_n = table[row + comptime!(K_N)] * (dist - r0)
                    + table[row + comptime!(GAMMA_N)] * v_n;

                // Shear: reproject onto the current normal plane, then integrate.
                let s_n = state[s] * nx + state[s + 1] * ny + state[s + 2] * nz;
                let sx = state[s] - s_n * nx + vt_x * dt * advance;
                let sy = state[s + 1] - s_n * ny + vt_y * dt * advance;
                let sz = state[s + 2] - s_n * nz + vt_z * dt * advance;
                let k_t = table[row + comptime!(K_T)];
                let gamma_t = table[row + comptime!(GAMMA_T)];
                let ft_x = k_t * sx + gamma_t * vt_x;
                let ft_y = k_t * sy + gamma_t * vt_y;
                let ft_z = k_t * sz + gamma_t * vt_z;

                // Rotation: split into twist along the axis and bending across it.
                let wr_x = wjx - wix;
                let wr_y = wjy - wiy;
                let wr_z = wjz - wiz;
                let wr_n = wr_x * nx + wr_y * ny + wr_z * nz;
                let th_x = state[s + 3] + wr_x * dt * advance;
                let th_y = state[s + 4] + wr_y * dt * advance;
                let th_z = state[s + 5] + wr_z * dt * advance;
                let th_n = th_x * nx + th_y * ny + th_z * nz;
                let k_tor = table[row + comptime!(K_TOR)];
                let gamma_tor = table[row + comptime!(GAMMA_TOR)];
                let k_bend = table[row + comptime!(K_BEND)];
                let gamma_bend = table[row + comptime!(GAMMA_BEND)];
                let mt_x = k_tor * (th_n * nx) + gamma_tor * (wr_n * nx);
                let mt_y = k_tor * (th_n * ny) + gamma_tor * (wr_n * ny);
                let mt_z = k_tor * (th_n * nz) + gamma_tor * (wr_n * nz);
                let mb_x = k_bend * (th_x - th_n * nx) + gamma_bend * (wr_x - wr_n * nx);
                let mb_y = k_bend * (th_y - th_n * ny) + gamma_bend * (wr_y - wr_n * ny);
                let mb_z = k_bend * (th_z - th_n * nz) + gamma_bend * (wr_z - wr_n * nz);

                // Combined stress. An unbreakable bond has infinite limits.
                let area = table[row + comptime!(AREA)];
                let r_b = table[row + comptime!(R_B)];
                let m_bend = F::sqrt(mb_x * mb_x + mb_y * mb_y + mb_z * mb_z);
                let m_tor = F::sqrt(mt_x * mt_x + mt_y * mt_y + mt_z * mt_z);
                let f_t = F::sqrt(ft_x * ft_x + ft_y * ft_y + ft_z * ft_z);
                let tensile =
                    F::max(f_n, zero) / area + r_b * m_bend / table[row + comptime!(IBEN)];
                let shear = f_t / area + r_b * m_tor / table[row + comptime!(JPOL)];
                let broken = tensile > table[row + comptime!(LIMIT_TENSILE)]
                    || shear > table[row + comptime!(LIMIT_SHEAR)];

                if broken {
                    if advancing {
                        partner[slot] = comptime!(SENTINEL);
                    }
                } else {
                    fx += f_n * nx + ft_x;
                    fy += f_n * ny + ft_y;
                    fz += f_n * nz + ft_z;
                    tx += r1y * ft_z - r1z * ft_y + mt_x + mb_x;
                    ty += r1z * ft_x - r1x * ft_z + mt_y + mb_y;
                    tz += r1x * ft_y - r1y * ft_x + mt_z + mb_z;
                    // The history is in place, so the entry prime must not write it back.
                    if advancing {
                        state[s] = sx;
                        state[s + 1] = sy;
                        state[s + 2] = sz;
                        state[s + 3] = th_x;
                        state[s + 4] = th_y;
                        state[s + 5] = th_z;
                    }
                }
            }
        }
    }

    force[3 * i] = fx;
    force[3 * i + 1] = fy;
    force[3 * i + 2] = fz;
    torque[3 * i] = tx;
    torque[3 * i + 1] = ty;
    torque[3 * i + 2] = tz;
}

/// A bond topology as flat columns: each atom's ends, in compressed rows.
#[derive(Clone, Debug, Default)]
pub struct BondTopology {
    /// Atom count + 1 offsets into the rest.
    pub start: Vec<u32>,
    /// The other atom's row at each end, or `SENTINEL` once broken.
    pub partner: Vec<u32>,
    /// `BOND_FIELDS` values per end.
    pub table: Vec<GpuReal>,
    /// `BOND_STATE` values per end.
    pub state: Vec<GpuReal>,
}

struct Live {
    n: usize,
    start: DeviceColumn,
    table: DeviceColumn,
    partner: DeviceColumn,
    state: DeviceColumn,
    ends: usize,
}

/// What the host keeps doing after the hook is handed to a [`GpuState`].
#[derive(Clone)]
pub struct BondColumns {
    gpu: GpuDevice,
    live: Rc<RefCell<Live>>,
}

fn live(gpu: &GpuDevice, topology: &BondTopology) -> Live {
    // A kernel argument cannot be empty.
    let or_one = |column: &[u32]| {
        if column.is_empty() {
            vec![SENTINEL]
        } else {
            column.to_vec()
        }
    };
    let or_zeros = |column: &[GpuReal], len| {
        if column.is_empty() {
            vec![0.0; len]
        } else {
            column.to_vec()
        }
    };
    Live {
        n: topology.start.len() - 1,
        start: DeviceColumn::from_slice(gpu, &topology.start),
        table: DeviceColumn::from_slice(gpu, &or_zeros(&topology.table, BOND_FIELDS)),
        partner: DeviceColumn::from_slice(gpu, &or_one(&topology.partner)),
        state: DeviceColumn::from_slice(gpu, &or_zeros(&topology.state, BOND_STATE)),
        ends: topology.partner.len(),
    }
}

impl BondColumns {
    /// Replaces the topology, for a reorder or a change in atoms.
    pub fn write(&self, topology: &BondTopology) {
        *self.live.borrow_mut() = live(&self.gpu, topology);
    }

    /// The partner column itself, which holds `SENTINEL` at broken ends. Contact reads it.
    pub fn alive(&self) -> DeviceColumn {
        self.live.borrow().partner.clone()
    }

    /// The partner column as last written, waiting for kernels already launched.
    pub fn read_partners(&self) -> Vec<u32> {
        let live = self.live.borrow();
        live.partner.read(&self.gpu, live.ends)
    }

    /// The history as last written, waiting for kernels already launched.
    pub fn read_state(&self) -> Vec<GpuReal> {
        let live = self.live.borrow();
        live.state.read(&self.gpu, live.ends * BOND_STATE)
    }
}

/// Beam bond hook.
pub struct BondForce {
    omega: usize,
    sizes: DeviceColumn,
    params: DeviceColumn,
    params_prime: DeviceColumn,
    columns: BondColumns,
}

impl BondForce {
    /// Build the hook. `omega` is the rotation aux index, `params` is
    /// `[dt, lx, ly, lz, tilt_xy]` with 0 lengths on non-periodic axes.
    pub fn new(
        state: &GpuState,
        omega: usize,
        topology: &BondTopology,
        params: [GpuReal; 5],
    ) -> Self {
        let gpu = state.device().clone();
        let with = |advance| {
            [
                params[0], advance, params[1], params[2], params[3], params[4],
            ]
        };
        Self {
            omega,
            sizes: DeviceColumn::from_slice(&gpu, &[(topology.start.len() - 1) as u32]),
            params: DeviceColumn::from_slice(&gpu, &with(1.0)),
            params_prime: DeviceColumn::from_slice(&gpu, &with(0.0)),
            columns: BondColumns {
                live: Rc::new(RefCell::new(live(&gpu, topology))),
                gpu,
            },
        }
    }

    /// A handle for the host's side of the columns.
    pub fn columns(&self) -> BondColumns {
        self.columns.clone()
    }

    fn launch_with(&mut self, c: &StepColumns<'_>, params: BufferArg) {
        let live = self.columns.live.borrow();
        if live.ends == 0 {
            return;
        }
        let aux = &c.aux[self.omega];
        let sizes = DeviceColumn::from_slice(c.gpu, &[live.n as u32]);
        self.sizes = sizes;
        // SAFETY: rows and partners are atoms of these columns, by how the topology is built.
        unsafe {
            bond_force::launch_unchecked::<GpuReal>(
                &c.gpu.client,
                CubeCount::Static((live.n as u32).div_ceil(64).max(1), 1, 1),
                CubeDim::new_1d(64),
                self.sizes.arg(),
                params,
                live.start.arg(),
                live.table.arg(),
                c.pos.arg(),
                c.vel.arg(),
                aux.state.arg(),
                c.force.arg(),
                aux.rate.arg(),
                live.partner.arg(),
                live.state.arg(),
            );
        }
    }
}

impl GpuForce for BondForce {
    fn launch(&mut self, c: &StepColumns<'_>) {
        let params = self.params.arg();
        self.launch_with(c, params);
    }

    fn launch_prime(&mut self, c: &StepColumns<'_>) {
        let params = self.params_prime.arg();
        self.launch_with(c, params);
    }
}
