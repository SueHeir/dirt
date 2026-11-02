//! DEM granular contact force as a [`GpuForce`] hook (DEM-specific; lives in dirt).
//!
//! Hertz normal and Mindlin tangential contact with persistent per-contact spring history and
//! torque. It is i-centric with no atomics: each atom owns its half of every spring, and the
//! two halves stay mirror images (the normal flips and so does the relative velocity), so no
//! canonical frame is needed on the device. Soil only sees a force hook and an aux DOF.

use std::cell::RefCell;
use std::rc::Rc;

use cubecl::prelude::*;
use grass_gpu::cubecl;
use soil_gpu::neighbor_slots::{slot_clear_from, slot_lookup, slot_write, MAX_CONTACTS};
use soil_gpu::{
    DeviceColumn, GpuDevice, GpuForce, GpuReal, GpuState, Grid, NeighborSlots, StepColumns,
};

/// Material and contact parameters (effective moduli already reduced for the pair).
#[derive(Clone, Copy, Debug)]
pub struct GranularConfig {
    pub e_eff: GpuReal,
    pub beta: GpuReal,
    pub g_eff: GpuReal,
    pub mu: GpuReal,
    pub dt: GpuReal,
    /// Periodic box (0 length means a non-periodic axis): the stencil wraps and the pair vector
    /// uses the triclinic minimum image. A y-image pair gets the Lees-Edwards x tilt and
    /// velocity offset. Defaults to an open box.
    pub lx: GpuReal,
    pub ly: GpuReal,
    pub lz: GpuReal,
    pub tilt_xy: GpuReal,
    pub dv_xy: GpuReal,
}

impl GranularConfig {
    /// An open-box config.
    pub fn new(e_eff: GpuReal, beta: GpuReal, g_eff: GpuReal, mu: GpuReal, dt: GpuReal) -> Self {
        Self {
            e_eff,
            beta,
            g_eff,
            mu,
            dt,
            lx: 0.0,
            ly: 0.0,
            lz: 0.0,
            tilt_xy: 0.0,
            dv_xy: 0.0,
        }
    }

    /// `advance` is 1 in the step loop and 0 in the entry prime, which reprojects the spring
    /// but does not integrate it.
    fn params(&self, advance: GpuReal) -> [GpuReal; 11] {
        [
            self.e_eff,
            self.beta,
            self.g_eff,
            self.mu,
            self.dt,
            advance,
            self.lx,
            self.ly,
            self.lz,
            self.tilt_xy,
            self.dv_xy,
        ]
    }
}

const TANGENTIAL_EPSILON: f32 = 1.0e-10;

#[cube(launch_unchecked)]
#[allow(clippy::too_many_arguments)]
fn contact_force<F: Float>(
    sizes: &[i32],
    params: &[F],
    pos: &[F],
    vel: &[F],
    omega: &[F],
    radius: &[F],
    inv_mass: &[F],
    tag: &[u32],
    atom_cell: &[u32],
    cell_start: &[u32],
    sorted_atoms: &[u32],
    old_partner: &[u32],
    old_payload: &[F],
    force: &mut [F],
    torque: &mut [F],
    new_partner: &mut [u32],
    new_payload: &mut [F],
) {
    let i = ABSOLUTE_POS;
    if (i as i32) >= sizes[0] {
        terminate!();
    }
    let zero = F::new(0.0f32);
    let (nx_cells, ny_cells, nz_cells) = (sizes[1], sizes[2], sizes[3]);
    let (e_eff, beta, g_eff, mu, dt, advance) = (
        params[0], params[1], params[2], params[3], params[4], params[5],
    );
    let (lx, ly, lz, tilt_xy, dv_xy) = (params[6], params[7], params[8], params[9], params[10]);
    let sqrt_5_6 = F::sqrt(F::new(5.0f32) / F::new(6.0f32));
    let epsilon = F::new(comptime!(TANGENTIAL_EPSILON));

    let (xi, yi, zi) = (pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
    let (vix, viy, viz) = (vel[3 * i], vel[3 * i + 1], vel[3 * i + 2]);
    let (wix, wiy, wiz) = (omega[3 * i], omega[3 * i + 1], omega[3 * i + 2]);
    let ri = radius[i];
    let imi = inv_mass[i];

    let c = i32::cast_from(atom_cell[i]);
    let cz0 = c % nz_cells;
    let cy0 = (c / nz_cells) % ny_cells;
    let cx0 = c / (ny_cells * nz_cells);

    // Start from the body force seed. Torque belongs to this hook, so it starts at zero.
    let mut fx = force[3 * i];
    let mut fy = force[3 * i + 1];
    let mut fz = force[3 * i + 2];
    let mut tx = zero;
    let mut ty = zero;
    let mut tz = zero;
    let mut nslot = 0u32;
    let i_u32 = i as u32;

    let mut dx_cell = 0i32;
    while dx_cell < 3i32 {
        let mut cx = cx0 + dx_cell - 1i32;
        dx_cell += 1i32;
        if lx > zero {
            cx = (cx + nx_cells) % nx_cells;
        }
        if cx >= 0i32 && cx < nx_cells {
            let mut dy_cell = 0i32;
            while dy_cell < 3i32 {
                let mut cy = cy0 + dy_cell - 1i32;
                dy_cell += 1i32;
                if ly > zero {
                    cy = (cy + ny_cells) % ny_cells;
                }
                if cy >= 0i32 && cy < ny_cells {
                    let mut dz_cell = 0i32;
                    while dz_cell < 3i32 {
                        let mut cz = cz0 + dz_cell - 1i32;
                        dz_cell += 1i32;
                        if lz > zero {
                            cz = (cz + nz_cells) % nz_cells;
                        }
                        if cz >= 0i32 && cz < nz_cells {
                            let cc = ((cx * ny_cells + cy) * nz_cells + cz) as usize;
                            let end = cell_start[cc + 1];
                            let mut m = cell_start[cc];
                            while m < end {
                                let j = sorted_atoms[m as usize] as usize;
                                m += 1u32;
                                if j != i {
                                    let mut ddx = pos[3 * j] - xi;
                                    let mut ddy = pos[3 * j + 1] - yi;
                                    let mut ddz = pos[3 * j + 2] - zi;
                                    // Triclinic minimum image. A y-image pair is shifted in x by the tilt,
                                    // and `image_y` carries the velocity offset below.
                                    let mut image_y = zero;
                                    if lz > zero {
                                        ddz -= lz * F::round(ddz / lz);
                                    }
                                    if ly > zero {
                                        image_y = F::round(ddy / ly);
                                        ddy -= ly * image_y;
                                        ddx -= tilt_xy * image_y;
                                    }
                                    if lx > zero {
                                        ddx -= lx * F::round(ddx / lx);
                                    }
                                    let dist_sq = ddx * ddx + ddy * ddy + ddz * ddz;
                                    let rj = radius[j];
                                    let sum_r = ri + rj;
                                    if dist_sq < sum_r * sum_r && dist_sq != zero {
                                        let distance = F::sqrt(dist_sq);
                                        let delta = F::min(
                                            sum_r - distance,
                                            F::new(0.5f32) * F::min(ri, rj),
                                        );
                                        if delta > zero {
                                            let inv_dist = F::new(1.0f32) / distance;
                                            let nx = ddx * inv_dist;
                                            let ny = ddy * inv_dist;
                                            let nz = ddz * inv_dist;

                                            let r_eff = (ri * rj) / sum_r;
                                            let sdr = F::sqrt(delta * r_eff);
                                            let s_n = F::new(2.0f32) * e_eff * sdr;
                                            let k_n = F::new(4.0f32) / F::new(3.0f32) * e_eff * sdr;
                                            let k_t = F::new(8.0f32) * g_eff * sdr;
                                            let m_r = F::new(1.0f32) / (imi + inv_mass[j]);

                                            let (r1x, r1y, r1z) = (ri * nx, ri * ny, ri * nz);
                                            let vc_ix = vix + (wiy * r1z - wiz * r1y);
                                            let vc_iy = viy + (wiz * r1x - wix * r1z);
                                            let vc_iz = viz + (wix * r1y - wiy * r1x);
                                            let (wjx, wjy, wjz) =
                                                (omega[3 * j], omega[3 * j + 1], omega[3 * j + 2]);
                                            let (r2x, r2y, r2z) = (rj * nx, rj * ny, rj * nz);
                                            let vc_jx = (vel[3 * j] - dv_xy * image_y)
                                                + (wjz * r2y - wjy * r2z);
                                            let vc_jy = vel[3 * j + 1] + (wjx * r2z - wjz * r2x);
                                            let vc_jz = vel[3 * j + 2] + (wjy * r2x - wjx * r2y);

                                            let vr_x = vc_jx - vc_ix;
                                            let vr_y = vc_jy - vc_iy;
                                            let vr_z = vc_jz - vc_iz;
                                            let v_n = vr_x * nx + vr_y * ny + vr_z * nz;

                                            let f_diss_n = F::new(2.0f32)
                                                * beta
                                                * sqrt_5_6
                                                * F::sqrt(s_n * m_r)
                                                * v_n;
                                            let f_n_mag = F::max(k_n * delta - f_diss_n, zero);
                                            fx -= f_n_mag * nx;
                                            fy -= f_n_mag * ny;
                                            fz -= f_n_mag * nz;

                                            let vt_x = vr_x - v_n * nx;
                                            let vt_y = vr_y - v_n * ny;
                                            let vt_z = vr_z - v_n * nz;

                                            let partner = tag[j];
                                            let old = slot_lookup::<F>(
                                                old_partner,
                                                old_payload,
                                                i_u32,
                                                partner,
                                            );
                                            let s_dot_n = old.0 * nx + old.1 * ny + old.2 * nz;
                                            let mut sx = old.0 - s_dot_n * nx + vt_x * dt * advance;
                                            let mut sy = old.1 - s_dot_n * ny + vt_y * dt * advance;
                                            let mut sz = old.2 - s_dot_n * nz + vt_z * dt * advance;

                                            let f_t_max = mu * F::abs(f_n_mag);
                                            let spring = k_t * F::sqrt(sx * sx + sy * sy + sz * sz);
                                            if spring > f_t_max && spring > epsilon {
                                                let scale = f_t_max / spring;
                                                sx *= scale;
                                                sy *= scale;
                                                sz *= scale;
                                            }

                                            let gamma_t = F::new(2.0f32)
                                                * sqrt_5_6
                                                * beta
                                                * F::sqrt(k_t * m_r);
                                            let mut ft_x = k_t * sx + gamma_t * vt_x;
                                            let mut ft_y = k_t * sy + gamma_t * vt_y;
                                            let mut ft_z = k_t * sz + gamma_t * vt_z;
                                            let f_t_mag =
                                                F::sqrt(ft_x * ft_x + ft_y * ft_y + ft_z * ft_z);
                                            if f_t_mag > f_t_max && f_t_mag > epsilon {
                                                let scale = f_t_max / f_t_mag;
                                                ft_x *= scale;
                                                ft_y *= scale;
                                                ft_z *= scale;
                                            }

                                            fx += ft_x;
                                            fy += ft_y;
                                            fz += ft_z;
                                            tx += r1y * ft_z - r1z * ft_y;
                                            ty += r1z * ft_x - r1x * ft_z;
                                            tz += r1x * ft_y - r1y * ft_x;

                                            if nslot < comptime!(MAX_CONTACTS as u32) {
                                                slot_write::<F>(
                                                    new_partner,
                                                    new_payload,
                                                    i_u32,
                                                    nslot,
                                                    partner,
                                                    (sx, sy, sz),
                                                );
                                                nslot += 1u32;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    slot_clear_from(new_partner, i_u32, nslot);

    force[3 * i] = fx;
    force[3 * i + 1] = fy;
    force[3 * i + 2] = fz;
    torque[3 * i] = tx;
    torque[3 * i + 1] = ty;
    torque[3 * i + 2] = tz;
}

/// Granular contact-force hook.
pub struct GranularForce {
    n_force: u32,
    omega: usize,
    sizes: DeviceColumn,
    params: DeviceColumn,
    params_prime: DeviceColumn,
    columns: GranularColumns,
}

/// What the host keeps writing after the hook is handed to a [`GpuState`].
#[derive(Clone)]
pub struct GranularColumns {
    gpu: GpuDevice,
    radius: Rc<RefCell<DeviceColumn>>,
    tag: Rc<RefCell<DeviceColumn>>,
    slots: Rc<RefCell<NeighborSlots>>,
}

impl GranularColumns {
    /// Replaces the radii, ghosts included.
    pub fn write_radius(&self, radius: &[GpuReal]) {
        *self.radius.borrow_mut() = DeviceColumn::from_slice(&self.gpu, radius);
    }

    /// Replaces the tags. Contact history is keyed by tag, so a spring keeps its partner when
    /// atoms reorder.
    pub fn write_tag(&self, tag: &[u32]) {
        *self.tag.borrow_mut() = DeviceColumn::from_slice(&self.gpu, tag);
    }

    /// The springs as last written, for the first `rows` atoms.
    pub fn read_slots(&self, rows: usize) -> (Vec<u32>, Vec<GpuReal>) {
        self.slots.borrow().read_latest(rows)
    }

    /// Replaces the springs, with slots for `rows` atoms.
    pub fn write_slots(&self, rows: usize, partner: &[u32], payload: &[GpuReal]) {
        let mut slots = self.slots.borrow_mut();
        if slots.len() != rows {
            *slots = NeighborSlots::new(self.gpu.clone(), rows);
        }
        slots.write_for_next(partner, payload);
    }
}

impl GranularForce {
    /// Build the hook. `omega` is the aux index registered on `state` for rotation (state =
    /// angular velocity, rate = torque). `radius` is per atom.
    pub fn new(
        state: &GpuState,
        grid: &Grid,
        omega: usize,
        radius: &[GpuReal],
        config: GranularConfig,
    ) -> Self {
        let gpu = state.device().clone();
        let n = state.n();
        assert_eq!(radius.len(), n);
        // Until the host says otherwise an atom's tag is its index, which is right for a
        // resident run where atoms never reorder.
        let identity: Vec<u32> = (0..n as u32).collect();
        let columns = GranularColumns {
            radius: Rc::new(RefCell::new(DeviceColumn::from_slice(&gpu, radius))),
            tag: Rc::new(RefCell::new(DeviceColumn::from_slice(&gpu, &identity))),
            slots: Rc::new(RefCell::new(NeighborSlots::new(gpu.clone(), n))),
            gpu: gpu.clone(),
        };
        Self {
            n_force: n as u32,
            omega,
            sizes: DeviceColumn::from_slice(&gpu, &[n as i32, grid.n[0], grid.n[1], grid.n[2]]),
            params: DeviceColumn::from_slice(&gpu, &config.params(1.0)),
            params_prime: DeviceColumn::from_slice(&gpu, &config.params(0.0)),
            columns,
        }
    }

    /// A handle for writing radii, tags and history after the hook is boxed.
    pub fn columns(&self) -> GranularColumns {
        self.columns.clone()
    }

    fn launch_with(&self, c: &StepColumns<'_>, params: BufferArg) {
        let slots = self.columns.slots.borrow();
        let (old_partner, old_payload, new_partner, new_payload) = slots.args();
        let aux = &c.aux[self.omega];
        // SAFETY: indices come from the cell list, whose size matches these columns, and from
        // the force count in `sizes`.
        unsafe {
            contact_force::launch_unchecked::<GpuReal>(
                &c.gpu.client,
                CubeCount::Static(self.n_force.div_ceil(64).max(1), 1, 1),
                CubeDim::new_1d(64),
                self.sizes.arg(),
                params,
                c.pos.arg(),
                c.vel.arg(),
                aux.state.arg(),
                self.columns.radius.borrow().arg(),
                c.inv_mass.arg(),
                self.columns.tag.borrow().arg(),
                c.cells.atom_cell(),
                c.cells.cell_start(),
                c.cells.sorted_atoms(),
                old_partner,
                old_payload,
                c.force.arg(),
                aux.rate.arg(),
                new_partner,
                new_payload,
            );
        }
    }
}

impl GpuForce for GranularForce {
    fn set_active(&mut self, _n_all: usize, n_force: usize, grid: &Grid) {
        self.n_force = n_force as u32;
        self.sizes = DeviceColumn::from_slice(
            &self.columns.gpu,
            &[n_force as i32, grid.n[0], grid.n[1], grid.n[2]],
        );
    }

    fn launch(&mut self, c: &StepColumns<'_>) {
        self.launch_with(c, self.params.arg());
        self.columns.slots.borrow().swap();
    }

    /// The force from the current spring, with no integration and no swap, so the old slots
    /// the next step reads are untouched. Idempotent across window boundaries.
    fn launch_prime(&mut self, c: &StepColumns<'_>) {
        self.launch_with(c, self.params_prime.arg());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pair_state(
        gpu: GpuDevice,
        r: GpuReal,
        p0: &[GpuReal; 6],
        grid: Grid,
        dt: GpuReal,
    ) -> (GpuState, usize) {
        let mass = 4.0 / 3.0 * std::f64::consts::PI as GpuReal * r * r * r;
        let inertia = 0.4 * mass * r * r;
        let mut state = GpuState::new(gpu, 2, grid);
        state.set_params(dt, [0.0; 3]);
        state.set_state(p0, &[0.0; 6], &[1.0 / mass; 2]);
        let omega = state.add_aux_dof();
        state.set_aux_inv_coeff(omega, &[1.0 / inertia; 2]);
        state.set_aux_state(omega, &[0.0; 6]);
        (state, omega)
    }

    #[test]
    fn two_spheres_repel_symmetrically() {
        let Ok(gpu) = GpuDevice::default_device() else {
            eprintln!("no GPU adapter; skipping");
            return;
        };
        let r = 0.5;
        let p0 = [-0.45, 0.0, 0.0, 0.45, 0.0, 0.0];
        let grid = Grid::from_positions(&p0, 2.0 * r);
        let dt = 1.0e-5;
        let (mut state, omega) = pair_state(gpu, r, &p0, grid, dt);
        let hook = GranularForce::new(
            &state,
            &grid,
            omega,
            &[r; 2],
            GranularConfig::new(1.0e6, 0.1, 4.0e5, 0.5, dt),
        );
        state.add_force_hook(Box::new(hook));
        state.run_steps(3000);
        let p = state.pos().to_vec();
        let v = state.vel().to_vec();
        assert!(p[0] < -0.45 && p[3] > 0.45, "did not separate: {p:?}");
        assert!((p[0] + p[3]).abs() < 1e-4, "asymmetric: {p:?}");
        assert!((v[0] + v[3]).abs() < 1e-4, "momentum: {v:?}");
        assert!(p[1].abs() < 1e-5 && p[2].abs() < 1e-5, "off-axis: {p:?}");
    }

    #[test]
    fn contact_across_a_periodic_boundary() {
        let Ok(gpu) = GpuDevice::default_device() else {
            eprintln!("no GPU adapter; skipping");
            return;
        };
        // Grains straddle the x boundary: raw distance 0.85, minimum image 0.15 < 0.2.
        let r = 0.1;
        let p0 = [0.05, 0.5, 0.5, 0.90, 0.5, 0.5];
        let grid = Grid::for_box([0.0; 3], [1.0; 3], [true; 3], &p0, 2.0 * r, 1000);
        let dt = 1.0e-5;
        let (mut state, omega) = pair_state(gpu, r, &p0, grid, dt);
        state.set_box([1.0; 3], [0.0; 3], 0.0, 0.0);
        let mut config = GranularConfig::new(1.0e6, 0.1, 4.0e5, 0.5, dt);
        (config.lx, config.ly, config.lz) = (1.0, 1.0, 1.0);
        let hook = GranularForce::new(&state, &grid, omega, &[r; 2], config);
        state.add_force_hook(Box::new(hook));
        // Few enough steps that they have not yet travelled around the box and met again.
        state.run_steps(20);
        let v = state.vel().to_vec();
        // Repelled across the boundary: atom 0 pushed +x, atom 1 pushed -x.
        assert!(v[0] > 0.0 && v[3] < 0.0, "{v:?}");
    }
}
