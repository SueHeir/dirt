//! DEM planar wall contact as a [`GpuForce`] hook: sphere-plane Hertz normal and Mindlin
//! tangential, with friction history per atom per wall. Each atom owns its own springs, so the
//! kernel writes them in place. Register after [`GranularForce`](crate::GranularForce), which is
//! what starts each torque at zero.

use std::cell::RefCell;
use std::rc::Rc;

use cubecl::prelude::*;
use grass_gpu::cubecl;
use soil_gpu::{Boundary, DeviceColumn, GpuDevice, GpuForce, GpuReal, GpuState, Grid, StepColumns};

/// The most walls a run can hold.
pub const MAX_WALLS: usize = 8;

const TANGENTIAL_EPSILON: f32 = 1.0e-10;

#[cube(launch_unchecked)]
#[allow(clippy::too_many_arguments)]
fn wall_force<F: Float>(
    sizes: &[u32],
    params: &[F],
    planes: &[F],
    pos: &[F],
    vel: &[F],
    omega: &[F],
    radius: &[F],
    inv_mass: &[F],
    force: &mut [F],
    torque: &mut [F],
    wall_spring: &mut [F],
) {
    let i = ABSOLUTE_POS;
    if (i as u32) >= sizes[0] {
        terminate!();
    }
    let zero = F::new(0.0f32);
    let n_walls = sizes[1];
    let (e_eff, beta, g_eff, mu, dt, advance) = (
        params[0], params[1], params[2], params[3], params[4], params[5],
    );
    let sqrt_5_6 = F::sqrt(F::new(5.0f32) / F::new(6.0f32));
    let epsilon = F::new(comptime!(TANGENTIAL_EPSILON));

    let (px, py, pz) = (pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
    let (vix, viy, viz) = (vel[3 * i], vel[3 * i + 1], vel[3 * i + 2]);
    let (wix, wiy, wiz) = (omega[3 * i], omega[3 * i + 1], omega[3 * i + 2]);
    let ri = radius[i];
    let imi = inv_mass[i];
    let mut m_i = zero;
    if imi > zero {
        m_i = F::new(1.0f32) / imi;
    }
    let mut fx = force[3 * i];
    let mut fy = force[3 * i + 1];
    let mut fz = force[3 * i + 2];
    let mut tx = torque[3 * i];
    let mut ty = torque[3 * i + 1];
    let mut tz = torque[3 * i + 2];

    let mut w = 0u32;
    while w < n_walls {
        let wall = w as usize;
        let sb = (3u32 * ((i as u32) * comptime!(MAX_WALLS as u32) + w)) as usize;
        w += 1u32;
        let (nx, ny, nz) = (
            planes[6 * wall + 3],
            planes[6 * wall + 4],
            planes[6 * wall + 5],
        );
        let signed = (px - planes[6 * wall]) * nx
            + (py - planes[6 * wall + 1]) * ny
            + (pz - planes[6 * wall + 2]) * nz;
        let overlap = ri - signed;
        if overlap <= zero {
            if advance != zero {
                wall_spring[sb] = zero;
                wall_spring[sb + 1] = zero;
                wall_spring[sb + 2] = zero;
            }
        } else {
            let delta = F::min(overlap, F::new(0.5f32) * ri);
            let sdr = F::sqrt(delta * ri);
            let s_n = F::new(2.0f32) * e_eff * sdr;
            let k_n = F::new(4.0f32) / F::new(3.0f32) * e_eff * sdr;
            let k_t = F::new(8.0f32) * g_eff * sdr;

            let vs_x = vix - ri * (wiy * nz - wiz * ny);
            let vs_y = viy - ri * (wiz * nx - wix * nz);
            let vs_z = viz - ri * (wix * ny - wiy * nx);
            let v_n = vs_x * nx + vs_y * ny + vs_z * nz;

            let f_diss_n = F::new(2.0f32) * beta * sqrt_5_6 * F::sqrt(s_n * m_i) * v_n;
            let f_n_mag = F::max(k_n * delta - f_diss_n, zero);
            fx += f_n_mag * nx;
            fy += f_n_mag * ny;
            fz += f_n_mag * nz;

            let vt_x = vs_x - v_n * nx;
            let vt_y = vs_y - v_n * ny;
            let vt_z = vs_z - v_n * nz;

            let (ox, oy, oz) = (wall_spring[sb], wall_spring[sb + 1], wall_spring[sb + 2]);
            let s_dot_n = ox * nx + oy * ny + oz * nz;
            let mut sx = ox - s_dot_n * nx + vt_x * dt * advance;
            let mut sy = oy - s_dot_n * ny + vt_y * dt * advance;
            let mut sz = oz - s_dot_n * nz + vt_z * dt * advance;

            let f_t_max = mu * F::abs(f_n_mag);
            let spring = k_t * F::sqrt(sx * sx + sy * sy + sz * sz);
            if spring > f_t_max && spring > epsilon {
                let scale = f_t_max / spring;
                sx *= scale;
                sy *= scale;
                sz *= scale;
            }

            let gamma_t = F::new(2.0f32) * sqrt_5_6 * beta * F::sqrt(k_t * m_i);
            let mut ft_x = -(k_t * sx + gamma_t * vt_x);
            let mut ft_y = -(k_t * sy + gamma_t * vt_y);
            let mut ft_z = -(k_t * sz + gamma_t * vt_z);
            let f_t_mag = F::sqrt(ft_x * ft_x + ft_y * ft_y + ft_z * ft_z);
            if f_t_mag > f_t_max && f_t_mag > epsilon {
                let scale = f_t_max / f_t_mag;
                ft_x *= scale;
                ft_y *= scale;
                ft_z *= scale;
            }

            fx += ft_x;
            fy += ft_y;
            fz += ft_z;
            tx -= ri * (ny * ft_z - nz * ft_y);
            ty -= ri * (nz * ft_x - nx * ft_z);
            tz -= ri * (nx * ft_y - ny * ft_x);

            // The spring is stored in place, so the entry prime must not write it back.
            if advance != zero {
                wall_spring[sb] = sx;
                wall_spring[sb + 1] = sy;
                wall_spring[sb + 2] = sz;
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

fn springs_len(rows: usize) -> usize {
    rows * MAX_WALLS * 3
}

/// DEM wall-contact hook.
pub struct WallForce {
    n: usize,
    n_walls: usize,
    omega: usize,
    sizes: DeviceColumn,
    params: DeviceColumn,
    params_prime: DeviceColumn,
    planes: DeviceColumn,
    columns: WallColumns,
}

/// What the host keeps writing after the hook is handed to a [`GpuState`].
#[derive(Clone)]
pub struct WallColumns {
    gpu: GpuDevice,
    radius: Rc<RefCell<DeviceColumn>>,
    spring: Rc<RefCell<DeviceColumn>>,
}

impl WallColumns {
    /// Replaces the radii.
    pub fn write_radius(&self, radius: &[GpuReal]) {
        *self.radius.borrow_mut() = DeviceColumn::from_slice(&self.gpu, radius);
    }

    /// Each of the first `rows` atoms' springs: `[rows][MAX_WALLS][3]`, flattened.
    pub fn read_springs(&self, rows: usize) -> Vec<GpuReal> {
        self.spring.borrow().read(&self.gpu, springs_len(rows))
    }

    /// Replaces the springs, in the layout of [`read_springs`](Self::read_springs).
    pub fn write_springs(&self, springs: &[GpuReal]) {
        *self.spring.borrow_mut() = DeviceColumn::from_slice(&self.gpu, springs);
    }
}

impl WallForce {
    /// Build the hook. `omega` is the rotation aux index, the same one passed to
    /// [`GranularForce`](crate::GranularForce).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        state: &GpuState,
        omega: usize,
        radius: &[GpuReal],
        boundary: &Boundary,
        e_eff: GpuReal,
        beta: GpuReal,
        g_eff: GpuReal,
        mu: GpuReal,
        dt: GpuReal,
    ) -> Self {
        let gpu = state.device().clone();
        let n = state.n();
        assert_eq!(radius.len(), n);
        assert!(
            boundary.len() <= MAX_WALLS,
            "at most {MAX_WALLS} walls supported"
        );
        let params = |advance| [e_eff, beta, g_eff, mu, dt, advance];
        let columns = WallColumns {
            radius: Rc::new(RefCell::new(DeviceColumn::from_slice(&gpu, radius))),
            spring: Rc::new(RefCell::new(DeviceColumn::from_slice(
                &gpu,
                &vec![0.0 as GpuReal; springs_len(n.max(1))],
            ))),
            gpu: gpu.clone(),
        };
        Self {
            n,
            n_walls: boundary.len(),
            omega,
            sizes: DeviceColumn::from_slice(&gpu, &[n as u32, boundary.len() as u32]),
            params: DeviceColumn::from_slice(&gpu, &params(1.0)),
            params_prime: DeviceColumn::from_slice(&gpu, &params(0.0)),
            planes: DeviceColumn::from_slice(&gpu, &boundary.to_gpu()),
            columns,
        }
    }

    /// A handle for writing radii and springs after the hook is boxed.
    pub fn columns(&self) -> WallColumns {
        self.columns.clone()
    }

    fn launch_with(&self, c: &StepColumns<'_>, params: BufferArg) {
        let aux = &c.aux[self.omega];
        // SAFETY: bounded by `sizes`, which matches these columns.
        unsafe {
            wall_force::launch_unchecked::<GpuReal>(
                &c.gpu.client,
                CubeCount::Static((self.n as u32).div_ceil(64).max(1), 1, 1),
                CubeDim::new_1d(64),
                self.sizes.arg(),
                params,
                self.planes.arg(),
                c.pos.arg(),
                c.vel.arg(),
                aux.state.arg(),
                self.columns.radius.borrow().arg(),
                c.inv_mass.arg(),
                c.force.arg(),
                aux.rate.arg(),
                self.columns.spring.borrow().arg(),
            );
        }
    }
}

impl GpuForce for WallForce {
    fn set_active(&mut self, _n_all: usize, n_force: usize, _grid: &Grid) {
        self.n = n_force;
        self.sizes =
            DeviceColumn::from_slice(&self.columns.gpu, &[n_force as u32, self.n_walls as u32]);
    }

    fn launch(&mut self, c: &StepColumns<'_>) {
        self.launch_with(c, self.params.arg());
    }

    /// The force from the current springs, without integrating or writing them back.
    fn launch_prime(&mut self, c: &StepColumns<'_>) {
        self.launch_with(c, self.params_prime.arg());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GranularConfig, GranularForce};
    use soil_gpu::Plane;

    #[test]
    fn a_particle_settles_on_the_floor() {
        let Ok(gpu) = GpuDevice::default_device() else {
            eprintln!("no GPU adapter; skipping");
            return;
        };
        let r: GpuReal = 0.5;
        let mass = 4.0 / 3.0 * std::f64::consts::PI as GpuReal * r * r * r;
        let inertia = 0.4 * mass * r * r;
        let p0 = [0.0, 0.0, 0.7];
        let grid = Grid::from_positions(&[0.0, 0.0, 0.0, 0.0, 0.0, 2.0], 2.0 * r);
        let mut state = GpuState::new(gpu, 1, grid);
        let dt = 1.0e-5;
        state.set_params(dt, [0.0, 0.0, -9.81]);
        state.set_state(&p0, &[0.0; 3], &[1.0 / mass]);
        let omega = state.add_aux_dof();
        state.set_aux_inv_coeff(omega, &[1.0 / inertia]);
        state.set_aux_state(omega, &[0.0; 3]);
        let config = GranularConfig::new(1.0e6, 0.3, 4.0e5, 0.5, dt);
        let contact = GranularForce::new(&state, &grid, omega, &[r], config);
        state.add_force_hook(Box::new(contact));
        let mut boundary = Boundary::new();
        boundary.push(Plane::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]));
        let wall = WallForce::new(&state, omega, &[r], &boundary, 1.0e6, 0.3, 4.0e5, 0.5, dt);
        state.add_force_hook(Box::new(wall));
        state.run_steps(5000);
        let p = state.pos().to_vec();
        assert!(
            p[2] > 0.3 && p[2] < 0.7,
            "did not stay above the floor: {p:?}"
        );
        assert!(p[0].abs() < 1e-3 && p[1].abs() < 1e-3, "drifted: {p:?}");
    }
}
