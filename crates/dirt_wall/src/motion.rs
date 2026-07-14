//! Plane-wall motion modes and update systems.

use crate::geometry::Walls;
use grass_scheduler::prelude::*;
use soil_core::Atom;

/// Wall motion mode at runtime (plane walls only).
pub enum WallMotion {
    /// Wall is stationary.
    Static,
    /// Wall translates at a constant velocity each timestep.
    ConstantVelocity,
    /// Wall oscillates sinusoidally: `displacement = amplitude * sin(2π * frequency * t)`.
    Oscillate {
        /// Peak displacement (meters).
        amplitude: f64,
        /// Oscillation frequency (Hz).
        frequency: f64,
    },
    /// Proportional servo: velocity = clamp(gain × (target_force − measured_force), ±max_velocity).
    Servo {
        /// Desired total contact force (N).
        target_force: f64,
        /// Maximum wall speed (m/s).
        max_velocity: f64,
        /// Proportional gain (m/s per N).
        gain: f64,
    },
}

/// Runtime representation of an infinite plane wall.
///
/// Contact geometry for a plane wall:
/// ```text
///          ·  particle center (px, py, pz)
///         /|
///        / |  distance = dot(P - W, n̂)
///       /  |
///   n̂ ↑   |  delta (overlap) = radius - distance
///      |   |
///  ────W───┴──── wall surface
///      (point_x, point_y, point_z)
/// ```
///
/// The signed distance from the particle center to the wall plane is
/// `dot(pos - point, normal)`. Overlap is `delta = radius - distance`.

// ── Systems ─────────────────────────────────────────────────────────────────

/// Update wall positions and velocities according to their motion mode.
///
/// Runs in [`ParticleSimScheduleSet::PreInitialIntegration`] so walls are moved
/// *before* the integration step each timestep. Advances `walls.time` by `dt`.
pub fn wall_move(mut walls: ResMut<Walls>, atoms: Res<Atom>) {
    let dt = atoms.dt;
    let time = walls.time;

    let nplanes = walls.planes.len();
    for idx in 0..nplanes {
        if !walls.active[idx] {
            continue;
        }

        let wall = &mut walls.planes[idx];
        match wall.motion {
            WallMotion::Static => {}
            WallMotion::ConstantVelocity => {
                wall.point_x += wall.velocity[0] * dt;
                wall.point_y += wall.velocity[1] * dt;
                wall.point_z += wall.velocity[2] * dt;
            }
            WallMotion::Oscillate {
                amplitude,
                frequency,
            } => {
                let phase = 2.0 * std::f64::consts::PI * frequency * (time + dt);
                let disp = amplitude * phase.sin();
                wall.point_x = wall.origin[0] + disp * wall.normal_x;
                wall.point_y = wall.origin[1] + disp * wall.normal_y;
                wall.point_z = wall.origin[2] + disp * wall.normal_z;
                // Velocity = d(disp)/dt = amplitude * 2*pi*freq * cos(phase)
                let vel_mag = amplitude * 2.0 * std::f64::consts::PI * frequency * phase.cos();
                wall.velocity = [
                    vel_mag * wall.normal_x,
                    vel_mag * wall.normal_y,
                    vel_mag * wall.normal_z,
                ];
            }
            WallMotion::Servo {
                target_force,
                max_velocity,
                gain,
            } => {
                let error = target_force - wall.force_accumulator;
                let vel_mag = (gain * error).clamp(-max_velocity, max_velocity);
                wall.velocity = [
                    vel_mag * wall.normal_x,
                    vel_mag * wall.normal_y,
                    vel_mag * wall.normal_z,
                ];
                wall.point_x += wall.velocity[0] * dt;
                wall.point_y += wall.velocity[1] * dt;
                wall.point_z += wall.velocity[2] * dt;
            }
        }
    }

    walls.time += dt;
}

/// Zero all per-wall force accumulators before the force computation pass.
///
/// Runs in [`ParticleSimScheduleSet::PreForce`]. The accumulators are summed during
/// [`wall_contact_force`] and read by servo controllers in the next
/// [`wall_move`] call.
pub fn wall_zero_force_accumulators(mut walls: ResMut<Walls>) {
    for wall in &mut walls.planes {
        wall.force_accumulator = 0.0;
    }
    for wall in &mut walls.cylinders {
        wall.force_accumulator = 0.0;
    }
    for wall in &mut walls.spheres {
        wall.force_accumulator = 0.0;
    }
    for wall in &mut walls.regions {
        wall.force_accumulator = 0.0;
    }
}
