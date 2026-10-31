//! Planar (3-DOF) rigid-body dynamics: horizontal and vertical translation plus
//! pitch.
//!
//! Kernel types in SI. x is horizontal, y up; θ is the body axis angle from +x,
//! counter-clockwise, so a vertical rocket has θ = π/2. Mass properties are
//! supplied per evaluation because they change as propellant leaves; thrust and
//! other forces come in as [`Loads`].

/// Position, velocity, attitude and pitch rate.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PlanarState {
    /// Horizontal position, m.
    pub x: f64,
    /// Height, m.
    pub y: f64,
    /// Horizontal velocity, m/s.
    pub vx: f64,
    /// Vertical velocity, m/s.
    pub vy: f64,
    /// Body axis angle from +x, rad.
    pub theta: f64,
    /// Pitch rate, rad/s.
    pub omega: f64,
}

impl PlanarState {
    /// As an ODE state `[x, y, vx, vy, θ, ω]`.
    pub fn to_array(self) -> [f64; 6] {
        [self.x, self.y, self.vx, self.vy, self.theta, self.omega]
    }

    /// From an ODE state `[x, y, vx, vy, θ, ω]`.
    pub fn from_array(state: &[f64; 6]) -> Self {
        let [x, y, vx, vy, theta, omega] = *state;
        Self {
            x,
            y,
            vx,
            vy,
            theta,
            omega,
        }
    }

    /// Unit vector along the body axis, nose-wards.
    pub fn axis(self) -> [f64; 2] {
        let (sin, cos) = self.theta.sin_cos();
        [cos, sin]
    }
}

/// Mass, pitch inertia about the centre of gravity, and where the centre of
/// gravity is.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassProperties {
    /// kg, > 0.
    pub mass: f64,
    /// kg·m², > 0.
    pub inertia: f64,
    /// Centre of gravity, m aft of the nose tip.
    pub cg_from_nose: f64,
}

/// Forces through the centre of gravity and a pitching moment about it.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Loads {
    /// Force, N, world frame `[x, y]`.
    pub force: [f64; 2],
    /// Pitching moment, N·m, counter-clockwise positive.
    pub moment: f64,
}

impl Loads {
    /// Sum of two sets of loads.
    #[must_use]
    pub fn plus(self, other: Loads) -> Loads {
        Loads {
            force: [
                self.force[0] + other.force[0],
                self.force[1] + other.force[1],
            ],
            moment: self.moment + other.moment,
        }
    }
}

/// Time derivative of `[x, y, vx, vy, θ, ω]` under `loads` and gravity
/// `gravity` (m/s², acting along −y).
#[inline]
pub fn derivatives(
    state: &PlanarState,
    mass: &MassProperties,
    loads: &Loads,
    gravity: f64,
) -> [f64; 6] {
    [
        state.vx,
        state.vy,
        loads.force[0] / mass.mass,
        loads.force[1] / mass.mass - gravity,
        state.omega,
        loads.moment / mass.inertia,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_planar_state_array_round_trip_and_axis() -> anyhow::Result<()> {
        let state = PlanarState {
            x: 1.0,
            y: 2.0,
            vx: 3.0,
            vy: 4.0,
            theta: std::f64::consts::FRAC_PI_2,
            omega: 6.0,
        };
        let expected = [1.0, 2.0, 3.0, 4.0, std::f64::consts::FRAC_PI_2, 6.0];
        anyhow::ensure!(
            state
                .to_array()
                .iter()
                .zip(expected)
                .all(|(a, b)| a.to_bits() == b.to_bits()),
            "array"
        );
        anyhow::ensure!(
            PlanarState::from_array(&state.to_array()) == state,
            "round trip"
        );
        let [ax, ay] = state.axis();
        anyhow::ensure!(
            ax.abs() < 1e-15 && (ay - 1.0).abs() < 1e-15,
            "vertical axis {ax} {ay}"
        );
        Ok(())
    }

    #[test]
    fn test_derivatives_apply_forces_moment_and_gravity() -> anyhow::Result<()> {
        let state = PlanarState {
            x: 0.0,
            y: 10.0,
            vx: 2.0,
            vy: -1.0,
            theta: 0.3,
            omega: 0.5,
        };
        let mass = MassProperties {
            mass: 2.0,
            inertia: 0.5,
            cg_from_nose: 0.3,
        };
        let loads = Loads {
            force: [4.0, 6.0],
            moment: 1.0,
        };
        let d = derivatives(&state, &mass, &loads, 9.81);
        let expected = [2.0, -1.0, 2.0, 3.0 - 9.81, 0.5, 2.0];
        anyhow::ensure!(
            d.iter().zip(expected).all(|(a, b)| (a - b).abs() < 1e-12),
            "{d:?}"
        );
        let sum = loads.plus(Loads {
            force: [1.0, -1.0],
            moment: -0.5,
        });
        anyhow::ensure!(
            sum == Loads {
                force: [5.0, 5.0],
                moment: 0.5
            },
            "{sum:?}"
        );
        Ok(())
    }
}
