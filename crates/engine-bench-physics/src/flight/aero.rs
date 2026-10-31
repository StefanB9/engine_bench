//! Quasi-steady aerodynamic loads of a slender rocket in planar flight.
//!
//! Each Barrowman component produces a normal force from the air's velocity
//! normal to the body *at that component*: the rocket's own velocity plus the
//! rotation ω about the centre of gravity. That one formula yields both the
//! static restoring moment and the aerodynamic pitch damping −½ρ·V·A·ω·Σ
//! `C_Nα·d²` (d = component distance aft of the CG). Axial drag uses a drag
//! coefficient on the same reference area.

use anyhow::Context;

use super::{
    barrowman::Stability,
    rigid_body::{Loads, PlanarState},
};
use crate::PhysicsError;

/// Aerodynamic description of a rocket.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AeroModel {
    stability: Stability,
    axial_coefficient: f64,
}

impl AeroModel {
    /// Barrowman `stability` plus an axial drag coefficient (≥ 0) on its
    /// reference area.
    pub fn new(stability: Stability, axial_coefficient: f64) -> anyhow::Result<Self> {
        if !(axial_coefficient >= 0.0 && axial_coefficient.is_finite()) {
            return Err(PhysicsError::InvalidState(
                "axial drag coefficient must be non-negative",
            ))
            .with_context(|| format!("C_A = {axial_coefficient}"));
        }
        Ok(Self {
            stability,
            axial_coefficient,
        })
    }
}

/// Aerodynamic loads on a rocket in `state`, with its centre of gravity
/// `cg_from_nose` m aft of the nose, in air of `density` kg/m³ moving
/// horizontally at `wind` m/s.
#[inline]
pub fn aerodynamic_loads(
    state: &PlanarState,
    aero: &AeroModel,
    cg_from_nose: f64,
    density: f64,
    wind: f64,
) -> Loads {
    let [ex, ey] = state.axis();
    // Normal unit vector, 90° counter-clockwise from the axis.
    let (nx, ny) = (-ey, ex);
    let (rx, ry) = (state.vx - wind, state.vy);
    let speed = rx.hypot(ry);
    let half_rho_area = 0.5 * density * aero.stability.reference_area();
    let axial_velocity = rx * ex + ry * ey;
    let normal_velocity = rx * nx + ry * ny;

    // Forces proportional to V·u keep their magnitude bounded at any angle of
    // attack.
    let axial = -half_rho_area * aero.axial_coefficient * speed * axial_velocity;
    let mut normal = 0.0;
    let mut moment = 0.0;
    for component in &aero.stability.components {
        // Distance aft of the CG; rotation moves this point by −d·ω along the
        // normal.
        let lever = component.centre_of_pressure - cg_from_nose;
        let local_normal_velocity = normal_velocity - lever * state.omega;
        let force = -half_rho_area * component.normal_force_slope * speed * local_normal_velocity;
        normal += force;
        moment -= lever * force;
    }
    Loads {
        force: [axial * ex + normal * nx, axial * ey + normal * ny],
        moment,
    }
}

/// Jet damping moment, N·m: exhaust leaving at `nozzle_from_nose` with mass
/// flow `mass_flow` (≥ 0) carries away the sideways momentum the rotation gave
/// it.
#[inline]
pub fn jet_damping(mass_flow: f64, nozzle_from_nose: f64, cg_from_nose: f64, omega: f64) -> f64 {
    let lever = nozzle_from_nose - cg_from_nose;
    -mass_flow * lever * lever * omega
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::{
        flight::barrowman::{Fins, NoseShape, RocketGeometry, barrowman},
        units::{Length, length::meter},
    };

    fn m(value: f64) -> Length {
        Length::new::<meter>(value)
    }

    /// A 1 m rocket, 0.1 m diameter, with fins at the tail: CP ≈ 0.8 m.
    pub(crate) fn rocket() -> anyhow::Result<AeroModel> {
        let geometry = RocketGeometry {
            nose_shape: NoseShape::Ogive,
            nose_length: m(0.2),
            body_diameter: m(0.1),
            fins: Some(Fins {
                count: 4,
                root_chord: m(0.12),
                tip_chord: m(0.06),
                semi_span: m(0.08),
                sweep: m(0.06),
                position: m(0.88),
            }),
        };
        AeroModel::new(barrowman(&geometry)?, 0.4)
    }

    #[test]
    fn test_aero_model_validates_drag_coefficient() -> anyhow::Result<()> {
        let stability = barrowman(&RocketGeometry {
            nose_shape: NoseShape::Conical,
            nose_length: m(0.2),
            body_diameter: m(0.1),
            fins: None,
        })?;
        for bad in [-0.1, f64::NAN] {
            anyhow::ensure!(
                AeroModel::new(stability, bad).is_err(),
                "accepted C_A {bad}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_aero_loads_vanish_in_still_air_at_rest() -> anyhow::Result<()> {
        let loads = aerodynamic_loads(
            &PlanarState {
                theta: 1.0,
                ..Default::default()
            },
            &rocket()?,
            0.5,
            1.2,
            0.0,
        );
        anyhow::ensure!(loads == Loads::default(), "{loads:?}");
        Ok(())
    }

    #[test]
    fn test_aero_axial_drag_when_aligned() -> anyhow::Result<()> {
        // Flying straight up at 20 m/s with the nose up: pure drag, down.
        let state = PlanarState {
            vy: 20.0,
            theta: std::f64::consts::FRAC_PI_2,
            ..Default::default()
        };
        let loads = aerodynamic_loads(&state, &rocket()?, 0.5, 1.2, 0.0);
        let area = std::f64::consts::FRAC_PI_4 * 0.01;
        let drag = 0.5 * 1.2 * 400.0 * area * 0.4;
        anyhow::ensure!(
            loads.force[0].abs() < 1e-12 && (loads.force[1] + drag).abs() < 1e-12,
            "{loads:?}"
        );
        anyhow::ensure!(loads.moment.abs() < 1e-12, "moment {}", loads.moment);
        Ok(())
    }

    #[test]
    fn test_aero_restores_a_stable_rocket_and_upsets_an_unstable_one() -> anyhow::Result<()> {
        // Moving along +x with the nose pitched up 5°.
        let state = PlanarState {
            vx: 30.0,
            theta: 5f64.to_radians(),
            ..Default::default()
        };
        let aero = rocket()?;
        let stable = aerodynamic_loads(&state, &aero, 0.4, 1.2, 0.0);
        anyhow::ensure!(
            stable.moment < 0.0,
            "stable rocket not pushed nose-down: {}",
            stable.moment
        );
        let unstable = aerodynamic_loads(&state, &aero, 0.95, 1.2, 0.0);
        anyhow::ensure!(
            unstable.moment > 0.0,
            "unstable rocket not pushed nose-up: {}",
            unstable.moment
        );
        Ok(())
    }

    #[test]
    fn test_jet_damping_opposes_rotation() -> anyhow::Result<()> {
        let moment = jet_damping(0.5, 1.0, 0.6, 2.0);
        anyhow::ensure!((moment + 0.5 * 0.16 * 2.0).abs() < 1e-12, "{moment}");
        anyhow::ensure!(
            jet_damping(0.0, 1.0, 0.6, 2.0) == 0.0,
            "no flow, no damping"
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_aero_damping_opposes_pitch_rate(
            speed in 1.0..100.0f64, theta in -3.0..3.0f64, omega in -20.0..20.0f64,
            cg in 0.1..0.95f64,
        ) {
            // Velocity along the body axis: only the rotation makes a normal velocity.
            let state = PlanarState {
                vx: speed * theta.cos(), vy: speed * theta.sin(), theta, omega, ..Default::default()
            };
            let aero = rocket().map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let loads = aerodynamic_loads(&state, &aero, cg, 1.2, 0.0);
            prop_assert!(loads.moment * omega <= 1e-12, "moment {} with ω {omega}", loads.moment);
        }
    }
}
