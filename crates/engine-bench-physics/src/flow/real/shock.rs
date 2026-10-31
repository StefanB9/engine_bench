//! Normal shock with a real equation of state: the Rankine–Hugoniot jump
//! conditions (mass, momentum and energy conserved) solved for the downstream
//! density.

use anyhow::Context;

use super::{FlowState, RealFluid, isentrope::solve};
use crate::{
    PhysicsError,
    units::{
        MassDensity, Pressure, Velocity, mass_density::kilogram_per_cubic_meter, pressure::pascal,
        velocity::meter_per_second,
    },
};

/// Downstream state and velocity of a normal shock in `fluid` with upstream
/// static state `upstream` moving at supersonic `velocity`.
///
/// With mass flux m = ρ₁u₁, momentum gives p₂ = p₁ + m(u₁ − u₂) and energy
/// h₂ = h₁ + (u₁² − u₂²)/2, so the downstream density ρ₂ = m/u₂ alone must
/// make the equation of state's h(p₂, ρ₂) match. ρ₂ = ρ₁ (no shock) solves
/// that too; the bracket is seeded from the perfect-gas jump with γ from the
/// local speed of sound to stay on the compressive root.
pub fn normal_shock<F: RealFluid>(
    fluid: &F,
    upstream: &FlowState,
    velocity: Velocity,
) -> anyhow::Result<(FlowState, Velocity)> {
    let u1 = velocity.get::<meter_per_second>();
    let (p1, rho1, h1, w1) = (upstream.p(), upstream.rho(), upstream.h(), upstream.w());
    if !(u1 > w1 && u1.is_finite()) {
        return Err(PhysicsError::InvalidState(
            "a normal shock needs supersonic upstream flow",
        ))
        .with_context(|| format!("u = {u1} m/s, w = {w1} m/s"));
    }
    let m = rho1 * u1;
    let state = |rho2: f64| {
        let u2 = m / rho2;
        fluid.state_prho(
            Pressure::new::<pascal>(p1 + m * (u1 - u2)),
            MassDensity::new::<kilogram_per_cubic_meter>(rho2),
        )
    };
    let residual = |rho2: f64| -> anyhow::Result<f64> {
        let u2 = m / rho2;
        Ok(state(rho2)?.h() - (h1 + 0.5 * (u1 * u1 - u2 * u2)))
    };
    let gamma = w1 * w1 * rho1 / p1;
    let mach2 = (u1 / w1).powi(2);
    let guess = (gamma + 1.0) * mach2 / ((gamma - 1.0) * mach2 + 2.0);
    let (mut lo, mut hi) = (1.0 + 0.25 * (guess - 1.0), 1.0 + 2.0 * (guess - 1.0));
    let positive_lo = residual(lo * rho1)? > 0.0;
    let mut tries = 0;
    while (residual(hi * rho1)? > 0.0) == positive_lo {
        tries += 1;
        if tries > 60 {
            return Err(PhysicsError::InvalidState("no compressive shock solution").into());
        }
        // Alternate: push the upper end out, pull the lower one towards ρ₁.
        if tries % 2 == 1 {
            hi = 1.0 + 2.0 * (hi - 1.0);
        } else {
            lo = 1.0 + 0.25 * (lo - 1.0);
            if (residual(lo * rho1)? > 0.0) != positive_lo {
                hi = lo;
                lo = 1.0 + 0.25 * (lo - 1.0);
            }
        }
    }
    let rho2 = solve(residual, lo * rho1, hi * rho1, 1e-14 * rho1).context("shock density")?;
    Ok((state(rho2)?, Velocity::new::<meter_per_second>(m / rho2)))
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::flow::{
        real::perfect::tests::{air, state},
        shock,
    };

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_real_shock_matches_perfect_gas_relations_and_conserves(mach in 1.05f64..8.0, t in 100.0f64..1500.0) {
            let gas = air().map_err(|e| fail(&e))?;
            let up = state(&gas, 1e5, t).map_err(|e| fail(&e))?;
            let u1 = mach * up.w();
            let (down, u2) = normal_shock(&gas, &up, Velocity::new::<meter_per_second>(u1)).map_err(|e| fail(&e))?;
            let u2 = u2.get::<meter_per_second>();
            let jump = shock::normal_shock(1.4, mach);
            prop_assert!((down.p() / up.p() / jump.pressure_ratio - 1.0).abs() < 1e-9, "p2/p1 {} vs {}", down.p() / up.p(), jump.pressure_ratio);
            prop_assert!((u2 / down.w() / jump.downstream_mach - 1.0).abs() < 1e-9, "M2");
            // Mass, momentum and energy across the shock.
            let (m1, m2) = (up.rho() * u1, down.rho() * u2);
            prop_assert!((m2 / m1 - 1.0).abs() < 1e-12);
            prop_assert!(((down.p() + m2 * u2) / (up.p() + m1 * u1) - 1.0).abs() < 1e-12);
            prop_assert!(((down.h() + u2 * u2 / 2.0) / (up.h() + u1 * u1 / 2.0) - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_real_shock_needs_supersonic_upstream() -> anyhow::Result<()> {
        let gas = air()?;
        let up = state(&gas, 1e5, 300.0)?;
        let result = normal_shock(&gas, &up, Velocity::new::<meter_per_second>(0.9 * up.w()));
        anyhow::ensure!(result.is_err(), "subsonic shock accepted");
        Ok(())
    }
}
