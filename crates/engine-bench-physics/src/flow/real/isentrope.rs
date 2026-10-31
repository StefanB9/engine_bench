//! Isentropic expansion from a stagnation state: the static state, velocity
//! and mass flux at each pressure, the choked (sonic) point, and the pressure
//! where the flux takes a given value on either branch.

use anyhow::Context;

use super::{FlowState, RealFluid, Stagnation};
use crate::{
    PhysicsError,
    flow::isentropic::Branch,
    numerics::roots::{RootOptions, brent},
    units::{
        Pressure, SpecificHeatCapacity, pressure::pascal,
        specific_heat_capacity::joule_per_kilogram_kelvin,
    },
};

/// A point on an isentrope.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Point {
    pub(crate) state: FlowState,
    /// Flow velocity, m/s.
    pub(crate) u: f64,
    /// Mass flux ρ·u, kg/(m²·s).
    pub(crate) g: f64,
}

/// The isentrope through a stagnation state of `fluid`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Isentrope<'a, F: RealFluid> {
    pub(crate) fluid: &'a F,
    pub(crate) stagnation: Stagnation,
}

impl<F: RealFluid> Isentrope<'_, F> {
    /// The point at pressure `p` Pa (≤ stagnation pressure).
    pub(crate) fn at(&self, p: f64) -> anyhow::Result<Point> {
        let state = self.fluid.state_ps(
            Pressure::new::<pascal>(p),
            SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(self.stagnation.s),
        )?;
        let u = (2.0 * (self.stagnation.h - state.h())).max(0.0).sqrt();
        Ok(Point {
            state,
            u,
            g: state.rho() * u,
        })
    }

    /// The sonic point, where the mass flux is greatest: u = w. Inside the dome
    /// the equilibrium speed of sound jumps at the phase boundary, and the
    /// root sits on the jump, which is then the flux maximum.
    pub(crate) fn choked(&self) -> anyhow::Result<Point> {
        let p0 = self.stagnation.p;
        let excess = |p: f64| -> anyhow::Result<f64> {
            let x = self.at(p)?;
            Ok(x.u * x.u - x.state.w() * x.state.w())
        };
        // Expand until the flow is supersonic.
        let (mut hi, mut lo) = (p0, 0.5 * p0);
        while excess(lo)? <= 0.0 {
            if lo < MIN_PRESSURE_RATIO * p0 {
                return Err(PhysicsError::InvalidState(
                    "the isentrope never reaches sonic speed",
                ))
                .with_context(|| format!("expanding from {p0} Pa"));
            }
            (hi, lo) = (lo, 0.5 * lo);
        }
        let p = solve(excess, lo, hi, 1e-14 * p0).context("sonic point")?;
        self.at(p)
    }

    /// The point on `branch` where the mass flux is `g` (≤ the choked flux of
    /// `sonic`).
    pub(crate) fn with_flux(&self, g: f64, sonic: &Point, branch: Branch) -> anyhow::Result<Point> {
        if g >= sonic.g {
            return Ok(*sonic);
        }
        let p_star = sonic.state.p();
        let (lo, hi) = match branch {
            // G falls from G* at the throat to 0 at stagnation.
            Branch::Subsonic => (p_star, self.stagnation.p),
            Branch::Supersonic => {
                let mut lo = 0.5 * p_star;
                while self.at(lo)?.g > g {
                    if lo < MIN_PRESSURE_RATIO * self.stagnation.p {
                        return Err(PhysicsError::InvalidState(
                            "mass flux not reached on the supersonic branch",
                        )
                        .into());
                    }
                    lo *= 0.5;
                }
                (lo, p_star)
            }
        };
        let p = solve(|p| Ok(self.at(p)?.g - g), lo, hi, 1e-14 * self.stagnation.p)
            .with_context(|| format!("{branch:?} point with mass flux {g} kg/(m² s)"))?;
        self.at(p)
    }
}

/// Lowest pressure, relative to stagnation, an expansion is followed to.
const MIN_PRESSURE_RATIO: f64 = 1e-12;

/// Brent on [lo, hi] where evaluating f can fail.
pub(crate) fn solve(
    mut f: impl FnMut(f64) -> anyhow::Result<f64>,
    lo: f64,
    hi: f64,
    x_tol: f64,
) -> anyhow::Result<f64> {
    let mut failure = None;
    let x = brent(
        |x| {
            f(x).unwrap_or_else(|e| {
                failure.get_or_insert(e);
                f64::NAN
            })
        },
        lo,
        hi,
        RootOptions::new(x_tol, 300)?,
    );
    if let Some(e) = failure {
        return Err(e);
    }
    x
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::flow::{
        isentropic,
        real::perfect::tests::{air, state},
    };

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
    }

    #[test]
    fn test_isentrope_chokes_at_the_perfect_gas_critical_ratio() -> anyhow::Result<()> {
        let gas = air()?;
        let line = Isentrope {
            fluid: &gas,
            stagnation: Stagnation::at_rest(&state(&gas, 5e5, 300.0)?),
        };
        let sonic = line.choked()?;
        let ratio = sonic.state.p() / 5e5;
        let critical = isentropic::critical_pressure_ratio(1.4);
        anyhow::ensure!(
            (ratio / critical - 1.0).abs() < 1e-9,
            "p*/p0 = {ratio} vs {critical}"
        );
        anyhow::ensure!(
            (sonic.u / sonic.state.w() - 1.0).abs() < 1e-8,
            "M = {}",
            sonic.u / sonic.state.w()
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_isentrope_flux_matches_the_area_mach_relation(epsilon in 1.01f64..20.0, supersonic in any::<bool>()) {
            let gas = air().map_err(|e| fail(&e))?;
            let line = Isentrope { fluid: &gas, stagnation: Stagnation::at_rest(&state(&gas, 5e5, 300.0).map_err(|e| fail(&e))?) };
            let sonic = line.choked().map_err(|e| fail(&e))?;
            let branch = if supersonic { Branch::Supersonic } else { Branch::Subsonic };
            let point = line.with_flux(sonic.g / epsilon, &sonic, branch).map_err(|e| fail(&e))?;
            let mach = isentropic::mach_from_area_ratio(1.4, epsilon, branch).map_err(|e| fail(&e))?;
            let ratio = isentropic::pressure_ratio(1.4, mach);
            prop_assert!((point.state.p() / 5e5 / ratio - 1.0).abs() < 1e-8, "p/p0 {} vs {}", point.state.p() / 5e5, ratio);
            prop_assert!((point.u / point.state.w() / mach - 1.0).abs() < 1e-7, "M {} vs {}", point.u / point.state.w(), mach);
        }
    }
}
