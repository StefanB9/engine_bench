//! Nozzle operating point with a real equation of state, for every
//! back-pressure regime: the same classification as the perfect-gas nozzle,
//! along real-fluid isentropes and through real-gas shocks.

use anyhow::Context;

use super::{
    FlowState, RealFluid, Stagnation,
    isentrope::{Isentrope, Point, solve},
    normal_shock,
};
use crate::{
    PhysicsError,
    flow::{
        Nozzle, NozzleOperatingPoint, NozzleRegime,
        isentropic::Branch,
        nozzle::{BOUNDARY_TOLERANCE, SEPARATION_PRESSURE_RATIO},
    },
    units::{
        Force, MassRate, Pressure, Velocity, force::newton, mass_rate::kilogram_per_second,
        pressure::pascal, velocity::meter_per_second,
    },
};

impl Nozzle {
    /// Steady operating point for `fluid` from `chamber` into `back_pressure`
    /// (≥ 0; 0 is vacuum).
    pub fn real_operating_point<F: RealFluid>(
        &self,
        fluid: &F,
        chamber: &Stagnation,
        back_pressure: Pressure,
    ) -> anyhow::Result<NozzleOperatingPoint> {
        let pb = back_pressure.get::<pascal>();
        if !(pb >= 0.0 && pb.is_finite()) {
            return Err(PhysicsError::InvalidState(
                "back pressure must be non-negative",
            ))
            .with_context(|| format!("p_b = {pb} Pa"));
        }
        let p0 = chamber.p;
        let line = Isentrope {
            fluid,
            stagnation: *chamber,
        };
        if pb >= p0 {
            let at_rest = line.at(p0)?;
            return Ok(Nozzle::no_flow(p0, at_rest.state.t()));
        }
        let sonic = line.choked().context("nozzle throat")?;
        let epsilon = self.expansion_ratio();
        let g_exit = sonic.g / epsilon;
        let subsonic = line.with_flux(g_exit, &sonic, Branch::Subsonic)?;
        let supersonic = line.with_flux(g_exit, &sonic, Branch::Supersonic)?;
        // A shock at the exit of a nozzle without divergence is no shock.
        let exit_shock = if supersonic.u > supersonic.state.w() * (1.0 + 1e-9) {
            Some(normal_shock(
                fluid,
                &supersonic.state,
                velocity(supersonic.u),
            )?)
        } else {
            None
        };
        let p_shock_exit = exit_shock.map_or(supersonic.state.p(), |(down, _)| down.p());
        let ratio = pb / p0;
        let near = |p: f64| (ratio - p / p0).abs() <= BOUNDARY_TOLERANCE;
        let choked = self.throat_area * sonic.g;
        let (regime, exit, u_exit, mass_flow) = if pb > subsonic.state.p() {
            let exit = line.at(pb)?;
            (
                NozzleRegime::Subsonic,
                exit.state,
                exit.u,
                self.exit_area * exit.g,
            )
        } else if near(supersonic.state.p()) {
            (NozzleRegime::Ideal, supersonic.state, supersonic.u, choked)
        } else if let Some((down, u2)) = exit_shock.filter(|_| near(p_shock_exit)) {
            (
                NozzleRegime::ShockAtExit,
                down,
                u2.get::<meter_per_second>(),
                choked,
            )
        } else if pb > p_shock_exit {
            let (shock_area_ratio, exit, u) = shock_in_nozzle(&line, &sonic, epsilon, pb)?;
            (
                NozzleRegime::ShockInNozzle { shock_area_ratio },
                exit,
                u,
                choked,
            )
        } else if pb > supersonic.state.p() {
            (
                NozzleRegime::OverExpanded,
                supersonic.state,
                supersonic.u,
                choked,
            )
        } else {
            (
                NozzleRegime::UnderExpanded,
                supersonic.state,
                supersonic.u,
                choked,
            )
        };
        let exit_pressure = exit.p();
        let thrust = mass_flow * u_exit + (exit_pressure - pb) * self.exit_area;
        Ok(NozzleOperatingPoint {
            regime,
            mass_flow: MassRate::new::<kilogram_per_second>(mass_flow),
            exit_mach: u_exit / exit.w(),
            exit_pressure: exit.pressure,
            exit_temperature: exit.temperature,
            exit_velocity: velocity(u_exit),
            thrust: Force::new::<newton>(thrust),
            flow_separation: regime == NozzleRegime::OverExpanded
                && exit_pressure < SEPARATION_PRESSURE_RATIO * pb,
        })
    }
}

fn velocity(u: f64) -> Velocity {
    Velocity::new::<meter_per_second>(u)
}

/// Shock position (area ratio) in the divergent section for which the
/// subsonic flow behind it reaches the exit at `pb`; the exit state and
/// velocity. The exit pressure falls monotonically from the subsonic design
/// exit (shock at the throat) to the exit-shock pressure (shock at the exit).
fn shock_in_nozzle<F: RealFluid>(
    line: &Isentrope<'_, F>,
    sonic: &Point,
    epsilon: f64,
    pb: f64,
) -> anyhow::Result<(f64, FlowState, f64)> {
    let g_exit = sonic.g / epsilon;
    let exit_for = |area_ratio: f64| -> anyhow::Result<(FlowState, f64)> {
        let up = line.with_flux(sonic.g / area_ratio, sonic, Branch::Supersonic)?;
        let (down, u2) = normal_shock(line.fluid, &up.state, velocity(up.u))?;
        subsonic_exit(line, &down, u2.get::<meter_per_second>(), g_exit)
    };
    // Just past the throat the shock is vanishingly weak and the exit is the
    // subsonic one; at the exit plane it is the exit shock.
    let lo = 1.0 + 1e-9 * epsilon;
    let a = solve(
        |r| Ok(exit_for(r)?.0.p() - pb),
        lo,
        epsilon,
        1e-13 * epsilon,
    )
    .context("shock position")?;
    let (exit, u) = exit_for(a)?;
    Ok((a, exit, u))
}

/// Subsonic expansion behind a shock: along the post-shock isentrope (same
/// stagnation enthalpy, higher entropy) to where the mass flux is `g_exit`.
fn subsonic_exit<F: RealFluid>(
    line: &Isentrope<'_, F>,
    shocked: &FlowState,
    u_shocked: f64,
    g_exit: f64,
) -> anyhow::Result<(FlowState, f64)> {
    let behind = Isentrope {
        fluid: line.fluid,
        stagnation: Stagnation {
            p: f64::NAN,
            h: line.stagnation.h,
            s: shocked.s(),
        },
    };
    // The flux falls from the shock towards the post-shock stagnation
    // pressure. With the shock (to round-off) in the exit plane there is
    // nothing left to expand.
    let p_shock = shocked.p();
    if behind.at(p_shock)?.g <= g_exit {
        return Ok((*shocked, u_shocked));
    }
    let mut hi = p_shock;
    let mut step = 0.01 * p_shock;
    loop {
        hi += step;
        step *= 2.0;
        if behind.at(hi)?.g <= g_exit {
            break;
        }
        if step > 1e3 * p_shock {
            return Err(PhysicsError::InvalidState("no subsonic exit behind the shock").into());
        }
    }
    let p = solve(|p| Ok(behind.at(p)?.g - g_exit), p_shock, hi, 1e-14 * hi)
        .context("subsonic exit behind the shock")?;
    let exit = behind.at(p)?;
    Ok((exit.state, exit.u))
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::{
        flow::{
            Reservoir,
            real::perfect::tests::{air, state},
        },
        units::{
            Area, ThermodynamicTemperature, area::square_meter, thermodynamic_temperature::kelvin,
        },
    };

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
    }

    fn same_regime(a: NozzleRegime, b: NozzleRegime) -> bool {
        std::mem::discriminant(&a) == std::mem::discriminant(&b)
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_real_nozzle_equals_the_perfect_gas_nozzle(epsilon in 1.0f64..12.0, ratio in 0.0f64..1.05) {
            let gas = air().map_err(|e| fail(&e))?;
            let nozzle = Nozzle::new(Area::new::<square_meter>(1e-4), Area::new::<square_meter>(epsilon * 1e-4)).map_err(|e| fail(&e))?;
            let (p0, t0) = (2e6, 900.0);
            let chamber = Stagnation::at_rest(&state(&gas, p0, t0).map_err(|e| fail(&e))?);
            let pb = Pressure::new::<pascal>(ratio * p0);
            let real = nozzle.real_operating_point(&gas, &chamber, pb).map_err(|e| fail(&e))?;
            let reservoir = Reservoir::new(Pressure::new::<pascal>(p0), ThermodynamicTemperature::new::<kelvin>(t0), gas).map_err(|e| fail(&e))?;
            let perfect = nozzle.operating_point(&reservoir, pb).map_err(|e| fail(&e))?;
            prop_assert!(same_regime(real.regime, perfect.regime), "{:?} vs {:?}", real.regime, perfect.regime);
            let close = |a: f64, b: f64, tol: f64| (a - b).abs() <= tol * b.abs().max(1e-9);
            prop_assert!(close(real.mass_flow.get::<kilogram_per_second>(), perfect.mass_flow.get::<kilogram_per_second>(), 1e-7));
            prop_assert!(close(real.thrust.get::<newton>(), perfect.thrust.get::<newton>(), 1e-6), "F {:?} vs {:?}", real.thrust, perfect.thrust);
            prop_assert!(close(real.exit_pressure.get::<pascal>(), perfect.exit_pressure.get::<pascal>(), 1e-6));
            prop_assert!(close(real.exit_mach, perfect.exit_mach, 1e-6), "M {} vs {}", real.exit_mach, perfect.exit_mach);
        }
    }
}
