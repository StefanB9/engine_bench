//! Orifice flow with a real equation of state: isentropic from the upstream
//! stagnation state to the downstream pressure (or the choked throat), with
//! the same C¹ quadratic as the perfect-gas orifice near pressure equality.

use anyhow::Context;

use super::{RealFluid, Stagnation, isentrope::Isentrope};
use crate::{
    PhysicsError,
    flow::{Orifice, orifice::LINEAR_BAND},
    units::{MassRate, Pressure, mass_rate::kilogram_per_second, pressure::pascal},
};

impl Orifice {
    /// Mass flow out of a reservoir of `fluid` at `upstream` into `p_down`
    /// (≤ upstream pressure), along the real-fluid isentrope.
    pub fn real_mass_flow<F: RealFluid>(
        &self,
        fluid: &F,
        upstream: &Stagnation,
        p_down: Pressure,
    ) -> anyhow::Result<MassRate> {
        let (p0, pd) = (upstream.p, p_down.get::<pascal>());
        if !(pd >= 0.0 && pd <= p0) {
            return Err(PhysicsError::InvalidState(
                "downstream pressure must lie between 0 and the upstream pressure",
            ))
            .with_context(|| format!("p_down = {pd} Pa, p_up = {p0} Pa"));
        }
        let line = Isentrope {
            fluid,
            stagnation: *upstream,
        };
        let sonic = line.choked().context("orifice choking")?;
        let area = self.effective_area;
        let x = 1.0 - pd / p0;
        let flow = if x >= LINEAR_BAND {
            if pd <= sonic.state.p() {
                area * sonic.g
            } else {
                area * line.at(pd)?.g
            }
        } else {
            // f(x) = a·x + b·x² matching value and slope at the band edge; on
            // an isentrope dG/dp = u/w² − 1/u.
            let edge = line.at(p0 * (1.0 - LINEAR_BAND))?;
            let value = area * edge.g;
            let w2 = edge.state.w().powi(2);
            let slope = area * p0 * (1.0 / edge.u - edge.u / w2);
            let b = (slope * LINEAR_BAND - value) / (LINEAR_BAND * LINEAR_BAND);
            let a = slope - 2.0 * b * LINEAR_BAND;
            x * (a + b * x)
        };
        Ok(MassRate::new::<kilogram_per_second>(flow))
    }
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

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_real_orifice_equals_the_perfect_gas_orifice(ratio in 0.0f64..1.0, t in 150.0f64..900.0) {
            // Unchoked, choked and inside the near-equality band.
            let gas = air().map_err(|e| fail(&e))?;
            let orifice = Orifice::new(Area::new::<square_meter>(1e-4), 0.8).map_err(|e| fail(&e))?;
            let p0 = 5e5;
            let up = Stagnation::at_rest(&state(&gas, p0, t).map_err(|e| fail(&e))?);
            let real = orifice
                .real_mass_flow(&gas, &up, Pressure::new::<pascal>(ratio * p0))
                .map_err(|e| fail(&e))?
                .get::<kilogram_per_second>();
            let a = Reservoir::new(Pressure::new::<pascal>(p0), ThermodynamicTemperature::new::<kelvin>(t), gas).map_err(|e| fail(&e))?;
            let b = Reservoir::new(Pressure::new::<pascal>((ratio * p0).max(1.0)), ThermodynamicTemperature::new::<kelvin>(t), gas).map_err(|e| fail(&e))?;
            let perfect = orifice.mass_flow(&a, &b).get::<kilogram_per_second>();
            prop_assert!((real - perfect).abs() <= 1e-8 * perfect.abs().max(1e-6), "{} vs {} at p/p0 = {}", real, perfect, ratio);
        }
    }

    #[test]
    fn test_real_orifice_has_no_flow_at_equal_pressure() -> anyhow::Result<()> {
        let gas = air()?;
        let orifice = Orifice::new(Area::new::<square_meter>(1e-4), 0.8)?;
        let up = Stagnation::at_rest(&state(&gas, 5e5, 300.0)?);
        let flow = orifice.real_mass_flow(&gas, &up, Pressure::new::<pascal>(5e5))?;
        anyhow::ensure!(flow.get::<kilogram_per_second>().abs() < 1e-15, "{flow:?}");
        Ok(())
    }
}
