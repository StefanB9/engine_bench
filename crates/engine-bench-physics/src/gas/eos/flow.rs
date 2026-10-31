//! A volume's equation of state and composition as a [`RealFluid`], so its
//! contents can expand through real-fluid orifices and nozzles.

use super::EquationOfState;
use crate::{
    flow::real::{FlowState, RealFluid},
    thermo::{Mixture, SpeciesSet},
    units::{
        MassDensity, Pressure, SpecificHeatCapacity, mass_density::kilogram_per_cubic_meter,
        pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
    },
};

/// `eos` applied to `mixture`, as a [`RealFluid`].
#[derive(Debug, Clone, Copy)]
pub struct EosFluid<'a, S: SpeciesSet<N>, const N: usize, E: EquationOfState<S, N>> {
    /// The equation of state.
    pub eos: &'a E,
    /// The composition flowing.
    pub mixture: &'a Mixture<S, N>,
}

impl<S: SpeciesSet<N>, const N: usize, E: EquationOfState<S, N>> RealFluid
    for EosFluid<'_, S, N, E>
{
    fn state_ps(&self, p: Pressure, s: SpecificHeatCapacity) -> anyhow::Result<FlowState> {
        self.eos.flow_state_ps(
            self.mixture,
            p.get::<pascal>(),
            s.get::<joule_per_kilogram_kelvin>(),
        )
    }

    fn state_prho(&self, p: Pressure, rho: MassDensity) -> anyhow::Result<FlowState> {
        self.eos.flow_state_prho(
            self.mixture,
            p.get::<pascal>(),
            rho.get::<kilogram_per_cubic_meter>(),
        )
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::{
        fluid::PureFluid,
        gas::eos::{IdealGas, PengRobinson, Reference},
        thermo::{Air3, ColdGas},
        units::{thermodynamic_temperature::kelvin, velocity::meter_per_second},
    };

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
    }

    /// Round trip (p, T) → s → `state_ps` and → ρ → `state_prho`; the speed of
    /// sound against √(Δp/Δρ) along the isentrope.
    fn check<S: SpeciesSet<N>, const N: usize, E: EquationOfState<S, N>>(
        eos: &E,
        mixture: &Mixture<S, N>,
        p: f64,
        t: f64,
    ) -> Result<(), TestCaseError> {
        let fluid = EosFluid { eos, mixture };
        let state = eos.state_pt(mixture, p, t).map_err(|e| fail(&e))?;
        let s = eos.entropy(mixture, &state);
        let pressure = Pressure::new::<pascal>(p);
        let by_s = fluid
            .state_ps(
                pressure,
                SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(s),
            )
            .map_err(|e| fail(&e))?;
        let by_rho = fluid
            .state_prho(
                pressure,
                MassDensity::new::<kilogram_per_cubic_meter>(1.0 / state.v),
            )
            .map_err(|e| fail(&e))?;
        for x in [by_s, by_rho] {
            prop_assert!(
                (x.temperature.get::<kelvin>() / t - 1.0).abs() < 1e-9,
                "T {:?} vs {}",
                x.temperature,
                t
            );
        }
        let dp = 1e-5 * p;
        let up = fluid
            .state_ps(Pressure::new::<pascal>(p + dp), by_s.entropy)
            .map_err(|e| fail(&e))?;
        let down = fluid
            .state_ps(Pressure::new::<pascal>(p - dp), by_s.entropy)
            .map_err(|e| fail(&e))?;
        let w = (2.0 * dp / (up.density - down.density).get::<kilogram_per_cubic_meter>()).sqrt();
        let got = by_s.speed_of_sound.get::<meter_per_second>();
        prop_assert!((got / w - 1.0).abs() < 1e-5, "w {} vs {}", got, w);
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_eos_flow_states_are_consistent(t in 150.0f64..900.0, log_p in 4.0f64..7.4) {
            let p = 10f64.powf(log_p);
            let air = Mixture::<Air3, 3>::dry_air().map_err(|e| fail(&e))?;
            check(&IdealGas, &air, p, t)?;
            check(&PengRobinson::<Air3, 3>::new().map_err(|e| fail(&e))?, &air, p, t)?;
            let n2 = Mixture::<ColdGas, 4>::from_mass_fractions([1.0, 0.0, 0.0, 0.0]).map_err(|e| fail(&e))?;
            check(&Reference::<ColdGas, 4>::new(PureFluid::Nitrogen).map_err(|e| fail(&e))?, &n2, p, t)?;
        }
    }
}
