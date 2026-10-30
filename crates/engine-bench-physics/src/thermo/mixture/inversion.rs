//! Temperature from internal energy or enthalpy.
//!
//! Newton iteration warm-started from a guess (usually the previous step's
//! temperature) and clamped to `[TEMPERATURE_FLOOR, set maximum]`. u(T) and
//! h(T) increase monotonically (cv, cp > 0), so the clamped iteration
//! terminates without bracketing evaluations. A target outside the range
//! returns the nearest bound; callers validate the result against the
//! fit range.

use super::Mixture;
use crate::{
    thermo::{SpeciesSet, TEMPERATURE_FLOOR},
    units::{
        AvailableEnergy, ThermodynamicTemperature, available_energy::joule_per_kilogram,
        thermodynamic_temperature::kelvin,
    },
};

impl<S: SpeciesSet<N>, const N: usize> Mixture<S, N> {
    /// Temperature at which the specific internal energy equals `u`.
    pub fn temperature_from_internal_energy(
        &self,
        u: AvailableEnergy,
        guess: ThermodynamicTemperature,
    ) -> ThermodynamicTemperature {
        ThermodynamicTemperature::new::<kelvin>(self.temperature_from_internal_energy_si(
            u.get::<joule_per_kilogram>(),
            guess.get::<kelvin>(),
        ))
    }

    /// Temperature at which the specific enthalpy equals `h`.
    pub fn temperature_from_enthalpy(
        &self,
        h: AvailableEnergy,
        guess: ThermodynamicTemperature,
    ) -> ThermodynamicTemperature {
        ThermodynamicTemperature::new::<kelvin>(
            self.temperature_from_enthalpy_si(h.get::<joule_per_kilogram>(), guess.get::<kelvin>()),
        )
    }

    /// Kernel form of
    /// [`temperature_from_internal_energy`](Self::temperature_from_internal_energy).
    #[inline]
    pub(crate) fn temperature_from_internal_energy_si(&self, u: f64, guess: f64) -> f64 {
        newton(u, guess, S::TEMPERATURE_RANGE.1, |t| {
            self.internal_energy_and_cv_si(t)
        })
    }

    /// Kernel form of
    /// [`temperature_from_enthalpy`](Self::temperature_from_enthalpy).
    #[inline]
    pub(crate) fn temperature_from_enthalpy_si(&self, h: f64, guess: f64) -> f64 {
        newton(h, guess, S::TEMPERATURE_RANGE.1, |t| {
            self.enthalpy_and_cp_si(t)
        })
    }
}

/// Generous cap: from a warm start Newton needs 2-3 iterations, from a cold one
/// about 10.
const MAX_ITERATIONS: u32 = 50;
/// Newton converges quadratically and u(T), h(T) are only mildly curved
/// (|u''/u'| < 1e-3 per K), so after a step below 1e-7·T the new iterate is
/// already accurate to ~1e-15 relative and no confirming evaluation is needed.
const RELATIVE_STEP: f64 = 1e-7;

/// Solves `energy(t) = target` for a monotone increasing `energy` whose
/// derivative `value_and_slope(t).1` is positive, clamping iterates to
/// `[TEMPERATURE_FLOOR, t_max]`.
#[inline]
fn newton(target: f64, guess: f64, t_max: f64, value_and_slope: impl Fn(f64) -> (f64, f64)) -> f64 {
    let start = if guess.is_nan() { 300.0 } else { guess };
    let mut t = start.clamp(TEMPERATURE_FLOOR, t_max);
    for _ in 0..MAX_ITERATIONS {
        let (value, slope) = value_and_slope(t);
        let next = (t - (value - target) / slope).clamp(TEMPERATURE_FLOOR, t_max);
        if (next - t).abs() <= RELATIVE_STEP * t {
            return next;
        }
        t = next;
    }
    t
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::thermo::{Air3, Major7};

    fn major7(raw: [f64; 7]) -> Result<Mixture<Major7, 7>, TestCaseError> {
        let total: f64 = raw.iter().sum();
        Mixture::from_mass_fractions(raw.map(|v| v / total))
            .map_err(|e| TestCaseError::fail(format!("{e:#}")))
    }

    #[test]
    fn test_temperature_inversion_uom_api_round_trips() -> anyhow::Result<()> {
        let air = Mixture::dry_air()?;
        let t = ThermodynamicTemperature::new::<kelvin>(812.5);
        let guess = ThermodynamicTemperature::new::<kelvin>(300.0);
        let from_u = air.temperature_from_internal_energy(air.internal_energy(t), guess);
        let from_h = air.temperature_from_enthalpy(air.enthalpy(t), guess);
        for got in [from_u, from_h] {
            let k = got.get::<kelvin>();
            anyhow::ensure!((k - 812.5).abs() < 1e-8, "T = {k}");
        }
        Ok(())
    }

    #[test]
    fn test_temperature_inversion_clamps_out_of_range_targets() -> anyhow::Result<()> {
        let m = Mixture::<Major7, 7>::from_mass_fractions([0.7, 0.2, 0.0, 0.1, 0.0, 0.0, 0.0])?;
        let (_, t_max) = Major7::temperature_range_kelvin();
        let too_hot = m.internal_energy_si(t_max * 1.5);
        let t = m.temperature_from_internal_energy_si(too_hot, 1000.0);
        anyhow::ensure!((t - t_max).abs() < 1e-9, "hot target gave {t}");
        let too_cold = m.internal_energy_si(TEMPERATURE_FLOOR) - 1e6;
        let t = m.temperature_from_internal_energy_si(too_cold, 1000.0);
        anyhow::ensure!((t - TEMPERATURE_FLOOR).abs() < 1e-9, "cold target gave {t}");
        Ok(())
    }

    #[test]
    fn test_temperature_inversion_below_fit_range() -> anyhow::Result<()> {
        let air = Mixture::<Air3, 3>::dry_air()?;
        let t = air.temperature_from_internal_energy_si(air.internal_energy_si(80.0), 300.0);
        anyhow::ensure!((t - 80.0).abs() < 1e-8, "T = {t}");
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_temperature_inversion_round_trips(
            raw in proptest::array::uniform7(0.01..1.0f64),
            t in 200.0..6000.0f64,
            guess in 50.0..6000.0f64,
        ) {
            let m = major7(raw)?;
            let from_u = m.temperature_from_internal_energy_si(m.internal_energy_si(t), guess);
            let from_h = m.temperature_from_enthalpy_si(m.enthalpy_si(t), guess);
            prop_assert!((from_u / t - 1.0).abs() < 1e-10, "T(u) = {from_u}, T = {t}");
            prop_assert!((from_h / t - 1.0).abs() < 1e-10, "T(h) = {from_h}, T = {t}");
        }
    }
}
