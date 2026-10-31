//! Temperature from (p, h) and (p, s), exactly consistent with the forward
//! equations: the IF97 backward equations give the starting point, and a
//! bracketed Newton solve on the forward equations removes their millikelvin
//! error. Inside the two-phase dome the temperature is the saturation
//! temperature.

use anyhow::Context;

use super::{
    P_MAX, P_MAX_HIGH_TEMPERATURE, T_MAX, T_MAX_VAPOUR, T_MIN, T_REGION3_MIN, backward,
    gibbs::{self, RawProperties},
    region::{Region, region_si},
    saturation,
};
use crate::{
    PhysicsError,
    numerics::roots::{RootOptions, newton_bracketed},
    units::{
        AvailableEnergy, Pressure, SpecificHeatCapacity, ThermodynamicTemperature,
        available_energy::joule_per_kilogram, pressure::pascal,
        specific_heat_capacity::joule_per_kilogram_kelvin, thermodynamic_temperature::kelvin,
    },
};

/// Temperature of water or steam at pressure `p` with specific enthalpy `h`.
///
/// Returns the saturation temperature inside the two-phase dome. Fails outside
/// IF97's range, and for now in region 3.
pub fn temperature_from_enthalpy(
    p: Pressure,
    h: AvailableEnergy,
) -> anyhow::Result<ThermodynamicTemperature> {
    let (p, h) = (p.get::<pascal>(), h.get::<joule_per_kilogram>());
    let (t, _) = temperature_from_enthalpy_si(p, h)
        .with_context(|| format!("steam temperature at {p} Pa, h = {h} J/kg"))?;
    Ok(ThermodynamicTemperature::new::<kelvin>(t))
}

/// Temperature of water or steam at pressure `p` with specific entropy `s`.
///
/// Returns the saturation temperature inside the two-phase dome. Fails outside
/// IF97's range, and for now in region 3.
pub fn temperature_from_entropy(
    p: Pressure,
    s: SpecificHeatCapacity,
) -> anyhow::Result<ThermodynamicTemperature> {
    let (p, s) = (p.get::<pascal>(), s.get::<joule_per_kilogram_kelvin>());
    let (t, _) = temperature_from_entropy_si(p, s)
        .with_context(|| format!("steam temperature at {p} Pa, s = {s} J/(kg·K)"))?;
    Ok(ThermodynamicTemperature::new::<kelvin>(t))
}

/// Temperature, K, and region at p Pa and h J/kg.
pub(crate) fn temperature_from_enthalpy_si(p: f64, h: f64) -> anyhow::Result<(f64, Region)> {
    temperature_si(p, h, Given::Enthalpy)
}

/// Temperature, K, and region at p Pa and s J/(kg·K).
pub(crate) fn temperature_from_entropy_si(p: f64, s: f64) -> anyhow::Result<(f64, Region)> {
    temperature_si(p, s, Given::Entropy)
}

/// The property that is given alongside the pressure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Given {
    Enthalpy,
    Entropy,
}

impl Given {
    fn name(self) -> &'static str {
        match self {
            Self::Enthalpy => "enthalpy",
            Self::Entropy => "entropy",
        }
    }

    fn of(self, x: &RawProperties) -> f64 {
        match self {
            Self::Enthalpy => x.h,
            Self::Entropy => x.s,
        }
    }

    /// ∂(given)/∂T at constant p: cp, or cp / T.
    fn slope(self, x: &RawProperties, t: f64) -> f64 {
        match self {
            Self::Enthalpy => x.cp,
            Self::Entropy => x.cp / t,
        }
    }
}

type Kernel = fn(f64, f64) -> RawProperties;

fn temperature_si(p: f64, x: f64, given: Given) -> anyhow::Result<(f64, Region)> {
    if !(p > 0.0 && p <= P_MAX) {
        return Err(PhysicsError::PressureOutOfRange {
            pressure: p,
            min: 0.0,
            max: P_MAX,
        }
        .into());
    }
    if !x.is_finite() {
        return Err(PhysicsError::InvalidState("non-finite enthalpy or entropy").into());
    }
    let value = |kernel: Kernel, t: f64| given.of(&kernel(p, t));
    let out_of_range = || {
        Err(PhysicsError::InvalidState(
            "enthalpy or entropy outside IF97 at this pressure",
        ))
        .with_context(|| format!("{} {x} at {p} Pa", given.name()))
    };

    let lowest: Kernel = if region_si(p, T_MIN) == Region::CompressedLiquid {
        gibbs::region1
    } else {
        gibbs::region2
    };
    if x < value(lowest, T_MIN) {
        return out_of_range();
    }

    // Below 16.53 MPa liquid meets vapour on the saturation line; above it,
    // region 3 separates them.
    let vapour_from = if p <= saturation::saturation_pressure_si(T_REGION3_MIN) {
        // The same test region selection uses at T_MIN, so both agree on
        // whether liquid exists at this pressure.
        if p < saturation::saturation_pressure_si(T_MIN) {
            T_MIN
        } else {
            let t_sat = saturation::saturation_temperature_si(p);
            if x <= value(gibbs::region1, t_sat) {
                return solve(p, x, given, Region::CompressedLiquid, T_MIN, t_sat);
            }
            if x < value(gibbs::region2, t_sat) {
                return Ok((t_sat, Region::Saturation));
            }
            t_sat
        }
    } else {
        if x <= value(gibbs::region1, T_REGION3_MIN) {
            return solve(p, x, given, Region::CompressedLiquid, T_MIN, T_REGION3_MIN);
        }
        let t_b23 = saturation::b23_temperature_si(p);
        if x < value(gibbs::region2, t_b23) {
            return Err(PhysicsError::InvalidState(
                "IF97 region 3 (near-critical) is not supported yet",
            ))
            .with_context(|| format!("{} {x} at {p} Pa", given.name()));
        }
        t_b23
    };

    if x <= value(gibbs::region2, T_MAX_VAPOUR) {
        return solve(p, x, given, Region::Vapour, vapour_from, T_MAX_VAPOUR);
    }
    if p > P_MAX_HIGH_TEMPERATURE || x > value(gibbs::region5, T_MAX) {
        return out_of_range();
    }
    // Regions 2 and 5 differ by ~1e-5 at 1073.15 K; a value in that gap is the
    // boundary temperature.
    if x <= value(gibbs::region5, T_MAX_VAPOUR) {
        return Ok((T_MAX_VAPOUR, Region::Vapour));
    }
    solve(p, x, given, Region::HighTemperature, T_MAX_VAPOUR, T_MAX)
}

/// Solves given(p, T) = x for T in [lo, hi] on the forward equations of
/// `region`, starting from its backward equation where IF97 has one.
fn solve(
    p: f64,
    x: f64,
    given: Given,
    region: Region,
    lo: f64,
    hi: f64,
) -> anyhow::Result<(f64, Region)> {
    let (kernel, guess): (Kernel, f64) = match (region, given) {
        (Region::CompressedLiquid, Given::Enthalpy) => {
            (gibbs::region1, backward::region1_t_ph(p, x))
        }
        (Region::CompressedLiquid, Given::Entropy) => {
            (gibbs::region1, backward::region1_t_ps(p, x))
        }
        (Region::Vapour, Given::Enthalpy) => (gibbs::region2, backward::region2_t_ph(p, x)),
        (Region::Vapour, Given::Entropy) => (gibbs::region2, backward::region2_t_ps(p, x)),
        _ => {
            // Region 5 has no backward equation: one Newton step from the
            // lower end.
            let start = gibbs::region5(p, lo);
            (
                gibbs::region5,
                lo + (x - given.of(&start)) / given.slope(&start, lo),
            )
        }
    };
    // A guess outside [lo, hi] (or NaN from a backward equation evaluated at
    // the edge of its range) is clamped or replaced by the midpoint.
    let guess = if guess.is_finite() {
        guess
    } else {
        f64::midpoint(lo, hi)
    };
    let t = newton_bracketed(
        |t| {
            let state = kernel(p, t);
            (given.of(&state) - x, given.slope(&state, t))
        },
        lo,
        hi,
        guess,
        RootOptions::new(1e-10, 60)?,
    )
    .with_context(|| format!("{region:?} temperature from {} {x} at {p} Pa", given.name()))?;
    Ok((t, region))
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::steam::{
        properties::properties_si,
        region::{region_si, validate_pt},
    };

    /// Forward then inverse at (p, T): returns (T back from h, T back from s,
    /// region).
    fn round_trip(p: f64, t: f64) -> anyhow::Result<(f64, f64, Region)> {
        validate_pt(p, t)?;
        let region = region_si(p, t);
        let x = properties_si(region, p, t)?;
        let (t_h, r_h) = temperature_from_enthalpy_si(p, x.h)?;
        let (t_s, r_s) = temperature_from_entropy_si(p, x.s)?;
        anyhow::ensure!(
            r_h == region && r_s == region,
            "{region:?} → {r_h:?}, {r_s:?}"
        );
        Ok((t_h, t_s, region))
    }

    #[test]
    fn test_inverse_recovers_the_iapws_table_states() -> anyhow::Result<()> {
        for (p_mpa, t) in [
            (3.0, 300.0),
            (80.0, 300.0),
            (3.0, 500.0),
            (0.0035, 300.0),
            (0.0035, 700.0),
            (30.0, 700.0),
            (0.5, 1500.0),
            (30.0, 2000.0),
        ] {
            let (t_h, t_s, region) = round_trip(p_mpa * 1e6, t)?;
            for back in [t_h, t_s] {
                anyhow::ensure!(
                    (back - t).abs() < 1e-9 * t,
                    "{region:?} {p_mpa} MPa, {t} K → {back} K"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn test_inverse_returns_saturation_temperature_inside_the_dome() -> anyhow::Result<()> {
        let p = 1e6;
        let t_sat = saturation::saturation_temperature_si(p);
        let liquid = gibbs::region1(p, t_sat);
        let vapour = gibbs::region2(p, t_sat);
        let (t, region) = temperature_from_enthalpy_si(p, f64::midpoint(liquid.h, vapour.h))?;
        anyhow::ensure!(
            region == Region::Saturation && t.to_bits() == t_sat.to_bits(),
            "{region:?}, {t} K"
        );
        let (t, region) = temperature_from_entropy_si(p, f64::midpoint(liquid.s, vapour.s))?;
        anyhow::ensure!(
            region == Region::Saturation && t.to_bits() == t_sat.to_bits(),
            "{region:?}, {t} K"
        );
        Ok(())
    }

    #[test]
    fn test_inverse_reports_region_3_and_out_of_range_inputs() -> anyhow::Result<()> {
        // h at 25 MPa between region 1 at 623.15 K and region 2 on B23.
        let error = temperature_from_enthalpy_si(25e6, 2.0e6)
            .err()
            .context("region 3 accepted")?;
        anyhow::ensure!(format!("{error:#}").contains("region 3"), "{error:#}");
        for (p, h) in [
            (1e6, -1e5),
            (1e6, 1e8),
            (2e8, 1e6),
            (0.0, 1e6),
            (1e6, f64::NAN),
        ] {
            anyhow::ensure!(
                temperature_from_enthalpy_si(p, h).is_err(),
                "{p} Pa, h = {h} accepted"
            );
        }
        for (p, s) in [(1e6, -100.0), (1e6, 2e4), (2e8, 1e3)] {
            anyhow::ensure!(
                temperature_from_entropy_si(p, s).is_err(),
                "{p} Pa, s = {s} accepted"
            );
        }
        Ok(())
    }

    #[test]
    fn test_temperature_from_enthalpy_and_entropy_use_si_units() -> anyhow::Result<()> {
        let p = Pressure::new::<pascal>(3e6);
        let x = gibbs::region1(3e6, 300.0);
        let t_h = temperature_from_enthalpy(p, AvailableEnergy::new::<joule_per_kilogram>(x.h))?;
        let t_s = temperature_from_entropy(
            p,
            SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(x.s),
        )?;
        for t in [t_h, t_s] {
            let t = t.get::<kelvin>();
            anyhow::ensure!((t - 300.0).abs() < 1e-9 * 300.0, "{t} K");
        }
        Ok(())
    }

    fn single_phase_state() -> impl Strategy<Value = (f64, f64)> {
        // log-uniform pressure 700 Pa – 100 MPa, uniform temperature over IF97.
        (2.845f64..8.0, 273.15f64..2273.15).prop_map(|(log_p, t)| (10f64.powf(log_p), t))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_inverse_round_trips_every_supported_region((p, t) in single_phase_state()) {
            prop_assume!(validate_pt(p, t).is_ok() && region_si(p, t) != Region::NearCritical);
            let (t_h, t_s, region) = round_trip(p, t).map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            prop_assert!((t_h - t).abs() < 1e-8 * t, "{:?} {} Pa, {} K: T(p, h) = {}", region, p, t, t_h);
            prop_assert!((t_s - t).abs() < 1e-8 * t, "{:?} {} Pa, {} K: T(p, s) = {}", region, p, t, t_s);
        }
    }
}
