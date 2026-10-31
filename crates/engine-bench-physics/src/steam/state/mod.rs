//! The equilibrium state from specific internal energy and volume, which is
//! what a closed volume conserves: single-phase in any region, or liquid and
//! vapour in equilibrium inside the dome.

use anyhow::Context;

mod single_phase;
mod two_phase;

use two_phase::Dome;

use super::{RHO_CRIT, region::Region};
use crate::{
    PhysicsError,
    units::{
        AvailableEnergy, Pressure, SpecificHeatCapacity, SpecificVolume, ThermodynamicTemperature,
        available_energy::joule_per_kilogram, pressure::pascal,
        specific_heat_capacity::joule_per_kilogram_kelvin,
        specific_volume::cubic_meter_per_kilogram, thermodynamic_temperature::kelvin,
    },
};

/// Equilibrium state of water at given specific internal energy and volume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SteamState {
    /// IF97 region; [`Region::Saturation`] inside the two-phase dome.
    pub region: Region,
    /// Pressure.
    pub pressure: Pressure,
    /// Temperature.
    pub temperature: ThermodynamicTemperature,
    /// Vapour mass fraction inside the dome, `None` for a single phase.
    pub vapour_quality: Option<f64>,
    /// Specific volume.
    pub specific_volume: SpecificVolume,
    /// Specific internal energy.
    pub internal_energy: AvailableEnergy,
    /// Specific enthalpy.
    pub enthalpy: AvailableEnergy,
    /// Specific entropy.
    pub entropy: SpecificHeatCapacity,
}

/// Raw SI state: p Pa, T K, quality (NaN for a single phase), v, u, h, s.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct RawState {
    pub(crate) region: Region,
    pub(crate) p: f64,
    pub(crate) t: f64,
    pub(crate) quality: f64,
    pub(crate) v: f64,
    pub(crate) u: f64,
    pub(crate) h: f64,
    pub(crate) s: f64,
    /// The coexisting phases inside the dome; NaN for a single phase.
    pub(crate) phases: Phases,
}

/// Specific volume, m³/kg, and enthalpy, J/kg, of saturated liquid and
/// vapour.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Phases {
    pub(crate) v_liquid: f64,
    pub(crate) v_vapour: f64,
    pub(crate) h_liquid: f64,
    pub(crate) h_vapour: f64,
}

impl Phases {
    /// A single phase has none.
    pub(crate) const NONE: Self = Self {
        v_liquid: f64::NAN,
        v_vapour: f64::NAN,
        h_liquid: f64::NAN,
        h_vapour: f64::NAN,
    };
}

impl RawState {
    pub(crate) fn to_state(self) -> SteamState {
        let energy = AvailableEnergy::new::<joule_per_kilogram>;
        SteamState {
            region: self.region,
            pressure: Pressure::new::<pascal>(self.p),
            temperature: ThermodynamicTemperature::new::<kelvin>(self.t),
            vapour_quality: (self.region == Region::Saturation).then_some(self.quality),
            specific_volume: SpecificVolume::new::<cubic_meter_per_kilogram>(self.v),
            internal_energy: energy(self.u),
            enthalpy: energy(self.h),
            entropy: SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(self.s),
        }
    }
}

/// Equilibrium state of water with specific internal energy `u` and specific
/// volume `v`.
pub fn state_from_uv(u: AvailableEnergy, v: SpecificVolume) -> anyhow::Result<SteamState> {
    let (u, v) = (
        u.get::<joule_per_kilogram>(),
        v.get::<cubic_meter_per_kilogram>(),
    );
    let state = state_from_uv_si(u, v, None)
        .with_context(|| format!("steam state at u = {u} J/kg, v = {v} m³/kg"))?;
    Ok(state.to_state())
}

/// Equilibrium state at u J/kg and v m³/kg. `hint` is a nearby earlier state
/// (the previous step of a simulation); trying its region first from its p and
/// T makes the usual case a few Newton iterations.
pub(crate) fn state_from_uv_si(
    u: f64,
    v: f64,
    hint: Option<&RawState>,
) -> anyhow::Result<RawState> {
    if !(u.is_finite() && v.is_finite() && v > 0.0) {
        return Err(PhysicsError::InvalidState(
            "internal energy and volume must be finite, volume positive",
        )
        .into());
    }
    if let Some(hint) = hint {
        let near = match hint.region {
            Region::Saturation => two_phase::from_uv_near(u, v, hint.t),
            region => single_phase::solve(region, u, v, Some((hint.p, hint.t))),
        };
        if let Some(state) = near {
            return Ok(state);
        }
    }
    let edge = match two_phase::from_uv(u, v)? {
        Dome::Inside(state) => return Ok(state),
        Dome::Outside(edge) => edge,
    };
    // A stable single-phase solution is the equilibrium state, so the first
    // region that has one wins; liquid-like volumes try the liquid first.
    let order = if v < 1.0 / RHO_CRIT {
        [
            Region::CompressedLiquid,
            Region::NearCritical,
            Region::Vapour,
            Region::HighTemperature,
        ]
    } else {
        [
            Region::Vapour,
            Region::NearCritical,
            Region::HighTemperature,
            Region::CompressedLiquid,
        ]
    };
    order
        .into_iter()
        .find_map(|region| {
            // The dome's edge is a good start only for the regions that border
            // it.
            let start =
                edge.filter(|_| matches!(region, Region::CompressedLiquid | Region::Vapour));
            single_phase::solve(region, u, v, start)
        })
        .ok_or_else(|| {
            PhysicsError::InvalidState("no IF97 state has this internal energy and volume").into()
        })
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::steam::{
        T_CRIT, T_MAX_VAPOUR, T_MIN, T_REGION3_MIN,
        properties::properties_si,
        region::{region_si, validate_pt},
        region3::saturated_si,
        saturation::b23_temperature_si,
    };

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
    }

    /// Single-phase (u, v) → state, checked against the (p, T) it came from.
    fn round_trip(p: f64, t: f64, hint: Option<&RawState>) -> anyhow::Result<RawState> {
        let region = region_si(p, t);
        let x = properties_si(region, p, t)?;
        let state = state_from_uv_si(x.u, x.v, hint)?;
        anyhow::ensure!(
            state.region == region,
            "{p} Pa, {t} K: {:?} instead of {region:?}",
            state.region
        );
        anyhow::ensure!(
            (state.t / t - 1.0).abs() < 1e-9,
            "{p} Pa, {t} K: T = {}",
            state.t
        );
        anyhow::ensure!(
            (state.p / p - 1.0).abs() < 1e-7,
            "{p} Pa, {t} K: p = {}",
            state.p
        );
        anyhow::ensure!(
            (state.h / x.h - 1.0).abs() < 1e-9,
            "{p} Pa, {t} K: h = {}",
            state.h
        );
        Ok(state)
    }

    #[test]
    fn test_state_from_uv_recovers_a_state_in_every_single_phase_region() -> anyhow::Result<()> {
        for (p, t) in [
            (3e6, 300.0),
            (80e6, 500.0),
            (3500.0, 300.0),
            (1e6, 700.0),
            (25.583_701_8e6, 650.0),
            (50e6, 750.0),
            (30e6, 2000.0),
        ] {
            round_trip(p, t, None)?;
        }
        Ok(())
    }

    #[test]
    #[expect(
        clippy::many_single_char_names,
        reason = "thermodynamic symbols p, T, u, v, h, s, x"
    )]
    fn test_state_from_uv_returns_a_two_phase_mixture_inside_the_dome() -> anyhow::Result<()> {
        for t in [300.0, 373.15, 500.0, 620.0, 640.0, 646.0] {
            let (p, liquid, vapour) = saturated_si(t)?;
            let x = 0.3;
            let v = liquid.v + x * (vapour.v - liquid.v);
            let u = liquid.u + x * (vapour.u - liquid.u);
            let state = state_from_uv_si(u, v, None)?;
            anyhow::ensure!(
                state.region == Region::Saturation,
                "{t} K: {:?}",
                state.region
            );
            anyhow::ensure!((state.t - t).abs() < 1e-8 * t, "{t} K: T = {}", state.t);
            anyhow::ensure!((state.p / p - 1.0).abs() < 1e-7, "{t} K: p = {}", state.p);
            anyhow::ensure!(
                (state.quality - x).abs() < 1e-8,
                "{t} K: x = {}",
                state.quality
            );
            let h = liquid.h + x * (vapour.h - liquid.h);
            anyhow::ensure!(
                (state.h / h - 1.0).abs() < 1e-8,
                "{t} K: h = {} vs {h}",
                state.h
            );
        }
        Ok(())
    }

    #[test]
    fn test_state_from_uv_stays_within_if97_consistency_at_623_k() -> anyhow::Result<()> {
        for t in [
            T_REGION3_MIN - 0.005,
            T_REGION3_MIN,
            T_REGION3_MIN + 1e-4,
            T_REGION3_MIN + 0.005,
        ] {
            for x in [0.1, 0.6, 0.95] {
                let (_, liquid, vapour) = saturated_si(t)?;
                let v = liquid.v + x * (vapour.v - liquid.v);
                let u = liquid.u + x * (vapour.u - liquid.u);
                let state = state_from_uv_si(u, v, None)?;
                anyhow::ensure!(
                    (state.t - t).abs() < 0.01,
                    "{t} K, x = {x}: T = {}",
                    state.t
                );
                anyhow::ensure!(
                    (state.quality - x).abs() < 1e-3,
                    "{t} K, x = {x}: x = {}",
                    state.quality
                );
            }
        }
        Ok(())
    }

    #[test]
    fn test_state_from_uv_wraps_the_state_in_si_units() -> anyhow::Result<()> {
        let x = properties_si(Region::Vapour, 1e5, 400.0)?;
        let state = state_from_uv(
            AvailableEnergy::new::<joule_per_kilogram>(x.u),
            SpecificVolume::new::<cubic_meter_per_kilogram>(x.v),
        )?;
        anyhow::ensure!(state.vapour_quality.is_none(), "{state:?}");
        let t = state.temperature.get::<kelvin>();
        anyhow::ensure!((t - 400.0).abs() < 1e-9 * 400.0, "{t} K");
        Ok(())
    }

    #[test]
    fn test_state_from_uv_rejects_states_outside_if97() -> anyhow::Result<()> {
        let (_, liquid, _) = saturated_si(T_MIN)?;
        for (u, v) in [
            (1e6, 0.0),
            (1e6, -1.0),
            (f64::NAN, 1.0),
            (1e6, f64::INFINITY),
            // Colder than the triple point: ice.
            (liquid.u - 5e4, 0.1),
            // Hotter than 2273.15 K.
            (2e7, 1.0),
        ] {
            anyhow::ensure!(
                state_from_uv_si(u, v, None).is_err(),
                "u = {u}, v = {v} accepted"
            );
        }
        Ok(())
    }

    /// Near a region boundary IF97's neighbouring equations disagree slightly,
    /// so a (u, v) pair can be valid in both; those bands are left out.
    fn away_from_boundaries(p: f64, t: f64) -> bool {
        let band = 0.2;
        let near = |boundary: f64| (t - boundary).abs() < band;
        !(near(T_MAX_VAPOUR)
            || (p > 16.5e6
                && (near(T_REGION3_MIN)
                    || (p >= 16.5e6 && t < 870.0 && near(b23_temperature_si(p))))))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_state_from_uv_round_trips_single_phase_states(
            log_p in 2.845f64..8.0,
            t in T_MIN..2273.15,
        ) {
            let p = 10f64.powf(log_p);
            prop_assume!(validate_pt(p, t).is_ok() && away_from_boundaries(p, t));
            // Exclude a 1e-4 band around saturation, where the liquid and
            // vapour sides meet the dome.
            prop_assume!(t >= T_CRIT || p < 611.0 || {
                let p_sat = crate::steam::saturation::saturation_pressure_si(t.min(T_CRIT));
                (p / p_sat - 1.0).abs() > 1e-4
            });
            round_trip(p, t, None).map_err(|e| fail(&e))?;
        }

        #[test]
        fn test_state_from_uv_round_trips_two_phase_mixtures(t in 274.0f64..646.9, x in 0.001f64..0.999) {
            // At 623.15 K the saturated phases switch from regions 1 and 2 to
            // region 3, whose equations differ slightly; checked separately.
            prop_assume!((t - T_REGION3_MIN).abs() > 0.01);
            let (p, liquid, vapour) = saturated_si(t).map_err(|e| fail(&e))?;
            let v = liquid.v + x * (vapour.v - liquid.v);
            let u = liquid.u + x * (vapour.u - liquid.u);
            let state = state_from_uv_si(u, v, None).map_err(|e| fail(&e))?;
            prop_assert!(state.region == Region::Saturation, "{} K, x = {}: {:?}", t, x, state.region);
            prop_assert!((state.t - t).abs() < 1e-7 * t, "{} K, x = {}: T = {}", t, x, state.t);
            prop_assert!((state.quality - x).abs() < 1e-6, "{} K, x = {}: x = {}", t, x, state.quality);
            prop_assert!((state.p / p - 1.0).abs() < 1e-6);
        }

        #[test]
        fn test_state_from_uv_with_a_nearby_hint_agrees_with_the_cold_solve(
            t in 280.0f64..640.0,
            x in 0.0f64..1.2,
            dt in -2.0f64..2.0,
        ) {
            // x > 1 continues as superheated vapour at the same pressure.
            let (p, liquid, vapour) = saturated_si(t).map_err(|e| fail(&e))?;
            let (u, v) = if x <= 1.0 {
                (liquid.u + x * (vapour.u - liquid.u), liquid.v + x * (vapour.v - liquid.v))
            } else {
                let y = properties_si(Region::Vapour, p, t + 50.0 * (x - 1.0)).map_err(|e| fail(&e))?;
                (y.u, y.v)
            };
            let cold = state_from_uv_si(u, v, None).map_err(|e| fail(&e))?;
            let hint_state = state_from_uv_si(u * (1.0 + 1e-4 * dt), v * (1.0 - 1e-4 * dt), None)
                .map_err(|e| fail(&e))?;
            let warm = state_from_uv_si(u, v, Some(&hint_state)).map_err(|e| fail(&e))?;
            prop_assert!(warm.region == cold.region, "{:?} vs {:?}", warm.region, cold.region);
            prop_assert!((warm.t - cold.t).abs() < 1e-8 * cold.t, "T {} vs {}", warm.t, cold.t);
            prop_assert!((warm.p / cold.p - 1.0).abs() < 1e-7, "p {} vs {}", warm.p, cold.p);
        }
    }
}
