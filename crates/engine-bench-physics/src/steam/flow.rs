//! Water and steam (IAPWS-IF97) as a [`RealFluid`]: isentropic states for
//! nozzle and orifice flow, with the homogeneous equilibrium model inside the
//! dome.

use anyhow::Context;

use super::{
    P_MAX_HIGH_TEMPERATURE, T_MAX, T_MAX_VAPOUR, T_MIN,
    inverse::temperature_from_entropy_si,
    properties::properties_si,
    region::{Region, region_si, validate_pt},
    region3::saturated_si,
};
use crate::{
    PhysicsError,
    flow::real::{FlowState, RealFluid},
    numerics::roots::{RootOptions, brent},
    units::{
        MassDensity, Pressure, SpecificHeatCapacity, mass_density::kilogram_per_cubic_meter,
        pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
    },
};

/// Water and steam from IAPWS-IF97.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Water;

impl RealFluid for Water {
    #[expect(
        clippy::many_single_char_names,
        reason = "the thermodynamic symbols p, s, T, h, w"
    )]
    fn state_ps(&self, p: Pressure, s: SpecificHeatCapacity) -> anyhow::Result<FlowState> {
        let (p, s) = (p.get::<pascal>(), s.get::<joule_per_kilogram_kelvin>());
        let (t, region) = temperature_from_entropy_si(p, s).context("isentropic water state")?;
        if region != Region::Saturation {
            let x = properties_si(region, p, t)?;
            return Ok(FlowState::from_si(p, t, 1.0 / x.v, x.h, x.s, x.w, f64::NAN));
        }
        let (rho, h, quality) = mixture(t, s)?;
        // Equilibrium speed of sound √(∂p/∂ρ)_s from the mixture density on
        // the same isentrope a little above and below p.
        let density_at = |p: f64| -> anyhow::Result<f64> {
            let (t, _) = temperature_from_entropy_si(p, s)?;
            Ok(mixture(t, s)?.0)
        };
        let dp = 1e-6 * p;
        let w = (2.0 * dp / (density_at(p + dp)? - density_at(p - dp)?)).sqrt();
        Ok(FlowState::from_si(p, t, rho, h, s, w, quality))
    }

    fn state_prho(&self, p: Pressure, rho: MassDensity) -> anyhow::Result<FlowState> {
        let (p, rho) = (p.get::<pascal>(), rho.get::<kilogram_per_cubic_meter>());
        if !(rho > 0.0 && rho.is_finite()) {
            return Err(PhysicsError::InvalidState("density must be positive"))
                .with_context(|| format!("ρ = {rho} kg/m³"));
        }
        let t_max = if p > P_MAX_HIGH_TEMPERATURE {
            T_MAX_VAPOUR
        } else {
            T_MAX
        };
        let volume = |t: f64| -> anyhow::Result<f64> {
            validate_pt(p, t)?;
            Ok(properties_si(region_si(p, t), p, t)?.v)
        };
        // Along an isobar v rises with T, jumping across saturation, except in
        // liquid colder than its density maximum near 4 °C. A density that the
        // liquid has twice there is given its warmer temperature.
        let lo = if volume(T_MIN)? * rho > 1.0 {
            T_DENSEST
        } else {
            T_MIN
        };
        let mut failure = None;
        let t = brent(
            |t| {
                volume(t).map_or_else(
                    |e| {
                        failure.get_or_insert(e);
                        f64::NAN
                    },
                    |v| v * rho - 1.0,
                )
            },
            lo,
            t_max,
            RootOptions::new(1e-12 * t_max, 200)?,
        );
        if let Some(e) = failure {
            return Err(e).context("water temperature at (p, ρ)");
        }
        let t = t.with_context(|| format!("water at {p} Pa, {rho} kg/m³"))?;
        let region = region_si(p, t);
        let x = properties_si(region, p, t)?;
        if (x.v * rho - 1.0).abs() > 1e-9 {
            return Err(PhysicsError::InvalidState(
                "density only a two-phase mixture has",
            ))
            .with_context(|| format!("water at {p} Pa, {rho} kg/m³"));
        }
        Ok(FlowState::from_si(p, t, rho, x.h, x.s, x.w, f64::NAN))
    }
}

/// Liquid water is densest here, K.
const T_DENSEST: f64 = 277.13;

/// Saturated liquid and vapour at `t` K mixed to entropy `s`: density, kg/m³,
/// enthalpy, J/kg, and vapour quality.
#[expect(
    clippy::many_single_char_names,
    reason = "IF97 symbols: l and v for liquid and vapour"
)]
fn mixture(t: f64, s: f64) -> anyhow::Result<(f64, f64, f64)> {
    let (_, l, v) = saturated_si(t)?;
    let x = (s - l.s) / (v.s - l.s);
    Ok((
        1.0 / x.mul_add(v.v - l.v, l.v),
        x.mul_add(v.h - l.h, l.h),
        x,
    ))
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::{
        steam::{
            properties::properties_si,
            region::{Region, region_si, validate_pt},
            region3::saturated_si,
        },
        units::{
            mass_density::kilogram_per_cubic_meter, pressure::pascal,
            specific_heat_capacity::joule_per_kilogram_kelvin, thermodynamic_temperature::kelvin,
            velocity::meter_per_second,
        },
    };

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
    }

    #[expect(
        clippy::many_single_char_names,
        reason = "IF97 symbols: l and v for liquid and vapour"
    )]
    #[test]
    fn test_water_state_ps_inside_the_dome_is_an_equilibrium_mixture() -> anyhow::Result<()> {
        let (p, l, v) = saturated_si(400.0)?;
        let s = l.s + 0.5 * (v.s - l.s);
        let x = Water.state_ps(
            Pressure::new::<pascal>(p),
            SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(s),
        )?;
        let quality = x
            .vapour_quality
            .ok_or_else(|| anyhow::anyhow!("single phase"))?;
        anyhow::ensure!((quality - 0.5).abs() < 1e-9, "x = {quality}");
        anyhow::ensure!(
            (x.temperature.get::<kelvin>() - 400.0).abs() < 1e-8,
            "{x:?}"
        );
        let rho = 1.0 / (0.5 * l.v + 0.5 * v.v);
        anyhow::ensure!(
            (x.density.get::<kilogram_per_cubic_meter>() / rho - 1.0).abs() < 1e-9,
            "ρ"
        );
        let w = x.speed_of_sound.get::<meter_per_second>();
        anyhow::ensure!(w > 0.0 && w < v.w, "w = {w} vs vapour {}", v.w);
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_water_flow_states_round_trip_single_phase(log_p in 3.0f64..7.6, t in 280.0f64..1500.0) {
            let p = 10f64.powf(log_p);
            prop_assume!(validate_pt(p, t).is_ok());
            // Away from IF97's region seams (623.15 K, B23, 1073.15 K).
            prop_assume!((t - 623.15).abs() > 0.5 && (t - 1073.15).abs() > 0.5);
            let region = region_si(p, t);
            prop_assume!(region != Region::NearCritical || t > 650.0);
            let x = properties_si(region, p, t).map_err(|e| fail(&e))?;
            let pressure = Pressure::new::<pascal>(p);
            let by_s = Water.state_ps(pressure, SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(x.s)).map_err(|e| fail(&e))?;
            prop_assert!((by_s.temperature.get::<kelvin>() / t - 1.0).abs() < 1e-8, "{:?} {} Pa: T {:?} vs {}", region, p, by_s.temperature, t);
            prop_assert!((by_s.speed_of_sound.get::<meter_per_second>() / x.w - 1.0).abs() < 1e-7);
            let by_rho = Water.state_prho(pressure, MassDensity::new::<kilogram_per_cubic_meter>(1.0 / x.v)).map_err(|e| fail(&e))?;
            prop_assert!((by_rho.temperature.get::<kelvin>() / t - 1.0).abs() < 1e-8, "{:?} {} Pa: T {:?} vs {}", region, p, by_rho.temperature, t);
        }
    }
}
