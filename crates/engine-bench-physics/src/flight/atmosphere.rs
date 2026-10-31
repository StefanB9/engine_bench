//! International Standard Atmosphere (ISO 2533) from −500 m to 20 km: the
//! troposphere with a constant lapse rate and the isothermal lower
//! stratosphere. Altitude is geopotential.

use anyhow::Context;

use crate::{
    PhysicsError,
    units::{
        Length, MassDensity, Pressure, ThermodynamicTemperature, Velocity, length::meter,
        mass_density::kilogram_per_cubic_meter, pressure::pascal,
        thermodynamic_temperature::kelvin, velocity::meter_per_second,
    },
};

/// Standard gravity, m/s².
pub const STANDARD_GRAVITY: f64 = 9.806_65;
/// Specific gas constant of ISA air, J/(kg·K).
const R_AIR: f64 = 287.052_87;
/// Ratio of specific heats of ISA air.
const GAMMA_AIR: f64 = 1.4;
/// Sea-level temperature, K, and pressure, Pa.
const T0: f64 = 288.15;
const P0: f64 = 101_325.0;
/// Tropospheric lapse rate, K/m.
const LAPSE: f64 = 0.0065;
/// Tropopause altitude, m.
const TROPOPAUSE: f64 = 11_000.0;
/// Lowest and highest altitudes covered, m.
const MIN_ALTITUDE: f64 = -500.0;
const MAX_ALTITUDE: f64 = 20_000.0;

/// Air properties at one altitude.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtmosphereState {
    /// Static temperature.
    pub temperature: ThermodynamicTemperature,
    /// Static pressure.
    pub pressure: Pressure,
    /// Density.
    pub density: MassDensity,
    /// Speed of sound.
    pub speed_of_sound: Velocity,
}

/// The standard atmosphere at `altitude`, which must lie in −500 m ..= 20 km.
pub fn standard_atmosphere(altitude: Length) -> anyhow::Result<AtmosphereState> {
    let h = altitude.get::<meter>();
    if !(MIN_ALTITUDE..=MAX_ALTITUDE).contains(&h) {
        return Err(PhysicsError::InvalidState(
            "standard atmosphere covers -500 m to 20 km",
        ))
        .with_context(|| format!("altitude {h} m"));
    }
    Ok(standard_atmosphere_clamped(altitude))
}

/// The standard atmosphere at `altitude` clamped into −500 m ..= 20 km. Kernel
/// for ODE derivatives: never fails.
#[inline]
pub fn standard_atmosphere_clamped(altitude: Length) -> AtmosphereState {
    let h = altitude.get::<meter>().clamp(MIN_ALTITUDE, MAX_ALTITUDE);
    // Hydrostatic exponent g₀/(R·L) of the constant-lapse layer.
    let exponent = STANDARD_GRAVITY / (R_AIR * LAPSE);
    let t_tropopause = T0 - LAPSE * TROPOPAUSE;
    let (t, p) = if h <= TROPOPAUSE {
        let t = T0 - LAPSE * h;
        (t, P0 * (t / T0).powf(exponent))
    } else {
        let p_tropopause = P0 * (t_tropopause / T0).powf(exponent);
        let decay = -STANDARD_GRAVITY * (h - TROPOPAUSE) / (R_AIR * t_tropopause);
        (t_tropopause, p_tropopause * decay.exp())
    };
    AtmosphereState {
        temperature: ThermodynamicTemperature::new::<kelvin>(t),
        pressure: Pressure::new::<pascal>(p),
        density: MassDensity::new::<kilogram_per_cubic_meter>(p / (R_AIR * t)),
        speed_of_sound: Velocity::new::<meter_per_second>((GAMMA_AIR * R_AIR * t).sqrt()),
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    fn at(altitude: f64) -> anyhow::Result<(f64, f64, f64, f64)> {
        let s = standard_atmosphere(Length::new::<meter>(altitude))?;
        Ok((
            s.temperature.get::<kelvin>(),
            s.pressure.get::<pascal>(),
            s.density.get::<kilogram_per_cubic_meter>(),
            s.speed_of_sound.get::<meter_per_second>(),
        ))
    }

    fn close(what: &str, got: f64, want: f64, rel: f64) -> anyhow::Result<()> {
        anyhow::ensure!((got / want - 1.0).abs() <= rel, "{what}: {got} vs {want}");
        Ok(())
    }

    /// ISA table values (ISO 2533 / US Standard Atmosphere 1976).
    #[test]
    fn test_standard_atmosphere_matches_tables() -> anyhow::Result<()> {
        for (altitude, t, p, rho, a) in [
            (0.0, 288.15, 101_325.0, 1.2250, 340.29),
            (-500.0, 291.40, 107_477.0, 1.2849, 342.21),
            (1000.0, 281.65, 89_874.6, 1.1117, 336.43),
            (5000.0, 255.65, 54_019.9, 0.73612, 320.53),
            (11_000.0, 216.65, 22_632.1, 0.36392, 295.07),
            (20_000.0, 216.65, 5_474.89, 0.088_035, 295.07),
        ] {
            let (gt, gp, grho, ga) = at(altitude)?;
            close(&format!("T({altitude})"), gt, t, 1e-5)?;
            close(&format!("p({altitude})"), gp, p, 2e-4)?;
            close(&format!("rho({altitude})"), grho, rho, 3e-4)?;
            close(&format!("a({altitude})"), ga, a, 1e-4)?;
        }
        Ok(())
    }

    #[test]
    fn test_standard_atmosphere_rejects_out_of_range_altitudes() -> anyhow::Result<()> {
        for bad in [-501.0, 20_001.0, f64::NAN] {
            anyhow::ensure!(at(bad).is_err(), "accepted {bad} m");
        }
        Ok(())
    }

    #[test]
    fn test_standard_atmosphere_clamped_holds_the_end_values() -> anyhow::Result<()> {
        let top = standard_atmosphere(Length::new::<meter>(MAX_ALTITUDE))?;
        let above = standard_atmosphere_clamped(Length::new::<meter>(30_000.0));
        anyhow::ensure!(above == top, "above range: {above:?}");
        let bottom = standard_atmosphere(Length::new::<meter>(MIN_ALTITUDE))?;
        let below = standard_atmosphere_clamped(Length::new::<meter>(-2000.0));
        anyhow::ensure!(below == bottom, "below range: {below:?}");
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_standard_atmosphere_is_hydrostatic_and_ideal(altitude in -490.0..19_990.0f64) {
            let state = |h: f64| standard_atmosphere_clamped(Length::new::<meter>(h));
            let s = state(altitude);
            let (t, p, rho) = (
                s.temperature.get::<kelvin>(),
                s.pressure.get::<pascal>(),
                s.density.get::<kilogram_per_cubic_meter>(),
            );
            // Ideal gas.
            prop_assert!((p / (rho * R_AIR * t) - 1.0).abs() < 1e-12);
            // Hydrostatic balance dp/dh = −ρ·g.
            let dh = 1.0;
            let dp = state(altitude + dh).pressure.get::<pascal>()
                - state(altitude - dh).pressure.get::<pascal>();
            let gradient = dp / (2.0 * dh);
            prop_assert!((gradient / (-rho * STANDARD_GRAVITY) - 1.0).abs() < 1e-6,
                "dp/dh {gradient} vs {}", -rho * STANDARD_GRAVITY);
        }
    }
}
