//! Region 4 (the saturation line) and the B23 and B2bc boundary curves.
//!
//! The `*_si` kernels are infallible on raw SI values; the unit-typed functions
//! validate the range.

use anyhow::Context;

use super::{
    P_CRIT, P_SAT_MIN, T_CRIT, T_MIN,
    coefficients::saturation::{B2BC, B23, SATURATION},
};
use crate::{
    PhysicsError,
    units::{
        Pressure, ThermodynamicTemperature, pressure::pascal, thermodynamic_temperature::kelvin,
    },
};

/// Saturation pressure of water at temperature `t` (273.15 K up to the critical
/// point).
pub fn saturation_pressure(t: ThermodynamicTemperature) -> anyhow::Result<Pressure> {
    let t = t.get::<kelvin>();
    if !(T_MIN..=T_CRIT).contains(&t) {
        return Err(PhysicsError::TemperatureOutOfRange {
            temperature: t,
            min: T_MIN,
            max: T_CRIT,
        })
        .context("saturation pressure");
    }
    Ok(Pressure::new::<pascal>(saturation_pressure_si(t)))
}

/// Saturation temperature of water at pressure `p` (611.213 Pa up to the
/// critical point).
pub fn saturation_temperature(p: Pressure) -> anyhow::Result<ThermodynamicTemperature> {
    let p = p.get::<pascal>();
    if !(P_SAT_MIN..=P_CRIT).contains(&p) {
        return Err(PhysicsError::PressureOutOfRange {
            pressure: p,
            min: P_SAT_MIN,
            max: P_CRIT,
        })
        .context("saturation temperature");
    }
    Ok(ThermodynamicTemperature::new::<kelvin>(
        saturation_temperature_si(p),
    ))
}

/// Saturation pressure, Pa, at `t` K (273.15 K up to the critical temperature).
#[expect(
    clippy::many_single_char_names,
    reason = "symbols of IAPWS R7-97 Eq. 30"
)]
pub(crate) fn saturation_pressure_si(t: f64) -> f64 {
    let n = &SATURATION;
    // IAPWS R7-97(2012) Eqs. 29a and 30.
    let theta = t + n[8] / (t - n[9]);
    let a = (theta + n[0]).mul_add(theta, n[1]);
    let b = n[2].mul_add(theta, n[3]).mul_add(theta, n[4]);
    let c = n[5].mul_add(theta, n[6]).mul_add(theta, n[7]);
    (2.0 * c / (-b + b.mul_add(b, -4.0 * a * c).sqrt())).powi(4) * 1e6
}

/// Saturation temperature, K, at `p` Pa (611.213 Pa up to the critical
/// pressure).
#[expect(
    clippy::many_single_char_names,
    reason = "symbols of IAPWS R7-97 Eq. 31"
)]
pub(crate) fn saturation_temperature_si(p: f64) -> f64 {
    let n = &SATURATION;
    // IAPWS R7-97(2012) Eqs. 29b and 31.
    let beta = (p * 1e-6).sqrt().sqrt();
    let e = (beta + n[2]).mul_add(beta, n[5]);
    let f = n[0].mul_add(beta, n[3]).mul_add(beta, n[6]);
    let g = n[1].mul_add(beta, n[4]).mul_add(beta, n[7]);
    let d = 2.0 * g / (-f - f.mul_add(f, -4.0 * e * g).sqrt());
    let n10_d = n[9] + d;
    0.5 * (n10_d - n10_d.mul_add(n10_d, -4.0 * n[9].mul_add(d, n[8])).sqrt())
}

/// Pressure, Pa, on the region 2–3 boundary at `t` K (623.15–863.15 K).
pub(crate) fn b23_pressure_si(t: f64) -> f64 {
    let n = &B23;
    n[2].mul_add(t, n[1]).mul_add(t, n[0]) * 1e6
}

/// Temperature, K, on the region 2–3 boundary at `p` Pa (16.53–100 MPa).
pub(crate) fn b23_temperature_si(p: f64) -> f64 {
    let n = &B23;
    n[3] + ((p * 1e-6 - n[4]) / n[2]).sqrt()
}

/// Enthalpy, J/kg, on the boundary between backward subregions 2b and 2c at
/// `p` Pa (6.547–100 MPa).
pub(crate) fn b2bc_enthalpy_si(p: f64) -> f64 {
    let n = &B2BC;
    (n[3] + ((p * 1e-6 - n[4]) / n[2]).sqrt()) * 1e3
}

#[cfg(test)]
#[expect(
    clippy::unreadable_literal,
    reason = "verification values exactly as IAPWS R7-97(2012) prints them"
)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    fn ensure_close(got: f64, want: f64, tolerance: f64, what: &str) -> anyhow::Result<()> {
        let error = (got / want - 1.0).abs();
        anyhow::ensure!(error < tolerance, "{what}: {got} vs {want} (rel {error:e})");
        Ok(())
    }

    #[test]
    fn test_saturation_pressure_matches_iapws_table_35() -> anyhow::Result<()> {
        for (t, p_mpa) in [
            (300.0, 0.353658941e-2),
            (500.0, 0.263889776e1),
            (600.0, 0.123443146e2),
        ] {
            ensure_close(
                saturation_pressure_si(t),
                p_mpa * 1e6,
                5e-9,
                &format!("p_sat({t} K)"),
            )?;
        }
        Ok(())
    }

    #[test]
    fn test_saturation_temperature_matches_iapws_table_36() -> anyhow::Result<()> {
        for (p_mpa, t) in [
            (0.1, 0.372755919e3),
            (1.0, 0.453035632e3),
            (10.0, 0.584149488e3),
        ] {
            ensure_close(
                saturation_temperature_si(p_mpa * 1e6),
                t,
                5e-9,
                &format!("T_sat({p_mpa} MPa)"),
            )?;
        }
        Ok(())
    }

    #[test]
    fn test_saturation_line_ends_at_the_critical_point() -> anyhow::Result<()> {
        ensure_close(saturation_pressure_si(T_CRIT), P_CRIT, 1e-8, "p_sat(T_c)")?;
        ensure_close(
            saturation_temperature_si(P_CRIT),
            T_CRIT,
            1e-8,
            "T_sat(p_c)",
        )
    }

    #[test]
    fn test_b23_matches_the_iapws_check_point() -> anyhow::Result<()> {
        // IAPWS R7-97(2012), after Eq. 6: T = 623.15 K ↔ p = 16.5291643 MPa.
        ensure_close(b23_pressure_si(623.15), 16.529_164_3e6, 5e-9, "p_B23")?;
        ensure_close(b23_temperature_si(16.529_164_3e6), 623.15, 5e-9, "T_B23")
    }

    #[test]
    fn test_b2bc_matches_the_iapws_check_point() -> anyhow::Result<()> {
        // IAPWS R7-97(2012), after Eq. 21: p = 100 MPa ↔ h = 3516.004323 kJ/kg.
        ensure_close(b2bc_enthalpy_si(100e6), 3_516_004.323, 5e-10, "h_2bc")
    }

    #[test]
    fn test_saturation_pressure_and_temperature_use_si_units() -> anyhow::Result<()> {
        let p = saturation_pressure(ThermodynamicTemperature::new::<kelvin>(500.0))?;
        anyhow::ensure!(
            p.get::<pascal>().to_bits() == saturation_pressure_si(500.0).to_bits(),
            "{p:?}"
        );
        let t = saturation_temperature(Pressure::new::<pascal>(1e6))?;
        anyhow::ensure!(
            t.get::<kelvin>().to_bits() == saturation_temperature_si(1e6).to_bits(),
            "{t:?}"
        );
        Ok(())
    }

    #[test]
    fn test_saturation_rejects_states_off_the_saturation_line() -> anyhow::Result<()> {
        for t in [273.0, 650.0, f64::NAN] {
            let result = saturation_pressure(ThermodynamicTemperature::new::<kelvin>(t));
            anyhow::ensure!(result.is_err(), "T = {t} K accepted");
        }
        for p in [600.0, 23e6, f64::NAN] {
            let result = saturation_temperature(Pressure::new::<pascal>(p));
            anyhow::ensure!(result.is_err(), "p = {p} Pa accepted");
        }
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_saturation_temperature_inverts_saturation_pressure(t in T_MIN..T_CRIT) {
            let back = saturation_temperature_si(saturation_pressure_si(t));
            prop_assert!((back - t).abs() < 1e-9 * t, "T {} → {}", t, back);
        }

        #[test]
        fn test_saturation_pressure_rises_with_temperature(a in T_MIN..T_CRIT, b in T_MIN..T_CRIT) {
            let (lo, hi) = (a.min(b), a.max(b));
            prop_assume!(hi - lo > 1e-6);
            prop_assert!(saturation_pressure_si(lo) < saturation_pressure_si(hi));
        }

        #[test]
        fn test_b23_temperature_inverts_b23_pressure(t in 623.15f64..863.15) {
            let back = b23_temperature_si(b23_pressure_si(t));
            prop_assert!((back - t).abs() < 1e-9 * t, "T {} → {}", t, back);
        }
    }
}
