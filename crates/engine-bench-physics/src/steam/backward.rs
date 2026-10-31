//! IF97 backward equations T(p, h) and T(p, s) for regions 1 and 2.
//!
//! They agree with the forward equations to a few millikelvin, so
//! [`super::inverse`] uses them only as starting points for a Newton solve on
//! the forward equations. Infallible kernels on raw SI values; callers have
//! picked the region.

use super::{
    coefficients::{
        Term,
        backward::{
            REGION1_PH, REGION1_PS, REGION2A_PH, REGION2A_PS, REGION2B_PH, REGION2B_PS,
            REGION2C_PH, REGION2C_PS,
        },
    },
    saturation::b2bc_enthalpy_si,
};

/// Upper pressure of backward subregion 2a, Pa.
const P_2A_MAX: f64 = 4e6;
/// Pressure below which subregion 2b extends to the saturation line, Pa (the
/// low end of the B2bc curve).
const P_2BC_MIN: f64 = 6.5467e6;
/// Entropy separating subregions 2b and 2c above 4 MPa, J/(kg·K).
const S_2BC: f64 = 5.85e3;

/// Region 1 temperature, K, from p Pa and h J/kg.
pub(crate) fn region1_t_ph(p: f64, h: f64) -> f64 {
    polynomial(&REGION1_PH, p * 1e-6, h / 2500e3 + 1.0)
}

/// Region 1 temperature, K, from p Pa and s J/(kg·K).
pub(crate) fn region1_t_ps(p: f64, s: f64) -> f64 {
    polynomial(&REGION1_PS, p * 1e-6, s * 1e-3 + 2.0)
}

/// Region 2 temperature, K, from p Pa and h J/kg (subregions 2a, 2b and 2c).
pub(crate) fn region2_t_ph(p: f64, h: f64) -> f64 {
    let (pi, eta) = (p * 1e-6, h / 2000e3);
    if p <= P_2A_MAX {
        polynomial(&REGION2A_PH, pi, eta - 2.1)
    } else if p <= P_2BC_MIN || h >= b2bc_enthalpy_si(p) {
        polynomial(&REGION2B_PH, pi - 2.0, eta - 2.6)
    } else {
        polynomial(&REGION2C_PH, pi + 25.0, eta - 1.8)
    }
}

/// Region 2 temperature, K, from p Pa and s J/(kg·K) (subregions 2a, 2b and
/// 2c).
pub(crate) fn region2_t_ps(p: f64, s: f64) -> f64 {
    let pi = p * 1e-6;
    if p <= P_2A_MAX {
        let y = s / 2e3 - 2.0;
        REGION2A_PS
            .iter()
            .map(|term| term.n * pi.powf(term.i) * y.powi(term.j))
            .sum()
    } else if s >= S_2BC {
        polynomial(&REGION2B_PS, pi, 10.0 - s / 785.3)
    } else {
        polynomial(&REGION2C_PS, pi, 2.0 - s / 2925.1)
    }
}

/// Σ n x^I y^J.
fn polynomial(terms: &[Term<i32>], x: f64, y: f64) -> f64 {
    terms
        .iter()
        .map(|term| term.n * x.powi(term.i) * y.powi(term.j))
        .sum()
}

#[cfg(test)]
#[expect(
    clippy::unreadable_literal,
    reason = "verification values exactly as IAPWS R7-97(2012) prints them"
)]
mod tests {
    use super::*;

    /// The tables give 9 significant digits (600.48404 K: 8).
    const TOLERANCE: f64 = 5e-8;

    fn check(
        rows: &[(f64, f64, f64)],
        kernel: fn(f64, f64) -> f64,
        x_unit: f64,
    ) -> anyhow::Result<()> {
        for &(p_mpa, x, t) in rows {
            let got = kernel(p_mpa * 1e6, x * x_unit);
            let error = (got / t - 1.0).abs();
            anyhow::ensure!(
                error < TOLERANCE,
                "{p_mpa} MPa, {x}: {got} K vs {t} K (rel {error:e})"
            );
        }
        Ok(())
    }

    #[test]
    fn test_region1_t_ph_matches_iapws_table_7() -> anyhow::Result<()> {
        let rows = [
            (3.0, 500.0, 391.798509),
            (80.0, 500.0, 378.108626),
            (80.0, 1500.0, 611.041229),
        ];
        check(&rows, region1_t_ph, 1e3)
    }

    #[test]
    fn test_region1_t_ps_matches_iapws_table_9() -> anyhow::Result<()> {
        let rows = [
            (3.0, 0.5, 307.842258),
            (80.0, 0.5, 309.979785),
            (80.0, 3.0, 565.899909),
        ];
        check(&rows, region1_t_ps, 1e3)
    }

    #[test]
    fn test_region2_t_ph_matches_iapws_table_24() -> anyhow::Result<()> {
        let rows = [
            (0.001, 3000.0, 534.433241),
            (3.0, 3000.0, 575.373370),
            (3.0, 4000.0, 1010.77577),
            (5.0, 3500.0, 801.299102),
            (5.0, 4000.0, 1015.31583),
            (25.0, 3500.0, 875.279054),
            (40.0, 2700.0, 743.056411),
            (60.0, 2700.0, 791.137067),
            (60.0, 3200.0, 882.756860),
        ];
        check(&rows, region2_t_ph, 1e3)
    }

    #[test]
    fn test_region2_t_ps_matches_iapws_table_29() -> anyhow::Result<()> {
        let rows = [
            (0.1, 7.5, 399.517097),
            (0.1, 8.0, 514.127081),
            (2.5, 8.0, 1039.84917),
            (8.0, 6.0, 600.48404),
            (8.0, 7.5, 1064.95556),
            (90.0, 6.0, 1038.01126),
            (20.0, 5.75, 697.992849),
            (80.0, 5.25, 854.011484),
            (80.0, 5.75, 949.017998),
        ];
        check(&rows, region2_t_ps, 1e3)
    }
}
