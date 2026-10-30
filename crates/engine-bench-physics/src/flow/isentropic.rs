//! Isentropic flow relations of a calorically perfect gas as functions of Mach
//! number M and γ. Kernels: no validation (γ > 1 and M ≥ 0 are the caller's
//! responsibility).

use anyhow::Context;

use crate::{
    PhysicsError,
    numerics::roots::{RootOptions, newton_bracketed},
    units::{Angle, angle::radian},
};

/// Which solution of A/A* = ε to return.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Branch {
    /// M ≤ 1.
    Subsonic,
    /// M ≥ 1.
    Supersonic,
}

/// 1 + (γ−1)/2·M² = T₀/T.
#[inline]
fn stagnation_factor(gamma: f64, mach: f64) -> f64 {
    1.0 + 0.5 * (gamma - 1.0) * mach * mach
}

/// (γ+1)/(2(γ−1)), the exponent in the area and mass-flux relations.
#[inline]
fn choking_exponent(gamma: f64) -> f64 {
    (gamma + 1.0) / (2.0 * (gamma - 1.0))
}

/// T/T₀.
#[inline]
pub fn temperature_ratio(gamma: f64, mach: f64) -> f64 {
    stagnation_factor(gamma, mach).recip()
}

/// p/p₀.
#[inline]
pub fn pressure_ratio(gamma: f64, mach: f64) -> f64 {
    stagnation_factor(gamma, mach).powf(-gamma / (gamma - 1.0))
}

/// ρ/ρ₀.
#[inline]
pub fn density_ratio(gamma: f64, mach: f64) -> f64 {
    stagnation_factor(gamma, mach).powf(-1.0 / (gamma - 1.0))
}

/// A/A*: area relative to the sonic throat area for the same mass flow.
#[inline]
pub fn area_ratio(gamma: f64, mach: f64) -> f64 {
    let base = 2.0 / (gamma + 1.0) * stagnation_factor(gamma, mach);
    base.powf(choking_exponent(gamma)) / mach
}

/// p*/p₀: the pressure ratio at which flow through a throat chokes.
#[inline]
pub fn critical_pressure_ratio(gamma: f64) -> f64 {
    (2.0 / (gamma + 1.0)).powf(gamma / (gamma - 1.0))
}

/// Mass flux parameter ṁ·√(R·T₀) / (A·p₀) at Mach `mach`; its maximum is at M =
/// 1.
#[inline]
pub fn mass_flux_parameter(gamma: f64, mach: f64) -> f64 {
    gamma.sqrt() * mach * stagnation_factor(gamma, mach).powf(-choking_exponent(gamma))
}

/// Prandtl–Meyer angle ν(M) for M ≥ 1 (zero below).
pub fn prandtl_meyer(gamma: f64, mach: f64) -> Angle {
    let g = ((gamma + 1.0) / (gamma - 1.0)).sqrt();
    let b = (mach * mach - 1.0).max(0.0).sqrt();
    Angle::new::<radian>(g * (b / g).atan() - b.atan())
}

/// Mach number at which A/A* equals `area_ratio` (≥ 1) on the given branch.
///
/// Solves on ln(A/A*), which stays finite where A/A* itself overflows (γ near
/// 1, large M).
pub fn mach_from_area_ratio(gamma: f64, area_ratio: f64, branch: Branch) -> anyhow::Result<f64> {
    if !(area_ratio >= 1.0 && area_ratio.is_finite()) {
        return Err(PhysicsError::InvalidState("area ratio must be at least 1"))
            .with_context(|| format!("A/A* = {area_ratio}"));
    }
    let exponent = choking_exponent(gamma);
    let target = area_ratio.ln();
    let ln_area = |m: f64| {
        let factor = stagnation_factor(gamma, m);
        let value = exponent * (2.0 / (gamma + 1.0) * factor).ln() - m.ln() - target;
        (value, (m * m - 1.0) / (m * factor))
    };
    // Asymptotes bound the root: ln(A/A*) → −ln M + c₀ as M → 0 and grows like
    // (2/(γ−1))·ln M + c∞ as M → ∞; the true curve lies above both.
    let (lo, hi, guess) = match branch {
        Branch::Subsonic => {
            let asymptote = (exponent * (2.0 / (gamma + 1.0)).ln() - target).exp();
            (0.5 * asymptote, 1.0, asymptote.min(1.0))
        }
        Branch::Supersonic => {
            let c_inf = exponent * ((gamma - 1.0) / (gamma + 1.0)).ln();
            let asymptote = ((target - c_inf) * (gamma - 1.0) / 2.0).exp();
            (1.0, 2.0 * asymptote.max(1.0), asymptote.max(1.0))
        }
    };
    let options = RootOptions::new(1e-13 * hi, 200).context("root options")?;
    newton_bracketed(ln_area, lo, hi, guess, options)
        .with_context(|| format!("mach from A/A* = {area_ratio}, gamma = {gamma}, {branch:?}"))
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::units::angle::degree;

    const AIR: f64 = 1.4;

    fn close(got: f64, want: f64, tol: f64) -> anyhow::Result<()> {
        anyhow::ensure!((got - want).abs() <= tol, "got {got}, want {want} ± {tol}");
        Ok(())
    }

    /// Anderson, Modern Compressible Flow, Table A.1 (γ = 1.4).
    #[test]
    fn test_isentropic_ratios_match_tables_at_mach_2() -> anyhow::Result<()> {
        close(temperature_ratio(AIR, 2.0), 0.5556, 1e-4).context("T/T0")?;
        close(pressure_ratio(AIR, 2.0), 0.1278, 1e-4).context("p/p0")?;
        close(density_ratio(AIR, 2.0), 0.2300, 1e-4).context("rho/rho0")?;
        close(area_ratio(AIR, 2.0), 1.6875, 1e-4).context("A/A*")?;
        close(area_ratio(AIR, 1.0), 1.0, 1e-15).context("A/A* at M = 1")
    }

    #[test]
    fn test_isentropic_critical_pressure_ratio() -> anyhow::Result<()> {
        close(critical_pressure_ratio(AIR), 0.5283, 1e-4)?;
        close(
            critical_pressure_ratio(AIR),
            pressure_ratio(AIR, 1.0),
            1e-15,
        )
    }

    #[test]
    fn test_isentropic_mass_flux_parameter_peaks_at_sonic() -> anyhow::Result<()> {
        // √γ·(2/(γ+1))^((γ+1)/(2(γ−1))) = 0.6847 for γ = 1.4.
        close(mass_flux_parameter(AIR, 1.0), 0.6847, 1e-4)?;
        for m in [0.5, 0.9, 1.1, 2.0] {
            anyhow::ensure!(
                mass_flux_parameter(AIR, m) < mass_flux_parameter(AIR, 1.0),
                "not maximal at M = 1 (M = {m})"
            );
        }
        Ok(())
    }

    #[test]
    fn test_isentropic_prandtl_meyer() -> anyhow::Result<()> {
        close(prandtl_meyer(AIR, 1.0).get::<degree>(), 0.0, 1e-12)?;
        close(prandtl_meyer(AIR, 2.0).get::<degree>(), 26.38, 1e-2)
    }

    #[test]
    fn test_isentropic_mach_from_area_ratio_both_branches() -> anyhow::Result<()> {
        close(
            mach_from_area_ratio(AIR, 1.6875, Branch::Supersonic)?,
            2.0,
            1e-4,
        )?;
        // Table A.1: A/A* = 2.0 at M = 0.3059 and 2.1972.
        close(
            mach_from_area_ratio(AIR, 2.0, Branch::Subsonic)?,
            0.3059,
            1e-4,
        )?;
        close(
            mach_from_area_ratio(AIR, 2.0, Branch::Supersonic)?,
            2.1972,
            1e-4,
        )?;
        close(mach_from_area_ratio(AIR, 1.0, Branch::Subsonic)?, 1.0, 1e-9)
    }

    #[test]
    fn test_isentropic_mach_from_area_ratio_rejects_below_one() -> anyhow::Result<()> {
        for bad in [0.5, f64::NAN] {
            anyhow::ensure!(
                mach_from_area_ratio(AIR, bad, Branch::Subsonic).is_err(),
                "accepted A/A* = {bad}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_isentropic_mach_from_area_ratio_survives_gamma_near_one() -> anyhow::Result<()> {
        // Exponent (γ+1)/(2(γ−1)) = 100.5: direct A/A* would overflow at large
        // M.
        let m = mach_from_area_ratio(1.01, 50.0, Branch::Supersonic)?;
        close(area_ratio(1.01, m) / 50.0, 1.0, 1e-9)
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_isentropic_area_ratio_round_trips(
            gamma in 1.05..1.67f64, mach in 0.01..8.0f64,
        ) {
            let branch = if mach < 1.0 { Branch::Subsonic } else { Branch::Supersonic };
            let got = mach_from_area_ratio(gamma, area_ratio(gamma, mach), branch)
                .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            prop_assert!((got / mach - 1.0).abs() < 1e-9, "M {mach} → {got}");
        }

        #[test]
        fn test_isentropic_ratios_are_consistent(gamma in 1.05..1.67f64, mach in 0.0..8.0f64) {
            // p/p0 = (ρ/ρ0)·(T/T0) and p/p0 = (T/T0)^(γ/(γ−1)).
            let (t, p, rho) = (
                temperature_ratio(gamma, mach),
                pressure_ratio(gamma, mach),
                density_ratio(gamma, mach),
            );
            prop_assert!((p - rho * t).abs() < 1e-12, "p {p} vs ρT {}", rho * t);
            prop_assert!((p - t.powf(gamma / (gamma - 1.0))).abs() < 1e-12);
        }
    }
}
