//! Region 3 at (p, T): the density solve on the Helmholtz equation, and the
//! saturated liquid and vapour states of the whole saturation line.

use anyhow::Context;

use super::{
    RHO_CRIT, T_CRIT, T_MIN, T_REGION3_MIN,
    gibbs::{self, RawProperties},
    helmholtz::region3,
    region3_backward::{Side, specific_volume_estimate},
    saturation::saturation_pressure_si,
};
use crate::{
    PhysicsError,
    numerics::roots::{RootOptions, newton_bracketed},
};

/// Search limits for region 3 densities, kg/m³ (the region spans about
/// 100–760 kg/m³).
const RHO_MIN: f64 = 1.0;
const RHO_MAX: f64 = 1200.0;

/// Density, kg/m³, of the region 3 state at p Pa and T K on `side` of the
/// saturation line (`Side::Either` off it).
///
/// Below the critical temperature p(ρ) has a van der Waals loop and three
/// roots near saturation, a few kg/m³ apart close to the critical point. The
/// search is kept on the requested phase's side of the critical density, a
/// bracket is grown from the SR5-05 estimate, and bracketed Newton finishes on
/// p(ρ, T) = p.
pub(crate) fn density_si(p: f64, t: f64, side: Side) -> anyhow::Result<f64> {
    let liquid_like = match side {
        Side::Liquid => true,
        Side::Vapour => false,
        Side::Either => p >= saturation_pressure_si(t.min(T_CRIT)),
    };
    let (floor, ceiling) = if t >= T_CRIT && side == Side::Either {
        (RHO_MIN, RHO_MAX)
    } else if liquid_like {
        (RHO_CRIT, RHO_MAX)
    } else {
        (RHO_MIN, RHO_CRIT)
    };
    let f = |rho: f64| region3(rho, t).p - p;
    let start = specific_volume_estimate(p, t, side)
        .map_or(f64::midpoint(floor, ceiling), |v| 1.0 / v)
        .clamp(floor, ceiling);
    let (lo, hi) = grow_bracket(&f, start, floor, ceiling)
        .with_context(|| format!("no region 3 density for {p} Pa at {t} K"))?;
    newton_bracketed(
        |rho| {
            let state = region3(rho, t);
            (state.p - p, state.dp_drho)
        },
        lo,
        hi,
        start,
        RootOptions::new(1e-11 * start, 100)?,
    )
    .with_context(|| format!("region 3 density at {p} Pa, {t} K"))
}

/// Grows [lo, hi] around `start` in doubling steps, within [floor, ceiling],
/// until f(lo) < 0 < f(hi) (p(ρ) rises with ρ on a stable branch).
fn grow_bracket(
    f: &impl Fn(f64) -> f64,
    start: f64,
    floor: f64,
    ceiling: f64,
) -> anyhow::Result<(f64, f64)> {
    let (mut lo, mut hi) = (start, start);
    let mut step = 1e-4 * start;
    while f(lo) >= 0.0 {
        if lo <= floor {
            return Err(PhysicsError::InvalidState("pressure below the region 3 branch").into());
        }
        lo = (lo - step).max(floor);
        step *= 2.0;
    }
    step = 1e-4 * start;
    while f(hi) <= 0.0 {
        if hi >= ceiling {
            return Err(PhysicsError::InvalidState("pressure above the region 3 branch").into());
        }
        hi = (hi + step).min(ceiling);
        step *= 2.0;
    }
    Ok((lo, hi))
}

/// Saturated liquid and vapour at `t` K (273.15 K up to the critical point),
/// with the saturation pressure, Pa.
pub(crate) fn saturated_si(t: f64) -> anyhow::Result<(f64, RawProperties, RawProperties)> {
    if !(T_MIN..=T_CRIT).contains(&t) {
        return Err(PhysicsError::TemperatureOutOfRange {
            temperature: t,
            min: T_MIN,
            max: T_CRIT,
        })
        .context("saturated states");
    }
    let p = saturation_pressure_si(t);
    if t <= T_REGION3_MIN {
        return Ok((p, gibbs::region1(p, t), gibbs::region2(p, t)));
    }
    let liquid = density_si(p, t, Side::Liquid).context("saturated liquid")?;
    let vapour = density_si(p, t, Side::Vapour).context("saturated vapour")?;
    Ok((
        p,
        region3(liquid, t).properties,
        region3(vapour, t).properties,
    ))
}

#[cfg(test)]
#[expect(
    clippy::unreadable_literal,
    reason = "verification values exactly as IAPWS R7-97(2012) prints them"
)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::steam::saturation::b23_pressure_si;

    #[test]
    fn test_region3_density_inverts_iapws_table_33() -> anyhow::Result<()> {
        for (t, p_mpa, rho) in [
            (650.0, 25.5837018, 500.0),
            (650.0, 22.2930643, 200.0),
            (750.0, 78.3095639, 500.0),
        ] {
            let got = density_si(p_mpa * 1e6, t, Side::Either)?;
            anyhow::ensure!(
                (got / rho - 1.0).abs() < 1e-7,
                "{t} K, {p_mpa} MPa: ρ = {got}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_saturated_states_join_regions_1_and_2_below_623_k() -> anyhow::Result<()> {
        let (p, liquid, vapour) = saturated_si(500.0)?;
        anyhow::ensure!(
            (p - saturation_pressure_si(500.0)).abs() < 1e-9 * p,
            "p = {p}"
        );
        anyhow::ensure!(
            liquid.h.to_bits() == gibbs::region1(p, 500.0).h.to_bits(),
            "liquid"
        );
        anyhow::ensure!(
            vapour.h.to_bits() == gibbs::region2(p, 500.0).h.to_bits(),
            "vapour"
        );
        Ok(())
    }

    #[test]
    fn test_saturated_states_meet_at_the_critical_point() -> anyhow::Result<()> {
        let (_, liquid, vapour) = saturated_si(T_CRIT - 1e-3)?;
        let (rho_l, rho_v) = (1.0 / liquid.v, 1.0 / vapour.v);
        anyhow::ensure!(
            rho_l > RHO_CRIT && rho_v < RHO_CRIT,
            "ρ' = {rho_l}, ρ'' = {rho_v}"
        );
        anyhow::ensure!(rho_l - rho_v < 40.0, "still {} kg/m³ apart", rho_l - rho_v);
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_region3_density_solves_the_pressure_equation(
            t in 623.15f64..863.15,
            fraction in 0.0f64..1.0,
        ) {
            // A pressure between the B23 boundary (or 16.53 MPa) and 100 MPa.
            let lo = b23_pressure_si(t).max(16.6e6);
            let p = lo + fraction * (100e6 - lo);
            let side_ok = t >= T_CRIT || (saturation_pressure_si(t) - p).abs() > 1e3;
            prop_assume!(side_ok);
            let rho = density_si(p, t, Side::Either).map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let state = region3(rho, t);
            prop_assert!((state.p / p - 1.0).abs() < 1e-11, "{} Pa, {} K: p(ρ) = {}", p, t, state.p);
            prop_assert!(state.dp_drho > 0.0, "unstable branch at {} Pa, {} K", p, t);
        }

        #[test]
        fn test_saturated_region3_phases_are_on_the_saturation_pressure(t in 623.2f64..647.0) {
            let (p, liquid, vapour) = saturated_si(t).map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            for x in [liquid, vapour] {
                let p3 = region3(1.0 / x.v, t).p;
                prop_assert!((p3 / p - 1.0).abs() < 1e-11, "{} K: p = {} vs {}", t, p3, p);
            }
            prop_assert!(liquid.v < 1.0 / RHO_CRIT && vapour.v > 1.0 / RHO_CRIT);
            // Equal Gibbs energy, as far as IF97's region 4 and 3 equations agree.
            let g = |x: &RawProperties| x.h - t * x.s;
            prop_assert!((g(&liquid) - g(&vapour)).abs() < 1e-4 * (vapour.h - liquid.h).max(1e4),
                "{} K: g' = {}, g'' = {}", t, g(&liquid), g(&vapour));
        }
    }
}
