//! Temperature from (p, h) or (p, s) in region 3, where every evaluation first
//! solves for the density.

use anyhow::Context;

use super::Given;
use crate::{
    numerics::roots::{RootOptions, newton_bracketed},
    steam::{
        P_CRIT, T_REGION3_MIN,
        helmholtz::region3,
        region::Region,
        region3::{density_si, saturated_si},
        region3_backward::Side,
        saturation,
    },
};

/// Region 3 at p Pa between 623.15 K and the B23 temperature `t_b23`: liquid,
/// two-phase and vapour below the critical pressure, one fluid above it.
pub(super) fn temperature(
    p: f64,
    x: f64,
    given: Given,
    t_b23: f64,
) -> anyhow::Result<(f64, Region)> {
    if p >= P_CRIT {
        return solve_region3(p, x, given, T_REGION3_MIN, t_b23, Side::Either);
    }
    let t_sat = saturation::saturation_temperature_si(p);
    let (_, liquid, vapour) = saturated_si(t_sat)?;
    if x <= given.of(&liquid) {
        solve_region3(p, x, given, T_REGION3_MIN, t_sat, Side::Liquid)
    } else if x < given.of(&vapour) {
        Ok((t_sat, Region::Saturation))
    } else {
        solve_region3(p, x, given, t_sat, t_b23, Side::Vapour)
    }
}

/// Solves given(p, T) = x for T in [lo, hi] in region 3, where each
/// evaluation first solves for the density.
fn solve_region3(
    p: f64,
    x: f64,
    given: Given,
    lo: f64,
    hi: f64,
    side: Side,
) -> anyhow::Result<(f64, Region)> {
    let eval = |t: f64| -> anyhow::Result<(f64, f64)> {
        let state = region3(density_si(p, t, side)?, t).properties;
        Ok((given.of(&state) - x, given.slope(&state, t)))
    };
    // IF97's regions disagree slightly at their boundaries; a value in such a
    // gap is the boundary temperature.
    if eval(lo)?.0 >= 0.0 {
        return Ok((lo, Region::NearCritical));
    }
    if eval(hi)?.0 <= 0.0 {
        return Ok((hi, Region::NearCritical));
    }
    let mut failure = None;
    let t = newton_bracketed(
        |t| {
            eval(t).unwrap_or_else(|e| {
                failure.get_or_insert(e);
                (f64::NAN, f64::NAN)
            })
        },
        lo,
        hi,
        f64::midpoint(lo, hi),
        RootOptions::new(1e-10, 100)?,
    );
    if let Some(e) = failure {
        return Err(e).context("region 3 temperature");
    }
    let t =
        t.with_context(|| format!("region 3 temperature from {} {x} at {p} Pa", given.name()))?;
    Ok((t, Region::NearCritical))
}
