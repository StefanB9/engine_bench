//! One phase: a 2D Newton solve on (ln p, T) for the Gibbs regions 1, 2 and 5,
//! and a 1D solve on T at fixed density for the Helmholtz region 3. A solution
//! counts only if it lies in the region it was solved in and is stable there.

use super::{Phases, RawState};
use crate::{
    numerics::roots::{RootOptions, newton_bracketed},
    steam::{
        P_MAX, R, RHO_CRIT, T_B23_MAX, T_CRIT, T_MAX, T_MAX_VAPOUR, T_MIN, T_REGION3_MIN,
        gibbs::{self, RawProperties},
        helmholtz::region3,
        region::{Region, region_si, validate_pt},
        saturation::saturation_pressure_si,
    },
};

/// Specific internal energy of saturated vapour near the triple point, J/kg,
/// for an ideal-gas starting guess.
const U_VAPOUR_TRIPLE: f64 = 2.3755e6;

/// The state in `region` with internal energy `u` and volume `v`, starting
/// from `start` = (p, T) if given.
pub(super) fn solve(region: Region, u: f64, v: f64, start: Option<(f64, f64)>) -> Option<RawState> {
    match region {
        Region::CompressedLiquid => {
            let start = start.unwrap_or_else(|| {
                let t = (T_MIN + u / 4180.0).clamp(T_MIN, T_REGION3_MIN);
                (1e7, t)
            });
            gibbs_region(region, gibbs::region1, u, v, start)
        }
        Region::Vapour => {
            // Near the dome the given start (its edge, or a hint) is close.
            // Dense, hot steam is reached from the hot end instead, where it
            // is nearly ideal: at fixed v, cooling from there stays in region
            // 2.
            let near = start.unwrap_or_else(|| {
                let t = (T_MIN + (u - U_VAPOUR_TRIPLE) / 1500.0).clamp(T_MIN, T_MAX_VAPOUR);
                (R * t / v, t)
            });
            let hot = ((R * T_MAX_VAPOUR / v).min(P_MAX), T_MAX_VAPOUR);
            gibbs_region(region, gibbs::region2, u, v, near)
                .or_else(|| gibbs_region(region, gibbs::region2, u, v, hot))
        }
        Region::HighTemperature => {
            let start = start.unwrap_or((R * 1500.0 / v, 1500.0));
            gibbs_region(region, gibbs::region5, u, v, start)
        }
        Region::NearCritical => region3_state(u, v, start.map(|(_, t)| t)),
        Region::Saturation => None,
    }
}

#[expect(
    clippy::many_single_char_names,
    reason = "thermodynamic symbols p, T, u, v, h, s, x"
)]
fn gibbs_region(
    region: Region,
    kernel: fn(f64, f64) -> RawProperties,
    u: f64,
    v: f64,
    (p0, t0): (f64, f64),
) -> Option<RawState> {
    const MAX_ITERATIONS: u32 = 80;
    // Vapour volume goes as 1/p, liquid volume nearly linearly in p: Newton
    // on ln p for the gases and on p for the liquid keeps the steps sensible.
    let logarithmic = region != Region::CompressedLiquid;
    let (mut q, mut t) = (if logarithmic { p0.ln() } else { p0 }, t0);
    for _ in 0..MAX_ITERATIONS {
        let p = if logarithmic { q.exp() } else { q };
        let x = kernel(p, t);
        let (r_v, r_u) = (x.v - v, x.u - u);
        // ∂v/∂q, ∂v/∂T and the same for u = h − pv.
        let dp_dq = if logarithmic { p } else { 1.0 };
        let v_q = dp_dq * x.dv_dp;
        let v_t = x.dv_dt;
        let u_q = -dp_dq * (t * x.dv_dt + p * x.dv_dp);
        let u_t = x.cp - p * x.dv_dt;
        let det = v_q * u_t - v_t * u_q;
        if !(det.is_finite() && det != 0.0) {
            return None;
        }
        let max_q_step = if logarithmic { 2.0 } else { 0.25 * P_MAX };
        let d_q = ((r_v * u_t - v_t * r_u) / det).clamp(-max_q_step, max_q_step);
        let d_t = ((v_q * r_u - r_v * u_q) / det).clamp(-100.0, 100.0);
        // A liquid's pressure step can overshoot below zero; shrink towards
        // zero instead.
        q = if logarithmic || q - d_q > 0.0 {
            q - d_q
        } else {
            0.1 * q
        };
        // Iterates may leave IF97's range a little; a solution on its edge
        // (273.15 K) would otherwise be clamped short of convergence.
        t = (t - d_t).clamp(T_MIN - 10.0, T_MAX + 10.0);
        let p = if logarithmic { q.exp() } else { q };
        // Iterates may pass through pressures above IF97's range; the result
        // is checked against it below.
        if !(p.is_finite() && p > 0.0 && p < 3.0 * P_MAX) {
            return None;
        }
        let relative_step = if logarithmic {
            d_q.abs()
        } else {
            (d_q / p).abs()
        };
        // A nearly incompressible liquid fixes p only to an ulp of v over
        // ∂v/∂p (~1e-10 at low pressure), so residuals at round-off also count
        // as converged.
        let at_round_off =
            r_v.abs() <= 4.0 * f64::EPSILON * v && r_u.abs() <= 1e-12 * (u.abs() + R * t);
        if at_round_off || (relative_step < 1e-12 && d_t.abs() < 1e-11 * t) {
            let t = snap_onto(t, T_MIN, T_MAX);
            let p = snap_onto(p, 0.0, P_MAX);
            if validate_pt(p, t).is_err() || region_si(p, t) != region {
                return None;
            }
            let x = kernel(p, t);
            return Some(RawState {
                region,
                p,
                t,
                quality: f64::NAN,
                v,
                u,
                h: x.h,
                s: x.s,
                phases: Phases::NONE,
            });
        }
    }
    None
}

/// `x` moved onto [lo, hi] if it is outside by no more than round-off.
fn snap_onto(x: f64, lo: f64, hi: f64) -> f64 {
    let tolerance = 1e-9 * x.abs();
    if x < lo && lo - x <= tolerance {
        lo
    } else if x > hi && x - hi <= tolerance {
        hi
    } else {
        x
    }
}

/// Region 3 at fixed density 1/v: u(ρ, T) rises with T (cv > 0).
#[expect(
    clippy::many_single_char_names,
    reason = "thermodynamic symbols p, T, u, v, h, s, x"
)]
fn region3_state(u: f64, v: f64, t_start: Option<f64>) -> Option<RawState> {
    let rho = 1.0 / v;
    let lo = T_REGION3_MIN;
    let hi = T_B23_MAX;
    let f = |t: f64| {
        let x = region3(rho, t).properties;
        (x.u - u, x.cv)
    };
    if f(lo).0 > 0.0 || f(hi).0 < 0.0 {
        return None;
    }
    let guess = t_start.unwrap_or_else(|| f64::midpoint(lo, hi));
    let t = newton_bracketed(f, lo, hi, guess, RootOptions::new(1e-11, 100).ok()?).ok()?;
    let state = region3(rho, t);
    let p = state.p;
    let stable = state.dp_drho > 0.0
        && (t >= T_CRIT || {
            let p_sat = saturation_pressure_si(t);
            if rho > RHO_CRIT {
                p >= p_sat
            } else {
                p <= p_sat
            }
        });
    if !stable || validate_pt(p, t).is_err() || region_si(p, t) != Region::NearCritical {
        return None;
    }
    Some(RawState {
        region: Region::NearCritical,
        p,
        t,
        quality: f64::NAN,
        v,
        u,
        h: state.properties.h,
        s: state.properties.s,
        phases: Phases::NONE,
    })
}
