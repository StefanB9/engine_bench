//! Liquid and vapour in equilibrium: the dome boundary at a given specific
//! volume, and the temperature at which the mixture has the given energy.

use std::sync::LazyLock;

use anyhow::Context;

use super::{Phases, RawState};
use crate::{
    PhysicsError,
    numerics::roots::{RootOptions, brent},
    steam::{T_CRIT, T_MIN, gibbs::RawProperties, region::Region, region3::saturated_si},
};

/// Highest temperature treated as two-phase. IF97's saturation-pressure
/// equation and its region 3 equation place the critical point ~0.2 kg/m³
/// apart, so the dome is closed 1 mK below the critical temperature; the
/// sliver above is single phase.
pub(super) const T_DOME_MAX: f64 = T_CRIT - 1e-3;
/// Liquid water is densest near 4 °C; v' rises monotonically above this.
const T_DENSEST: f64 = 277.13;

/// Saturated states at the fixed temperatures the boundary search starts from.
struct Anchors {
    triple: (RawProperties, RawProperties),
    densest: RawProperties,
    top: (RawProperties, RawProperties),
}

static ANCHORS: LazyLock<Option<Anchors>> = LazyLock::new(|| {
    let (_, triple_l, triple_v) = saturated_si(T_MIN).ok()?;
    let (_, densest, _) = saturated_si(T_DENSEST).ok()?;
    let (_, top_l, top_v) = saturated_si(T_DOME_MAX).ok()?;
    Some(Anchors {
        triple: (triple_l, triple_v),
        densest,
        top: (top_l, top_v),
    })
});

/// Mixture at `t` K with specific volume `v`: quality and the mixed state.
fn mixture(t: f64, v: f64) -> anyhow::Result<RawState> {
    let (p, liquid, vapour) = saturated_si(t)?;
    let x = (v - liquid.v) / (vapour.v - liquid.v);
    let mix = |a: f64, b: f64| x.mul_add(b - a, a);
    Ok(RawState {
        region: Region::Saturation,
        p,
        t,
        quality: x,
        v,
        u: mix(liquid.u, vapour.u),
        h: mix(liquid.h, vapour.h),
        s: mix(liquid.s, vapour.s),
        phases: Phases {
            v_liquid: liquid.v,
            v_vapour: vapour.v,
            h_liquid: liquid.h,
            h_vapour: vapour.h,
        },
    })
}

/// Solves f(T) = 0 on [lo, hi] with Brent, where evaluating f can fail.
fn solve(mut f: impl FnMut(f64) -> anyhow::Result<f64>, lo: f64, hi: f64) -> anyhow::Result<f64> {
    let mut failure = None;
    let t = brent(
        |t| {
            f(t).unwrap_or_else(|e| {
                failure.get_or_insert(e);
                f64::NAN
            })
        },
        lo,
        hi,
        RootOptions::new(1e-11, 200)?,
    );
    if let Some(e) = failure {
        return Err(e);
    }
    t
}

/// Temperature, K, at which the dome's edge has specific volume `v`, or
/// `None` if no two-phase state has this volume.
fn dome_edge(v: f64) -> anyhow::Result<Option<f64>> {
    let anchors = ANCHORS
        .as_ref()
        .ok_or(PhysicsError::InvalidState("saturated anchor states failed"))?;
    let (top_l, top_v) = &anchors.top;
    if v > top_l.v && v < top_v.v {
        return Ok(Some(T_DOME_MAX));
    }
    if v >= top_v.v {
        if v > anchors.triple.1.v {
            return Ok(None);
        }
        let t = solve(|t| Ok(saturated_si(t)?.2.v - v), T_MIN, T_DOME_MAX)
            .context("dome edge on the vapour side")?;
        return Ok(Some(t));
    }
    if v < anchors.densest.v {
        return Ok(None);
    }
    let t = solve(|t| Ok(saturated_si(t)?.1.v - v), T_DENSEST, T_DOME_MAX)
        .context("dome edge on the liquid side")?;
    Ok(Some(t))
}

/// Where (u, v) lies relative to the two-phase dome.
pub(super) enum Dome {
    /// Liquid and vapour in equilibrium.
    Inside(RawState),
    /// Single phase. The saturated state on the dome's edge with the same
    /// volume, (p, T), when there is one: a valid state close to the solution
    /// for a Newton solve to start from.
    Outside(Option<(f64, f64)>),
}

/// The two-phase state with internal energy `u` and volume `v`, if (u, v) lies
/// inside the dome. At fixed v the mixture energy rises with temperature up
/// to the dome's edge, so the state is two-phase exactly when u is below the
/// edge's energy.
pub(super) fn from_uv(u: f64, v: f64) -> anyhow::Result<Dome> {
    let Some(t_edge) = dome_edge(v)? else {
        return Ok(Dome::Outside(None));
    };
    let edge = mixture(t_edge, v)?;
    if u >= edge.u {
        return Ok(Dome::Outside(Some((edge.p, t_edge))));
    }
    // Below the mixture energy at 273.15 K is either ice (outside IF97, so no
    // region will find a state either) or, between 0 and 4 °C where liquid
    // water expands on cooling, compressed liquid.
    if u < mixture(T_MIN, v)?.u {
        return Ok(Dome::Outside(Some((edge.p, t_edge))));
    }
    let t = solve(|t| Ok(mixture(t, v)?.u - u), T_MIN, t_edge).context("two-phase temperature")?;
    let state = mixture(t, v)?;
    Ok(if (0.0..=1.0).contains(&state.quality) {
        Dome::Inside(state)
    } else {
        Dome::Outside(None)
    })
}

/// Warm start from a nearby two-phase temperature: secant iterations on the
/// mixture energy. `None` unless it converges to a quality in [0, 1].
pub(super) fn from_uv_near(u: f64, v: f64, t_hint: f64) -> Option<RawState> {
    const MAX_ITERATIONS: u32 = 30;
    let energy = |t: f64| mixture(t, v).ok().map(|m| m.u - u);
    let (mut t0, mut t1) = (t_hint, (t_hint + 0.01).min(T_DOME_MAX));
    if (t1 - t0).abs() < 1e-6 {
        t1 = t0 - 0.01;
    }
    let (mut f0, mut f1) = (energy(t0)?, energy(t1)?);
    for _ in 0..MAX_ITERATIONS {
        let slope = (f1 - f0) / (t1 - t0);
        // A NaN or non-positive slope means the secant has lost the root.
        if slope.is_nan() || slope <= 0.0 {
            return None;
        }
        let t2 = (t1 - f1 / slope).clamp(T_MIN, T_DOME_MAX);
        (t0, f0) = (t1, f1);
        t1 = t2;
        f1 = energy(t1)?;
        // Converged steps alone are not enough: a secant pinned at a clamp
        // stops moving without matching the energy.
        if (t1 - t0).abs() <= 1e-11 * t1 && f1.abs() <= 1e-9 * u.abs().max(1e5) {
            let state = mixture(t1, v).ok()?;
            return (state.quality >= 0.0 && state.quality <= 1.0).then_some(state);
        }
    }
    None
}
