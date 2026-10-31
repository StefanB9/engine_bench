//! States of a pure fluid from (p, T) and from what a closed volume conserves,
//! (u, ρ): one phase, or liquid and vapour in equilibrium inside the dome.

use anyhow::Context;

use super::{
    Fluid,
    helmholtz::{Molar, molar},
    saturation::{ancillary, dome_top, saturation_si},
};
use crate::{
    PhysicsError,
    numerics::roots::{RootOptions, brent, newton_bracketed},
};

/// One saturated phase: density, mol/m³, and its properties.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Phase {
    pub(crate) rho: f64,
    pub(crate) properties: Molar,
}

/// Equilibrium state, molar.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct RawState {
    /// K.
    pub(crate) t: f64,
    /// Bulk density, mol/m³.
    pub(crate) rho: f64,
    /// Vapour mass fraction; NaN for a single phase.
    pub(crate) quality: f64,
    /// The single phase, or the mixture's p, u, h and s (`c_v`, `c_p` and w are
    /// NaN for a mixture).
    pub(crate) bulk: Molar,
    /// Saturated liquid and vapour inside the dome.
    pub(crate) phases: Option<(Phase, Phase)>,
}

/// Density, mol/m³, of the stable phase at `p` Pa and `t` K: liquid above the
/// saturation pressure, vapour below it, one fluid above the critical point.
pub(crate) fn density_si(fluid: &Fluid, p: f64, t: f64) -> anyhow::Result<f64> {
    let data = fluid.data;
    check_temperature(fluid, t)?;
    if !(p > 0.0 && p <= data.p_max) {
        return Err(PhysicsError::PressureOutOfRange {
            pressure: p,
            min: 0.0,
            max: data.p_max,
        })
        .with_context(|| format!("{} at {t} K", data.species));
    }
    let ideal = p / (data.gas_constant * t);
    let f = |rho: f64| molar(fluid, rho, t).p - p;
    let (lo, hi, start) = if t < dome_top(fluid) {
        let (p_sat, rho_l, rho_v) = saturation_si(fluid, t)?;
        // On the saturation line itself (to round-off) the bracket would have
        // no sign change: the saturated phase is the answer.
        if p >= p_sat {
            if f(rho_l) >= 0.0 {
                return Ok(rho_l);
            }
            let hi = grow_up(&f, rho_l, max_density(fluid))?;
            (rho_l, hi, rho_l)
        } else {
            if f(rho_v) <= 0.0 {
                return Ok(rho_v);
            }
            (rho_v * 1e-12, rho_v, ideal.min(rho_v))
        }
    } else {
        let hi = grow_up(&f, ideal.min(data.rho_crit), max_density(fluid))?;
        (ideal.min(hi) * 1e-6, hi, ideal.min(hi))
    };
    newton_bracketed(
        |rho| {
            let x = molar(fluid, rho, t);
            (x.p - p, x.dp_drho)
        },
        lo,
        hi,
        start,
        RootOptions::new(1e-14 * hi, 200)?,
    )
    .with_context(|| format!("{} density at {p} Pa, {t} K", data.species))
}

/// Densest state searched, mol/m³: three times the liquid at the triple point
/// (cold, dense supercritical helium reaches ~1.5 times).
fn max_density(fluid: &Fluid) -> f64 {
    3.0 * ancillary(&fluid.data.liquid_density, fluid.data.t_triple)
}

/// Smallest density above `start` (doubling the step) where f turns positive.
fn grow_up(f: &impl Fn(f64) -> f64, start: f64, ceiling: f64) -> anyhow::Result<f64> {
    let mut rho = start;
    let mut step = 1e-3 * start;
    while f(rho) <= 0.0 {
        if rho >= ceiling {
            return Err(PhysicsError::InvalidState("pressure above the equation's range").into());
        }
        rho = (rho + step).min(ceiling);
        step *= 2.0;
    }
    Ok(rho)
}

fn check_temperature(fluid: &Fluid, t: f64) -> anyhow::Result<()> {
    let data = fluid.data;
    if !(data.t_triple..=data.t_max).contains(&t) {
        // Below the triple point the fluid is solid, which is not modelled.
        return Err(PhysicsError::TemperatureOutOfRange {
            temperature: t,
            min: data.t_triple,
            max: data.t_max,
        })
        .with_context(|| data.species.to_string());
    }
    Ok(())
}

/// Equilibrium state with molar internal energy `u` J/mol and density `rho`
/// mol/m³. `hint` is a nearby earlier state that makes the usual update a few
/// iterations.
pub(crate) fn state_uv_si(
    fluid: &Fluid,
    u: f64,
    rho: f64,
    hint: Option<&RawState>,
) -> anyhow::Result<RawState> {
    if !(u.is_finite() && rho.is_finite() && rho > 0.0) {
        return Err(PhysicsError::InvalidState(
            "internal energy and density must be finite, density positive",
        )
        .into());
    }
    if let Some(hint) = hint {
        let near = if hint.phases.is_some() {
            mixture_near(fluid, u, rho, hint.t)
        } else {
            single_phase(fluid, u, rho, Some(hint.t), None)
                .ok()
                .filter(|s| outside_dome(fluid, s.t, rho))
        };
        if let Some(state) = near {
            return Ok(state);
        }
    }
    match inside_dome(fluid, u, rho)? {
        Dome::Inside(state) => Ok(state),
        Dome::Outside(edge) => single_phase(fluid, u, rho, None, edge),
    }
}

/// Where (u, ρ) lies relative to the two-phase dome.
#[expect(
    clippy::large_enum_variant,
    reason = "returned once per cold solve; boxing would allocate on the simulation path"
)]
enum Dome {
    Inside(RawState),
    /// Single phase; the dome's edge temperature at this density, K, when the
    /// density lies within the dome's range. A single phase there is hotter
    /// than the edge: below it the homogeneous equation is unstable and need
    /// not even have `c_v` > 0, so the search must not go there.
    Outside(Option<f64>),
}

/// Single phase at fixed density: u(ρ, T) rises with T (`c_v` > 0).
fn single_phase(
    fluid: &Fluid,
    u: f64,
    rho: f64,
    start: Option<f64>,
    edge: Option<f64>,
) -> anyhow::Result<RawState> {
    let data = fluid.data;
    let (lo, hi) = (edge.unwrap_or(data.t_triple), data.t_max);
    let f = |t: f64| {
        let x = molar(fluid, rho, t);
        (x.u - u, x.cv)
    };
    if f(lo).0 > 0.0 {
        return Err(PhysicsError::InvalidState(
            "internal energy below the triple point (solid is not modelled)",
        ))
        .with_context(|| format!("{} at {rho} mol/m³", data.species));
    }
    if f(hi).0 < 0.0 {
        return Err(PhysicsError::InvalidState(
            "internal energy above the equation's range",
        ))
        .with_context(|| format!("{} at {rho} mol/m³", data.species));
    }
    let start = start.unwrap_or_else(|| f64::midpoint(lo, hi)).clamp(lo, hi);
    let t = newton_bracketed(f, lo, hi, start, RootOptions::new(1e-12, 200)?)
        .with_context(|| format!("{} temperature at {rho} mol/m³", data.species))?;
    Ok(RawState {
        t,
        rho,
        quality: f64::NAN,
        bulk: molar(fluid, rho, t),
        phases: None,
    })
}

/// Whether a single phase at (T, ρ) is stable: outside the dome.
fn outside_dome(fluid: &Fluid, t: f64, rho: f64) -> bool {
    let data = fluid.data;
    if t >= dome_top(fluid) {
        return true;
    }
    // Clearly outside by the ancillaries; otherwise ask the equation.
    let (l, v) = (
        ancillary(&data.liquid_density, t),
        ancillary(&data.vapour_density, t),
    );
    if rho > 1.05 * l || rho < 0.95 * v {
        return true;
    }
    saturation_si(fluid, t).is_ok_and(|(_, rho_l, rho_v)| rho >= rho_l || rho <= rho_v)
}

/// Liquid and vapour at `t` K mixed to bulk density `rho`.
#[expect(
    clippy::many_single_char_names,
    reason = "thermodynamic symbols p, T, u, v, h, s, x"
)]
fn mixture(fluid: &Fluid, t: f64, rho: f64) -> anyhow::Result<RawState> {
    let (p, rho_l, rho_v) = saturation_si(fluid, t)?;
    let (l, v) = (molar(fluid, rho_l, t), molar(fluid, rho_v, t));
    let x = (1.0 / rho - 1.0 / rho_l) / (1.0 / rho_v - 1.0 / rho_l);
    let mix = |a: f64, b: f64| x.mul_add(b - a, a);
    Ok(RawState {
        t,
        rho,
        quality: x,
        bulk: Molar {
            p,
            dp_drho: 0.0,
            dp_dt: f64::NAN,
            u: mix(l.u, v.u),
            h: mix(l.h, v.h),
            s: mix(l.s, v.s),
            cv: f64::NAN,
            cp: f64::NAN,
            w: f64::NAN,
        },
        phases: Some((
            Phase {
                rho: rho_l,
                properties: l,
            },
            Phase {
                rho: rho_v,
                properties: v,
            },
        )),
    })
}

/// Temperature, K, of the dome's edge at density `rho` mol/m³: where the
/// saturated liquid or vapour has that density (the dome's top between them).
/// `None` if no saturated state is that dense or that dilute.
pub(crate) fn dome_edge(fluid: &Fluid, rho: f64) -> anyhow::Result<Option<f64>> {
    let data = fluid.data;
    let top = dome_top(fluid);
    let (_, top_l, top_v) = saturation_si(fluid, top)?;
    if rho > top_v && rho < top_l {
        return Ok(Some(top));
    }
    if rho >= top_l {
        let (_, triple_l, _) = saturation_si(fluid, data.t_triple)?;
        if rho > triple_l {
            return Ok(None);
        }
        return Ok(Some(solve(
            |t| Ok(saturation_si(fluid, t)?.1 - rho),
            data.t_triple,
            top,
        )?));
    }
    let (_, _, triple_v) = saturation_si(fluid, data.t_triple)?;
    if rho < triple_v {
        return Ok(None);
    }
    Ok(Some(solve(
        |t| Ok(saturation_si(fluid, t)?.2 - rho),
        data.t_triple,
        top,
    )?))
}

/// The two-phase state with energy `u` at density `rho`, if (u, ρ) lies in the
/// dome. At fixed ρ the mixture energy rises with T up to the dome's edge, so
/// the state is two-phase exactly when u is below the edge's energy.
fn inside_dome(fluid: &Fluid, u: f64, rho: f64) -> anyhow::Result<Dome> {
    let data = fluid.data;
    let Some(t_edge) = dome_edge(fluid, rho)? else {
        return Ok(Dome::Outside(None));
    };
    if u >= mixture(fluid, t_edge, rho)?.bulk.u {
        return Ok(Dome::Outside(Some(t_edge)));
    }
    if u < mixture(fluid, data.t_triple, rho)?.bulk.u {
        return Err(PhysicsError::InvalidState(
            "internal energy below the triple point (solid is not modelled)",
        ))
        .with_context(|| format!("{} at {rho} mol/m³", data.species));
    }
    let t = solve(
        |t| Ok(mixture(fluid, t, rho)?.bulk.u - u),
        data.t_triple,
        t_edge,
    )?;
    let state = mixture(fluid, t, rho)?;
    Ok(if (0.0..=1.0).contains(&state.quality) {
        Dome::Inside(state)
    } else {
        Dome::Outside(Some(t_edge))
    })
}

/// Warm start from a nearby two-phase temperature: secant iterations on the
/// mixture energy, kept only if they converge to a quality in [0, 1].
fn mixture_near(fluid: &Fluid, u: f64, rho: f64, t_hint: f64) -> Option<RawState> {
    const MAX_ITERATIONS: u32 = 30;
    let data = fluid.data;
    let top = dome_top(fluid);
    let energy = |t: f64| mixture(fluid, t, rho).ok().map(|m| m.bulk.u - u);
    // A hint from the previous step is close: a small first secant step.
    let (mut t0, mut t1) = (t_hint, (t_hint + 1e-6 * t_hint).min(top));
    if (t1 - t0).abs() < 1e-9 * t0 {
        t1 = t0 - 1e-6 * t0;
    }
    let (mut f0, mut f1) = (energy(t0)?, energy(t1)?);
    for _ in 0..MAX_ITERATIONS {
        let slope = (f1 - f0) / (t1 - t0);
        if slope.is_nan() || slope <= 0.0 {
            return None;
        }
        let t2 = (t1 - f1 / slope).clamp(data.t_triple, top);
        (t0, f0) = (t1, f1);
        t1 = t2;
        f1 = energy(t1)?;
        // Converged steps alone are not enough: a secant pinned at a clamp
        // stops moving without matching the energy.
        if (t1 - t0).abs() <= 1e-12 * t1 && f1.abs() <= 1e-9 * u.abs().max(data.gas_constant * t1) {
            let state = mixture(fluid, t1, rho).ok()?;
            return (0.0..=1.0).contains(&state.quality).then_some(state);
        }
    }
    None
}

/// Brent on [lo, hi] where evaluating f can fail.
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
        RootOptions::new(1e-12 * hi, 200)?,
    );
    if let Some(e) = failure {
        return Err(e);
    }
    t
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::fluid::PureFluid;

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
    }

    fn fluids() -> impl Strategy<Value = PureFluid> {
        prop_oneof![
            Just(PureFluid::Nitrogen),
            Just(PureFluid::Argon),
            Just(PureFluid::Helium),
            Just(PureFluid::CarbonDioxide),
        ]
    }

    #[test]
    fn test_density_picks_liquid_above_and_vapour_below_saturation() -> anyhow::Result<()> {
        let co2 = PureFluid::CarbonDioxide.data();
        let (p_sat, rho_l, rho_v) = saturation_si(co2, 250.0)?;
        let liquid = density_si(co2, 1.01 * p_sat, 250.0)?;
        let vapour = density_si(co2, 0.99 * p_sat, 250.0)?;
        anyhow::ensure!(
            liquid > rho_l && vapour < rho_v,
            "ρ_l {liquid} vs {rho_l}, ρ_v {vapour} vs {rho_v}"
        );
        Ok(())
    }

    #[test]
    #[expect(
        clippy::many_single_char_names,
        reason = "thermodynamic symbols p, T, u, v, h, s, x"
    )]
    fn test_state_uv_returns_a_mixture_inside_the_dome() -> anyhow::Result<()> {
        let co2 = PureFluid::CarbonDioxide.data();
        let (p, rho_l, rho_v) = saturation_si(co2, 280.0)?;
        let (l, v) = (molar(co2, rho_l, 280.0), molar(co2, rho_v, 280.0));
        let x = 0.25;
        let volume = (1.0 - x) / rho_l + x / rho_v;
        let u = l.u + x * (v.u - l.u);
        let state = state_uv_si(co2, u, 1.0 / volume, None)?;
        anyhow::ensure!((state.t - 280.0).abs() < 1e-8, "T = {}", state.t);
        anyhow::ensure!((state.quality - x).abs() < 1e-9, "x = {}", state.quality);
        anyhow::ensure!(
            (state.bulk.p / p - 1.0).abs() < 1e-9,
            "p = {}",
            state.bulk.p
        );
        anyhow::ensure!(state.phases.is_some(), "no phases");
        Ok(())
    }

    #[test]
    fn test_state_uv_rejects_solid_and_invalid_input() -> anyhow::Result<()> {
        let co2 = PureFluid::CarbonDioxide.data();
        let (_, rho_l, _) = saturation_si(co2, 217.0)?;
        let cold = molar(co2, rho_l, 217.0).u - 2e4;
        for (u, rho) in [
            (cold, 0.5 * rho_l),
            (0.0, 0.0),
            (0.0, -1.0),
            (f64::NAN, 100.0),
        ] {
            anyhow::ensure!(
                state_uv_si(co2, u, rho, None).is_err(),
                "u {u}, ρ {rho} accepted"
            );
        }
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_state_uv_round_trips_single_phase(fluid in fluids(), tr in 0.0f64..1.0, log_p in 4.0f64..7.5) {
            let f = fluid.data();
            // Lowest temperatures whose melting pressure exceeds the 316 bar
            // tried here: the solid is not modelled.
            let t_low = match fluid {
                PureFluid::Nitrogen => 70.0,
                PureFluid::Argon => 95.0,
                PureFluid::Helium => 12.0,
                PureFluid::CarbonDioxide => 225.0,
            };
            let t = t_low + tr * (600.0 - t_low);
            let p = 10f64.powf(log_p);
            // Stay clear of the saturation line, where (u, ρ) sits on the dome's edge.
            if t < dome_top(f) {
                let (p_sat, _, _) = saturation_si(f, t).map_err(|e| fail(&e))?;
                prop_assume!((p / p_sat - 1.0).abs() > 1e-3);
            }
            let rho = density_si(f, p, t).map_err(|e| fail(&e))?;
            let x = molar(f, rho, t);
            prop_assert!((x.p / p - 1.0).abs() < 1e-9 && x.dp_drho > 0.0, "{:?}", x);
            let state = state_uv_si(f, x.u, rho, None).map_err(|e| fail(&e))?;
            prop_assert!(state.quality.is_nan(), "{:?} at {} K, {} Pa: x = {}", fluid, t, p, state.quality);
            prop_assert!((state.t / t - 1.0).abs() < 1e-9, "{:?}: T {} vs {}", fluid, state.t, t);
        }

        #[test]
        fn test_state_uv_round_trips_two_phase(fluid in fluids(), tr in 0.001f64..0.995, x in 0.001f64..0.999) {
            let f = fluid.data();
            let t = f.data.t_triple + tr * (f.data.t_crit - f.data.t_triple);
            let (p, rho_l, rho_v) = saturation_si(f, t).map_err(|e| fail(&e))?;
            let (l, v) = (molar(f, rho_l, t), molar(f, rho_v, t));
            let volume = (1.0 - x) / rho_l + x / rho_v;
            let u = l.u + x * (v.u - l.u);
            let state = state_uv_si(f, u, 1.0 / volume, None).map_err(|e| fail(&e))?;
            prop_assert!((state.t / t - 1.0).abs() < 1e-8, "{:?}: T {} vs {}", fluid, state.t, t);
            prop_assert!((state.quality - x).abs() < 1e-7, "{:?}: x {} vs {}", fluid, state.quality, x);
            prop_assert!((state.bulk.p / p - 1.0).abs() < 1e-7);
        }

        #[test]
        fn test_state_uv_warm_start_agrees_with_a_cold_solve(fluid in fluids(), tr in 0.01f64..0.98, x in 0.0f64..1.3, nudge in -1.0f64..1.0) {
            // x above 1 continues as vapour heated at constant volume.
            let f = fluid.data();
            let t = f.data.t_triple + tr * (f.data.t_crit - f.data.t_triple);
            let (_, rho_l, rho_v) = saturation_si(f, t).map_err(|e| fail(&e))?;
            let (l, v) = (molar(f, rho_l, t), molar(f, rho_v, t));
            let (u, rho) = if x <= 1.0 {
                (l.u + x * (v.u - l.u), 1.0 / ((1.0 - x) / rho_l + x / rho_v))
            } else {
                (molar(f, rho_v, t * x).u, rho_v)
            };
            let cold = state_uv_si(f, u, rho, None).map_err(|e| fail(&e))?;
            let near = state_uv_si(f, u * (1.0 + 1e-4 * nudge), rho * (1.0 - 1e-4 * nudge), None).map_err(|e| fail(&e))?;
            let warm = state_uv_si(f, u, rho, Some(&near)).map_err(|e| fail(&e))?;
            prop_assert!((warm.t / cold.t - 1.0).abs() < 1e-9, "{:?}: T {} vs {}", fluid, warm.t, cold.t);
            prop_assert!(warm.quality.is_nan() == cold.quality.is_nan());
        }
    }
}
