//! Vapour–liquid equilibrium on the reference equation: the pressure and the
//! liquid and vapour densities at which both phases have equal pressure and
//! equal Gibbs energy, solved by Newton from the ancillary equations.

use anyhow::Context;

use super::{
    Fluid,
    coefficients::{Ancillary, AncillaryForm},
    helmholtz::molar,
};
use crate::PhysicsError;

/// The dome is closed this fraction below `T_c`: the phases' densities merge
/// there and the equilibrium becomes ill-conditioned (0.5 mK for He, 30 mK
/// for CO₂).
pub(crate) const BELOW_CRITICAL: f64 = 1e-4;

/// Highest temperature treated as two-phase, K.
pub(crate) fn dome_top(fluid: &Fluid) -> f64 {
    fluid.data.t_crit * (1.0 - BELOW_CRITICAL)
}

/// Value of an ancillary equation at `t` K.
pub(crate) fn ancillary(a: &Ancillary, t: f64) -> f64 {
    let theta = (1.0 - t / a.t_reducing).max(0.0);
    let mut sum: f64 = a.terms.iter().map(|&(n, e)| n * theta.powf(e)).sum();
    if a.using_tau {
        sum *= a.t_reducing / t;
    }
    match a.form {
        AncillaryForm::VapourDensity => a.reducing_value * sum.exp(),
        AncillaryForm::LiquidDensity => a.reducing_value * (1.0 + sum),
    }
}

/// Saturation pressure, Pa, and liquid and vapour densities, mol/m³, at `t` K.
pub(crate) fn saturation_si(fluid: &Fluid, t: f64) -> anyhow::Result<(f64, f64, f64)> {
    const MAX_ITERATIONS: u32 = 100;
    let data = fluid.data;
    let t_max = dome_top(fluid);
    if !(data.t_triple..=t_max).contains(&t) {
        return Err(PhysicsError::TemperatureOutOfRange {
            temperature: t,
            min: data.t_triple,
            max: t_max,
        })
        .context("saturation");
    }
    let mut rho_l = ancillary(&data.liquid_density, t).max(1.001 * data.rho_crit);
    let mut rho_v = ancillary(&data.vapour_density, t).min(0.999 * data.rho_crit);
    let mut residual = f64::INFINITY;
    for _ in 0..MAX_ITERATIONS {
        let (l, v) = (molar(fluid, rho_l, t), molar(fluid, rho_v, t));
        // Equal p and equal g = h − Ts; ∂g/∂ρ at constant T is (∂p/∂ρ)/ρ.
        let r1 = l.p - v.p;
        let r2 = (l.h - t * l.s) - (v.h - t * v.s);
        // Equal pressure and Gibbs energy to round-off.
        if r1.abs() <= 1e-12 * l.p.abs() && r2.abs() <= 1e-12 * data.gas_constant * t {
            return Ok((l.p, rho_l, rho_v));
        }
        let (j11, j12) = (l.dp_drho, -v.dp_drho);
        let (j21, j22) = (l.dp_drho / rho_l, -v.dp_drho / rho_v);
        let det = j11 * j22 - j12 * j21;
        if !(det.is_finite() && det != 0.0) {
            break;
        }
        let d_l = (r1 * j22 - j12 * r2) / det;
        let d_v = (j11 * r2 - j21 * r1) / det;
        // Keep each phase on its side of the critical density.
        rho_l = (rho_l - d_l.clamp(-0.2 * rho_l, 0.2 * rho_l)).max(1.000_001 * data.rho_crit);
        rho_v = (rho_v - d_v.clamp(-0.5 * rho_v, 0.5 * rho_v)).min(0.999_999 * data.rho_crit);
        residual = (d_l / rho_l).abs().max((d_v / rho_v).abs());
        if residual < 1e-12 {
            return Ok((molar(fluid, rho_l, t).p, rho_l, rho_v));
        }
    }
    Err(PhysicsError::NoConvergence {
        iterations: MAX_ITERATIONS,
        residual,
    })
    .with_context(|| format!("saturation of {} at {t} K", data.species))
}
