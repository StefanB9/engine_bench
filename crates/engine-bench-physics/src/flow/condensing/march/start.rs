//! Where the march starts: the frozen (dry, isentropic) flow from the
//! reservoir through the sonic throat to just downstream of it, where the
//! condensing march takes over.

use anyhow::Context;

use super::state::{Station, saturation};
use crate::{
    PhysicsError,
    flow::{
        condensing::{CondensingFluid, VapourPoint, saturation_temperature},
        isentropic::{Branch, mach_from_area_ratio, pressure_ratio, temperature_ratio},
    },
};

/// The reservoir: stagnation enthalpy, J/kg, and entropy, J/(kg·K), with
/// the vapour's cp/cv there.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Reservoir {
    pub(super) p0: f64,
    pub(super) t0: f64,
    pub(super) h0: f64,
    pub(super) s0: f64,
    pub(super) gamma: f64,
}

impl Reservoir {
    pub(super) fn new<F: CondensingFluid>(fluid: &F, p0: f64, t0: f64) -> anyhow::Result<Self> {
        let a = fluid
            .vapour_variable(p0, t0)
            .context("the reservoir must hold vapour")?;
        let v = fluid.vapour(a, t0);
        Ok(Self {
            p0,
            t0,
            h0: v.h,
            s0: v.s,
            gamma: v.cp / v.cv,
        })
    }

    /// A perfect-gas estimate of the vapour's (p, T) at Mach `mach`.
    fn estimate(&self, mach: f64) -> (f64, f64) {
        let g = self.gamma;
        (
            self.p0 * pressure_ratio(g, mach),
            self.t0 * temperature_ratio(g, mach),
        )
    }
}

/// Newton's method in (a, T) on two residuals with a difference Jacobian;
/// `None` if a residual cannot be evaluated or it does not converge.
#[expect(
    clippy::many_single_char_names,
    reason = "the symbols of the equations: p, T, a, u, h, s"
)]
fn newton2(
    mut f: impl FnMut(f64, f64) -> Option<[f64; 2]>,
    mut a: f64,
    mut t: f64,
) -> Option<(f64, f64)> {
    for _ in 0..100 {
        let r = f(a, t)?;
        let (ha, ht) = (1e-7 * a, 1e-7 * t);
        let ra = f(a + ha, t)?;
        let rt = f(a, t + ht)?;
        let j = [
            [(ra[0] - r[0]) / ha, (rt[0] - r[0]) / ht],
            [(ra[1] - r[1]) / ha, (rt[1] - r[1]) / ht],
        ];
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if !det.is_normal() {
            return None;
        }
        let da = (r[0] * j[1][1] - r[1] * j[0][1]) / det;
        let dt = (j[0][0] * r[1] - j[1][0] * r[0]) / det;
        let shrink = (0.2 * a / da.abs()).min(0.05 * t / dt.abs()).min(1.0);
        a -= shrink * da;
        t -= shrink * dt;
        if shrink >= 1.0 && da.abs() <= 1e-12 * a && dt.abs() <= 1e-12 * t {
            return Some((a, t));
        }
    }
    None
}

/// Speed along the isentrope, m/s: √(2(h₀ − h)), `None` above h₀.
fn speed(r: &Reservoir, v: &VapourPoint) -> Option<f64> {
    let k = 2.0 * (r.h0 - v.h);
    (k > 0.0).then(|| k.sqrt())
}

/// The frozen sonic throat: on the reservoir's isentrope, where the flow
/// speed equals the vapour's speed of sound. Returns the vapour there and
/// the mass flux, kg/(m²·s).
#[expect(
    clippy::many_single_char_names,
    reason = "the symbols of the equations: p, T, a, u, h, s"
)]
pub(super) fn sonic_throat<F: CondensingFluid>(
    fluid: &F,
    r: &Reservoir,
) -> anyhow::Result<(VapourPoint, f64)> {
    let (p, t) = r.estimate(1.0);
    let a = fluid
        .vapour_variable(p, t)
        .context("vapour at the estimated throat")?;
    let (a, t) = newton2(
        |a, t| {
            let v = fluid.vapour(a, t);
            let u = speed(r, &v)?;
            Some([(v.s - r.s0) / v.cp, (u * u - v.w * v.w) / (v.w * v.w)])
        },
        a,
        t,
    )
    .ok_or(PhysicsError::NoConvergence {
        iterations: 100,
        residual: f64::NAN,
    })
    .context("frozen sonic throat")?;
    let v = fluid.vapour(a, t);
    Ok((v, v.rho * v.w))
}

/// The frozen supersonic station where the area is `area_ratio` times the
/// throat's, carrying mass flux `throat_flux` through the throat.
#[expect(
    clippy::many_single_char_names,
    reason = "the symbols of the equations: p, T, a, u, h, s"
)]
pub(super) fn supersonic_station<F: CondensingFluid>(
    fluid: &F,
    r: &Reservoir,
    throat_flux: f64,
    area_ratio: f64,
    gas_constant: f64,
) -> anyhow::Result<Station> {
    let mach = mach_from_area_ratio(r.gamma, area_ratio, Branch::Supersonic)
        .context("start of the march")?;
    let (p, t) = r.estimate(mach);
    let a = fluid
        .vapour_variable(p, t)
        .context("vapour at the start of the march")?;
    let flux = throat_flux / area_ratio;
    let (a, t) = newton2(
        |a, t| {
            let v = fluid.vapour(a, t);
            let u = speed(r, &v)?;
            Some([(v.s - r.s0) / v.cp, (v.rho * u - flux) / flux])
        },
        a,
        t,
    )
    .ok_or(PhysicsError::NoConvergence {
        iterations: 100,
        residual: f64::NAN,
    })
    .context("frozen state at the start of the march")?;
    let vapour = fluid.vapour(a, t);
    let u = speed(r, &vapour).ok_or(PhysicsError::InvalidState("no flow speed"))?;
    if u < vapour.w {
        return Err(PhysicsError::InvalidState("start of the march is not supersonic").into());
    }
    let t_sat = saturation_temperature(fluid, vapour.p, t, gas_constant);
    let droplets = if t_sat.is_finite() {
        saturation(fluid, t_sat).unwrap_or_default()
    } else {
        super::state::Saturation::default()
    };
    Ok(Station {
        a,
        t,
        vapour,
        u,
        rho: vapour.rho,
        y: 0.0,
        t_sat,
        droplets,
    })
}
