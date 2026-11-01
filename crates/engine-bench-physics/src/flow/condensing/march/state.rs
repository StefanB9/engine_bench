//! The flow state at one station: from the conserved quantities (mass flow,
//! stagnation enthalpy, impulse function F = pA + ṁu) and the liquid
//! fraction, the vapour's state by Newton's method on its natural variables.
//!
//! The droplets are taken at the saturation temperature of the local
//! pressure (Gyarmathy: all but the smallest are within a fraction of a
//! kelvin of it), so the liquid's enthalpy is the saturated vapour's less
//! the latent heat, with the saturated vapour's enthalpy reached from the
//! local vapour by its cp.

use crate::{
    flow::condensing::{CondensingFluid, VapourPoint, saturation_temperature},
    units::{
        ThermodynamicTemperature, available_energy::joule_per_kilogram,
        mass_density::kilogram_per_cubic_meter, pressure::pascal,
        surface_tension::newton_per_meter, thermodynamic_temperature::kelvin,
    },
};

/// Newton iterations before a station is given up.
const MAX_ITERATIONS: usize = 60;
/// Relative change in a and T at which Newton has converged.
const TOLERANCE: f64 = 1e-13;

/// Liquid–vapour properties at one temperature, raw SI.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(super) struct Saturation {
    /// Pa.
    pub(super) p: f64,
    /// kg/m³.
    pub(super) rho_liquid: f64,
    /// J/kg.
    pub(super) latent: f64,
    /// N/m.
    pub(super) sigma: f64,
}

/// The saturation properties at `t` K, `None` outside the fluid's range.
pub(super) fn saturation<F: CondensingFluid>(fluid: &F, t: f64) -> Option<Saturation> {
    let x = fluid
        .condensation(ThermodynamicTemperature::new::<kelvin>(t))
        .ok()?;
    Some(Saturation {
        p: x.vapour_pressure.get::<pascal>(),
        rho_liquid: x.liquid_density.get::<kilogram_per_cubic_meter>(),
        latent: x.latent_heat.get::<joule_per_kilogram>(),
        sigma: x.surface_tension.get::<newton_per_meter>(),
    })
}

/// What stays fixed along the march: mass flow, kg/s, stagnation enthalpy,
/// J/kg, and the vapour's gas constant, J/(kg·K).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Conserved {
    pub(super) mass_flow: f64,
    pub(super) h0: f64,
    pub(super) gas_constant: f64,
}

/// The flow at one station.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Station {
    /// The vapour's natural variable and temperature, K.
    pub(super) a: f64,
    pub(super) t: f64,
    pub(super) vapour: VapourPoint,
    /// Velocity, m/s.
    pub(super) u: f64,
    /// Mixture density, kg/m³.
    pub(super) rho: f64,
    /// Liquid mass fraction.
    pub(super) y: f64,
    /// Saturation temperature at the local pressure, K (NaN above the dome).
    pub(super) t_sat: f64,
    /// Saturation properties at `t_sat`, where the droplets are.
    pub(super) droplets: Saturation,
}

impl Station {
    /// Mixture enthalpy, J/kg.
    pub(super) fn enthalpy(&self) -> f64 {
        self.vapour.h + self.y * self.liquid_offset()
    }

    /// hₗ − hᵥ: cp·(Tₛ − T) − L(Tₛ); zero without liquid.
    fn liquid_offset(&self) -> f64 {
        if self.y == 0.0 {
            0.0
        } else {
            self.vapour.cp * (self.t_sat - self.t) - self.droplets.latent
        }
    }

    /// Mixture entropy, J/(kg·K): the liquid's is the saturated vapour's less
    /// L/Tₛ.
    pub(super) fn entropy(&self) -> f64 {
        if self.y == 0.0 {
            return self.vapour.s;
        }
        let liquid = self.vapour.s + self.vapour.cp * (self.t_sat / self.t).ln()
            - self.droplets.latent / self.t_sat;
        (1.0 - self.y) * self.vapour.s + self.y * liquid
    }

    /// Frozen speed of sound of the mixture, m/s: only the vapour is
    /// compressed, (1 − y) of the mass.
    pub(super) fn frozen_sound_speed(&self) -> f64 {
        self.vapour.w * (1.0 - self.y).sqrt()
    }
}

/// The station with impulse function `impulse` N at area `area` m² and
/// liquid fraction `y`, by Newton's method from `guess` (a station nearby
/// on the same branch); `None` if it does not converge.
pub(super) fn solve<F: CondensingFluid>(
    fluid: &F,
    c: &Conserved,
    area: f64,
    impulse: f64,
    y: f64,
    guess: &Station,
) -> Option<Station> {
    let mut s = *guess;
    s.y = y;
    // The droplets' temperature follows the pressure; it enters only
    // through the liquid, so a few outer passes settle it.
    for _ in 0..6 {
        s = newton(fluid, c, area, impulse, &s)?;
        let t_sat = if s.vapour.p.is_finite() {
            let start = if s.t_sat.is_finite() { s.t_sat } else { s.t };
            saturation_temperature(fluid, s.vapour.p, start, c.gas_constant)
        } else {
            f64::NAN
        };
        let changed =
            !(t_sat.is_finite() && s.t_sat.is_finite()) || (t_sat - s.t_sat).abs() > 1e-10 * t_sat;
        s.t_sat = t_sat;
        if t_sat.is_finite() {
            s.droplets = saturation(fluid, t_sat)?;
        } else if y > 0.0 {
            return None;
        }
        if !changed || y == 0.0 {
            return Some(s);
        }
    }
    Some(s)
}

/// Newton's method on mass and energy with the droplet temperature held.
#[expect(
    clippy::many_single_char_names,
    clippy::similar_names,
    reason = "the symbols of the equations: a, T, u, v and their derivatives"
)]
fn newton<F: CondensingFluid>(
    fluid: &F,
    c: &Conserved,
    area: f64,
    impulse: f64,
    start: &Station,
) -> Option<Station> {
    let (mut a, mut t) = (start.a, start.t);
    let y = start.y;
    for _ in 0..MAX_ITERATIONS {
        let v = fluid.vapour(a, t);
        let mut s = Station {
            a,
            t,
            vapour: v,
            u: (impulse - v.p * area) / c.mass_flow,
            rho: 0.0,
            ..*start
        };
        let liquid_volume = if y > 0.0 {
            y / s.droplets.rho_liquid
        } else {
            0.0
        };
        let volume = (1.0 - y) / v.rho + liquid_volume;
        s.rho = 1.0 / volume;
        // Mass: u·A = ṁ·v; energy: h + u²/2 = h₀; both scaled to O(1).
        let mass = (c.mass_flow * volume - s.u * area) / (s.u * area).abs().max(1e-300);
        let energy = (s.enthalpy() + 0.5 * s.u * s.u - c.h0) / (s.u * s.u).max(1e-300);
        let du_da = -area * v.dp_da / c.mass_flow;
        let du_dt = -area * v.dp_dt / c.mass_flow;
        let dv_da = -(1.0 - y) * v.drho_da / (v.rho * v.rho);
        let dv_dt = -(1.0 - y) * v.drho_dt / (v.rho * v.rho);
        let mass_scale = (s.u * area).abs().max(1e-300);
        let energy_scale = (s.u * s.u).max(1e-300);
        let j11 = (c.mass_flow * dv_da - area * du_da) / mass_scale;
        let j12 = (c.mass_flow * dv_dt - area * du_dt) / mass_scale;
        let cooling = if y > 0.0 { y * v.cp } else { 0.0 };
        let j21 = (v.dh_da + s.u * du_da) / energy_scale;
        let j22 = (v.dh_dt - cooling + s.u * du_dt) / energy_scale;
        let det = j11 * j22 - j12 * j21;
        if !det.is_normal() || !mass.is_finite() || !energy.is_finite() {
            return None;
        }
        let da = (mass * j22 - energy * j12) / det;
        let dt = (j11 * energy - j21 * mass) / det;
        // Bounded steps keep the iterate on the branch it started on.
        let shrink = (0.2 * a.abs() / da.abs()).min(0.05 * t / dt.abs()).min(1.0);
        a -= shrink * da;
        t -= shrink * dt;
        if shrink >= 1.0 && da.abs() <= TOLERANCE * a.abs() && dt.abs() <= TOLERANCE * t {
            let v = fluid.vapour(a, t);
            s.a = a;
            s.t = t;
            s.vapour = v;
            s.u = (impulse - v.p * area) / c.mass_flow;
            s.rho = 1.0 / ((1.0 - y) / v.rho + liquid_volume);
            return (s.u > 0.0 && s.rho.is_finite()).then_some(s);
        }
    }
    None
}
