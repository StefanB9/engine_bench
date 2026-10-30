//! One NASA Glenn 9-coefficient interval (McBride et al.,
//! `NASA/TP-2002-211556`, eqs. 1-3):
//!
//! - cp/R  = a1 T⁻² + a2 T⁻¹ + a3 + a4 T + a5 T² + a6 T³ + a7 T⁴
//! - h/RT  = −a1 T⁻² + a2 ln(T)/T + a3 + a4 T/2 + a5 T²/3 + a6 T³/4 + a7 T⁴/5 +
//!   b1/T
//! - s°/R  = −a1 T⁻²/2 − a2 T⁻¹ + a3 ln T + a4 T + a5 T²/2 + a6 T³/3 + a7 T⁴/4
//!   + b2

use crate::units::{ThermodynamicTemperature, thermodynamic_temperature::kelvin};

/// Powers of one temperature, shared by every species evaluated at it.
///
/// Requires T > 0; kernels do not check (validation happens at the edges).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TemperaturePowers {
    pub(super) t: f64,
    pub(super) inv_t: f64,
    pub(super) inv_t2: f64,
    pub(super) t2: f64,
    pub(super) t3: f64,
    pub(super) t4: f64,
    pub(super) ln_t: f64,
}

impl TemperaturePowers {
    /// Precomputes the powers of `t`.
    pub fn new(t: ThermodynamicTemperature) -> Self {
        Self::from_kelvin(t.get::<kelvin>())
    }

    /// Precomputes the powers of `t` given in kelvin.
    #[inline]
    pub(crate) fn from_kelvin(t: f64) -> Self {
        let inv_t = t.recip();
        let t2 = t * t;
        Self {
            t,
            inv_t,
            inv_t2: inv_t * inv_t,
            t2,
            t3: t2 * t,
            t4: t2 * t2,
            ln_t: t.ln(),
        }
    }

    /// The temperature in kelvin.
    pub fn kelvin(&self) -> f64 {
        self.t
    }
}

/// cp/R of one interval.
#[inline]
pub(super) fn cp_over_r(a: &[f64; 7], p: &TemperaturePowers) -> f64 {
    a[0] * p.inv_t2 + a[1] * p.inv_t + a[2] + a[3] * p.t + a[4] * p.t2 + a[5] * p.t3 + a[6] * p.t4
}

/// h/(RT) of one interval.
#[inline]
pub(super) fn h_over_rt(a: &[f64; 7], b: &[f64; 2], p: &TemperaturePowers) -> f64 {
    -a[0] * p.inv_t2
        + a[1] * p.ln_t * p.inv_t
        + a[2]
        + a[3] * p.t / 2.0
        + a[4] * p.t2 / 3.0
        + a[5] * p.t3 / 4.0
        + a[6] * p.t4 / 5.0
        + b[0] * p.inv_t
}

/// s°/R of one interval at the standard pressure.
#[inline]
pub(super) fn s_over_r(a: &[f64; 7], b: &[f64; 2], p: &TemperaturePowers) -> f64 {
    -a[0] * p.inv_t2 / 2.0 - a[1] * p.inv_t
        + a[2] * p.ln_t
        + a[3] * p.t
        + a[4] * p.t2 / 2.0
        + a[5] * p.t3 / 3.0
        + a[6] * p.t4 / 4.0
        + b[1]
}
