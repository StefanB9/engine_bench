//! Water and steam from IAPWS-IF97, the industrial formulation of the
//! thermodynamic properties of water (IAPWS R7-97(2012)).
//!
//! Regions 1 (liquid), 2 (vapour), 4 (saturation) and 5 (high-temperature
//! vapour) are implemented; region 3 (near-critical) and states from (u, v)
//! follow. Temperatures from (p, h) and (p, s) are solved on the forward
//! equations, so they are exactly consistent with [`properties`].

mod backward;
#[allow(
    clippy::unreadable_literal,
    clippy::excessive_precision,
    reason = "generated tables reproduce the IF97 coefficients digit for digit"
)]
mod coefficients;
mod gibbs;
mod inverse;
mod properties;
mod region;
mod saturation;

pub use inverse::{temperature_from_enthalpy, temperature_from_entropy};
pub use properties::{SteamProperties, properties};
pub use region::{Region, region};
pub use saturation::{saturation_pressure, saturation_temperature};

/// Specific gas constant of water in IF97, J/(kg·K).
const R: f64 = 461.526;
/// Critical temperature, K.
const T_CRIT: f64 = 647.096;
/// Critical pressure, Pa.
const P_CRIT: f64 = 22.064e6;
/// Lowest temperature of IF97, K.
const T_MIN: f64 = 273.15;
/// Upper temperature of region 2, K.
const T_MAX_VAPOUR: f64 = 1073.15;
/// Upper temperature of region 5, K.
const T_MAX: f64 = 2273.15;
/// Upper pressure up to 1073.15 K, Pa.
const P_MAX: f64 = 100e6;
/// Upper pressure of region 5, Pa.
const P_MAX_HIGH_TEMPERATURE: f64 = 50e6;
/// Lower temperature of region 3 (and of the B23 boundary), K.
const T_REGION3_MIN: f64 = 623.15;
/// Upper temperature of the B23 boundary, K.
const T_B23_MAX: f64 = 863.15;
/// Saturation pressure at 273.15 K, the low end of the saturation equations,
/// Pa.
const P_SAT_MIN: f64 = 611.213;
