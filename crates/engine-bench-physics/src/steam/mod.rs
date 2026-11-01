//! Water and steam from IAPWS-IF97, the industrial formulation of the
//! thermodynamic properties of water (IAPWS R7-97(2012)).
//!
//! All five regions: 1 (liquid), 2 (vapour), 3 (near-critical, Helmholtz
//! form, densities from the SR5-05 backward equations refined by Newton),
//! 4 (saturation) and 5 (high-temperature vapour).
//!
//! - [`properties`] at (p, T); [`saturated`] liquid and vapour at T.
//! - [`temperature_from_enthalpy`] / [`temperature_from_entropy`], solved on
//!   the forward equations so they are exactly consistent with [`properties`].
//! - [`state_from_uv`], the equilibrium state from what a closed volume
//!   conserves, liquid and vapour inside the dome; [`SteamVolume`] tracks one.
//!
//! IF97's neighbouring regions disagree by up to ~1e-4 at their boundaries
//! (623.15 K, B23, 1073.15 K); inverse solves there are as exact as that.

mod backward;
#[allow(
    clippy::unreadable_literal,
    clippy::excessive_precision,
    reason = "generated tables reproduce the IF97 coefficients digit for digit"
)]
mod coefficients;
mod condensation;
mod flow;
mod gibbs;
mod helmholtz;
mod inverse;
mod properties;
mod region;
mod region3;
mod region3_backward;
mod saturation;
mod state;
mod transport;
mod vapour;
mod volume;

pub use flow::Water;
pub use inverse::{temperature_from_enthalpy, temperature_from_entropy};
pub use properties::{Saturated, SteamProperties, properties, saturated};
pub use region::{Region, region};
pub use saturation::{saturation_pressure, saturation_temperature};
pub use state::{SteamState, state_from_uv};
pub use volume::SteamVolume;

/// Specific gas constant of water in IF97, J/(kg·K).
const R: f64 = 461.526;
/// Critical temperature, K.
const T_CRIT: f64 = 647.096;
/// Critical pressure, Pa.
const P_CRIT: f64 = 22.064e6;
/// Critical density, kg/m³.
const RHO_CRIT: f64 = 322.0;
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
