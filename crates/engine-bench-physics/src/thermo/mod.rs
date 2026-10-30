//! Ideal-gas thermodynamic properties from NASA Glenn 9-coefficient
//! polynomials.
//!
//! Per-species kernels work in dimensionless form (cp/R, h/RT, s/R) on
//! precomputed [`TemperaturePowers`], so a mixture evaluates `ln T` and the
//! powers of T once for all its species. Mixtures and their uom-typed API build
//! on top.

#[allow(
    clippy::unreadable_literal,
    clippy::excessive_precision,
    reason = "generated tables reproduce the NASA coefficients digit for digit"
)]
pub mod gases;
mod mixture;
mod nasa9;
mod species;
mod species_set;

pub use mixture::Mixture;
pub use nasa9::TemperaturePowers;
pub use species::{Element, Species};
pub use species_set::{Air3, Major7, SpeciesSet};

/// Lowest temperature, in kelvin, that temperature inversion will return; below
/// the fit range cp is extrapolated as constant.
pub const TEMPERATURE_FLOOR: f64 = 10.0;

/// Molar gas constant R in J/(mol·K) (exact since the 2019 SI redefinition).
pub const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24;

/// Standard-state pressure p° of the NASA data, in Pa (1 bar).
pub const STANDARD_PRESSURE: f64 = 100_000.0;
