//! What a condensing flow needs to know of its fluid.

use crate::units::{
    AvailableEnergy, Mass, MassDensity, Pressure, SurfaceTension, ThermodynamicTemperature,
};

/// Liquid–vapour properties at one temperature, the liquid supercooled below
/// the triple point (the metastable liquid nucleates first: Ostwald's rule).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CondensationProperties {
    /// Vapour pressure over a flat liquid surface.
    pub vapour_pressure: Pressure,
    /// Density of the liquid.
    pub liquid_density: MassDensity,
    /// Latent heat of vaporisation.
    pub latent_heat: AvailableEnergy,
    /// Liquid–vapour surface tension.
    pub surface_tension: SurfaceTension,
}

/// A fluid that can condense.
pub trait Condensable {
    /// Properties at temperature `t`, up to (not including) the critical
    /// point.
    fn condensation(&self, t: ThermodynamicTemperature) -> anyhow::Result<CondensationProperties>;

    /// Mass of one molecule.
    fn molecular_mass(&self) -> Mass;
}
