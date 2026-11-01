//! Condensation properties of water: IF97 for the saturated phases, the IAPWS
//! release on the surface tension of ordinary water substance (R1-76(2014)).

use anyhow::Context;

use super::{T_CRIT, T_MIN, Water, properties::saturated};
use crate::{
    PhysicsError,
    flow::condensing::{Condensable, CondensationProperties},
    units::{
        Mass, SurfaceTension, ThermodynamicTemperature, mass::kilogram,
        surface_tension::newton_per_meter, thermodynamic_temperature::kelvin,
    },
};

/// Molar mass of water in IF97, kg/mol.
const MOLAR_MASS: f64 = 0.018_015_268;
/// Avogadro's number, 1/mol.
const AVOGADRO: f64 = 6.022_140_76e23;

/// IAPWS R1-76(2014): σ = B·τ^μ·(1 + b·τ), τ = 1 − T/Tc, in N/m.
fn surface_tension(t: f64) -> f64 {
    const B: f64 = 235.8e-3;
    const LITTLE_B: f64 = -0.625;
    const MU: f64 = 1.256;
    let tau = (1.0 - t / T_CRIT).max(0.0);
    B * tau.powf(MU) * LITTLE_B.mul_add(tau, 1.0)
}

impl Condensable for Water {
    /// From 273.15 K (IF97's lowest temperature: supercooled water is not
    /// covered) up to the critical point.
    fn condensation(&self, t: ThermodynamicTemperature) -> anyhow::Result<CondensationProperties> {
        let kelvin_t = t.get::<kelvin>();
        if !(T_MIN..T_CRIT).contains(&kelvin_t) {
            return Err(PhysicsError::TemperatureOutOfRange {
                temperature: kelvin_t,
                min: T_MIN,
                max: T_CRIT,
            })
            .context("water condensation properties");
        }
        let sat = saturated(t).context("water condensation properties")?;
        Ok(CondensationProperties {
            vapour_pressure: sat.pressure,
            liquid_density: sat.liquid.density(),
            latent_heat: sat.vapour.enthalpy - sat.liquid.enthalpy,
            surface_tension: SurfaceTension::new::<newton_per_meter>(surface_tension(kelvin_t)),
        })
    }

    fn molecular_mass(&self) -> Mass {
        Mass::new::<kilogram>(MOLAR_MASS / AVOGADRO)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::{
        available_energy::kilojoule_per_kilogram, mass_density::kilogram_per_cubic_meter,
    };

    fn at(celsius: f64) -> anyhow::Result<CondensationProperties> {
        Water.condensation(ThermodynamicTemperature::new::<kelvin>(celsius + 273.15))
    }

    #[test]
    fn test_water_surface_tension_matches_the_iapws_table() -> anyhow::Result<()> {
        // IAPWS R1-76(2014), Table 1, mN/m.
        for (celsius, millinewton) in [(0.01, 75.65), (25.0, 71.97), (100.0, 58.91)] {
            let sigma = at(celsius)?.surface_tension.get::<newton_per_meter>() * 1e3;
            anyhow::ensure!(
                (sigma - millinewton).abs() < 0.006,
                "{celsius} °C: {sigma} mN/m"
            );
        }
        Ok(())
    }

    #[test]
    fn test_water_latent_heat_and_liquid_at_the_boiling_point() -> anyhow::Result<()> {
        // At 100 °C: about 2256.5 kJ/kg and 958.35 kg/m³ (IF97).
        let x = at(99.974)?;
        let latent = x.latent_heat.get::<kilojoule_per_kilogram>();
        anyhow::ensure!((latent / 2256.5 - 1.0).abs() < 5e-4, "latent heat {latent}");
        let rho = x.liquid_density.get::<kilogram_per_cubic_meter>();
        anyhow::ensure!((rho / 958.35 - 1.0).abs() < 5e-4, "liquid density {rho}");
        Ok(())
    }

    #[test]
    fn test_water_condensation_needs_if97_temperatures() {
        assert!(at(-5.0).is_err());
        assert!(at(374.0).is_err());
        let m = Water.molecular_mass().get::<kilogram>();
        assert!((m / 2.991_5e-26 - 1.0).abs() < 1e-4, "{m}");
    }
}
