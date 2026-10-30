//! Calorically perfect gas: constant γ and specific gas constant.

use anyhow::Context;

use crate::{
    PhysicsError,
    thermo::{Mixture, SpeciesSet},
    units::{
        SpecificHeatCapacity, ThermodynamicTemperature,
        specific_heat_capacity::joule_per_kilogram_kelvin, thermodynamic_temperature::kelvin,
    },
};

/// A gas with constant γ = cp/cv and specific gas constant R.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PerfectGas {
    pub(crate) gamma: f64,
    /// J/(kg·K).
    pub(crate) gas_constant: f64,
}

impl PerfectGas {
    /// A gas with ratio of specific heats `gamma` (> 1) and specific gas
    /// constant `gas_constant` (> 0).
    pub fn new(gamma: f64, gas_constant: SpecificHeatCapacity) -> anyhow::Result<Self> {
        if !(gamma > 1.0 && gamma.is_finite()) {
            return Err(PhysicsError::InvalidState("gamma must exceed 1"))
                .with_context(|| format!("gamma = {gamma}"));
        }
        let gas_constant = gas_constant.get::<joule_per_kilogram_kelvin>();
        if !(gas_constant > 0.0 && gas_constant.is_finite()) {
            return Err(PhysicsError::InvalidState("gas constant must be positive"))
                .with_context(|| format!("R = {gas_constant} J/(kg K)"));
        }
        Ok(Self {
            gamma,
            gas_constant,
        })
    }

    /// Freezes `mixture` at temperature `t`: γ and R evaluated there.
    pub fn frozen<S: SpeciesSet<N>, const N: usize>(
        mixture: &Mixture<S, N>,
        t: ThermodynamicTemperature,
    ) -> anyhow::Result<Self> {
        let kelvin_value = t.get::<kelvin>();
        if !(kelvin_value > 0.0 && kelvin_value.is_finite()) {
            return Err(PhysicsError::InvalidState("temperature must be positive"))
                .with_context(|| format!("T = {kelvin_value} K"));
        }
        Self::new(mixture.gamma(t), mixture.gas_constant()).context("freezing mixture")
    }

    /// Ratio of specific heats.
    pub fn gamma(&self) -> f64 {
        self.gamma
    }

    /// Specific gas constant.
    pub fn gas_constant(&self) -> SpecificHeatCapacity {
        SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(self.gas_constant)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(value: f64) -> SpecificHeatCapacity {
        SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(value)
    }

    #[test]
    fn test_perfect_gas_validates_inputs() -> anyhow::Result<()> {
        for (gamma, gas_constant) in [(1.0, 287.0), (0.9, 287.0), (f64::NAN, 287.0), (1.4, 0.0)] {
            anyhow::ensure!(
                PerfectGas::new(gamma, r(gas_constant)).is_err(),
                "accepted γ {gamma}, R {gas_constant}"
            );
        }
        let air = PerfectGas::new(1.4, r(287.0))?;
        anyhow::ensure!((air.gamma() - 1.4).abs() < 1e-15, "γ {}", air.gamma());
        let got = air.gas_constant().get::<joule_per_kilogram_kelvin>();
        anyhow::ensure!((got - 287.0).abs() < 1e-12, "R {got}");
        Ok(())
    }

    #[test]
    fn test_perfect_gas_frozen_from_mixture() -> anyhow::Result<()> {
        let air = Mixture::dry_air()?;
        let t = ThermodynamicTemperature::new::<kelvin>(300.0);
        let gas = PerfectGas::frozen(&air, t)?;
        anyhow::ensure!(
            (gas.gamma() - air.gamma(t)).abs() < 1e-15,
            "γ {}",
            gas.gamma()
        );
        anyhow::ensure!(gas.gas_constant() == air.gas_constant(), "R differs");
        Ok(())
    }
}
