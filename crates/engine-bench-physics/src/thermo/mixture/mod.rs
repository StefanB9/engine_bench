//! Ideal-gas mixtures over a compile-time species set.
//!
//! `mod.rs` owns the type, constructors and composition; `properties` evaluates
//! cp, h, u, s and γ; `inversion` recovers temperature from u or h.

mod inversion;
mod properties;

use std::marker::PhantomData;

use anyhow::Context;

use super::{MOLAR_GAS_CONSTANT, SpeciesSet};
use crate::{
    PhysicsError,
    units::{
        MolarMass, SpecificHeatCapacity, molar_mass::kilogram_per_mole,
        specific_heat_capacity::joule_per_kilogram_kelvin,
    },
};

/// Composition of an ideal-gas mixture of the species in `S`, as mass
/// fractions.
///
/// Thermodynamic properties are per unit mass. The composition is fixed
/// (frozen); chemistry changes it by building a new mixture.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mixture<S: SpeciesSet<N>, const N: usize> {
    mass_fractions: [f64; N],
    /// yᵢ/Mᵢ in mol/kg: the weight of each species' molar property in a
    /// per-mass property.
    molar_weights: [f64; N],
    /// Σ yᵢ/Mᵢ in mol/kg, the inverse mean molar mass.
    moles_per_kg: f64,
    set: PhantomData<S>,
}

/// Fractions must sum to one within this tolerance; they are then renormalised
/// exactly.
const SUM_TOLERANCE: f64 = 1e-6;

impl<S: SpeciesSet<N>, const N: usize> Mixture<S, N> {
    /// Builds a mixture from mass fractions in set order.
    pub fn from_mass_fractions(mass_fractions: [f64; N]) -> anyhow::Result<Self> {
        let y = normalised(mass_fractions).context("mass fractions")?;
        Ok(Self::from_normalised_mass_fractions(y))
    }

    /// Builds a mixture from mole fractions in set order.
    pub fn from_mole_fractions(mole_fractions: [f64; N]) -> anyhow::Result<Self> {
        let x = normalised(mole_fractions).context("mole fractions")?;
        let mut masses = [0.0; N];
        for ((m, xi), species) in masses.iter_mut().zip(x).zip(S::SPECIES) {
            *m = xi * species.molar_mass_si();
        }
        Ok(Self::from_masses_unchecked(&masses))
    }

    /// Builds a mixture from species masses (any positive total), for kernels
    /// that track per-species mass. Inputs must be finite, non-negative and not
    /// all zero; this is not checked (validation happens in the solver step).
    #[inline]
    pub(crate) fn from_masses_unchecked(masses: &[f64; N]) -> Self {
        let inverse_total = masses.iter().sum::<f64>().recip();
        Self::from_normalised_mass_fractions(masses.map(|m| m * inverse_total))
    }

    #[inline]
    fn from_normalised_mass_fractions(mass_fractions: [f64; N]) -> Self {
        let mut molar_weights = [0.0; N];
        for ((w, y), species) in molar_weights.iter_mut().zip(mass_fractions).zip(S::SPECIES) {
            *w = y / species.molar_mass_si();
        }
        Self {
            mass_fractions,
            molar_weights,
            moles_per_kg: molar_weights.iter().sum(),
            set: PhantomData,
        }
    }

    /// Mass fractions in set order.
    pub fn mass_fractions(&self) -> &[f64; N] {
        &self.mass_fractions
    }

    /// Mole fractions in set order.
    pub fn mole_fractions(&self) -> [f64; N] {
        let inverse = self.moles_per_kg.recip();
        self.molar_weights.map(|w| w * inverse)
    }

    /// Mean molar mass.
    pub fn molar_mass(&self) -> MolarMass {
        MolarMass::new::<kilogram_per_mole>(self.moles_per_kg.recip())
    }

    /// Specific gas constant R / M.
    pub fn gas_constant(&self) -> SpecificHeatCapacity {
        SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(self.gas_constant_si())
    }

    /// R / M in J/(kg·K).
    #[inline]
    pub(crate) fn gas_constant_si(&self) -> f64 {
        MOLAR_GAS_CONSTANT * self.moles_per_kg
    }
}

/// Checks fractions are finite, non-negative and sum to one within
/// [`SUM_TOLERANCE`], then rescales them to sum to one exactly.
fn normalised<const N: usize>(fractions: [f64; N]) -> anyhow::Result<[f64; N]> {
    if fractions.iter().any(|f| !f.is_finite() || *f < 0.0) {
        return Err(PhysicsError::InvalidComposition(
            "fractions must be finite and non-negative",
        ))
        .with_context(|| format!("{fractions:?}"));
    }
    let sum: f64 = fractions.iter().sum();
    if (sum - 1.0).abs() > SUM_TOLERANCE {
        return Err(PhysicsError::InvalidComposition("fractions must sum to 1"))
            .with_context(|| format!("sum = {sum}"));
    }
    Ok(fractions.map(|f| f / sum))
}

impl Mixture<super::Air3, 3> {
    /// Dry air (CO₂ omitted): mole fractions N₂ 0.78084, O₂ 0.209476, Ar
    /// 0.00934, normalised.
    pub fn dry_air() -> anyhow::Result<Self> {
        let x: [f64; 3] = [0.780_84, 0.209_476, 0.009_34];
        // Leaving out CO₂ (~0.04 %) means these sum to 0.999656: rescale
        // deliberately.
        let sum: f64 = x.iter().sum();
        Self::from_mole_fractions(x.map(|xi| xi / sum)).context("dry air")
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::thermo::{Air3, Major7};

    fn composition_error(result: anyhow::Result<Mixture<Air3, 3>>) -> bool {
        result.is_err_and(|e| {
            matches!(
                e.downcast_ref::<PhysicsError>(),
                Some(PhysicsError::InvalidComposition(_))
            )
        })
    }

    #[test]
    fn test_mixture_rejects_invalid_fractions() -> anyhow::Result<()> {
        for bad in [
            [0.5, 0.6, -0.1],
            [f64::NAN, 0.5, 0.5],
            [0.5, 0.2, 0.2],
            [0.0, 0.0, 0.0],
            [f64::INFINITY, 0.0, 0.0],
        ] {
            anyhow::ensure!(
                composition_error(Mixture::<Air3, 3>::from_mass_fractions(bad)),
                "mass fractions {bad:?} accepted"
            );
            anyhow::ensure!(
                composition_error(Mixture::<Air3, 3>::from_mole_fractions(bad)),
                "mole fractions {bad:?} accepted"
            );
        }
        Ok(())
    }

    #[test]
    fn test_mixture_renormalises_rounding() -> anyhow::Result<()> {
        let m = Mixture::<Air3, 3>::from_mass_fractions([0.7, 0.2, 0.100_000_1])?;
        let sum: f64 = m.mass_fractions().iter().sum();
        anyhow::ensure!((sum - 1.0).abs() < 1e-15, "sum {sum}");
        Ok(())
    }

    #[test]
    fn test_mixture_dry_air_molar_mass_and_gas_constant() -> anyhow::Result<()> {
        let air = Mixture::dry_air()?;
        let m = air.molar_mass().get::<kilogram_per_mole>();
        anyhow::ensure!((m - 0.028_960).abs() < 2e-6, "M = {m}");
        let r = air.gas_constant().get::<joule_per_kilogram_kelvin>();
        anyhow::ensure!((r - 287.10).abs() < 0.05, "R = {r}");
        Ok(())
    }

    #[test]
    fn test_mixture_from_masses_matches_mass_fractions() -> anyhow::Result<()> {
        let a = Mixture::<Major7, 7>::from_masses_unchecked(&[7.0, 2.0, 0.1, 0.5, 0.3, 0.1, 0.0]);
        let b = Mixture::<Major7, 7>::from_mass_fractions([0.7, 0.2, 0.01, 0.05, 0.03, 0.01, 0.0])?;
        for (x, y) in a.mass_fractions().iter().zip(b.mass_fractions()) {
            anyhow::ensure!((x - y).abs() < 1e-15, "{a:?} vs {b:?}");
        }
        anyhow::ensure!(
            (a.moles_per_kg - b.moles_per_kg).abs() < 1e-12,
            "cached moles differ"
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_mixture_mole_mass_fraction_round_trip(
            raw in proptest::array::uniform7(0.0..1.0f64),
        ) {
            let total: f64 = raw.iter().sum();
            prop_assume!(total > 1e-3);
            let x = raw.map(|v| v / total);
            let m = Mixture::<Major7, 7>::from_mole_fractions(x)
                .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            for (got, want) in m.mole_fractions().iter().zip(x) {
                prop_assert!((got - want).abs() < 1e-12, "{got} vs {want}");
            }
        }
    }
}
