//! The condensation model's settings: nucleation corrections and the growth
//! law's empirical parameters.

use anyhow::Context;

use crate::PhysicsError;

/// Which theory gives the nucleation rate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NucleationTheory {
    /// Classical (Volmer–Frenkel–Zeldovich): capillarity down to the
    /// critical cluster.
    Classical,
    /// Kalikmanov's mean-field kinetic theory: microscopic surface tension
    /// for small clusters, the classical limit for large ones.
    MeanFieldKinetic,
}

/// Nucleation (with Kantrowitz's non-isothermal correction) and Young's
/// (1982) droplet growth law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CondensationModel {
    pub(crate) theory: NucleationTheory,
    /// Multiplies the nucleation rate: 1 for the theories as published,
    /// fitted to measured onsets where they fall short.
    pub(crate) rate_factor: f64,
    /// Young's α, which relates the condensation and evaporation
    /// coefficients.
    pub(crate) alpha: f64,
    /// Young's β, the Langmuir layer thickness in mean free paths.
    pub(crate) beta: f64,
    /// Condensation coefficient `q_c` (1: every molecule that hits sticks).
    pub(crate) condensation_coefficient: f64,
    /// Kantrowitz's correction for the latent heat that warms nucleating
    /// clusters.
    pub(crate) kantrowitz: bool,
    /// Off only for tests of the frozen limit.
    pub(crate) nucleation: bool,
    /// Multiplies the growth rate: 1 except in tests of the equilibrium
    /// limit.
    pub(crate) growth_factor: f64,
}

impl Default for CondensationModel {
    /// Young's own α = 9, β = 0, `q_c` = 1, with the Kantrowitz correction.
    fn default() -> Self {
        Self {
            theory: NucleationTheory::Classical,
            rate_factor: 1.0,
            alpha: 9.0,
            beta: 0.0,
            condensation_coefficient: 1.0,
            kantrowitz: true,
            nucleation: true,
            growth_factor: 1.0,
        }
    }
}

impl CondensationModel {
    /// Young's growth law with parameters `alpha` (0 to 20: the Wet Steam
    /// Modeling Project's participants used 0 to 14) and `beta` (0 to 2).
    pub fn young(alpha: f64, beta: f64) -> anyhow::Result<Self> {
        if !(0.0..=20.0).contains(&alpha) || !(0.0..=2.0).contains(&beta) {
            return Err(PhysicsError::InvalidSetting(
                "growth parameters outside α in [0, 20], β in [0, 2]",
            ))
            .with_context(|| format!("α = {alpha}, β = {beta}"));
        }
        Ok(Self {
            alpha,
            beta,
            ..Self::default()
        })
    }

    /// The same model with condensation coefficient `q` in (0, 1].
    pub fn with_condensation_coefficient(self, q: f64) -> anyhow::Result<Self> {
        if !(q > 0.0 && q <= 1.0) {
            return Err(PhysicsError::InvalidSetting(
                "condensation coefficient outside (0, 1]",
            ))
            .with_context(|| format!("q_c = {q}"));
        }
        Ok(Self {
            condensation_coefficient: q,
            ..self
        })
    }

    /// No nucleation at all: the frozen limit, vapour expanding as far
    /// past saturation as the nozzle takes it.
    #[must_use]
    pub fn frozen() -> Self {
        Self {
            nucleation: false,
            ..Self::default()
        }
    }

    /// The same model with nucleation by `theory`.
    #[must_use]
    pub fn with_theory(self, theory: NucleationTheory) -> Self {
        Self { theory, ..self }
    }

    /// The same model with its nucleation rate multiplied by `factor`
    /// (positive): a calibration to measured onsets.
    pub fn with_rate_factor(self, factor: f64) -> anyhow::Result<Self> {
        if !(factor > 0.0 && factor.is_finite()) {
            return Err(PhysicsError::InvalidSetting(
                "rate factor must be positive and finite",
            ))
            .with_context(|| format!("factor = {factor}"));
        }
        Ok(Self {
            rate_factor: factor,
            ..self
        })
    }

    /// The same model without the Kantrowitz correction (isothermal rate).
    #[must_use]
    pub fn without_kantrowitz(self) -> Self {
        Self {
            kantrowitz: false,
            ..self
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_condensation_model_validates_its_parameters() -> anyhow::Result<()> {
        let m = CondensationModel::young(11.0, 0.0)?;
        anyhow::ensure!((m.alpha - 11.0).abs() < 1e-12 && m.kantrowitz, "{m:?}");
        for (a, b) in [(-1.0, 0.0), (21.0, 0.0), (9.0, 3.0), (f64::NAN, 0.0)] {
            anyhow::ensure!(CondensationModel::young(a, b).is_err(), "accepted {a}, {b}");
        }
        let m = CondensationModel::default().with_condensation_coefficient(0.5)?;
        anyhow::ensure!((m.condensation_coefficient - 0.5).abs() < 1e-12, "q_c");
        anyhow::ensure!(m.with_condensation_coefficient(0.0).is_err(), "q_c = 0");
        anyhow::ensure!(!m.without_kantrowitz().kantrowitz, "kantrowitz");
        anyhow::ensure!(!CondensationModel::frozen().nucleation, "frozen");
        let m = m.with_theory(NucleationTheory::MeanFieldKinetic);
        anyhow::ensure!(m.theory == NucleationTheory::MeanFieldKinetic, "theory");
        anyhow::ensure!(
            m.with_rate_factor(0.0).is_err() && m.with_rate_factor(f64::NAN).is_err(),
            "factor"
        );
        anyhow::ensure!(
            (m.with_rate_factor(1e5)?.rate_factor - 1e5).abs() < 1e-6,
            "factor"
        );
        Ok(())
    }
}
