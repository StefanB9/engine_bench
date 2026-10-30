//! Per-mass thermodynamic properties of a mixture.
//!
//! The uom-typed methods are the public API; the `*_si` kernels take and return
//! raw SI values for solvers inside this crate, evaluating T powers once for
//! all species.

use super::Mixture;
use crate::{
    thermo::{MOLAR_GAS_CONSTANT, STANDARD_PRESSURE, SpeciesSet, TemperaturePowers},
    units::{
        AvailableEnergy, Pressure, SpecificHeatCapacity, ThermodynamicTemperature,
        available_energy::joule_per_kilogram, pressure::pascal,
        specific_heat_capacity::joule_per_kilogram_kelvin, thermodynamic_temperature::kelvin,
    },
};

impl<S: SpeciesSet<N>, const N: usize> Mixture<S, N> {
    /// Specific heat at constant pressure.
    pub fn cp(&self, t: ThermodynamicTemperature) -> SpecificHeatCapacity {
        SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(self.cp_si(t.get::<kelvin>()))
    }

    /// Specific heat at constant volume.
    pub fn cv(&self, t: ThermodynamicTemperature) -> SpecificHeatCapacity {
        SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(self.cv_si(t.get::<kelvin>()))
    }

    /// Ratio of specific heats cp/cv.
    pub fn gamma(&self, t: ThermodynamicTemperature) -> f64 {
        let t = t.get::<kelvin>();
        self.cp_si(t) / self.cv_si(t)
    }

    /// Specific enthalpy, including enthalpies of formation.
    pub fn enthalpy(&self, t: ThermodynamicTemperature) -> AvailableEnergy {
        AvailableEnergy::new::<joule_per_kilogram>(self.enthalpy_si(t.get::<kelvin>()))
    }

    /// Specific internal energy u = h − RT.
    pub fn internal_energy(&self, t: ThermodynamicTemperature) -> AvailableEnergy {
        AvailableEnergy::new::<joule_per_kilogram>(self.internal_energy_si(t.get::<kelvin>()))
    }

    /// Specific entropy at temperature `t` and pressure `p`, including entropy
    /// of mixing.
    pub fn entropy(&self, t: ThermodynamicTemperature, p: Pressure) -> SpecificHeatCapacity {
        SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(
            self.entropy_si(t.get::<kelvin>(), p.get::<pascal>()),
        )
    }

    /// cp in J/(kg·K).
    #[inline]
    pub(crate) fn cp_si(&self, t: f64) -> f64 {
        let p = TemperaturePowers::from_kelvin(t);
        let mut sum = 0.0;
        for (w, species) in self.molar_weights.iter().zip(S::SPECIES) {
            sum += w * species.cp_over_r(&p);
        }
        MOLAR_GAS_CONSTANT * sum
    }

    /// cv in J/(kg·K).
    #[inline]
    pub(crate) fn cv_si(&self, t: f64) -> f64 {
        self.cp_si(t) - self.gas_constant_si()
    }

    /// h in J/kg.
    #[inline]
    pub(crate) fn enthalpy_si(&self, t: f64) -> f64 {
        self.enthalpy_and_cp_si(t).0
    }

    /// u in J/kg.
    #[inline]
    pub(crate) fn internal_energy_si(&self, t: f64) -> f64 {
        self.enthalpy_si(t) - self.gas_constant_si() * t
    }

    /// (h in J/kg, cp in J/(kg·K)) from one evaluation of the temperature
    /// powers.
    #[inline]
    pub(crate) fn enthalpy_and_cp_si(&self, t: f64) -> (f64, f64) {
        let p = TemperaturePowers::from_kelvin(t);
        let (mut h, mut cp) = (0.0, 0.0);
        for (w, species) in self.molar_weights.iter().zip(S::SPECIES) {
            h += w * species.h_over_rt(&p);
            cp += w * species.cp_over_r(&p);
        }
        (MOLAR_GAS_CONSTANT * t * h, MOLAR_GAS_CONSTANT * cp)
    }

    /// (u in J/kg, cv in J/(kg·K)) from one evaluation of the temperature
    /// powers.
    #[inline]
    pub(crate) fn internal_energy_and_cv_si(&self, t: f64) -> (f64, f64) {
        let (h, cp) = self.enthalpy_and_cp_si(t);
        let r = self.gas_constant_si();
        (h - r * t, cp - r)
    }

    /// s in J/(kg·K) at `p` in Pa: Σ (yᵢ/Mᵢ)·R·(s°ᵢ/R − ln(xᵢ p/p°)).
    #[inline]
    pub(crate) fn entropy_si(&self, t: f64, p: f64) -> f64 {
        let powers = TemperaturePowers::from_kelvin(t);
        let inverse_moles = self.moles_per_kg.recip();
        let mut sum = 0.0;
        for (w, species) in self.molar_weights.iter().zip(S::SPECIES) {
            // Absent species contribute nothing (x ln x → 0).
            if *w > 0.0 {
                sum += w * (species.s_over_r(&powers) - (w * inverse_moles).ln());
            }
        }
        MOLAR_GAS_CONSTANT * (sum - self.moles_per_kg * (p / STANDARD_PRESSURE).ln())
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::thermo::{Air3, Major7, gases};

    fn kelvin_t(t: f64) -> ThermodynamicTemperature {
        ThermodynamicTemperature::new::<kelvin>(t)
    }

    fn major7(raw: [f64; 7]) -> Result<Mixture<Major7, 7>, TestCaseError> {
        let total: f64 = raw.iter().sum();
        Mixture::from_mass_fractions(raw.map(|v| v / total))
            .map_err(|e| TestCaseError::fail(format!("{e:#}")))
    }

    #[test]
    fn test_mixture_air_cp_and_gamma_at_300_k() -> anyhow::Result<()> {
        let air = Mixture::dry_air()?;
        let cp = air.cp(kelvin_t(300.0)).get::<joule_per_kilogram_kelvin>();
        anyhow::ensure!((cp / 1005.0 - 1.0).abs() < 5e-3, "cp = {cp}");
        let gamma = air.gamma(kelvin_t(300.0));
        anyhow::ensure!((gamma - 1.400).abs() < 2e-3, "γ = {gamma}");
        // Vibrational modes lower γ when hot.
        let hot = air.gamma(kelvin_t(2000.0));
        anyhow::ensure!((1.28..1.32).contains(&hot), "γ(2000 K) = {hot}");
        Ok(())
    }

    #[test]
    fn test_mixture_pure_species_matches_species_kernel() -> anyhow::Result<()> {
        let water = Mixture::<Major7, 7>::from_mass_fractions([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0])?;
        let t = 1500.0;
        let p = TemperaturePowers::from_kelvin(t);
        let per_kg = MOLAR_GAS_CONSTANT / gases::H2O.molar_mass_si();
        let h = gases::H2O.h_over_rt(&p) * per_kg * t;
        let got = water.enthalpy(kelvin_t(t)).get::<joule_per_kilogram>();
        anyhow::ensure!((got - h).abs() < 1e-9 * h.abs(), "h {got} vs {h}");
        let s = gases::H2O.s_over_r(&p) * per_kg;
        let got = water
            .entropy(kelvin_t(t), Pressure::new::<pascal>(STANDARD_PRESSURE))
            .get::<joule_per_kilogram_kelvin>();
        anyhow::ensure!((got - s).abs() < 1e-9 * s.abs(), "s {got} vs {s}");
        Ok(())
    }

    #[test]
    fn test_mixture_entropy_of_mixing_is_r_ln_2_per_mole() -> anyhow::Result<()> {
        let half = Mixture::<Air3, 3>::from_mole_fractions([0.5, 0.5, 0.0])?;
        let n2 = Mixture::<Air3, 3>::from_mole_fractions([1.0, 0.0, 0.0])?;
        let o2 = Mixture::<Air3, 3>::from_mole_fractions([0.0, 1.0, 0.0])?;
        let (t, p) = (500.0, 2e5);
        let per_mole = |m: &Mixture<Air3, 3>| m.entropy_si(t, p) / m.moles_per_kg;
        let mixing = per_mole(&half) - f64::midpoint(per_mole(&n2), per_mole(&o2));
        let expected = MOLAR_GAS_CONSTANT * std::f64::consts::LN_2;
        anyhow::ensure!(
            (mixing - expected).abs() < 1e-9,
            "Δs_mix {mixing} vs {expected}"
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_mixture_enthalpy_is_internal_energy_plus_rt(
            raw in proptest::array::uniform7(0.01..1.0f64), t in 100.0..6000.0f64,
        ) {
            let m = major7(raw)?;
            let r = m.gas_constant().get::<joule_per_kilogram_kelvin>();
            let (h, u) = (m.enthalpy_si(t), m.internal_energy_si(t));
            prop_assert!((h - u - r * t).abs() <= 1e-9 * h.abs().max(r * t), "h {h}, u {u}");
            prop_assert!((m.cp_si(t) - m.cv_si(t) - r).abs() < 1e-9 * r, "cp − cv ≠ R");
        }

        #[test]
        fn test_mixture_enthalpy_derivative_is_cp(
            raw in proptest::array::uniform7(0.01..1.0f64), t in 250.0..5900.0f64,
        ) {
            let m = major7(raw)?;
            let dt = 1e-2;
            prop_assume!((t - 1000.0).abs() > 2.0 * dt);
            let dh_dt = (m.enthalpy_si(t + dt) - m.enthalpy_si(t - dt)) / (2.0 * dt);
            prop_assert!((dh_dt / m.cp_si(t) - 1.0).abs() < 1e-6, "dh/dT {dh_dt}");
        }

        #[test]
        fn test_mixture_entropy_pressure_dependence(
            raw in proptest::array::uniform7(0.01..1.0f64), t in 200.0..6000.0f64,
            p in 1e3..1e7f64,
        ) {
            let m = major7(raw)?;
            let r = m.gas_constant().get::<joule_per_kilogram_kelvin>();
            let ds = m.entropy_si(t, 2.0 * p) - m.entropy_si(t, p);
            prop_assert!((ds + r * std::f64::consts::LN_2).abs() < 1e-9 * r, "Δs {ds}");
        }
    }
}
