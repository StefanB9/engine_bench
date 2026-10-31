//! A gas-phase species and its NASA Glenn coefficient intervals.

use super::nasa9::{self, TemperaturePowers};
use crate::units::{
    MolarEnergy, MolarMass, ThermodynamicTemperature, molar_energy::joule_per_mole,
    molar_mass::kilogram_per_mole, thermodynamic_temperature::kelvin,
};

/// Chemical elements present in the vendored species data.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Element {
    /// Hydrogen.
    H,
    /// Carbon.
    C,
    /// Nitrogen.
    N,
    /// Oxygen.
    O,
    /// Argon.
    Ar,
    /// Helium.
    He,
}

/// One temperature interval `[t_min, t_max]` of a species fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    pub(super) t_min: f64,
    pub(super) t_max: f64,
    pub(super) a: [f64; 7],
    pub(super) b: [f64; 2],
}

/// An ideal-gas species with its thermodynamic fit.
///
/// Constants for every vendored species live in [`gases`](super::gases).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Species {
    pub(super) name: &'static str,
    /// kg/mol.
    pub(super) molar_mass: f64,
    pub(super) elements: &'static [(Element, f64)],
    /// J/mol at 298.15 K.
    pub(super) enthalpy_of_formation: f64,
    /// Contiguous, ascending, non-empty (checked by `build.rs`).
    pub(super) intervals: &'static [Interval],
}

impl Species {
    /// Formula name as in the NASA data, e.g. `"CO2"`.
    pub fn name(&self) -> &'static str {
        self.name
    }

    /// Molar mass.
    pub fn molar_mass(&self) -> MolarMass {
        MolarMass::new::<kilogram_per_mole>(self.molar_mass)
    }

    /// Molar mass in kg/mol, for kernels.
    #[inline]
    pub(crate) fn molar_mass_si(&self) -> f64 {
        self.molar_mass
    }

    /// Atoms per molecule, e.g. CO2 → `[(C, 1.0), (O, 2.0)]`.
    pub fn elements(&self) -> &'static [(Element, f64)] {
        self.elements
    }

    /// Standard enthalpy of formation at 298.15 K.
    pub fn enthalpy_of_formation(&self) -> MolarEnergy {
        MolarEnergy::new::<joule_per_mole>(self.enthalpy_of_formation)
    }

    /// Temperature range covered by the fit. Outside it, values are
    /// extrapolated (constant cp below, last polynomial above); callers
    /// validate against this range.
    pub fn temperature_range(&self) -> (ThermodynamicTemperature, ThermodynamicTemperature) {
        let (t_min, t_max) = self.range_kelvin();
        (
            ThermodynamicTemperature::new::<kelvin>(t_min),
            ThermodynamicTemperature::new::<kelvin>(t_max),
        )
    }

    /// Fit range in kelvin, for kernels and compile-time set ranges.
    #[inline]
    pub(crate) const fn range_kelvin(&self) -> (f64, f64) {
        (
            self.intervals[0].t_min,
            self.intervals[self.intervals.len() - 1].t_max,
        )
    }

    /// cp/R at the given temperature.
    #[inline]
    pub fn cp_over_r(&self, p: &TemperaturePowers) -> f64 {
        if let Some((first, at_min)) = self.below_range(p) {
            nasa9::cp_over_r(&first.a, &at_min)
        } else {
            nasa9::cp_over_r(&self.interval(p.t).a, p)
        }
    }

    /// h/(RT) at the given temperature; h includes the enthalpy of formation.
    #[inline]
    pub fn h_over_rt(&self, p: &TemperaturePowers) -> f64 {
        if let Some((first, at_min)) = self.below_range(p) {
            let h_min = nasa9::h_over_rt(&first.a, &first.b, &at_min) * first.t_min;
            let cp = nasa9::cp_over_r(&first.a, &at_min);
            (h_min + cp * (p.t - first.t_min)) * p.inv_t
        } else {
            let interval = self.interval(p.t);
            nasa9::h_over_rt(&interval.a, &interval.b, p)
        }
    }

    /// s°/R at the given temperature and the standard pressure (1 bar).
    #[inline]
    pub fn s_over_r(&self, p: &TemperaturePowers) -> f64 {
        if let Some((first, at_min)) = self.below_range(p) {
            let s_min = nasa9::s_over_r(&first.a, &first.b, &at_min);
            let cp = nasa9::cp_over_r(&first.a, &at_min);
            s_min + cp * (p.ln_t - at_min.ln_t)
        } else {
            let interval = self.interval(p.t);
            nasa9::s_over_r(&interval.a, &interval.b, p)
        }
    }

    fn first(&self) -> &Interval {
        &self.intervals[0]
    }

    /// The interval containing `t`; the last one above the range (polynomial
    /// extrapolation).
    #[inline]
    fn interval(&self, t: f64) -> &Interval {
        let last = self.intervals.len() - 1;
        let mut index = 0;
        while index < last && t > self.intervals[index].t_max {
            index += 1;
        }
        &self.intervals[index]
    }

    /// Below the fit range cp is held constant at its `t_min` value: returns
    /// the first interval and the powers at `t_min`.
    #[inline]
    fn below_range(&self, p: &TemperaturePowers) -> Option<(&Interval, TemperaturePowers)> {
        let first = self.first();
        (p.t < first.t_min).then(|| (first, TemperaturePowers::from_kelvin(first.t_min)))
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::thermo::{MOLAR_GAS_CONSTANT as R, gases};

    fn powers(t: f64) -> TemperaturePowers {
        TemperaturePowers::from_kelvin(t)
    }

    #[test]
    fn test_gases_table_parses_n2_record() -> anyhow::Result<()> {
        let n2 = gases::N2;
        anyhow::ensure!(n2.name() == "N2", "name {}", n2.name());
        anyhow::ensure!(
            (n2.molar_mass - 0.028_013_4).abs() < 1e-12,
            "M {}",
            n2.molar_mass
        );
        anyhow::ensure!(n2.intervals.len() == 3, "{} intervals", n2.intervals.len());
        let low = n2.intervals[0];
        anyhow::ensure!(
            (low.t_min, low.t_max) == (200.0, 1000.0),
            "range {}..{}",
            low.t_min,
            low.t_max
        );
        anyhow::ensure!(
            (low.a[0] - 2.210_371_497e4).abs() < 1e-8,
            "a1 = {}",
            low.a[0]
        );
        anyhow::ensure!(
            (low.a[6] - 2.519_705_809e-12).abs() < 1e-24,
            "a7 = {}",
            low.a[6]
        );
        anyhow::ensure!(
            (low.b[0] - 7.108_460_860e2).abs() < 1e-9,
            "b1 = {}",
            low.b[0]
        );
        anyhow::ensure!(
            (low.b[1] + 1.076_003_744e1).abs() < 1e-11,
            "b2 = {}",
            low.b[1]
        );
        Ok(())
    }

    #[test]
    fn test_gases_all_has_thirteen_unique_species() -> anyhow::Result<()> {
        let mut names: Vec<&str> = gases::ALL.iter().map(|s| s.name()).collect();
        names.sort_unstable();
        names.dedup();
        anyhow::ensure!(names.len() == 13, "species {names:?}");
        Ok(())
    }

    #[test]
    fn test_species_elements_match_formula() -> anyhow::Result<()> {
        let co2 = gases::CO2.elements();
        anyhow::ensure!(co2 == [(Element::C, 1.0), (Element::O, 2.0)], "CO2 {co2:?}");
        let ar = gases::AR.elements();
        anyhow::ensure!(ar == [(Element::Ar, 1.0)], "Ar {ar:?}");
        let h2o = gases::H2O.elements();
        anyhow::ensure!(h2o == [(Element::H, 2.0), (Element::O, 1.0)], "H2O {h2o:?}");
        Ok(())
    }

    #[test]
    fn test_species_reports_uom_quantities() -> anyhow::Result<()> {
        let m = gases::H2O.molar_mass().get::<kilogram_per_mole>();
        anyhow::ensure!((m - 0.018_015_28).abs() < 1e-12, "M = {m}");
        let hf = gases::H2O.enthalpy_of_formation().get::<joule_per_mole>();
        anyhow::ensure!((hf + 241_826.0).abs() < 1e-9, "Hf = {hf}");
        Ok(())
    }

    #[test]
    fn test_species_temperature_range() -> anyhow::Result<()> {
        for (species, lo, hi) in [(gases::H2O, 200.0, 6000.0), (gases::N2, 200.0, 20_000.0)] {
            let (t_min, t_max) = species.temperature_range();
            let (t_min, t_max) = (t_min.get::<kelvin>(), t_max.get::<kelvin>());
            anyhow::ensure!(
                (t_min - lo).abs() < 1e-9 && (t_max - hi).abs() < 1e-9,
                "{}: {t_min}..{t_max}",
                species.name()
            );
        }
        Ok(())
    }

    #[test]
    fn test_species_enthalpy_at_298_equals_heat_of_formation() -> anyhow::Result<()> {
        let t = 298.15;
        for species in gases::ALL {
            let h = species.h_over_rt(&powers(t)) * R * t;
            let hf = species.enthalpy_of_formation;
            anyhow::ensure!(
                (h - hf).abs() < 1e-5 * hf.abs() + 1.0,
                "{}: h(298.15) = {h} J/mol, Hf = {hf}",
                species.name()
            );
        }
        Ok(())
    }

    #[test]
    fn test_species_argon_cp_is_five_halves_r() -> anyhow::Result<()> {
        // Monatomic translation only; electronic excitation adds ~1e-6 by 6000
        // K.
        for t in [200.0, 298.15, 1000.0, 3000.0, 6000.0] {
            let cp = gases::AR.cp_over_r(&powers(t));
            anyhow::ensure!((cp - 2.5).abs() < 1e-5, "cp/R({t}) = {cp}");
        }
        Ok(())
    }

    /// JANAF (Chase 1998): cp and S° at 298.15 K in J/(mol·K).
    #[test]
    fn test_species_matches_janaf_at_298() -> anyhow::Result<()> {
        let reference = [
            (gases::N2, 29.124, 191.609),
            (gases::O2, 29.376, 205.147),
            (gases::AR, 20.786, 154.845),
            (gases::H2O, 33.587, 188.834),
            (gases::CO2, 37.135, 213.795),
            (gases::CO, 29.142, 197.653),
            (gases::H2, 28.836, 130.680),
            (gases::OH, 29.886, 183.737),
            (gases::H, 20.786, 114.716),
            (gases::O, 21.911, 161.058),
            (gases::NO, 29.845, 210.758),
            (gases::N, 20.786, 153.300),
        ];
        let p = powers(298.15);
        for (species, cp_ref, s_ref) in reference {
            let cp = species.cp_over_r(&p) * R;
            let s = species.s_over_r(&p) * R;
            anyhow::ensure!(
                (cp / cp_ref - 1.0).abs() < 2e-3 && (s - s_ref).abs() < 0.1,
                "{}: cp {cp} vs {cp_ref}, S {s} vs {s_ref}",
                species.name()
            );
        }
        Ok(())
    }

    /// JANAF (Chase 1998): cp at 1000 K in J/(mol·K).
    #[test]
    fn test_species_matches_janaf_cp_at_1000() -> anyhow::Result<()> {
        let p = powers(1000.0);
        for (species, cp_ref) in [
            (gases::O2, 34.870),
            (gases::H2O, 41.268),
            (gases::CO2, 54.308),
            (gases::H2, 30.205),
        ] {
            let cp = species.cp_over_r(&p) * R;
            anyhow::ensure!(
                (cp / cp_ref - 1.0).abs() < 2e-3,
                "{}: cp {cp} vs {cp_ref}",
                species.name()
            );
        }
        Ok(())
    }

    #[test]
    fn test_species_extrapolates_constant_cp_below_range() -> anyhow::Result<()> {
        for species in gases::ALL {
            let t_min = species.intervals[0].t_min;
            let at_min = powers(t_min);
            let below = powers(t_min - 100.0);
            let cp_min = species.cp_over_r(&at_min);
            anyhow::ensure!(
                (species.cp_over_r(&below) - cp_min).abs() < 1e-12,
                "{}: cp not constant below range",
                species.name()
            );
            // h continuous at t_min, then linear: H(T) = H(t_min) + cp·(T −
            // t_min).
            let h_min = species.h_over_rt(&at_min) * t_min;
            let h_below = species.h_over_rt(&below) * (t_min - 100.0);
            let expected = h_min - cp_min * 100.0;
            anyhow::ensure!(
                (h_below - expected).abs() < 1e-9 * h_min.abs().max(1.0),
                "{}: H/R below range {h_below}, expected {expected}",
                species.name()
            );
        }
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_species_is_continuous_at_interval_boundaries(index in 0..12usize) {
            let species = gases::ALL[index];
            for pair in species.intervals.windows(2) {
                let t = pair[0].t_max;
                let (lo, hi) = (powers(t * (1.0 - 1e-12)), powers(t * (1.0 + 1e-12)));
                for (what, a, b) in [
                    ("cp", species.cp_over_r(&lo), species.cp_over_r(&hi)),
                    ("h", species.h_over_rt(&lo), species.h_over_rt(&hi)),
                    ("s", species.s_over_r(&lo), species.s_over_r(&hi)),
                ] {
                    prop_assert!((a - b).abs() <= 1e-4 * a.abs().max(1.0),
                        "{} {what} jumps at {t}: {a} → {b}", species.name());
                }
            }
        }

        #[test]
        fn test_species_enthalpy_and_entropy_derivatives_match_cp(
            index in 0..12usize, t in 100.0..6000.0f64,
        ) {
            let species = gases::ALL[index];
            let dt = 1e-3 * t;
            // Stay inside one interval so the finite difference sees one polynomial.
            let crosses = species.intervals.iter().any(|i| (i.t_max - t).abs() < 2.0 * dt)
                || (species.intervals[0].t_min - t).abs() < 2.0 * dt;
            prop_assume!(!crosses);
            let (lo, mid, hi) = (powers(t - dt), powers(t), powers(t + dt));
            let cp = species.cp_over_r(&mid);
            let dh_dt = (species.h_over_rt(&hi) * (t + dt) - species.h_over_rt(&lo) * (t - dt))
                / (2.0 * dt);
            let ds_dt = (species.s_over_r(&hi) - species.s_over_r(&lo)) / (2.0 * dt);
            prop_assert!((dh_dt - cp).abs() < 1e-5 * cp, "{}: dh/dT {dh_dt} vs cp {cp}",
                species.name());
            prop_assert!((ds_dt * t - cp).abs() < 1e-5 * cp, "{}: T ds/dT {} vs cp {cp}",
                species.name(), ds_dt * t);
        }
    }
}
