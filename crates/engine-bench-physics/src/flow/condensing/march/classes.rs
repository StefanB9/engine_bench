//! Droplet classes: the droplets nucleated within one march step share a
//! birth place and so, in steady flow, a size ever after. Each class keeps
//! its number per unit mass (fixed) and its liquid mass fraction (grown by
//! the march); the radius follows from the two.

use std::f64::consts::PI;

/// Classes kept at most; beyond, the two most alike are merged.
pub(super) const CLASSES: usize = 64;

/// The closed classes, in birth order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Classes {
    /// Droplets per kilogram of mixture.
    numbers: [f64; CLASSES],
    count: usize,
}

impl Default for Classes {
    fn default() -> Self {
        Self {
            numbers: [0.0; CLASSES],
            count: 0,
        }
    }
}

/// Radius, m, of droplets of liquid density `rho_liquid` when `number` of
/// them per kilogram hold liquid fraction `y`; zero when empty.
pub(super) fn radius(number: f64, y: f64, rho_liquid: f64) -> f64 {
    if number > 0.0 && y > 0.0 {
        (3.0 * y / (4.0 * PI * rho_liquid * number)).cbrt()
    } else {
        0.0
    }
}

impl Classes {
    pub(super) fn count(&self) -> usize {
        self.count
    }

    /// Droplets per kilogram in class `k`.
    pub(super) fn number(&self, k: usize) -> f64 {
        self.numbers.get(k).copied().unwrap_or(0.0)
    }

    /// Droplets per kilogram in all classes.
    pub(super) fn total(&self) -> f64 {
        self.numbers.iter().take(self.count).sum()
    }

    /// Adds a class of `number` droplets per kilogram holding liquid fraction
    /// `y`, whose fractions are kept in `fractions` alongside; when every slot
    /// is taken, first merges the adjacent pair of most alike radii (number
    /// and liquid summed, so both are conserved).
    pub(super) fn push(
        &mut self,
        number: f64,
        y: f64,
        fractions: &mut [f64; CLASSES],
        rho_liquid: f64,
    ) {
        if self.count == CLASSES {
            self.merge_closest(fractions, rho_liquid);
        }
        let k = self.count;
        if let (Some(n), Some(f)) = (self.numbers.get_mut(k), fractions.get_mut(k)) {
            *n = number;
            *f = y;
            self.count += 1;
        }
    }

    /// Drops the classes that have evaporated (no liquid left): their mass is
    /// back in the vapour already.
    pub(super) fn drop_empty(&mut self, fractions: &mut [f64; CLASSES]) {
        let mut kept = 0;
        for k in 0..self.count {
            let (n, f) = (self.number(k), fractions.get(k).copied().unwrap_or(0.0));
            if f > 0.0 && n > 0.0 {
                if let (Some(slot_n), Some(slot_f)) =
                    (self.numbers.get_mut(kept), fractions.get_mut(kept))
                {
                    *slot_n = n;
                    *slot_f = f;
                }
                kept += 1;
            }
        }
        for k in kept..self.count {
            if let (Some(n), Some(f)) = (self.numbers.get_mut(k), fractions.get_mut(k)) {
                *n = 0.0;
                *f = 0.0;
            }
        }
        self.count = kept;
    }

    fn merge_closest(&mut self, fractions: &mut [f64; CLASSES], rho_liquid: f64) {
        let ln_r = |k: usize| {
            radius(
                self.number(k),
                fractions.get(k).copied().unwrap_or(0.0),
                rho_liquid,
            )
            .max(f64::MIN_POSITIVE)
            .ln()
        };
        let Some(best) = (0..self.count.saturating_sub(1)).min_by(|&i, &j| {
            (ln_r(i) - ln_r(i + 1))
                .abs()
                .total_cmp(&(ln_r(j) - ln_r(j + 1)).abs())
        }) else {
            return;
        };
        let (n_next, f_next) = (
            self.number(best + 1),
            fractions.get(best + 1).copied().unwrap_or(0.0),
        );
        if let (Some(n), Some(f)) = (self.numbers.get_mut(best), fractions.get_mut(best)) {
            *n += n_next;
            *f += f_next;
        }
        for k in best + 1..self.count - 1 {
            let (n, f) = (
                self.number(k + 1),
                fractions.get(k + 1).copied().unwrap_or(0.0),
            );
            if let (Some(slot_n), Some(slot_f)) = (self.numbers.get_mut(k), fractions.get_mut(k)) {
                *slot_n = n;
                *slot_f = f;
            }
        }
        let last = self.count - 1;
        if let (Some(n), Some(f)) = (self.numbers.get_mut(last), fractions.get_mut(last)) {
            *n = 0.0;
            *f = 0.0;
        }
        self.count -= 1;
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    #[test]
    fn test_classes_drop_evaporated_droplets_keeping_order() {
        let mut c = Classes::default();
        let mut f = [0.0; CLASSES];
        for (n, y) in [(1e15, 1e-3), (2e15, 0.0), (3e15, 2e-3)] {
            c.push(n, y, &mut f, 1000.0);
        }
        c.drop_empty(&mut f);
        assert_eq!(c.count(), 2);
        assert!((c.number(1) - 3e15).abs() < 1.0 && (f[1] - 2e-3).abs() < 1e-18);
        assert!(f[2].abs() < f64::MIN_POSITIVE && c.number(2).abs() < f64::MIN_POSITIVE);
    }

    #[test]
    fn test_classes_radius_holds_the_liquid() {
        // 1e18 droplets of 10 nm water per kg hold 4.19e-3 kg/kg.
        let r = radius(1e18, 4.0 / 3.0 * PI * 1e-24 * 1000.0 * 1e18, 1000.0);
        assert!((r / 1e-8 - 1.0).abs() < 1e-12, "{r}");
        assert!(radius(0.0, 1e-3, 1000.0).abs() < f64::MIN_POSITIVE);
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        /// Merging keeps the droplet count and the liquid exactly, and joins
        /// neighbours of similar size.
        #[test]
        fn test_classes_merging_conserves_droplets_and_liquid(
            seeds in prop::collection::vec((1e10f64..1e18, 1e-9f64..1e-2), CLASSES + 1..3 * CLASSES),
        ) {
            let mut c = Classes::default();
            let mut f = [0.0; CLASSES];
            let (mut number, mut liquid) = (0.0, 0.0);
            for &(n, y) in &seeds {
                c.push(n, y, &mut f, 1000.0);
                number += n;
                liquid += y;
                prop_assert!(c.count() <= CLASSES);
            }
            let kept_liquid: f64 = f.iter().sum();
            prop_assert!((c.total() / number - 1.0).abs() < 1e-12);
            prop_assert!((kept_liquid / liquid - 1.0).abs() < 1e-12);
        }
    }
}
