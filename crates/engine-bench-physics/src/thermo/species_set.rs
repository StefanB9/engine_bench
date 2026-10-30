//! Compile-time species sets: each model picks the species it tracks, and
//! mixtures are fixed-size arrays over that set.

use super::{Species, gases};

/// A fixed, ordered set of `N` species. Mixture arrays are indexed in this
/// order.
pub trait SpeciesSet<const N: usize>: Copy + 'static {
    /// The species, in array order.
    const SPECIES: [&'static Species; N];

    /// Temperature range in kelvin where every species has fit data: the
    /// highest lower bound and the lowest upper bound. Evaluated at compile
    /// time.
    const TEMPERATURE_RANGE: (f64, f64) = intersect_ranges(&Self::SPECIES);

    /// [`TEMPERATURE_RANGE`](Self::TEMPERATURE_RANGE) as a function.
    fn temperature_range_kelvin() -> (f64, f64) {
        Self::TEMPERATURE_RANGE
    }

    /// Index of `species` in this set, if present.
    fn index_of(species: &Species) -> Option<usize> {
        Self::SPECIES
            .iter()
            .position(|s| s.name() == species.name())
    }
}

const fn intersect_ranges(species: &[&Species]) -> (f64, f64) {
    let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
    let mut i = 0;
    while i < species.len() {
        let (t_min, t_max) = species[i].range_kelvin();
        lo = lo.max(t_min);
        hi = hi.min(t_max);
        i += 1;
    }
    (lo, hi)
}

/// Dry air without CO₂: N₂, O₂, Ar.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Air3;

impl SpeciesSet<3> for Air3 {
    const SPECIES: [&'static Species; 3] = [&gases::N2, &gases::O2, &gases::AR];
}

/// Air plus major hydrocarbon/hydrogen combustion products: N₂, O₂, Ar, H₂O,
/// CO₂, CO, H₂.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Major7;

impl SpeciesSet<7> for Major7 {
    const SPECIES: [&'static Species; 7] = [
        &gases::N2,
        &gases::O2,
        &gases::AR,
        &gases::H2O,
        &gases::CO2,
        &gases::CO,
        &gases::H2,
    ];
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_species_set_temperature_range_is_intersection() -> anyhow::Result<()> {
        let air = Air3::temperature_range_kelvin();
        anyhow::ensure!(air == (200.0, 20_000.0), "Air3 {air:?}");
        // H2O data stops at 6000 K.
        let major = Major7::temperature_range_kelvin();
        anyhow::ensure!(major == (200.0, 6000.0), "Major7 {major:?}");
        Ok(())
    }

    #[test]
    fn test_species_set_index_of() -> anyhow::Result<()> {
        anyhow::ensure!(Major7::index_of(&gases::CO2) == Some(4), "CO2 index");
        anyhow::ensure!(Air3::index_of(&gases::AR) == Some(2), "Ar index");
        anyhow::ensure!(Air3::index_of(&gases::H2O).is_none(), "H2O in Air3");
        Ok(())
    }
}
