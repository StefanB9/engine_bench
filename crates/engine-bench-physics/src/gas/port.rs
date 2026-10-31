//! Gas ports and orifice transfers between them.

use anyhow::Context;

use super::control_volume::check_temperature;
use crate::{
    PhysicsError,
    flow::{Orifice, PerfectGas, Reservoir},
    thermo::{Mixture, SpeciesSet},
    units::{
        AvailableEnergy, Pressure, ThermodynamicTemperature, available_energy::joule_per_kilogram,
        pressure::pascal, thermodynamic_temperature::kelvin,
    },
};

/// Something gas can flow out of or into through an orifice.
pub trait GasPort<S: SpeciesSet<N>, const N: usize> {
    /// Stagnation state and frozen perfect-gas properties.
    fn reservoir(&self) -> Reservoir;
    /// Composition of gas leaving through the port.
    fn mass_fractions(&self) -> &[f64; N];
    /// Specific (stagnation) enthalpy carried by gas leaving through the port.
    fn specific_enthalpy(&self) -> AvailableEnergy;
}

/// A fixed reservoir too large to change: the atmosphere, a test-stand supply,
/// vacuum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Environment<S: SpeciesSet<N>, const N: usize> {
    mixture: Mixture<S, N>,
    reservoir: Reservoir,
    /// J/kg.
    specific_enthalpy: f64,
}

impl<S: SpeciesSet<N>, const N: usize> Environment<S, N> {
    /// Gas of composition `mixture` at pressure `p` (≥ 0; 0 is vacuum) and
    /// temperature `t`.
    pub fn new(
        mixture: Mixture<S, N>,
        p: Pressure,
        t: ThermodynamicTemperature,
    ) -> anyhow::Result<Self> {
        let (pressure, temperature) = (p.get::<pascal>(), t.get::<kelvin>());
        if !(pressure >= 0.0 && pressure.is_finite()) {
            return Err(PhysicsError::InvalidState("pressure must be non-negative"))
                .with_context(|| format!("environment p = {pressure} Pa"));
        }
        check_temperature::<S, N>(temperature).context("environment")?;
        let r = mixture.gas_constant_si();
        let (specific_enthalpy, cp) = mixture.enthalpy_and_cp_si(temperature);
        Ok(Self {
            mixture,
            reservoir: Reservoir {
                pressure,
                temperature,
                gas: PerfectGas {
                    gamma: cp / (cp - r),
                    gas_constant: r,
                },
            },
            specific_enthalpy,
        })
    }
}

impl<S: SpeciesSet<N>, const N: usize> Environment<S, N> {
    /// Composition of the surroundings.
    pub(super) fn mixture(&self) -> &Mixture<S, N> {
        &self.mixture
    }
}

impl<S: SpeciesSet<N>, const N: usize> GasPort<S, N> for Environment<S, N> {
    fn reservoir(&self) -> Reservoir {
        self.reservoir
    }

    fn mass_fractions(&self) -> &[f64; N] {
        self.mixture.mass_fractions()
    }

    fn specific_enthalpy(&self) -> AvailableEnergy {
        AvailableEnergy::new::<joule_per_kilogram>(self.specific_enthalpy)
    }
}

/// What flows from port `a` to port `b` (negative values flow from `b` to `a`).
/// Subtract it from `a`'s rates and add it to `b`'s.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transfer<const N: usize> {
    /// Mass flow per species in kg/s, in set order.
    pub species: [f64; N],
    /// Enthalpy flow in W.
    pub energy: f64,
}

impl<const N: usize> Transfer<N> {
    /// Total mass flow in kg/s.
    pub fn mass_flow(&self) -> f64 {
        self.species.iter().sum()
    }
}

/// Flow through `orifice` from `a` to `b`, carrying the upstream port's
/// composition and enthalpy. Kernel: never fails.
#[inline]
pub fn orifice_transfer<S, const N: usize>(
    orifice: Orifice,
    a: &impl GasPort<S, N>,
    b: &impl GasPort<S, N>,
) -> Transfer<N>
where
    S: SpeciesSet<N>,
{
    let mass_flow = orifice.mass_flow_si(&a.reservoir(), &b.reservoir());
    let (fractions, enthalpy) = if mass_flow >= 0.0 {
        (a.mass_fractions(), a.specific_enthalpy())
    } else {
        (b.mass_fractions(), b.specific_enthalpy())
    };
    Transfer {
        species: fractions.map(|y| mass_flow * y),
        energy: mass_flow * enthalpy.get::<joule_per_kilogram>(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        gas::ControlVolume,
        thermo::Air3,
        units::{Area, Volume, area::square_meter, volume::cubic_meter},
    };

    fn tank(p: f64, t: f64) -> anyhow::Result<ControlVolume<Air3, 3>> {
        ControlVolume::new(
            Mixture::dry_air()?,
            Pressure::new::<pascal>(p),
            ThermodynamicTemperature::new::<kelvin>(t),
            Volume::new::<cubic_meter>(0.01),
        )
    }

    fn orifice() -> anyhow::Result<Orifice> {
        Orifice::new(Area::new::<square_meter>(1e-5), 0.8)
    }

    #[test]
    fn test_environment_allows_vacuum_but_not_negative_pressure() -> anyhow::Result<()> {
        let air = Mixture::dry_air()?;
        let t = ThermodynamicTemperature::new::<kelvin>(288.15);
        anyhow::ensure!(
            Environment::new(air, Pressure::new::<pascal>(0.0), t).is_ok(),
            "vacuum"
        );
        for p in [-1.0, f64::NAN] {
            anyhow::ensure!(
                Environment::new(air, Pressure::new::<pascal>(p), t).is_err(),
                "accepted p {p}"
            );
        }
        anyhow::ensure!(
            Environment::new(
                air,
                Pressure::new::<pascal>(1e5),
                ThermodynamicTemperature::new::<kelvin>(0.0)
            )
            .is_err(),
            "accepted T = 0"
        );
        Ok(())
    }

    #[test]
    fn test_orifice_transfer_carries_upstream_enthalpy_and_composition() -> anyhow::Result<()> {
        let (hot, cold) = (tank(5e5, 600.0)?, tank(1e5, 300.0)?);
        let t = orifice_transfer(orifice()?, &hot, &cold);
        let mdot = orifice()?.mass_flow(&hot.reservoir(), &cold.reservoir());
        let mdot = mdot.get::<crate::units::mass_rate::kilogram_per_second>();
        anyhow::ensure!(
            (t.mass_flow() / mdot - 1.0).abs() < 1e-14,
            "ṁ {} vs {mdot}",
            t.mass_flow()
        );
        let h = hot.specific_enthalpy().get::<joule_per_kilogram>();
        anyhow::ensure!(
            (t.energy / (mdot * h) - 1.0).abs() < 1e-14,
            "energy {}",
            t.energy
        );
        for (species, y) in t.species.iter().zip(hot.mass_fractions()) {
            anyhow::ensure!((species - mdot * y).abs() < 1e-15, "species flow {species}");
        }
        Ok(())
    }

    #[test]
    fn test_orifice_transfer_is_antisymmetric() -> anyhow::Result<()> {
        let (a, b) = (tank(5e5, 600.0)?, tank(1e5, 300.0)?);
        let forward = orifice_transfer(orifice()?, &a, &b);
        let backward = orifice_transfer(orifice()?, &b, &a);
        // Exact negation, so coupled volumes conserve to roundoff.
        anyhow::ensure!(
            forward.energy.to_bits() == (-backward.energy).to_bits(),
            "energy"
        );
        anyhow::ensure!(
            forward
                .species
                .iter()
                .zip(backward.species)
                .all(|(f, b)| f.to_bits() == (-b).to_bits()),
            "species"
        );
        Ok(())
    }

    #[test]
    fn test_orifice_transfer_into_vacuum_is_choked() -> anyhow::Result<()> {
        let air = Mixture::dry_air()?;
        let vacuum = Environment::new(
            air,
            Pressure::new::<pascal>(0.0),
            ThermodynamicTemperature::new::<kelvin>(300.0),
        )?;
        let tank = tank(5e5, 300.0)?;
        let t = orifice_transfer(orifice()?, &tank, &vacuum);
        anyhow::ensure!(
            t.mass_flow() > 0.0 && t.mass_flow().is_finite(),
            "ṁ {}",
            t.mass_flow()
        );
        let back = orifice_transfer(orifice()?, &vacuum, &tank);
        anyhow::ensure!(
            (back.mass_flow() + t.mass_flow()).abs() < 1e-18,
            "reverse {}",
            back.mass_flow()
        );
        Ok(())
    }
}
