//! A well-mixed, ideal-gas control volume.

use anyhow::Context;

use super::GasPort;
use crate::{
    PhysicsError,
    flow::{PerfectGas, Reservoir},
    thermo::{Mixture, SpeciesSet},
    units::{
        AvailableEnergy, Energy, Mass, Pressure, SpecificHeatCapacity, ThermodynamicTemperature,
        Volume, available_energy::joule_per_kilogram, energy::joule, mass::kilogram,
        pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
        thermodynamic_temperature::kelvin, volume::cubic_meter,
    },
};

/// Gas in a volume, tracked by its conserved quantities.
///
/// The ODE state of a model holds the species masses and the internal energy;
/// after every state change call [`set_conserved`](Self::set_conserved), which
/// recomputes temperature (warm-started from the previous one), pressure and
/// the cached properties.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ControlVolume<S: SpeciesSet<N>, const N: usize> {
    /// kg per species, in set order.
    masses: [f64; N],
    /// J.
    internal_energy: f64,
    /// m³.
    volume: f64,
    // Derived by `refresh`.
    mixture: Mixture<S, N>,
    mass: f64,
    temperature: f64,
    pressure: f64,
    specific_enthalpy: f64,
    gas: PerfectGas,
}

impl<S: SpeciesSet<N>, const N: usize> ControlVolume<S, N> {
    /// Gas of composition `mixture` at pressure `p` and temperature `t` filling
    /// `volume`.
    pub fn new(
        mixture: Mixture<S, N>,
        p: Pressure,
        t: ThermodynamicTemperature,
        volume: Volume,
    ) -> anyhow::Result<Self> {
        let (pressure, temperature) = (p.get::<pascal>(), t.get::<kelvin>());
        let volume = volume.get::<cubic_meter>();
        if !(pressure > 0.0 && pressure.is_finite()) {
            return Err(PhysicsError::InvalidState("pressure must be positive"))
                .with_context(|| format!("p = {pressure} Pa"));
        }
        if !(volume > 0.0 && volume.is_finite()) {
            return Err(PhysicsError::InvalidGeometry("volume must be positive"))
                .with_context(|| format!("V = {volume} m³"));
        }
        check_temperature::<S, N>(temperature).context("control volume")?;
        let mass = pressure * volume / (mixture.gas_constant_si() * temperature);
        let mut cv = Self {
            masses: mixture.mass_fractions().map(|y| y * mass),
            internal_energy: mass * mixture.internal_energy_si(temperature),
            volume,
            mixture,
            mass,
            temperature,
            pressure,
            specific_enthalpy: 0.0,
            gas: PerfectGas {
                gamma: 1.4,
                gas_constant: 1.0,
            },
        };
        cv.refresh();
        Ok(cv)
    }

    /// Replaces the conserved state: species masses in kg (set order), internal
    /// energy in J and volume in m³. Kernel: never fails; check the result with
    /// [`validate`](Self::validate) once per solver step.
    #[inline]
    pub fn set_conserved(&mut self, masses: &[f64; N], internal_energy: f64, volume: f64) {
        self.masses = *masses;
        self.internal_energy = internal_energy;
        self.volume = volume;
        self.refresh();
    }

    /// Recomputes the derived state, warm-starting the temperature inversion
    /// from the previous temperature.
    #[inline]
    fn refresh(&mut self) {
        self.mass = self.masses.iter().sum();
        self.mixture = Mixture::from_masses_unchecked(&self.masses);
        let r = self.mixture.gas_constant_si();
        self.temperature = self.mixture.temperature_from_internal_energy_si(
            self.internal_energy / self.mass,
            self.temperature,
        );
        self.pressure = self.mass * r * self.temperature / self.volume;
        let (h, cp) = self.mixture.enthalpy_and_cp_si(self.temperature);
        self.specific_enthalpy = h;
        self.gas = PerfectGas {
            gamma: cp / (cp - r),
            gas_constant: r,
        };
    }

    /// Checks that the derived state is physical: finite, non-negative masses
    /// with a positive total, positive volume and pressure, temperature inside
    /// the species data range.
    ///
    /// Species masses down to −[`NEGATIVE_MASS_TOLERANCE`]·total pass:
    /// high-order Runge–Kutta weights include negative ones, so a species that
    /// is absent on one side of a reversing flow can undershoot zero by
    /// integration-error-sized amounts. Clamping them would break exact
    /// conservation; rejecting them would stop valid simulations.
    pub fn validate(&self) -> anyhow::Result<()> {
        let floor = -NEGATIVE_MASS_TOLERANCE * self.mass;
        if self.masses.iter().any(|m| !(m.is_finite() && *m >= floor)) {
            return Err(PhysicsError::InvalidState(
                "species masses must be finite and non-negative",
            ))
            .with_context(|| format!("masses {:?} kg", self.masses));
        }
        if self.mass.is_nan() || self.mass <= 0.0 {
            return Err(PhysicsError::InvalidState("control volume holds no gas").into());
        }
        if !self.internal_energy.is_finite() {
            return Err(PhysicsError::InvalidState("internal energy must be finite"))
                .with_context(|| format!("U = {} J", self.internal_energy));
        }
        if !(self.volume > 0.0 && self.volume.is_finite()) {
            return Err(PhysicsError::InvalidGeometry("volume must be positive"))
                .with_context(|| format!("V = {} m³", self.volume));
        }
        check_temperature::<S, N>(self.temperature).context("control volume")?;
        if !(self.pressure > 0.0 && self.pressure.is_finite()) {
            return Err(PhysicsError::InvalidState("pressure must be positive"))
                .with_context(|| format!("p = {} Pa", self.pressure));
        }
        Ok(())
    }

    /// Species masses in kg, in set order (the ODE state).
    pub fn species_masses(&self) -> &[f64; N] {
        &self.masses
    }

    /// Total internal energy (the ODE state).
    pub fn internal_energy(&self) -> Energy {
        Energy::new::<joule>(self.internal_energy)
    }

    /// Volume.
    pub fn volume(&self) -> Volume {
        Volume::new::<cubic_meter>(self.volume)
    }

    /// Total mass.
    pub fn mass(&self) -> Mass {
        Mass::new::<kilogram>(self.mass)
    }

    /// Current composition.
    pub fn mixture(&self) -> &Mixture<S, N> {
        &self.mixture
    }

    /// Temperature.
    pub fn temperature(&self) -> ThermodynamicTemperature {
        ThermodynamicTemperature::new::<kelvin>(self.temperature)
    }

    /// Pressure.
    pub fn pressure(&self) -> Pressure {
        Pressure::new::<pascal>(self.pressure)
    }

    /// Specific entropy of the contents.
    pub fn specific_entropy(&self) -> SpecificHeatCapacity {
        SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(
            self.mixture.entropy_si(self.temperature, self.pressure),
        )
    }

    /// Rate of internal-energy change from boundary work, −p·dV/dt, in W, for a
    /// volume changing at `volume_rate` m³/s.
    #[inline]
    pub fn boundary_work_rate(&self, volume_rate: f64) -> f64 {
        -self.pressure * volume_rate
    }
}

/// Relative undershoot below zero a species mass may have before
/// [`ControlVolume::validate`] rejects it. The undershoot scales with the
/// integrator's local error, not machine epsilon (1e-9 relative at rtol 1e-9),
/// and real-time fixed-step loops run loose; 1 ppm of the volume's mass is
/// physically negligible yet well above typical integration error.
pub const NEGATIVE_MASS_TOLERANCE: f64 = 1e-6;

/// Fails unless `temperature` (K) lies inside the species set's data range.
pub(super) fn check_temperature<S: SpeciesSet<N>, const N: usize>(
    temperature: f64,
) -> anyhow::Result<()> {
    let (min, max) = S::TEMPERATURE_RANGE;
    if (min..=max).contains(&temperature) {
        Ok(())
    } else {
        Err(PhysicsError::TemperatureOutOfRange {
            temperature,
            min,
            max,
        }
        .into())
    }
}

impl<S: SpeciesSet<N>, const N: usize> GasPort<S, N> for ControlVolume<S, N> {
    fn reservoir(&self) -> Reservoir {
        Reservoir {
            pressure: self.pressure,
            temperature: self.temperature,
            gas: self.gas,
        }
    }

    fn mass_fractions(&self) -> &[f64; N] {
        self.mixture.mass_fractions()
    }

    fn specific_enthalpy(&self) -> AvailableEnergy {
        AvailableEnergy::new::<joule_per_kilogram>(self.specific_enthalpy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermo::{Air3, MOLAR_GAS_CONSTANT};

    fn tank(p: f64, t: f64, v: f64) -> anyhow::Result<ControlVolume<Air3, 3>> {
        ControlVolume::new(
            Mixture::dry_air()?,
            Pressure::new::<pascal>(p),
            ThermodynamicTemperature::new::<kelvin>(t),
            Volume::new::<cubic_meter>(v),
        )
    }

    fn physics_error(result: anyhow::Result<()>) -> Option<PhysicsError> {
        result.err()?.downcast_ref::<PhysicsError>().copied()
    }

    #[test]
    fn test_control_volume_new_reproduces_state() -> anyhow::Result<()> {
        let cv = tank(5e5, 350.0, 0.02)?;
        let p = cv.pressure().get::<pascal>();
        let t = cv.temperature().get::<kelvin>();
        anyhow::ensure!(
            (p / 5e5 - 1.0).abs() < 1e-12 && (t / 350.0 - 1.0).abs() < 1e-12,
            "{p} {t}"
        );
        let r = cv
            .mixture()
            .gas_constant()
            .get::<joule_per_kilogram_kelvin>();
        let m = cv.mass().get::<kilogram>();
        anyhow::ensure!(
            (m / (5e5 * 0.02 / (r * 350.0)) - 1.0).abs() < 1e-12,
            "m = {m}"
        );
        let sum: f64 = cv.species_masses().iter().sum();
        anyhow::ensure!((sum / m - 1.0).abs() < 1e-15, "species sum {sum}");
        cv.validate()
    }

    #[test]
    fn test_control_volume_new_rejects_invalid_state() -> anyhow::Result<()> {
        for (p, t, v) in [
            (0.0, 300.0, 1.0),
            (1e5, -1.0, 1.0),
            (1e5, 300.0, 0.0),
            (1e5, 150.0, 1.0),
        ] {
            anyhow::ensure!(tank(p, t, v).is_err(), "accepted p {p}, T {t}, V {v}");
        }
        Ok(())
    }

    #[test]
    fn test_control_volume_set_conserved_round_trips() -> anyhow::Result<()> {
        let original = tank(2e5, 800.0, 0.5)?;
        let mut copy = tank(1e5, 300.0, 0.1)?;
        let u = original.internal_energy().get::<joule>();
        copy.set_conserved(original.species_masses(), u, 0.5);
        let (p, t) = (
            copy.pressure().get::<pascal>(),
            copy.temperature().get::<kelvin>(),
        );
        anyhow::ensure!(
            (p / 2e5 - 1.0).abs() < 1e-10 && (t / 800.0 - 1.0).abs() < 1e-10,
            "{p} {t}"
        );
        copy.validate()
    }

    #[test]
    fn test_control_volume_validate_flags_bad_states() -> anyhow::Result<()> {
        let mut cv = tank(1e5, 300.0, 0.1)?;
        let masses = *cv.species_masses();
        let u = cv.internal_energy().get::<joule>();

        // Cooled far below the 200 K data range.
        let cold_u = u - cv.mass().get::<kilogram>() * 1000.0 * 200.0;
        cv.set_conserved(&masses, cold_u, 0.1);
        anyhow::ensure!(
            matches!(
                physics_error(cv.validate()),
                Some(PhysicsError::TemperatureOutOfRange { .. })
            ),
            "cold state accepted at T = {:?}",
            cv.temperature()
        );

        let total = cv.mass().get::<kilogram>();
        let mut roundoff = masses;
        roundoff[2] = -1e-8 * total;
        cv.set_conserved(&roundoff, u, 0.1);
        cv.validate()
            .context("roundoff-sized undershoot rejected")?;

        let mut negative = masses;
        negative[1] = -1e-3;
        cv.set_conserved(&negative, u, 0.1);
        anyhow::ensure!(cv.validate().is_err(), "negative species mass accepted");

        cv.set_conserved(&masses, f64::NAN, 0.1);
        anyhow::ensure!(cv.validate().is_err(), "NaN energy accepted");

        cv.set_conserved(&[0.0; 3], 0.0, 0.1);
        anyhow::ensure!(cv.validate().is_err(), "empty volume accepted");
        Ok(())
    }

    #[test]
    fn test_control_volume_boundary_work_rate() -> anyhow::Result<()> {
        let cv = tank(3e5, 300.0, 0.1)?;
        // Compressing at 1 L/s against 3 bar adds 300 W.
        anyhow::ensure!((cv.boundary_work_rate(-1e-3) - 300.0).abs() < 1e-9, "work");
        Ok(())
    }

    #[test]
    fn test_control_volume_entropy_matches_mixture() -> anyhow::Result<()> {
        let cv = tank(3e5, 450.0, 0.1)?;
        let want = cv.mixture().entropy(cv.temperature(), cv.pressure());
        anyhow::ensure!(cv.specific_entropy() == want, "entropy differs");
        let _ = MOLAR_GAS_CONSTANT;
        Ok(())
    }
}
