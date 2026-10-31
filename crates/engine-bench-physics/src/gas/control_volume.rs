//! A well-mixed, ideal-gas control volume.

use anyhow::Context;

use super::{
    GasPort,
    eos::{EquationOfState, FluidState, IdealGas},
};
use crate::{
    PhysicsError,
    flow::Reservoir,
    thermo::{Mixture, SpeciesSet},
    units::{
        AvailableEnergy, Energy, Mass, Pressure, Ratio, SpecificHeatCapacity,
        ThermodynamicTemperature, Volume, available_energy::joule_per_kilogram, energy::joule,
        mass::kilogram, pressure::pascal, ratio::ratio,
        specific_heat_capacity::joule_per_kilogram_kelvin, thermodynamic_temperature::kelvin,
        volume::cubic_meter,
    },
};

/// Gas in a volume, tracked by its conserved quantities, with the equation of
/// state `E` (ideal gas unless chosen otherwise).
///
/// The ODE state of a model holds the species masses and the internal energy;
/// after every state change call [`set_conserved`](Self::set_conserved), which
/// recomputes temperature (warm-started from the previous one), pressure and
/// the cached properties.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ControlVolume<S: SpeciesSet<N>, const N: usize, E: EquationOfState<S, N> = IdealGas> {
    /// kg per species, in set order.
    masses: [f64; N],
    /// J.
    internal_energy: f64,
    /// m³.
    volume: f64,
    eos: E,
    // Derived by `refresh`.
    mixture: Mixture<S, N>,
    mass: f64,
    state: FluidState,
    /// The equation of state found a stable state in the last refresh.
    solved: bool,
}

impl<S: SpeciesSet<N>, const N: usize> ControlVolume<S, N> {
    /// Ideal gas of composition `mixture` at pressure `p` and temperature `t`
    /// filling `volume`.
    pub fn new(
        mixture: Mixture<S, N>,
        p: Pressure,
        t: ThermodynamicTemperature,
        volume: Volume,
    ) -> anyhow::Result<Self> {
        Self::with_eos(IdealGas, mixture, p, t, volume)
    }
}

impl<S: SpeciesSet<N>, const N: usize, E: EquationOfState<S, N>> ControlVolume<S, N, E> {
    /// Gas of composition `mixture` at pressure `p` and temperature `t` filling
    /// `volume`, described by the equation of state `eos`.
    pub fn with_eos(
        eos: E,
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
        let state = eos
            .state_pt(&mixture, pressure, temperature)
            .context("control volume")?;
        let mass = volume / state.v;
        let mut cv = Self {
            masses: mixture.mass_fractions().map(|y| y * mass),
            internal_energy: mass * state.u,
            volume,
            eos,
            mixture,
            mass,
            state,
            solved: true,
        };
        cv.refresh();
        Ok(cv)
    }

    /// Liquid and vapour in equilibrium at temperature `t`, the liquid filling
    /// `liquid_fraction` of `volume`. Needs an equation of state with a
    /// two-phase region, such as [`Reference`](super::eos::Reference).
    pub fn with_saturated(
        eos: E,
        mixture: Mixture<S, N>,
        t: ThermodynamicTemperature,
        liquid_fraction: Ratio,
        volume: Volume,
    ) -> anyhow::Result<Self> {
        let volume = volume.get::<cubic_meter>();
        if !(volume > 0.0 && volume.is_finite()) {
            return Err(PhysicsError::InvalidGeometry("volume must be positive"))
                .with_context(|| format!("V = {volume} m³"));
        }
        let alpha = liquid_fraction.get::<ratio>();
        if !(0.0..=1.0).contains(&alpha) {
            return Err(PhysicsError::InvalidState(
                "liquid fraction must be within 0..1",
            ))
            .with_context(|| format!("liquid fraction {alpha}"));
        }
        let state = eos
            .state_saturated(&mixture, t.get::<kelvin>(), alpha)
            .context("saturated control volume")?;
        let mass = volume / state.v;
        let mut cv = Self {
            masses: mixture.mass_fractions().map(|y| y * mass),
            internal_energy: mass * state.u,
            volume,
            eos,
            mixture,
            mass,
            state,
            solved: true,
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

    /// Recomputes the derived state, warm-starting the equation of state from
    /// the previous temperature.
    #[inline]
    fn refresh(&mut self) {
        self.mass = self.masses.iter().sum();
        self.mixture = Mixture::from_masses_unchecked(&self.masses);
        let state = self.eos.state_uv(
            &self.mixture,
            self.internal_energy / self.mass,
            self.volume / self.mass,
            &self.state,
        );
        // Keep the last good temperature as the next warm start.
        self.state = if state.t.is_nan() {
            FluidState {
                t: self.state.t,
                ..state
            }
        } else {
            state
        };
        self.solved = !state.t.is_nan();
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
        if !self.solved {
            return Err(PhysicsError::InvalidState(
                "the equation of state has no stable state here",
            ))
            .with_context(|| {
                format!(
                    "U = {} J, V = {} m³, m = {} kg",
                    self.internal_energy, self.volume, self.mass
                )
            });
        }
        check_temperature::<S, N>(self.state.t).context("control volume")?;
        if !(self.state.p > 0.0 && self.state.p.is_finite()) {
            return Err(PhysicsError::InvalidState("pressure must be positive"))
                .with_context(|| format!("p = {} Pa", self.state.p));
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

    /// Vapour mass fraction inside the two-phase dome, `None` for one phase.
    pub fn vapour_quality(&self) -> Option<f64> {
        (!self.state.quality.is_nan()).then_some(self.state.quality)
    }

    /// Fraction of the volume filled with liquid.
    pub fn liquid_fraction(&self) -> Ratio {
        Ratio::new::<ratio>(self.state.liquid_fraction)
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
        ThermodynamicTemperature::new::<kelvin>(self.state.t)
    }

    /// Pressure.
    pub fn pressure(&self) -> Pressure {
        Pressure::new::<pascal>(self.state.p)
    }

    /// Specific entropy of the contents.
    pub fn specific_entropy(&self) -> SpecificHeatCapacity {
        SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(
            self.eos.entropy(&self.mixture, &self.state),
        )
    }

    /// Rate of internal-energy change from boundary work, −p·dV/dt, in W, for a
    /// volume changing at `volume_rate` m³/s.
    #[inline]
    pub fn boundary_work_rate(&self, volume_rate: f64) -> f64 {
        -self.state.p * volume_rate
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

impl<S: SpeciesSet<N>, const N: usize, E: EquationOfState<S, N>> GasPort<S, N>
    for ControlVolume<S, N, E>
{
    fn reservoir(&self) -> Reservoir {
        Reservoir {
            pressure: self.state.p,
            temperature: self.state.t,
            gas: self.state.perfect_gas(),
        }
    }

    fn mass_fractions(&self) -> &[f64; N] {
        self.mixture.mass_fractions()
    }

    fn specific_enthalpy(&self) -> AvailableEnergy {
        AvailableEnergy::new::<joule_per_kilogram>(self.state.h)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        fluid::PureFluid,
        gas::eos::{IdealGas, PengRobinson, Reference},
        thermo::{Air3, ColdGas, MOLAR_GAS_CONSTANT},
    };

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

    fn nitrogen_tank(
        p: f64,
        t: f64,
    ) -> anyhow::Result<ControlVolume<ColdGas, 4, PengRobinson<ColdGas, 4>>> {
        ControlVolume::with_eos(
            PengRobinson::new()?,
            Mixture::from_mass_fractions([1.0, 0.0, 0.0, 0.0])?,
            Pressure::new::<pascal>(p),
            ThermodynamicTemperature::new::<kelvin>(t),
            Volume::new::<cubic_meter>(0.01),
        )
    }

    #[test]
    fn test_control_volume_with_ideal_gas_eos_matches_new() -> anyhow::Result<()> {
        let via_new = tank(4e5, 500.0, 0.2)?;
        let via_eos = ControlVolume::with_eos(
            IdealGas,
            Mixture::dry_air()?,
            Pressure::new::<pascal>(4e5),
            ThermodynamicTemperature::new::<kelvin>(500.0),
            Volume::new::<cubic_meter>(0.2),
        )?;
        // Debug output, not ==: the NaN quality of a single phase never
        // compares equal.
        let (a, b) = (format!("{via_new:?}"), format!("{via_eos:?}"));
        anyhow::ensure!(a == b, "{a} vs {b}");
        Ok(())
    }

    #[test]
    fn test_control_volume_peng_robinson_tank_holds_more_than_ideal() -> anyhow::Result<()> {
        // 100 bar of N2 at 200 K is ~16 % denser than an ideal gas (Z = 0.845).
        let cv = nitrogen_tank(1e7, 200.0)?;
        cv.validate()?;
        let (p, t) = (
            cv.pressure().get::<pascal>(),
            cv.temperature().get::<kelvin>(),
        );
        anyhow::ensure!(
            (p / 1e7 - 1.0).abs() < 1e-9 && (t / 200.0 - 1.0).abs() < 1e-10,
            "{p} Pa, {t} K"
        );
        let r = cv.mixture().gas_constant_si();
        let ideal_mass = 1e7 * 0.01 / (r * 200.0);
        anyhow::ensure!(
            cv.mass().get::<kilogram>() > 1.15 * ideal_mass,
            "mass {:?}",
            cv.mass()
        );
        Ok(())
    }

    #[test]
    fn test_control_volume_peng_robinson_round_trips_and_heats() -> anyhow::Result<()> {
        let mut cv = nitrogen_tank(1e7, 250.0)?;
        let masses = *cv.species_masses();
        let u = cv.internal_energy().get::<joule>();
        cv.set_conserved(&masses, u, 0.01);
        cv.validate()?;
        let t = cv.temperature().get::<kelvin>();
        anyhow::ensure!((t / 250.0 - 1.0).abs() < 1e-10, "T = {t}");
        let mut last = t;
        for step in 1..20 {
            cv.set_conserved(&masses, u + f64::from(step) * 1e4, 0.01);
            cv.validate()?;
            let t = cv.temperature().get::<kelvin>();
            anyhow::ensure!(t > last, "T fell to {t}");
            last = t;
        }
        Ok(())
    }

    #[test]
    fn test_control_volume_peng_robinson_rejects_a_state_inside_the_spinodal() -> anyhow::Result<()>
    {
        let co2 = Mixture::<ColdGas, 4>::from_mass_fractions([0.0, 0.0, 0.0, 1.0])?;
        let mut cv = ControlVolume::with_eos(
            PengRobinson::new()?,
            co2,
            Pressure::new::<pascal>(30e5),
            ThermodynamicTemperature::new::<kelvin>(250.0),
            Volume::new::<cubic_meter>(0.01),
        )?;
        cv.validate()?;
        // Same mass and energy spread over four times the volume.
        let masses = *cv.species_masses();
        cv.set_conserved(&masses, cv.internal_energy().get::<joule>(), 0.04);
        anyhow::ensure!(
            cv.validate().is_err(),
            "spinodal state accepted at {:?}",
            cv.pressure()
        );
        Ok(())
    }

    fn co2_bottle(
        t: f64,
        liquid: f64,
    ) -> anyhow::Result<ControlVolume<ColdGas, 4, Reference<ColdGas, 4>>> {
        ControlVolume::with_saturated(
            Reference::new(PureFluid::CarbonDioxide)?,
            Mixture::from_mass_fractions([0.0, 0.0, 0.0, 1.0])?,
            ThermodynamicTemperature::new::<kelvin>(t),
            Ratio::new::<ratio>(liquid),
            Volume::new::<cubic_meter>(1e-3),
        )
    }

    #[test]
    fn test_control_volume_reference_eos_needs_its_species_in_the_set() -> anyhow::Result<()> {
        for fluid in PureFluid::ALL {
            Reference::<ColdGas, 4>::new(fluid)?;
        }
        anyhow::ensure!(
            Reference::<Air3, 3>::new(PureFluid::Helium).is_err(),
            "He in Air3"
        );
        Ok(())
    }

    #[test]
    fn test_control_volume_saturated_co2_bottle_sits_on_the_saturation_line() -> anyhow::Result<()>
    {
        let cv = co2_bottle(290.0, 0.6)?;
        cv.validate()?;
        let sat =
            PureFluid::CarbonDioxide.saturation(ThermodynamicTemperature::new::<kelvin>(290.0))?;
        let p = cv.pressure().get::<pascal>();
        anyhow::ensure!(
            (p / sat.pressure.get::<pascal>() - 1.0).abs() < 1e-9,
            "p = {p}"
        );
        let alpha = cv.liquid_fraction().get::<ratio>();
        anyhow::ensure!((alpha - 0.6).abs() < 1e-9, "liquid fraction {alpha}");
        let x = cv.vapour_quality().context("single phase")?;
        anyhow::ensure!(x > 0.0 && x < 0.2, "quality {x}");
        // What leaves through the top is saturated vapour.
        let vapour = PureFluid::CarbonDioxide.homogeneous(
            ThermodynamicTemperature::new::<kelvin>(290.0),
            sat.vapour_density,
        )?;
        let h = cv.specific_enthalpy().get::<joule_per_kilogram>();
        anyhow::ensure!(
            (h - vapour.enthalpy.get::<joule_per_kilogram>()).abs() < 1e-6,
            "h = {h}"
        );
        Ok(())
    }

    #[test]
    fn test_control_volume_venting_vapour_flash_cools_a_co2_bottle() -> anyhow::Result<()> {
        let mut cv = co2_bottle(290.0, 0.6)?;
        let mut masses = *cv.species_masses();
        let mut u = cv.internal_energy().get::<joule>();
        let mut last = cv.temperature().get::<kelvin>();
        for _ in 0..40 {
            let dm = 2e-3;
            u -= dm * cv.specific_enthalpy().get::<joule_per_kilogram>();
            masses[3] -= dm;
            cv.set_conserved(&masses, u, 1e-3);
            cv.validate()?;
            let t = cv.temperature().get::<kelvin>();
            anyhow::ensure!(t < last, "T rose from {last} to {t}");
            anyhow::ensure!(cv.vapour_quality().is_some(), "left the dome at {t} K");
            last = t;
        }
        Ok(())
    }

    #[test]
    fn test_control_volume_reference_nitrogen_matches_the_fluid() -> anyhow::Result<()> {
        let cv = ControlVolume::with_eos(
            Reference::new(PureFluid::Nitrogen)?,
            Mixture::<ColdGas, 4>::from_mass_fractions([1.0, 0.0, 0.0, 0.0])?,
            Pressure::new::<pascal>(3e7),
            ThermodynamicTemperature::new::<kelvin>(300.0),
            Volume::new::<cubic_meter>(0.01),
        )?;
        cv.validate()?;
        let rho = PureFluid::Nitrogen
            .properties(
                Pressure::new::<pascal>(3e7),
                ThermodynamicTemperature::new::<kelvin>(300.0),
            )?
            .density
            .get::<crate::units::mass_density::kilogram_per_cubic_meter>();
        let m = cv.mass().get::<kilogram>();
        anyhow::ensure!((m / (rho * 0.01) - 1.0).abs() < 1e-12, "mass {m}");
        anyhow::ensure!(
            cv.vapour_quality().is_none() && cv.liquid_fraction().get::<ratio>() < 1.0,
            "phase"
        );
        Ok(())
    }

    #[test]
    fn test_control_volume_reference_eos_rejects_other_species() -> anyhow::Result<()> {
        let mut cv = co2_bottle(290.0, 0.6)?;
        let mut masses = *cv.species_masses();
        masses[0] = 1e-3; // a little nitrogen
        cv.set_conserved(&masses, cv.internal_energy().get::<joule>(), 1e-3);
        anyhow::ensure!(cv.validate().is_err(), "mixture accepted");
        let ideal = ControlVolume::with_saturated(
            IdealGas,
            Mixture::<Air3, 3>::dry_air()?,
            ThermodynamicTemperature::new::<kelvin>(100.0),
            Ratio::new::<ratio>(0.5),
            Volume::new::<cubic_meter>(1e-3),
        );
        anyhow::ensure!(ideal.is_err(), "ideal gas saturated");
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
