//! Pure fluids described by reference multiparameter Helmholtz equations of
//! state: Span et al. (2000) for N₂, Tegeler et al. (1999) for Ar, Ortiz-Vega
//! et al. (2023) for He and Span and Wagner (1996) for CO₂, as implemented in
//! `CoolProp`.
//!
//! Each equation's own reference state is replaced by NASA's: at 298.15 K the
//! ideal-gas enthalpy and entropy at 1 bar equal those of the NASA species
//! data, so energy carried between these and NASA-based volumes is consistent.

#[allow(
    clippy::unreadable_literal,
    clippy::excessive_precision,
    reason = "generated tables reproduce the published coefficients digit for digit"
)]
mod coefficients;
mod helmholtz;
mod saturation;
mod state;

use std::sync::LazyLock;

use anyhow::Context;
use coefficients::FluidData;
pub(crate) use helmholtz::{Molar, molar};
pub(crate) use saturation::saturation_si;
pub(crate) use state::{Phase, RawState, density_si, state_uv_si};

use crate::{
    PhysicsError,
    thermo::{MOLAR_GAS_CONSTANT, STANDARD_PRESSURE, Species, TemperaturePowers, gases},
    units::{
        AvailableEnergy, MassDensity, Pressure, SpecificHeatCapacity, SpecificVolume,
        ThermodynamicTemperature, Velocity, available_energy::joule_per_kilogram,
        mass_density::kilogram_per_cubic_meter, pressure::pascal,
        specific_heat_capacity::joule_per_kilogram_kelvin,
        specific_volume::cubic_meter_per_kilogram, thermodynamic_temperature::kelvin,
        velocity::meter_per_second,
    },
};

/// A pure fluid with a reference equation of state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PureFluid {
    /// N₂, Span et al. (2000).
    Nitrogen,
    /// Ar, Tegeler et al. (1999).
    Argon,
    /// He-4, Ortiz-Vega et al. (2023).
    Helium,
    /// CO₂, Span and Wagner (1996).
    CarbonDioxide,
}

/// A fluid's equation with its reference state aligned to NASA's: α⁰ gains
/// c₁ + c₂τ.
#[derive(Debug)]
pub(crate) struct Fluid {
    pub(crate) data: &'static FluidData,
    pub(crate) c1: f64,
    pub(crate) c2: f64,
}

/// Properties of a single-phase state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluidProperties {
    /// Pressure.
    pub pressure: Pressure,
    /// Temperature.
    pub temperature: ThermodynamicTemperature,
    /// Density.
    pub density: MassDensity,
    /// Specific internal energy (NASA reference state).
    pub internal_energy: AvailableEnergy,
    /// Specific enthalpy (NASA reference state).
    pub enthalpy: AvailableEnergy,
    /// Specific entropy (NASA reference state).
    pub entropy: SpecificHeatCapacity,
    /// Isochoric heat capacity.
    pub cv: SpecificHeatCapacity,
    /// Isobaric heat capacity.
    pub cp: SpecificHeatCapacity,
    /// Speed of sound.
    pub speed_of_sound: Velocity,
}

/// Saturated liquid and vapour at one temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Saturation {
    /// Saturation pressure.
    pub pressure: Pressure,
    /// Saturated liquid density.
    pub liquid_density: MassDensity,
    /// Saturated vapour density.
    pub vapour_density: MassDensity,
}

/// Equilibrium state of a pure fluid from its internal energy and volume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PureState {
    /// Pressure.
    pub pressure: Pressure,
    /// Temperature.
    pub temperature: ThermodynamicTemperature,
    /// Bulk density.
    pub density: MassDensity,
    /// Vapour mass fraction inside the two-phase dome, `None` for one phase.
    pub vapour_quality: Option<f64>,
    /// Specific internal energy (NASA reference state).
    pub internal_energy: AvailableEnergy,
    /// Specific enthalpy (NASA reference state).
    pub enthalpy: AvailableEnergy,
    /// Specific entropy (NASA reference state).
    pub entropy: SpecificHeatCapacity,
}

impl Fluid {
    /// `data` with c₁ and c₂ chosen so that at 298.15 K the ideal-gas enthalpy,
    /// and the entropy at 1 bar, equal those of the NASA `species`.
    fn aligned(data: &'static FluidData, species: &Species) -> Self {
        let t = 298.15;
        let tau = data.t_crit / t;
        let delta = STANDARD_PRESSURE / (data.gas_constant * t) / data.rho_crit;
        let powers = TemperaturePowers::from_kelvin(t);
        let h = species.h_over_rt(&powers) * MOLAR_GAS_CONSTANT * t;
        let s = species.s_over_r(&powers) * MOLAR_GAS_CONSTANT;
        let (a0, t0, _) = helmholtz::ideal_shape(data, tau);
        // h/RT = 1 + τα⁰_τ and s/R = τα⁰_τ − α⁰ for the ideal gas.
        let c2 = (h / (data.gas_constant * t) - 1.0 - t0) / tau;
        let c1 = t0 - delta.ln() - a0 - s / data.gas_constant;
        Self { data, c1, c2 }
    }
}

static NITROGEN: LazyLock<Fluid> =
    LazyLock::new(|| Fluid::aligned(&coefficients::nitrogen::NITROGEN, &gases::N2));
static ARGON: LazyLock<Fluid> =
    LazyLock::new(|| Fluid::aligned(&coefficients::argon::ARGON, &gases::AR));
static HELIUM: LazyLock<Fluid> =
    LazyLock::new(|| Fluid::aligned(&coefficients::helium::HELIUM, &gases::HE));
static CARBON_DIOXIDE: LazyLock<Fluid> =
    LazyLock::new(|| Fluid::aligned(&coefficients::carbon_dioxide::CARBON_DIOXIDE, &gases::CO2));

impl PureFluid {
    /// All four fluids.
    pub const ALL: [Self; 4] = [
        Self::Nitrogen,
        Self::Argon,
        Self::Helium,
        Self::CarbonDioxide,
    ];

    /// The fluid's equation, reference state aligned with NASA's.
    pub(crate) fn data(self) -> &'static Fluid {
        match self {
            Self::Nitrogen => &NITROGEN,
            Self::Argon => &ARGON,
            Self::Helium => &HELIUM,
            Self::CarbonDioxide => &CARBON_DIOXIDE,
        }
    }

    /// The NASA species this fluid is.
    pub fn species(self) -> &'static Species {
        match self {
            Self::Nitrogen => &gases::N2,
            Self::Argon => &gases::AR,
            Self::Helium => &gases::HE,
            Self::CarbonDioxide => &gases::CO2,
        }
    }

    /// The equation evaluated as one homogeneous phase at temperature `t` and
    /// density `rho`, whether or not that phase is stable there.
    pub fn homogeneous(
        self,
        t: ThermodynamicTemperature,
        rho: MassDensity,
    ) -> anyhow::Result<FluidProperties> {
        let (t, rho) = (t.get::<kelvin>(), rho.get::<kilogram_per_cubic_meter>());
        if !(t > 0.0 && t.is_finite() && rho > 0.0 && rho.is_finite()) {
            return Err(PhysicsError::InvalidState(
                "temperature and density must be positive",
            ))
            .with_context(|| format!("{self:?} at {t} K, {rho} kg/m³"));
        }
        let fluid = self.data();
        let m = fluid.data.molar_mass;
        Ok(self.properties_from(&helmholtz::molar(fluid, rho / m, t), t, rho))
    }

    /// Saturated liquid and vapour at temperature `t`, from the triple point to
    /// 0.01 % below the critical temperature.
    pub fn saturation(self, t: ThermodynamicTemperature) -> anyhow::Result<Saturation> {
        let fluid = self.data();
        let m = fluid.data.molar_mass;
        let (p, rho_l, rho_v) = saturation::saturation_si(fluid, t.get::<kelvin>())?;
        let density = MassDensity::new::<kilogram_per_cubic_meter>;
        Ok(Saturation {
            pressure: Pressure::new::<pascal>(p),
            liquid_density: density(rho_l * m),
            vapour_density: density(rho_v * m),
        })
    }

    /// Properties of the stable phase at pressure `p` and temperature `t`:
    /// liquid above the saturation pressure, vapour below it.
    pub fn properties(
        self,
        p: Pressure,
        t: ThermodynamicTemperature,
    ) -> anyhow::Result<FluidProperties> {
        let (p, t) = (p.get::<pascal>(), t.get::<kelvin>());
        let fluid = self.data();
        let rho = state::density_si(fluid, p, t)?;
        let m = fluid.data.molar_mass;
        Ok(self.properties_from(&helmholtz::molar(fluid, rho, t), t, rho * m))
    }

    /// Equilibrium state with specific internal energy `u` and specific volume
    /// `v`: one phase, or liquid and vapour inside the dome.
    pub fn state_from_uv(self, u: AvailableEnergy, v: SpecificVolume) -> anyhow::Result<PureState> {
        let fluid = self.data();
        let m = fluid.data.molar_mass;
        let (u, v) = (
            u.get::<joule_per_kilogram>(),
            v.get::<cubic_meter_per_kilogram>(),
        );
        let state = state::state_uv_si(fluid, u * m, 1.0 / (v * m), None)
            .with_context(|| format!("{self:?} at u = {u} J/kg, v = {v} m³/kg"))?;
        let energy = |q: f64| AvailableEnergy::new::<joule_per_kilogram>(q / m);
        Ok(PureState {
            pressure: Pressure::new::<pascal>(state.bulk.p),
            temperature: ThermodynamicTemperature::new::<kelvin>(state.t),
            density: MassDensity::new::<kilogram_per_cubic_meter>(state.rho * m),
            vapour_quality: (!state.quality.is_nan()).then_some(state.quality),
            internal_energy: energy(state.bulk.u),
            enthalpy: energy(state.bulk.h),
            entropy: SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(state.bulk.s / m),
        })
    }

    fn properties_from(self, x: &helmholtz::Molar, t: f64, rho: f64) -> FluidProperties {
        let m = self.data().data.molar_mass;
        let energy = |q: f64| AvailableEnergy::new::<joule_per_kilogram>(q / m);
        let capacity = |q: f64| SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(q / m);
        FluidProperties {
            pressure: Pressure::new::<pascal>(x.p),
            temperature: ThermodynamicTemperature::new::<kelvin>(t),
            density: MassDensity::new::<kilogram_per_cubic_meter>(rho),
            internal_energy: energy(x.u),
            enthalpy: energy(x.h),
            entropy: capacity(x.s),
            cv: capacity(x.cv),
            cp: capacity(x.cp),
            speed_of_sound: Velocity::new::<meter_per_second>(x.w),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[expect(
        clippy::many_single_char_names,
        reason = "thermodynamic symbols p, T, u, v, h, s, x"
    )]
    fn test_pure_fluids_share_the_nasa_reference_state() -> anyhow::Result<()> {
        // Ideal-gas h and s at 298.15 K equal the NASA species data; at 10 Pa
        // the real-gas departure is negligible.
        let (t, p) = (298.15, 10.0);
        for fluid in PureFluid::ALL {
            let species = fluid.species();
            let m = fluid.data().data.molar_mass;
            let rho = p / (fluid.data().data.gas_constant * t) * m;
            let x = fluid.homogeneous(
                ThermodynamicTemperature::new::<crate::units::thermodynamic_temperature::kelvin>(t),
                MassDensity::new::<crate::units::mass_density::kilogram_per_cubic_meter>(rho),
            )?;
            let powers = TemperaturePowers::from_kelvin(t);
            let h_nasa = species.h_over_rt(&powers) * MOLAR_GAS_CONSTANT * t;
            let s_nasa = species.s_over_r(&powers) * MOLAR_GAS_CONSTANT
                - MOLAR_GAS_CONSTANT * (p / STANDARD_PRESSURE).ln();
            let h = x
                .enthalpy
                .get::<crate::units::available_energy::joule_per_kilogram>()
                * m;
            let s = x
                .entropy
                .get::<crate::units::specific_heat_capacity::joule_per_kilogram_kelvin>()
                * m;
            // The tiny real-gas part at 1 bar is a few J/mol at most.
            anyhow::ensure!(
                (h - h_nasa).abs() < 20.0,
                "{fluid:?}: h {h} vs NASA {h_nasa}"
            );
            anyhow::ensure!(
                (s - s_nasa).abs() < 0.1,
                "{fluid:?}: s {s} vs NASA {s_nasa}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_state_from_uv_reports_quality_inside_the_dome() -> anyhow::Result<()> {
        use crate::units::{
            specific_volume::cubic_meter_per_kilogram, thermodynamic_temperature::kelvin,
        };
        let co2 = PureFluid::CarbonDioxide;
        let t = ThermodynamicTemperature::new::<kelvin>(270.0);
        let sat = co2.saturation(t)?;
        let (liquid, vapour) = (
            co2.homogeneous(t, sat.liquid_density)?,
            co2.homogeneous(t, sat.vapour_density)?,
        );
        let x = 0.4;
        let v = (1.0 - x)
            / sat
                .liquid_density
                .get::<crate::units::mass_density::kilogram_per_cubic_meter>()
            + x / sat
                .vapour_density
                .get::<crate::units::mass_density::kilogram_per_cubic_meter>();
        let u = liquid.internal_energy + (vapour.internal_energy - liquid.internal_energy) * x;
        let state = co2.state_from_uv(u, SpecificVolume::new::<cubic_meter_per_kilogram>(v))?;
        let quality = state
            .vapour_quality
            .ok_or_else(|| anyhow::anyhow!("single phase"))?;
        anyhow::ensure!((quality - x).abs() < 1e-9, "x = {quality}");
        anyhow::ensure!(
            (state.temperature.get::<kelvin>() - 270.0).abs() < 1e-8,
            "{state:?}"
        );
        // One phase: properties at (p, T) round-trip through (u, v).
        let gas = co2.properties(Pressure::new::<crate::units::pressure::pascal>(1e5), t)?;
        let back = co2.state_from_uv(
            gas.internal_energy,
            SpecificVolume::new::<cubic_meter_per_kilogram>(
                1.0 / gas
                    .density
                    .get::<crate::units::mass_density::kilogram_per_cubic_meter>(),
            ),
        )?;
        anyhow::ensure!(back.vapour_quality.is_none(), "{back:?}");
        anyhow::ensure!(
            (back.pressure.get::<crate::units::pressure::pascal>() / 1e5 - 1.0).abs() < 1e-9,
            "{back:?}"
        );
        Ok(())
    }

    #[test]
    fn test_pure_fluids_name_their_species() -> anyhow::Result<()> {
        let names = PureFluid::ALL.map(|f| f.species().name());
        anyhow::ensure!(names == ["N2", "Ar", "He", "CO2"], "{names:?}");
        Ok(())
    }
}
