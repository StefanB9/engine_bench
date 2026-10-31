//! Real-fluid compressible flow: isentropic expansion along any equation of
//! state (the homogeneous equilibrium model inside the two-phase dome),
//! choking, orifices, nozzles and normal shocks with real-gas
//! Rankine–Hugoniot jumps.
//!
//! The same quasi-1D physics as the perfect-gas solvers, without γ: mass flux
//! along an isentrope is G(p) = ρ·√(2(h₀ − h)), maximal where the flow reaches
//! the equilibrium speed of sound.

mod isentrope;
mod nozzle;
mod orifice;
mod perfect;
mod shock;

pub use shock::normal_shock;

use crate::units::{
    AvailableEnergy, MassDensity, Pressure, SpecificHeatCapacity, ThermodynamicTemperature,
    Velocity, available_energy::joule_per_kilogram, mass_density::kilogram_per_cubic_meter,
    pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
    thermodynamic_temperature::kelvin, velocity::meter_per_second,
};

/// An equation of state that flow can expand through and shock in.
pub trait RealFluid {
    /// Equilibrium state at pressure `p` and specific entropy `s`: liquid and
    /// vapour in equilibrium inside the dome.
    fn state_ps(&self, p: Pressure, s: SpecificHeatCapacity) -> anyhow::Result<FlowState>;

    /// Single-phase state at pressure `p` and density `rho`.
    fn state_prho(&self, p: Pressure, rho: MassDensity) -> anyhow::Result<FlowState>;
}

/// A static state in a flow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlowState {
    /// Pressure.
    pub pressure: Pressure,
    /// Temperature.
    pub temperature: ThermodynamicTemperature,
    /// Density.
    pub density: MassDensity,
    /// Specific enthalpy.
    pub enthalpy: AvailableEnergy,
    /// Specific entropy.
    pub entropy: SpecificHeatCapacity,
    /// Equilibrium speed of sound, √(∂p/∂ρ) at constant entropy.
    pub speed_of_sound: Velocity,
    /// Vapour mass fraction inside the dome, `None` for one phase.
    pub vapour_quality: Option<f64>,
}

impl FlowState {
    /// A state from raw SI values: Pa, K, kg/m³, J/kg, J/(kg·K), m/s; quality
    /// NaN for one phase.
    #[expect(
        clippy::many_single_char_names,
        reason = "the thermodynamic symbols p, T, ρ, h, s, w"
    )]
    pub(crate) fn from_si(p: f64, t: f64, rho: f64, h: f64, s: f64, w: f64, quality: f64) -> Self {
        Self {
            pressure: Pressure::new::<pascal>(p),
            temperature: ThermodynamicTemperature::new::<kelvin>(t),
            density: MassDensity::new::<kilogram_per_cubic_meter>(rho),
            enthalpy: AvailableEnergy::new::<joule_per_kilogram>(h),
            entropy: SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(s),
            speed_of_sound: Velocity::new::<meter_per_second>(w),
            vapour_quality: (!quality.is_nan()).then_some(quality),
        }
    }

    pub(crate) fn p(&self) -> f64 {
        self.pressure.get::<pascal>()
    }

    pub(crate) fn t(&self) -> f64 {
        self.temperature.get::<kelvin>()
    }

    pub(crate) fn rho(&self) -> f64 {
        self.density.get::<kilogram_per_cubic_meter>()
    }

    pub(crate) fn h(&self) -> f64 {
        self.enthalpy.get::<joule_per_kilogram>()
    }

    pub(crate) fn s(&self) -> f64 {
        self.entropy.get::<joule_per_kilogram_kelvin>()
    }

    pub(crate) fn w(&self) -> f64 {
        self.speed_of_sound.get::<meter_per_second>()
    }
}

/// Stagnation (total) state of a flow: a reservoir at rest.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stagnation {
    /// Pa.
    pub(crate) p: f64,
    /// J/kg.
    pub(crate) h: f64,
    /// J/(kg·K).
    pub(crate) s: f64,
}

impl Stagnation {
    /// The stagnation state of fluid at rest in `state`.
    pub fn at_rest(state: &FlowState) -> Self {
        Self {
            p: state.p(),
            h: state.h(),
            s: state.s(),
        }
    }
}
