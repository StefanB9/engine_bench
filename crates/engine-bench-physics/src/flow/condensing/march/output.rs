//! What a condensing march returns: the flow along the nozzle and its
//! summary quantities.

use super::Local;
use crate::{
    flow::condensing::Contour,
    units::{
        Area, Force, Length, MassDensity, MassRate, Pressure, SpecificHeatCapacity,
        TemperatureInterval, ThermodynamicTemperature, Velocity, area::square_meter, force::newton,
        length::meter, mass_density::kilogram_per_cubic_meter, mass_rate::kilogram_per_second,
        pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin, temperature_interval,
        thermodynamic_temperature::kelvin, velocity::meter_per_second,
    },
};

/// The flow at one station of the march.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlowSample {
    /// Distance downstream of the throat.
    pub position: Length,
    /// Pressure.
    pub pressure: Pressure,
    /// Vapour temperature.
    pub temperature: ThermodynamicTemperature,
    /// Flow velocity (droplets and vapour together).
    pub velocity: Velocity,
    /// Mixture density.
    pub density: MassDensity,
    /// Mixture entropy.
    pub entropy: SpecificHeatCapacity,
    /// Velocity over the frozen speed of sound.
    pub frozen_mach: f64,
    /// Liquid mass fraction.
    pub liquid_fraction: f64,
    /// Droplets per kilogram of mixture.
    pub droplets_per_kilogram: f64,
    /// Sauter mean radius of the droplets (zero without).
    pub sauter_radius: Length,
    /// Nucleation rate, droplets per cubic metre per second.
    pub nucleation_rate: f64,
    /// Saturation temperature less vapour temperature (negative when
    /// superheated).
    pub subcooling: TemperatureInterval,
    /// Impulse function pA + ṁu.
    pub(crate) impulse: f64,
}

impl FlowSample {
    pub(super) fn from_local(x: f64, local: &Local, contour: &Contour) -> Self {
        let s = &local.station;
        let area = contour.area_si(x).0;
        Self {
            position: Length::new::<meter>(x),
            pressure: Pressure::new::<pascal>(s.vapour.p),
            temperature: ThermodynamicTemperature::new::<kelvin>(s.t),
            velocity: Velocity::new::<meter_per_second>(s.u),
            density: MassDensity::new::<kilogram_per_cubic_meter>(s.rho),
            entropy: SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(s.entropy()),
            frozen_mach: s.u / s.frozen_sound_speed(),
            liquid_fraction: s.y,
            droplets_per_kilogram: local.number,
            sauter_radius: Length::new::<meter>(local.sauter_radius()),
            nucleation_rate: local.rate,
            subcooling: TemperatureInterval::new::<temperature_interval::kelvin>(s.t_sat - s.t),
            impulse: s.vapour.p * area + s.rho * s.u * area * s.u,
        }
    }
}

/// A condensing nozzle flow: the stations along it and what it delivers.
#[derive(Debug, Clone, PartialEq)]
pub struct CondensingFlow {
    samples: Vec<FlowSample>,
    mass_flow: f64,
    exit_area: f64,
    shock: Option<f64>,
}

impl CondensingFlow {
    pub(super) fn new(
        samples: Vec<FlowSample>,
        mass_flow: f64,
        exit_area: f64,
        shock: Option<f64>,
    ) -> Self {
        Self {
            samples,
            mass_flow,
            exit_area,
            shock,
        }
    }

    /// The stations, from just past the throat to the exit.
    pub fn samples(&self) -> &[FlowSample] {
        &self.samples
    }

    /// The exit station.
    pub fn exit(&self) -> Option<&FlowSample> {
        self.samples.last()
    }

    /// Mass flow through the nozzle.
    pub fn mass_flow(&self) -> MassRate {
        MassRate::new::<kilogram_per_second>(self.mass_flow)
    }

    /// Exit area.
    pub fn exit_area(&self) -> Area {
        Area::new::<square_meter>(self.exit_area)
    }

    /// Where the normal shock stands, if there is one.
    pub fn shock(&self) -> Option<Length> {
        self.shock.map(Length::new::<meter>)
    }

    /// Thrust into `ambient` pressure pₐ: ṁ·uₑ + (pₑ − pₐ)·Aₑ.
    pub fn thrust(&self, ambient: Pressure) -> Option<Force> {
        let exit = self.exit()?;
        Some(Force::new::<newton>(
            exit.impulse - ambient.get::<pascal>() * self.exit_area,
        ))
    }

    /// The Wilson point: the greatest subcooling, where condensation takes
    /// over from the expansion; `None` if that is still rising at the exit
    /// (condensation has not set in) or there is no subcooling.
    pub fn wilson_point(&self) -> Option<Length> {
        let (index, peak) = self
            .samples
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.subcooling.value.total_cmp(&b.1.subcooling.value))?;
        (index + 1 < self.samples.len() && peak.subcooling.value > 0.0).then_some(peak.position)
    }

    /// Where the liquid fraction first reaches `fraction`.
    pub fn onset(&self, fraction: f64) -> Option<Length> {
        self.samples
            .iter()
            .find(|s| s.liquid_fraction >= fraction)
            .map(|s| s.position)
    }
}
