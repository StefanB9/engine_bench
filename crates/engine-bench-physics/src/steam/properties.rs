//! Single-phase properties at (p, T): the public, unit-typed API.

use anyhow::Context;

use super::{
    gibbs::{self, RawProperties},
    region::{Region, region_si, validate_pt},
};
use crate::{
    PhysicsError,
    units::{
        AvailableEnergy, MassDensity, Pressure, SpecificHeatCapacity, SpecificVolume,
        ThermodynamicTemperature, Velocity, available_energy::joule_per_kilogram,
        mass_density::kilogram_per_cubic_meter, pressure::pascal,
        specific_heat_capacity::joule_per_kilogram_kelvin,
        specific_volume::cubic_meter_per_kilogram, thermodynamic_temperature::kelvin,
        velocity::meter_per_second,
    },
};

/// Water or steam properties at one (p, T) state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SteamProperties {
    /// IF97 region the state lies in.
    pub region: Region,
    /// Specific volume.
    pub specific_volume: SpecificVolume,
    /// Specific enthalpy.
    pub enthalpy: AvailableEnergy,
    /// Specific internal energy.
    pub internal_energy: AvailableEnergy,
    /// Specific entropy.
    pub entropy: SpecificHeatCapacity,
    /// Isobaric heat capacity.
    pub cp: SpecificHeatCapacity,
    /// Isochoric heat capacity.
    pub cv: SpecificHeatCapacity,
    /// Speed of sound.
    pub speed_of_sound: Velocity,
}

impl SteamProperties {
    /// Density, 1 / v.
    pub fn density(&self) -> MassDensity {
        MassDensity::new::<kilogram_per_cubic_meter>(
            1.0 / self.specific_volume.get::<cubic_meter_per_kilogram>(),
        )
    }

    fn from_raw(region: Region, x: RawProperties) -> Self {
        let energy = AvailableEnergy::new::<joule_per_kilogram>;
        let capacity = SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>;
        Self {
            region,
            specific_volume: SpecificVolume::new::<cubic_meter_per_kilogram>(x.v),
            enthalpy: energy(x.h),
            internal_energy: energy(x.u),
            entropy: capacity(x.s),
            cp: capacity(x.cp),
            cv: capacity(x.cv),
            speed_of_sound: Velocity::new::<meter_per_second>(x.w),
        }
    }
}

/// Properties of water or steam at pressure `p` and temperature `t`.
///
/// Fails outside IF97's range, and for now in region 3 (near-critical), which
/// follows in the next step of milestone 3.
pub fn properties(p: Pressure, t: ThermodynamicTemperature) -> anyhow::Result<SteamProperties> {
    let (p, t) = (p.get::<pascal>(), t.get::<kelvin>());
    validate_pt(p, t).context("steam properties")?;
    let region = region_si(p, t);
    let raw = properties_si(region, p, t)
        .with_context(|| format!("steam properties at {p} Pa, {t} K"))?;
    Ok(SteamProperties::from_raw(region, raw))
}

/// Raw properties of a validated (p, T) state in `region`.
pub(crate) fn properties_si(region: Region, p: f64, t: f64) -> anyhow::Result<RawProperties> {
    match region {
        Region::CompressedLiquid => Ok(gibbs::region1(p, t)),
        Region::Vapour => Ok(gibbs::region2(p, t)),
        Region::HighTemperature => Ok(gibbs::region5(p, t)),
        Region::NearCritical => Err(PhysicsError::InvalidState(
            "IF97 region 3 (near-critical) is not supported yet",
        )
        .into()),
        Region::Saturation => Err(PhysicsError::InvalidState(
            "a saturated state needs a quality, not just (p, T)",
        )
        .into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn at(p_pa: f64, t: f64) -> anyhow::Result<SteamProperties> {
        properties(
            Pressure::new::<pascal>(p_pa),
            ThermodynamicTemperature::new::<kelvin>(t),
        )
    }

    #[test]
    fn test_properties_wraps_the_region_kernel_in_si_units() -> anyhow::Result<()> {
        let x = at(3e6, 300.0)?;
        anyhow::ensure!(x.region == Region::CompressedLiquid, "{:?}", x.region);
        let raw = gibbs::region1(3e6, 300.0);
        anyhow::ensure!(
            x.enthalpy.get::<joule_per_kilogram>().to_bits() == raw.h.to_bits(),
            "h"
        );
        anyhow::ensure!(
            x.speed_of_sound.get::<meter_per_second>().to_bits() == raw.w.to_bits(),
            "w"
        );
        let rho = x.density().get::<kilogram_per_cubic_meter>();
        anyhow::ensure!((rho * raw.v - 1.0).abs() < 1e-15, "ρ = {rho}");
        Ok(())
    }

    #[test]
    fn test_properties_uses_region_5_above_1073_k() -> anyhow::Result<()> {
        let x = at(30e6, 2000.0)?;
        anyhow::ensure!(x.region == Region::HighTemperature, "{:?}", x.region);
        anyhow::ensure!(
            x.cv.get::<joule_per_kilogram_kelvin>().to_bits()
                == gibbs::region5(30e6, 2000.0).cv.to_bits(),
            "cv"
        );
        Ok(())
    }

    #[test]
    fn test_properties_reports_region_3_as_not_supported_yet() -> anyhow::Result<()> {
        let error = at(25e6, 650.0).err().context("region 3 accepted")?;
        anyhow::ensure!(format!("{error:#}").contains("region 3"), "{error:#}");
        Ok(())
    }

    #[test]
    fn test_properties_rejects_states_outside_if97() -> anyhow::Result<()> {
        anyhow::ensure!(at(1e6, 250.0).is_err(), "250 K accepted");
        anyhow::ensure!(at(2e8, 300.0).is_err(), "200 MPa accepted");
        Ok(())
    }
}
