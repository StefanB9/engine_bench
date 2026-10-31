//! Single-phase properties at (p, T): the public, unit-typed API.

use anyhow::Context;

use super::{
    T_REGION3_MIN,
    gibbs::{self, RawProperties},
    helmholtz::region3,
    region::{Region, region_si, validate_pt},
    region3::{density_si, saturated_si},
    region3_backward::Side,
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
/// Fails outside IF97's range.
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
        Region::NearCritical => {
            let rho = density_si(p, t, Side::Either)?;
            Ok(region3(rho, t).properties)
        }
        Region::Saturation => Err(PhysicsError::InvalidState(
            "a saturated state needs a quality, not just (p, T)",
        )
        .into()),
    }
}

/// Saturated liquid and vapour at one temperature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Saturated {
    /// Saturation pressure.
    pub pressure: Pressure,
    /// Saturated liquid.
    pub liquid: SteamProperties,
    /// Saturated vapour.
    pub vapour: SteamProperties,
}

/// Saturated liquid and vapour at temperature `t`, from 273.15 K up to the
/// critical point (regions 1 and 2 up to 623.15 K, region 3 above).
pub fn saturated(t: ThermodynamicTemperature) -> anyhow::Result<Saturated> {
    let t = t.get::<kelvin>();
    let (p, liquid, vapour) = saturated_si(t).context("saturated liquid and vapour")?;
    let (liquid_region, vapour_region) = if t <= T_REGION3_MIN {
        (Region::CompressedLiquid, Region::Vapour)
    } else {
        (Region::NearCritical, Region::NearCritical)
    };
    Ok(Saturated {
        pressure: Pressure::new::<pascal>(p),
        liquid: SteamProperties::from_raw(liquid_region, liquid),
        vapour: SteamProperties::from_raw(vapour_region, vapour),
    })
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
    fn test_properties_solves_region_3_from_p_and_t() -> anyhow::Result<()> {
        // IAPWS R7-97(2012) Table 33, first row, entered as (p, T).
        let x = at(25.583_701_8e6, 650.0)?;
        anyhow::ensure!(x.region == Region::NearCritical, "{:?}", x.region);
        let rho = x.density().get::<kilogram_per_cubic_meter>();
        let h = x.enthalpy.get::<joule_per_kilogram>();
        anyhow::ensure!((rho / 500.0 - 1.0).abs() < 1e-8, "ρ = {rho}");
        anyhow::ensure!((h / 1_863_430.19 - 1.0).abs() < 5e-8, "h = {h}");
        Ok(())
    }

    #[test]
    fn test_saturated_marks_the_region_of_each_phase() -> anyhow::Result<()> {
        for (t, liquid, vapour) in [
            (500.0, Region::CompressedLiquid, Region::Vapour),
            (640.0, Region::NearCritical, Region::NearCritical),
        ] {
            let sat = saturated(ThermodynamicTemperature::new::<kelvin>(t))?;
            anyhow::ensure!(
                sat.liquid.region == liquid && sat.vapour.region == vapour,
                "{t} K: {sat:?}"
            );
            anyhow::ensure!(sat.liquid.density() > sat.vapour.density(), "{t} K");
        }
        anyhow::ensure!(
            saturated(ThermodynamicTemperature::new::<kelvin>(650.0)).is_err(),
            "650 K"
        );
        Ok(())
    }

    #[test]
    fn test_properties_rejects_states_outside_if97() -> anyhow::Result<()> {
        anyhow::ensure!(at(1e6, 250.0).is_err(), "250 K accepted");
        anyhow::ensure!(at(2e8, 300.0).is_err(), "200 MPa accepted");
        Ok(())
    }
}
