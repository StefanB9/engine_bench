//! IF97 regions and region selection from (p, T).

use anyhow::Context;

use super::{
    P_MAX, P_MAX_HIGH_TEMPERATURE, T_B23_MAX, T_MAX, T_MAX_VAPOUR, T_MIN, T_REGION3_MIN,
    saturation::{b23_pressure_si, saturation_pressure_si},
};
use crate::{
    PhysicsError,
    units::{
        Pressure, ThermodynamicTemperature, pressure::pascal, thermodynamic_temperature::kelvin,
    },
};

/// The five regions of IAPWS-IF97.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Region {
    /// Region 1: compressed liquid, 273.15–623.15 K, at or above the saturation
    /// pressure.
    CompressedLiquid,
    /// Region 2: vapour up to 1073.15 K, below the saturation pressure and the
    /// B23 boundary.
    Vapour,
    /// Region 3: dense fluid around the critical point, 623.15 K up to the B23
    /// boundary, above region 1 and 2 pressures.
    NearCritical,
    /// Region 4: the liquid–vapour saturation line.
    Saturation,
    /// Region 5: high-temperature vapour, 1073.15–2273.15 K, up to 50 MPa.
    HighTemperature,
}

impl Region {
    /// The region's number in IAPWS-IF97 (1 to 5).
    pub const fn number(self) -> u8 {
        match self {
            Self::CompressedLiquid => 1,
            Self::Vapour => 2,
            Self::NearCritical => 3,
            Self::Saturation => 4,
            Self::HighTemperature => 5,
        }
    }
}

/// Region of the single-phase state (p, T). A state exactly on the saturation
/// line counts as liquid, so this never returns [`Region::Saturation`].
pub fn region(p: Pressure, t: ThermodynamicTemperature) -> anyhow::Result<Region> {
    let (p, t) = (p.get::<pascal>(), t.get::<kelvin>());
    validate_pt(p, t).context("IF97 region")?;
    Ok(region_si(p, t))
}

/// Checks that (p, T) lies inside IF97's range of validity.
pub(crate) fn validate_pt(p: f64, t: f64) -> anyhow::Result<()> {
    if !(T_MIN..=T_MAX).contains(&t) {
        return Err(PhysicsError::TemperatureOutOfRange {
            temperature: t,
            min: T_MIN,
            max: T_MAX,
        }
        .into());
    }
    let p_max = if t > T_MAX_VAPOUR {
        P_MAX_HIGH_TEMPERATURE
    } else {
        P_MAX
    };
    if !(p > 0.0 && p <= p_max) {
        return Err(PhysicsError::PressureOutOfRange {
            pressure: p,
            min: 0.0,
            max: p_max,
        })
        .with_context(|| format!("at {t} K"));
    }
    Ok(())
}

/// Region of a validated (p, T) state.
pub(crate) fn region_si(p: f64, t: f64) -> Region {
    if t <= T_REGION3_MIN {
        if p >= saturation_pressure_si(t) {
            Region::CompressedLiquid
        } else {
            Region::Vapour
        }
    } else if t <= T_B23_MAX && p > b23_pressure_si(t) {
        Region::NearCritical
    } else if t <= T_MAX_VAPOUR {
        Region::Vapour
    } else {
        Region::HighTemperature
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region_at(p_mpa: f64, t: f64) -> anyhow::Result<Region> {
        region(
            Pressure::new::<pascal>(p_mpa * 1e6),
            ThermodynamicTemperature::new::<kelvin>(t),
        )
    }

    #[test]
    fn test_region_number_follows_iapws() -> anyhow::Result<()> {
        let numbers = [
            Region::CompressedLiquid,
            Region::Vapour,
            Region::NearCritical,
            Region::Saturation,
            Region::HighTemperature,
        ]
        .map(Region::number);
        anyhow::ensure!(numbers == [1, 2, 3, 4, 5], "{numbers:?}");
        Ok(())
    }

    #[test]
    fn test_region_selects_each_region_from_p_and_t() -> anyhow::Result<()> {
        let cases = [
            (3.0, 300.0, Region::CompressedLiquid),
            (100.0, 623.15, Region::CompressedLiquid),
            (0.0035, 300.0, Region::Vapour),
            (30.0, 700.0, Region::Vapour),
            (100.0, 900.0, Region::Vapour),
            (25.0, 650.0, Region::NearCritical),
            (100.0, 700.0, Region::NearCritical),
            (0.5, 1500.0, Region::HighTemperature),
            (50.0, 2273.15, Region::HighTemperature),
            // 1073.15 K belongs to region 2.
            (10.0, 1073.15, Region::Vapour),
        ];
        for (p, t, want) in cases {
            let got = region_at(p, t)?;
            anyhow::ensure!(got == want, "{p} MPa, {t} K: {got:?}, expected {want:?}");
        }
        Ok(())
    }

    #[test]
    fn test_region_puts_the_saturation_line_in_region_1() -> anyhow::Result<()> {
        let p_sat = crate::steam::saturation::saturation_pressure_si(373.15);
        let got = region(
            Pressure::new::<pascal>(p_sat),
            ThermodynamicTemperature::new::<kelvin>(373.15),
        )?;
        anyhow::ensure!(got == Region::CompressedLiquid, "{got:?}");
        Ok(())
    }

    #[test]
    fn test_region_rejects_states_outside_if97() -> anyhow::Result<()> {
        for (p, t) in [
            (1.0, 273.0),
            (1.0, 2300.0),
            (0.0, 300.0),
            (-1.0, 300.0),
            (101.0, 300.0),
            (51.0, 1500.0),
            (f64::NAN, 300.0),
            (1.0, f64::INFINITY),
        ] {
            anyhow::ensure!(region_at(p, t).is_err(), "{p} MPa, {t} K was accepted");
        }
        let error = region_at(101.0, 300.0).err().context("101 MPa accepted")?;
        anyhow::ensure!(
            matches!(
                error.downcast_ref::<PhysicsError>(),
                Some(PhysicsError::PressureOutOfRange { .. })
            ),
            "{error:#}"
        );
        let error = region_at(1.0, 273.0).err().context("273 K accepted")?;
        anyhow::ensure!(
            matches!(
                error.downcast_ref::<PhysicsError>(),
                Some(PhysicsError::TemperatureOutOfRange { .. })
            ),
            "{error:#}"
        );
        Ok(())
    }
}
