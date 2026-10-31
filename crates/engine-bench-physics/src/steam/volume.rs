//! A closed, well-mixed volume of water: liquid, vapour, or both in
//! equilibrium, tracked by its conserved mass and internal energy.

use anyhow::Context;

use super::{
    RHO_CRIT,
    properties::properties_si,
    region::{Region, region_si, validate_pt},
    region3::saturated_si,
    state::{Phases, RawState, SteamState, state_from_uv_si},
};
use crate::{
    PhysicsError,
    units::{
        AvailableEnergy, Energy, Mass, Pressure, Ratio, ThermodynamicTemperature, Volume,
        available_energy::joule_per_kilogram, energy::joule, mass::kilogram, pressure::pascal,
        ratio::ratio, thermodynamic_temperature::kelvin, volume::cubic_meter,
    },
};

/// Water in a volume of fixed or slowly varying size, tracked by its conserved
/// mass and internal energy.
///
/// The ODE state of a model holds the mass and internal energy; after every
/// change call [`set_conserved`](Self::set_conserved), which re-solves the
/// equilibrium state warm-started from the previous one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SteamVolume {
    /// kg.
    mass: f64,
    /// J.
    internal_energy: f64,
    /// m³.
    volume: f64,
    /// Equilibrium state; NaN fields after a failed solve.
    state: RawState,
    /// Last successfully solved state: the warm start for the next solve.
    hint: Option<RawState>,
    /// The last solve failed; `validate` reports it.
    failed: bool,
}

impl SteamVolume {
    /// Single-phase water at pressure `p` and temperature `t` filling `volume`.
    pub fn new(p: Pressure, t: ThermodynamicTemperature, volume: Volume) -> anyhow::Result<Self> {
        let (p, t, volume) = (
            p.get::<pascal>(),
            t.get::<kelvin>(),
            volume.get::<cubic_meter>(),
        );
        check_volume(volume)?;
        validate_pt(p, t).context("steam volume")?;
        let region = region_si(p, t);
        let x = properties_si(region, p, t).context("steam volume")?;
        let mass = volume / x.v;
        let state = RawState {
            region,
            p,
            t,
            quality: f64::NAN,
            v: x.v,
            u: x.u,
            h: x.h,
            s: x.s,
            phases: Phases::NONE,
        };
        Ok(Self {
            mass,
            internal_energy: mass * x.u,
            volume,
            state,
            hint: Some(state),
            failed: false,
        })
    }

    /// Saturated liquid and vapour at temperature `t`, the liquid taking
    /// `liquid_fraction` of `volume` (0 = all vapour, 1 = all liquid).
    pub fn saturated(
        t: ThermodynamicTemperature,
        liquid_fraction: Ratio,
        volume: Volume,
    ) -> anyhow::Result<Self> {
        let (t, alpha, volume) = (
            t.get::<kelvin>(),
            liquid_fraction.get::<ratio>(),
            volume.get::<cubic_meter>(),
        );
        check_volume(volume)?;
        if !(0.0..=1.0).contains(&alpha) {
            return Err(PhysicsError::InvalidState(
                "liquid fraction must be within 0..1",
            ))
            .with_context(|| format!("liquid fraction {alpha}"));
        }
        let (p, liquid, vapour) = saturated_si(t).context("saturated steam volume")?;
        let liquid_mass = alpha * volume / liquid.v;
        let vapour_mass = (1.0 - alpha) * volume / vapour.v;
        let mass = liquid_mass + vapour_mass;
        let x = vapour_mass / mass;
        let mix = |a: f64, b: f64| x.mul_add(b - a, a);
        let state = RawState {
            region: Region::Saturation,
            p,
            t,
            quality: x,
            v: volume / mass,
            u: mix(liquid.u, vapour.u),
            h: mix(liquid.h, vapour.h),
            s: mix(liquid.s, vapour.s),
            phases: Phases {
                v_liquid: liquid.v,
                v_vapour: vapour.v,
                h_liquid: liquid.h,
                h_vapour: vapour.h,
            },
        };
        Ok(Self {
            mass,
            internal_energy: mass * state.u,
            volume,
            state,
            hint: Some(state),
            failed: false,
        })
    }

    /// Replaces the conserved state: mass in kg, internal energy in J and
    /// volume in m³. Kernel: never fails; check the result with
    /// [`validate`](Self::validate) once per solver step.
    pub fn set_conserved(&mut self, mass: f64, internal_energy: f64, volume: f64) {
        self.mass = mass;
        self.internal_energy = internal_energy;
        self.volume = volume;
        let solved = (mass > 0.0 && mass.is_finite() && volume > 0.0 && volume.is_finite())
            .then(|| {
                state_from_uv_si(internal_energy / mass, volume / mass, self.hint.as_ref()).ok()
            })
            .flatten();
        if let Some(state) = solved {
            self.state = state;
            self.hint = Some(state);
            self.failed = false;
        } else {
            self.state = RawState {
                p: f64::NAN,
                t: f64::NAN,
                quality: f64::NAN,
                v: f64::NAN,
                u: f64::NAN,
                h: f64::NAN,
                s: f64::NAN,
                phases: Phases::NONE,
                ..self.state
            };
            self.failed = true;
        }
    }

    /// Checks that the last [`set_conserved`](Self::set_conserved) found an
    /// IF97 state for positive, finite mass and volume.
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.failed {
            return Err(PhysicsError::InvalidState(
                "no IF97 state for this mass, internal energy and volume",
            ))
            .with_context(|| {
                format!(
                    "m = {} kg, U = {} J, V = {} m³",
                    self.mass, self.internal_energy, self.volume
                )
            });
        }
        Ok(())
    }

    /// Mass of water.
    pub fn mass(&self) -> Mass {
        Mass::new::<kilogram>(self.mass)
    }

    /// Internal energy of the contents.
    pub fn internal_energy(&self) -> Energy {
        Energy::new::<joule>(self.internal_energy)
    }

    /// Size of the volume.
    pub fn volume(&self) -> Volume {
        Volume::new::<cubic_meter>(self.volume)
    }

    /// Pressure.
    pub fn pressure(&self) -> Pressure {
        Pressure::new::<pascal>(self.state.p)
    }

    /// Temperature.
    pub fn temperature(&self) -> ThermodynamicTemperature {
        ThermodynamicTemperature::new::<kelvin>(self.state.t)
    }

    /// The full equilibrium state of the contents.
    pub fn state(&self) -> SteamState {
        self.state.to_state()
    }

    /// Fraction of the volume taken by liquid: 1 for compressed liquid (and
    /// supercritical fluid denser than critical), 0 for vapour.
    pub fn liquid_fraction(&self) -> Ratio {
        let s = &self.state;
        let alpha = match s.region {
            Region::Saturation => (1.0 - s.quality) * s.phases.v_liquid / s.v,
            Region::CompressedLiquid => 1.0,
            Region::NearCritical if s.v < 1.0 / RHO_CRIT => 1.0,
            Region::NearCritical | Region::Vapour | Region::HighTemperature => 0.0,
        };
        Ratio::new::<ratio>(alpha)
    }

    /// Specific enthalpy of what leaves through a port near the top: saturated
    /// vapour inside the dome, the bulk otherwise.
    pub fn vapour_enthalpy(&self) -> AvailableEnergy {
        let h = if self.state.region == Region::Saturation {
            self.state.phases.h_vapour
        } else {
            self.state.h
        };
        AvailableEnergy::new::<joule_per_kilogram>(h)
    }

    /// Specific enthalpy of what leaves through a port near the bottom:
    /// saturated liquid inside the dome, the bulk otherwise.
    pub fn liquid_enthalpy(&self) -> AvailableEnergy {
        let h = if self.state.region == Region::Saturation {
            self.state.phases.h_liquid
        } else {
            self.state.h
        };
        AvailableEnergy::new::<joule_per_kilogram>(h)
    }
}

fn check_volume(volume: f64) -> anyhow::Result<()> {
    if !(volume > 0.0 && volume.is_finite()) {
        return Err(PhysicsError::InvalidGeometry("volume must be positive"))
            .with_context(|| format!("V = {volume} m³"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    fn kelvin(t: f64) -> ThermodynamicTemperature {
        ThermodynamicTemperature::new::<kelvin>(t)
    }

    fn litres(v: f64) -> Volume {
        Volume::new::<cubic_meter>(v * 1e-3)
    }

    fn boiler(t: f64, liquid_fraction: f64) -> anyhow::Result<SteamVolume> {
        SteamVolume::saturated(kelvin(t), Ratio::new::<ratio>(liquid_fraction), litres(1.0))
    }

    #[test]
    fn test_steam_volume_new_holds_the_given_single_phase_state() -> anyhow::Result<()> {
        let cv = SteamVolume::new(Pressure::new::<pascal>(1e5), kelvin(400.0), litres(10.0))?;
        cv.validate()?;
        let x = properties_si(Region::Vapour, 1e5, 400.0)?;
        anyhow::ensure!(
            (cv.mass().get::<kilogram>() * x.v / 0.01 - 1.0).abs() < 1e-12,
            "mass"
        );
        anyhow::ensure!(
            (cv.pressure().get::<pascal>() / 1e5 - 1.0).abs() < 1e-9,
            "p"
        );
        anyhow::ensure!((cv.temperature().get::<kelvin>() - 400.0).abs() < 1e-9, "T");
        anyhow::ensure!(
            cv.liquid_fraction().get::<ratio>().to_bits() == 0.0f64.to_bits(),
            "vapour"
        );
        Ok(())
    }

    #[test]
    fn test_steam_volume_saturated_fill_sits_on_the_saturation_line() -> anyhow::Result<()> {
        let cv = boiler(373.15, 0.4)?;
        cv.validate()?;
        let (p_sat, liquid, vapour) = saturated_si(373.15)?;
        let state = cv.state();
        anyhow::ensure!(state.region == Region::Saturation, "{:?}", state.region);
        anyhow::ensure!(
            (cv.pressure().get::<pascal>() / p_sat - 1.0).abs() < 1e-9,
            "p"
        );
        let alpha = cv.liquid_fraction().get::<ratio>();
        anyhow::ensure!((alpha - 0.4).abs() < 1e-9, "liquid fraction {alpha}");
        let mass = 1e-3 * (0.4 / liquid.v + 0.6 / vapour.v);
        anyhow::ensure!(
            (cv.mass().get::<kilogram>() / mass - 1.0).abs() < 1e-12,
            "mass"
        );
        anyhow::ensure!(
            (cv.vapour_enthalpy().get::<joule_per_kilogram>() - vapour.h).abs() < 1e-6,
            "h''"
        );
        anyhow::ensure!(
            (cv.liquid_enthalpy().get::<joule_per_kilogram>() - liquid.h).abs() < 1e-6,
            "h'"
        );
        Ok(())
    }

    #[test]
    #[expect(
        clippy::many_single_char_names,
        reason = "thermodynamic symbols p, T, u, v, h, s, x"
    )]
    fn test_steam_volume_heated_at_constant_volume_follows_the_saturation_line()
    -> anyhow::Result<()> {
        // A sealed vessel 70 % full of water at 100 °C: pressure follows
        // p_sat(T) and the liquid swells until it fills the vessel (near 590 K,
        // so as region 1 liquid; a fill below the critical density would dry
        // out).
        let mut cv = boiler(373.15, 0.7)?;
        let (m, v) = (cv.mass().get::<kilogram>(), 1e-3);
        let mut u = cv.internal_energy().get::<joule>();
        let mut last_t = cv.temperature().get::<kelvin>();
        let mut last_alpha = 0.7;
        let mut steps_after_dome = 0;
        for _ in 0..400 {
            u += 5e3;
            cv.set_conserved(m, u, v);
            cv.validate()?;
            let t = cv.temperature().get::<kelvin>();
            anyhow::ensure!(t > last_t, "T fell from {last_t} to {t}");
            last_t = t;
            if cv.state().region != Region::Saturation {
                steps_after_dome += 1;
                if steps_after_dome == 5 {
                    break;
                }
                continue;
            }
            let (p_sat, _, _) = saturated_si(t)?;
            let p = cv.pressure().get::<pascal>();
            anyhow::ensure!((p / p_sat - 1.0).abs() < 1e-8, "p off p_sat at {t} K");
            let alpha = cv.liquid_fraction().get::<ratio>();
            anyhow::ensure!(
                alpha >= last_alpha - 1e-12,
                "liquid shrank to {alpha} at {t} K"
            );
            last_alpha = alpha;
        }
        anyhow::ensure!(steps_after_dome == 5, "never left the dome");
        let state = cv.state();
        anyhow::ensure!(
            state.region == Region::CompressedLiquid,
            "{:?}",
            state.region
        );
        let full = cv.liquid_fraction().get::<ratio>();
        anyhow::ensure!(full.to_bits() == 1.0f64.to_bits(), "liquid fraction {full}");
        Ok(())
    }

    #[test]
    fn test_steam_volume_venting_vapour_flash_cools_the_boiler() -> anyhow::Result<()> {
        let mut cv = boiler(450.0, 0.3)?;
        let mut m = cv.mass().get::<kilogram>();
        let mut u = cv.internal_energy().get::<joule>();
        let mut last = cv.temperature().get::<kelvin>();
        for _ in 0..50 {
            // Saturated vapour leaves; it carries its enthalpy (flow work
            // included).
            let dm = 1e-4;
            u -= dm * cv.vapour_enthalpy().get::<joule_per_kilogram>();
            m -= dm;
            cv.set_conserved(m, u, 1e-3);
            cv.validate()?;
            let t = cv.temperature().get::<kelvin>();
            anyhow::ensure!(t < last, "T rose from {last} to {t}");
            last = t;
        }
        anyhow::ensure!(
            cv.state().region == Region::Saturation,
            "{:?}",
            cv.state().region
        );
        Ok(())
    }

    #[test]
    fn test_steam_volume_validate_reports_impossible_states() -> anyhow::Result<()> {
        let mut cv = boiler(373.15, 0.5)?;
        let (m, u) = (
            cv.mass().get::<kilogram>(),
            cv.internal_energy().get::<joule>(),
        );
        for (mass, energy, volume) in [
            (0.0, u, 1e-3),
            (-m, u, 1e-3),
            (m, f64::NAN, 1e-3),
            (m, u, 0.0),
            (m, -1e9, 1e-3),
        ] {
            cv.set_conserved(mass, energy, volume);
            anyhow::ensure!(
                cv.validate().is_err(),
                "m = {mass}, U = {energy}, V = {volume} accepted"
            );
        }
        // Recovers once the state is valid again.
        cv.set_conserved(m, u, 1e-3);
        cv.validate()?;
        Ok(())
    }

    #[test]
    fn test_steam_volume_rejects_bad_construction() -> anyhow::Result<()> {
        let bad_fraction = boiler(373.15, 1.5);
        anyhow::ensure!(bad_fraction.is_err(), "liquid fraction 1.5 accepted");
        let bad_t = boiler(700.0, 0.5);
        anyhow::ensure!(bad_t.is_err(), "700 K saturated accepted");
        let bad_v = SteamVolume::new(Pressure::new::<pascal>(1e5), kelvin(400.0), litres(0.0));
        anyhow::ensure!(bad_v.is_err(), "zero volume accepted");
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_steam_volume_warm_start_matches_a_fresh_solve(
            t in 280.0f64..640.0,
            fill in 0.0f64..1.0,
            heat in -0.05f64..0.05,
        ) {
            let fail = |e: anyhow::Error| TestCaseError::fail(format!("{e:#}"));
            let mut cv = boiler(t, fill).map_err(fail)?;
            let (m, u) = (cv.mass().get::<kilogram>(), cv.internal_energy().get::<joule>());
            // Move by up to ±5 % of U, then compare with a cold solve.
            let target = u + heat * u.abs().max(1e5);
            cv.set_conserved(m, target, 1e-3);
            prop_assume!(cv.validate().is_ok());
            let cold = state_from_uv_si(target / m, 1e-3 / m, None).map_err(fail)?;
            prop_assert!(cv.state.region == cold.region, "{:?} vs {:?}", cv.state.region, cold.region);
            prop_assert!((cv.state.t - cold.t).abs() < 1e-8 * cold.t, "T {} vs {}", cv.state.t, cold.t);
        }
    }
}
