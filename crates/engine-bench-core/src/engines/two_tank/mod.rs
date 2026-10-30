//! Two-tank demo: a high-pressure tank discharging through a valve into a
//! second tank with an optional heater. Not a catalogue engine; it exercises
//! the runtime end to end.
//!
//! `mod.rs` holds the model and its [`Engine`] implementation; `draw` turns a
//! channel sample into a scene.

mod draw;
mod dynamics;
mod spec;

use anyhow::Context;
use engine_bench_physics::{
    flow::Orifice,
    gas::{ControlVolume, orifice_transfer},
    numerics::ode::{AdaptiveOptions, DormandPrince45, NoEvents, rk4_step},
    thermo::{Air3, Mixture},
    units::{
        Area, Pressure, ThermodynamicTemperature, Time, Volume, area::square_meter, mass::kilogram,
        pressure::pascal, thermodynamic_temperature::kelvin, time::second, volume::cubic_meter,
    },
};

pub(crate) use self::{
    draw::draw,
    spec::{CHANNELS, CONTROLS, channel},
};
use self::{
    dynamics::Dynamics,
    spec::{HEATER, PRESSURE_A, PRESSURE_B, RESET, VALVE},
};
use crate::{
    CoreError,
    engine::{Engine, EngineDescriptor, EngineId},
    spec::{Accuracy, ChannelSpec, ControlId, ControlSpec, StepContext},
};

type Tank = ControlVolume<Air3, 3>;

const DESCRIPTOR: EngineDescriptor = EngineDescriptor {
    id: EngineId::TwoTankDemo,
    name: "Two-tank demo",
    catalogue_number: None,
    summary: "Compressed air flowing through a valve into a second tank with a heater.",
    schema_version: 1,
};

const VOLUME_A: f64 = 0.02;
const VOLUME_B: f64 = 0.05;
const INITIAL_TEMPERATURE: f64 = 300.0;
/// Fully open valve flow area, m², and its discharge coefficient.
const VALVE_AREA: f64 = 2e-5;
const DISCHARGE_COEFFICIENT: f64 = 0.8;
/// Largest real-time RK4 sub-step, s; well below the ~0.1 s fastest time
/// constant.
const REAL_TIME_STEP: f64 = 1e-3;

/// See the module documentation.
#[derive(Debug, Clone)]
pub struct TwoTankDemo {
    a: Tank,
    b: Tank,
    controls: [f64; 5],
    time: f64,
}

impl TwoTankDemo {
    /// The demo at its default controls.
    pub fn new() -> anyhow::Result<Self> {
        let controls = CONTROLS.map(|c| c.default_value());
        let (a, b) = initial_tanks(&controls)?;
        Ok(Self {
            a,
            b,
            controls,
            time: 0.0,
        })
    }

    /// The valve as an orifice, or `None` when closed.
    fn valve(&self) -> anyhow::Result<Option<Orifice>> {
        let opening = self.controls[VALVE];
        if opening <= 0.0 {
            return Ok(None);
        }
        let area = Area::new::<square_meter>(VALVE_AREA * opening);
        Ok(Some(
            Orifice::new(area, DISCHARGE_COEFFICIENT).context("valve")?,
        ))
    }

    pub(crate) fn descriptor_static() -> &'static EngineDescriptor {
        &DESCRIPTOR
    }
}

impl Engine for TwoTankDemo {
    fn descriptor(&self) -> &'static EngineDescriptor {
        &DESCRIPTOR
    }

    fn controls(&self) -> &'static [ControlSpec] {
        &CONTROLS
    }

    fn control_values(&self) -> &[f64] {
        &self.controls
    }

    fn channels(&self) -> &'static [ChannelSpec] {
        &CHANNELS
    }

    fn set_control(&mut self, id: ControlId, value: f64) -> anyhow::Result<()> {
        let index = id.index();
        let spec = CONTROLS
            .get(index)
            .ok_or(CoreError::UnknownControl(index))
            .context("two-tank demo")?;
        let value = spec.validate(value)?;
        if index == RESET {
            return self.reset();
        }
        self.controls[index] = value;
        Ok(())
    }

    /// On error the tanks may hold an intermediate state; reset before
    /// continuing.
    fn advance(&mut self, dt: Time, context: &StepContext) -> anyhow::Result<()> {
        let dt = dt.get::<second>();
        if !(dt > 0.0 && dt.is_finite()) {
            return Err(CoreError::InvalidSetting("time step must be positive"))
                .with_context(|| format!("dt = {dt} s"));
        }
        let valve = self.valve()?;
        let mut dynamics = Dynamics {
            a: &mut self.a,
            b: &mut self.b,
            valve,
            heater: self.controls[HEATER],
        };
        let mut y = dynamics.state();
        let mut t = 0.0;
        match context.accuracy() {
            Accuracy::RealTime => {
                let steps = (dt / REAL_TIME_STEP).ceil();
                let h = dt / steps;
                let mut taken = 0.0;
                while taken < steps {
                    rk4_step(&mut dynamics, &mut t, &mut y, h).context("two-tank RK4 step")?;
                    taken += 1.0;
                }
            }
            Accuracy::Offline { rtol } => {
                let options = AdaptiveOptions::new(rtol, 0.0).context("two-tank accuracy")?;
                DormandPrince45::new(options)
                    .integrate(&mut dynamics, &mut NoEvents, &mut t, &mut y, dt)
                    .context("two-tank adaptive integration")?;
            }
        }
        dynamics.unpack(&y);
        self.a.validate().context("tank A")?;
        self.b.validate().context("tank B")?;
        self.time += dt;
        Ok(())
    }

    fn time(&self) -> Time {
        Time::new::<second>(self.time)
    }

    fn sample(&self, out: &mut [f64]) -> anyhow::Result<()> {
        let actual = out.len();
        let out: &mut [f64; channel::COUNT] = out
            .try_into()
            .map_err(|_| CoreError::SampleLength {
                expected: channel::COUNT,
                actual,
            })
            .context("two-tank sample")?;
        let mass_flow = match self.valve()? {
            Some(valve) => orifice_transfer(valve, &self.a, &self.b).mass_flow(),
            None => 0.0,
        };
        out[channel::P_A] = self.a.pressure().get::<pascal>();
        out[channel::P_B] = self.b.pressure().get::<pascal>();
        out[channel::T_A] = self.a.temperature().get::<kelvin>();
        out[channel::T_B] = self.b.temperature().get::<kelvin>();
        out[channel::M_A] = self.a.mass().get::<kilogram>();
        out[channel::M_B] = self.b.mass().get::<kilogram>();
        out[channel::MASS_FLOW] = mass_flow;
        out[channel::V_A] = self.a.volume().get::<cubic_meter>();
        out[channel::V_B] = self.b.volume().get::<cubic_meter>();
        out[channel::VALVE] = self.controls[VALVE];
        Ok(())
    }

    fn reset(&mut self) -> anyhow::Result<()> {
        (self.a, self.b) = initial_tanks(&self.controls)?;
        self.time = 0.0;
        Ok(())
    }
}

/// Both tanks at the initial pressures in `controls` and the initial
/// temperature.
fn initial_tanks(controls: &[f64; 5]) -> anyhow::Result<(Tank, Tank)> {
    let tank = |p: f64, volume: f64| -> anyhow::Result<Tank> {
        ControlVolume::new(
            Mixture::dry_air()?,
            Pressure::new::<pascal>(p),
            ThermodynamicTemperature::new::<kelvin>(INITIAL_TEMPERATURE),
            Volume::new::<cubic_meter>(volume),
        )
    };
    Ok((
        tank(controls[PRESSURE_A], VOLUME_A).context("tank A")?,
        tank(controls[PRESSURE_B], VOLUME_B).context("tank B")?,
    ))
}

#[cfg(test)]
mod tests {
    use engine_bench_physics::units::energy::joule;

    use super::*;

    fn real_time() -> StepContext {
        StepContext::real_time()
    }

    fn seconds(value: f64) -> Time {
        Time::new::<second>(value)
    }

    fn sample(demo: &TwoTankDemo) -> anyhow::Result<[f64; channel::COUNT]> {
        let mut out = [0.0; channel::COUNT];
        demo.sample(&mut out)?;
        Ok(out)
    }

    /// (total mass in kg, total internal energy in J).
    fn totals(demo: &TwoTankDemo) -> (f64, f64) {
        (
            demo.a.mass().get::<kilogram>() + demo.b.mass().get::<kilogram>(),
            demo.a.internal_energy().get::<joule>() + demo.b.internal_energy().get::<joule>(),
        )
    }

    #[test]
    fn test_two_tank_starts_at_default_controls() -> anyhow::Result<()> {
        let demo = TwoTankDemo::new()?;
        let defaults: Vec<f64> = CONTROLS.iter().map(ControlSpec::default_value).collect();
        anyhow::ensure!(demo.control_values() == defaults.as_slice(), "controls");
        let s = sample(&demo)?;
        anyhow::ensure!(
            (s[channel::P_A] / 8e5 - 1.0).abs() < 1e-12,
            "p_a {}",
            s[channel::P_A]
        );
        anyhow::ensure!(
            (s[channel::P_B] / 1e5 - 1.0).abs() < 1e-12,
            "p_b {}",
            s[channel::P_B]
        );
        anyhow::ensure!(
            (s[channel::T_A] - 300.0).abs() < 1e-9,
            "T_a {}",
            s[channel::T_A]
        );
        anyhow::ensure!(s[channel::MASS_FLOW] > 0.0, "no flow with valve open");
        anyhow::ensure!(
            s[channel::V_A] == VOLUME_A && s[channel::V_B] == VOLUME_B,
            "volumes {s:?}"
        );
        anyhow::ensure!(demo.time().get::<second>() == 0.0, "time");
        Ok(())
    }

    #[test]
    fn test_two_tank_conserves_mass_and_energy_in_both_modes() -> anyhow::Result<()> {
        for context in [real_time(), StepContext::offline(1e-9)?] {
            let mut demo = TwoTankDemo::new()?;
            demo.set_control(ControlId::new(VALVE), 1.0)?;
            let (mass0, energy0) = totals(&demo);
            for _ in 0..4 {
                demo.advance(seconds(0.5), &context)?;
            }
            let (mass1, energy1) = totals(&demo);
            anyhow::ensure!(
                (mass1 / mass0 - 1.0).abs() < 1e-12,
                "{context:?}: mass {mass0} → {mass1}"
            );
            anyhow::ensure!(
                (energy1 / energy0 - 1.0).abs() < 1e-12,
                "{context:?}: energy {energy0} → {energy1}"
            );
            anyhow::ensure!((demo.time().get::<second>() - 2.0).abs() < 1e-12, "time");
        }
        Ok(())
    }

    #[test]
    fn test_two_tank_heater_adds_its_energy() -> anyhow::Result<()> {
        let mut demo = TwoTankDemo::new()?;
        demo.set_control(ControlId::new(HEATER), 1500.0)?;
        let (_, energy0) = totals(&demo);
        demo.advance(seconds(2.0), &real_time())?;
        let (_, energy1) = totals(&demo);
        anyhow::ensure!(
            ((energy1 - energy0) / 3000.0 - 1.0).abs() < 1e-9,
            "ΔU {}",
            energy1 - energy0
        );
        Ok(())
    }

    #[test]
    fn test_two_tank_closed_valve_stops_flow() -> anyhow::Result<()> {
        let mut demo = TwoTankDemo::new()?;
        demo.set_control(ControlId::new(VALVE), 0.0)?;
        let before = sample(&demo)?;
        demo.advance(seconds(1.0), &real_time())?;
        let after = sample(&demo)?;
        anyhow::ensure!(
            after[channel::MASS_FLOW] == 0.0,
            "flow {}",
            after[channel::MASS_FLOW]
        );
        // No flow means no change at all, bit for bit.
        anyhow::ensure!(
            after[channel::M_A].to_bits() == before[channel::M_A].to_bits(),
            "mass moved"
        );
        Ok(())
    }

    #[test]
    fn test_two_tank_pressures_equalise() -> anyhow::Result<()> {
        let mut demo = TwoTankDemo::new()?;
        demo.set_control(ControlId::new(VALVE), 1.0)?;
        demo.advance(seconds(120.0), &real_time())?;
        let s = sample(&demo)?;
        anyhow::ensure!(
            (s[channel::P_A] / s[channel::P_B] - 1.0).abs() < 1e-4,
            "{s:?}"
        );
        Ok(())
    }

    #[test]
    fn test_two_tank_real_time_agrees_with_offline() -> anyhow::Result<()> {
        let mut fast = TwoTankDemo::new()?;
        let mut accurate = TwoTankDemo::new()?;
        fast.advance(seconds(3.0), &real_time())?;
        accurate.advance(seconds(3.0), &StepContext::offline(1e-10)?)?;
        let (f, a) = (sample(&fast)?, sample(&accurate)?);
        for i in [channel::P_A, channel::P_B, channel::T_A, channel::T_B] {
            anyhow::ensure!(
                (f[i] / a[i] - 1.0).abs() < 1e-6,
                "channel {i}: {} vs {}",
                f[i],
                a[i]
            );
        }
        Ok(())
    }

    #[test]
    fn test_two_tank_rejects_invalid_controls_and_steps() -> anyhow::Result<()> {
        let mut demo = TwoTankDemo::new()?;
        anyhow::ensure!(
            demo.set_control(ControlId::new(VALVE), 1.5).is_err(),
            "valve 1.5"
        );
        let Err(e) = demo.set_control(ControlId::new(99), 0.0) else {
            anyhow::bail!("accepted control 99");
        };
        anyhow::ensure!(
            matches!(
                e.downcast_ref::<CoreError>(),
                Some(CoreError::UnknownControl(99))
            ),
            "{e:#}"
        );
        for dt in [0.0, -1.0, f64::NAN] {
            anyhow::ensure!(demo.advance(seconds(dt), &real_time()).is_err(), "dt {dt}");
        }
        anyhow::ensure!(demo.sample(&mut [0.0; 3]).is_err(), "short sample buffer");
        Ok(())
    }

    #[test]
    fn test_two_tank_reset_applies_initial_pressures() -> anyhow::Result<()> {
        let mut demo = TwoTankDemo::new()?;
        demo.advance(seconds(1.0), &real_time())?;
        demo.set_control(ControlId::new(PRESSURE_A), 1.5e6)?;
        // Initial pressures take effect on reset, not immediately.
        anyhow::ensure!(sample(&demo)?[channel::P_A] < 8e5, "applied early");
        demo.set_control(ControlId::new(RESET), 1.0)?;
        let s = sample(&demo)?;
        anyhow::ensure!(
            (s[channel::P_A] / 1.5e6 - 1.0).abs() < 1e-12,
            "p_a {}",
            s[channel::P_A]
        );
        anyhow::ensure!(demo.time().get::<second>() == 0.0, "time not reset");
        Ok(())
    }

    #[test]
    fn test_two_tank_survives_extreme_control_settings() -> anyhow::Result<()> {
        // Highest A, lowest B, valve open: tank A cools the most.
        let mut demo = TwoTankDemo::new()?;
        demo.set_control(ControlId::new(PRESSURE_A), 2e6)?;
        demo.set_control(ControlId::new(PRESSURE_B), 5e4)?;
        demo.set_control(ControlId::new(VALVE), 1.0)?;
        demo.reset()?;
        demo.advance(seconds(200.0), &real_time())?;
        let t_a = sample(&demo)?[channel::T_A];
        anyhow::ensure!(t_a > 200.0, "tank A cooled to {t_a} K");
        Ok(())
    }
}
