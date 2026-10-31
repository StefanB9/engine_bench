//! Cold-gas thruster (catalogue engine 2): a high-pressure tank feeding a
//! regulated plenum, a solenoid valve and a nozzle.
//!
//! Real fluids throughout: reference equations of state in every volume,
//! real-fluid orifices and frozen-vapour nozzle flow, tabulated in real time.
//!
//! `mod.rs` holds the model and its [`Engine`] implementation; `dynamics` the
//! state and flows, `regulator` and `valve` their components, `hardware` the
//! sizing, `draw` the scene.

mod draw;
mod dynamics;
mod fluids;
mod hardware;
mod regulator;
mod spec;
mod valve;

use anyhow::Context;
use engine_bench_physics::{
    flow::{Nozzle, NozzleRegime},
    fluid::PureFluid,
    gas::eos::Reference,
    numerics::ode::{AdaptiveOptions, DormandPrince45, NoEvents, rk4_step},
    thermo::Mixture,
    units::{
        Area, Pressure, ThermodynamicTemperature, Time, Volume as VolumeUnit, area::square_meter,
        force::newton, mass::kilogram, mass_rate::kilogram_per_second, pressure::pascal,
        thermodynamic_temperature::kelvin, time::second, volume::cubic_meter,
    },
};

pub(crate) use self::{
    draw::draw,
    spec::{CHANNELS, CONTROLS, channel},
};
use self::{
    dynamics::{Dynamics, Plant, STATES, Settings, VACUUM_FLOOR, Volume},
    fluids::Fluids,
    hardware::{OUTLET_VOLUME, PLENUM_VOLUME, Tank, circle_area},
    spec::{
        AMBIENT_PRESSURE, CONTROL_COUNT, EXIT_DIAMETER, FIRE, PROPELLANT, PROPELLANTS, PROPERTIES,
        PULSE_PERIOD, PULSE_WIDTH, RESET, SETPOINT, TANK_PRESSURE, TANK_VOLUME, TEMPERATURE,
        THROAT_DIAMETER, WALL_INSIDE, WALL_OUTSIDE,
    },
    valve::FireMode,
};
use crate::{
    CoreError,
    engine::{Engine, EngineDescriptor, EngineId},
    spec::{Accuracy, ChannelSpec, ControlId, ControlSpec, StepContext, choice_index},
};

const DESCRIPTOR: EngineDescriptor = EngineDescriptor {
    id: EngineId::ColdGasThruster,
    name: "Cold-gas thruster",
    catalogue_number: Some(2),
    summary: "Compressed gas from a tank through a regulator, a solenoid valve and a nozzle.",
    schema_version: 1,
};

/// Longest real-time RK4 step, s: the poppet's fastest mode decays at
/// ~6000/s (RK4 is stable to 2.8/6000 s).
const REAL_TIME_STEP: f64 = 2.5e-4;
/// Real-time steps are also at most this fraction of the valve outlet's
/// filling time: near equal pressures the valve's flow is very sensitive to
/// the outlet pressure, and a fast gas (helium) or a large valve makes the
/// small outlet stiff.
const OUTLET_STEP_FRACTION: f64 = 0.2;

/// See the module documentation.
#[derive(Debug, Clone)]
pub struct ColdGasThruster {
    controls: [f64; CONTROL_COUNT],
    /// Propellant the plant was filled with at the last reset.
    fluid: PureFluid,
    plant: Plant,
    /// The plant's tank as sized at the last reset.
    tank: Tank,
    /// State scales for the integrators, set at reset.
    scale: [f64; STATES],
    time: f64,
    /// The propellant's tables, built on the first real-time step.
    tables: Option<Fluids>,
}

impl ColdGasThruster {
    /// The thruster at its default controls.
    pub fn new() -> anyhow::Result<Self> {
        let controls = CONTROLS.map(|c| c.default_value());
        let fluid = propellant(&controls);
        let (plant, tank) = initial_plant(&controls, fluid)?;
        Ok(Self {
            controls,
            fluid,
            scale: plant.scales(),
            plant,
            tank,
            time: 0.0,
            tables: None,
        })
    }

    pub(crate) fn descriptor_static() -> &'static EngineDescriptor {
        &DESCRIPTOR
    }

    /// The settings the flows run under now.
    fn settings(&self) -> anyhow::Result<Settings> {
        let c = &self.controls;
        let throat = c[THROAT_DIAMETER];
        // An exit narrower than the throat is a plain convergent nozzle.
        let exit = c[EXIT_DIAMETER].max(throat);
        let nozzle = Nozzle::new(
            Area::new::<square_meter>(circle_area(throat)),
            Area::new::<square_meter>(circle_area(exit)),
        )
        .context("nozzle")?;
        Ok(Settings {
            tank: self.tank,
            setpoint: c[SETPOINT],
            ambient_pressure: c[AMBIENT_PRESSURE],
            ambient_temperature: c[TEMPERATURE],
            wall_inside: c[WALL_INSIDE],
            wall_outside: c[WALL_OUTSIDE],
            fire: FireMode::from_index(choice_index(c[FIRE], 3).unwrap_or(0)),
            pulse_width: c[PULSE_WIDTH],
            pulse_period: c[PULSE_PERIOD],
            nozzle,
            throat_area: circle_area(throat),
            exit_area: circle_area(exit),
        })
    }

    /// Whether the fluid-properties control asks for the exact equations.
    fn exact(&self) -> bool {
        choice_index(self.controls[PROPERTIES], 2) == Some(1)
    }

    /// The fluids to step with: the propellant's tables (built on first use),
    /// or its exact equations.
    fn fluids(&mut self) -> anyhow::Result<Fluids> {
        if self.exact() {
            return Ok(Fluids::exact(self.fluid));
        }
        if let Some(tables) = self.tables {
            return Ok(tables);
        }
        let tables = Fluids::tabulated(self.fluid).context("propellant tables")?;
        self.tables = Some(tables);
        Ok(tables)
    }
}

/// The real-time step for `plant`, s: at most [`REAL_TIME_STEP`], and a
/// fraction of the valve outlet's filling time through the fully open valve,
/// V / (0.56·Cd·A·√(RT)) (choked flux ≈ 0.56·ρ·√(RT) for γ from 1.3 to 5/3,
/// independent of pressure).
fn real_time_step(settings: &Settings, plant: &Plant) -> f64 {
    let r = plant
        .plenum
        .mixture()
        .gas_constant()
        .get::<engine_bench_physics::units::specific_heat_capacity::joule_per_kilogram_kelvin>(
    );
    let t = plant.plenum.temperature().get::<kelvin>();
    let valve = valve::DISCHARGE_COEFFICIENT * valve::AREA_PER_THROAT * settings.throat_area;
    let filling = OUTLET_VOLUME / (0.56 * valve * (r * t).sqrt());
    REAL_TIME_STEP.min(OUTLET_STEP_FRACTION * filling)
}

/// The propellant chosen in `controls`.
fn propellant(controls: &[f64; CONTROL_COUNT]) -> PureFluid {
    match choice_index(controls[PROPELLANT], PROPELLANTS.len()) {
        Some(1) => PureFluid::Argon,
        Some(2) => PureFluid::Helium,
        _ => PureFluid::Nitrogen,
    }
}

/// The plant as filled at reset: the tank at its pressure, the plenum locked
/// up at the setpoint (at tank pressure without a regulator), the valve outlet
/// at ambient pressure, all at the ambient temperature.
fn initial_plant(
    controls: &[f64; CONTROL_COUNT],
    fluid: PureFluid,
) -> anyhow::Result<(Plant, Tank)> {
    let (p_tank, t) = (controls[TANK_PRESSURE], controls[TEMPERATURE]);
    let tank = Tank::sized(controls[TANK_VOLUME], p_tank);
    let setpoint = controls[SETPOINT];
    let p_plenum = if setpoint > 0.0 {
        regulator::lock_up(p_tank, setpoint).min(p_tank)
    } else {
        p_tank
    };
    let p_outlet = controls[AMBIENT_PRESSURE].max(VACUUM_FLOOR);
    let volume = |p: f64, size: f64| -> anyhow::Result<Volume> {
        let mut fractions = [0.0; 4];
        let index = match fluid {
            PureFluid::Nitrogen => 0,
            PureFluid::Argon => 1,
            PureFluid::Helium => 2,
            PureFluid::CarbonDioxide => 3,
        };
        fractions[index] = 1.0;
        Volume::with_eos(
            Reference::new(fluid)?,
            Mixture::from_mass_fractions(fractions)?,
            Pressure::new::<pascal>(p),
            ThermodynamicTemperature::new::<kelvin>(t),
            VolumeUnit::new::<cubic_meter>(size),
        )
    };
    let plant = Plant {
        tank: volume(p_tank, tank.volume).context("tank")?,
        wall: t,
        plenum: volume(p_plenum, PLENUM_VOLUME).context("plenum")?,
        outlet: volume(p_outlet, OUTLET_VOLUME).context("valve outlet")?,
        lift: if setpoint > 0.0 {
            0.0
        } else {
            regulator::FULL_LIFT
        },
        lift_rate: 0.0,
        valve: 0.0,
        impulse: 0.0,
    };
    Ok((plant, tank))
}

/// Index of a nozzle regime, for the regime channel.
fn regime_index(regime: NozzleRegime) -> f64 {
    match regime {
        NozzleRegime::NoFlow => 0.0,
        NozzleRegime::Subsonic => 1.0,
        NozzleRegime::ShockInNozzle { .. } => 2.0,
        NozzleRegime::ShockAtExit => 3.0,
        NozzleRegime::OverExpanded => 4.0,
        NozzleRegime::Ideal => 5.0,
        NozzleRegime::UnderExpanded => 6.0,
    }
}

impl Engine for ColdGasThruster {
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
            .context("cold-gas thruster")?;
        let value = spec.validate(value)?;
        if index == RESET {
            return self.reset();
        }
        self.controls[index] = value;
        Ok(())
    }

    /// On error the plant may hold an intermediate state; reset before
    /// continuing.
    fn advance(&mut self, dt: Time, context: &StepContext) -> anyhow::Result<()> {
        let dt = dt.get::<second>();
        if !(dt > 0.0 && dt.is_finite()) {
            return Err(CoreError::InvalidSetting("time step must be positive"))
                .with_context(|| format!("dt = {dt} s"));
        }
        let settings = self.settings()?;
        let fluids = self.fluids()?;
        let start = self.time;
        let mut dynamics = Dynamics {
            plant: &mut self.plant,
            scale: &self.scale,
            settings: &settings,
            fluid: &fluids.equilibrium,
            frozen: &fluids.frozen,
            start,
            // Real time integrates the impulse by the trapezoidal rule, so the
            // stages need only the mass flow: no shock is solved four times
            // a step.
            impulse_in_ode: matches!(context.accuracy(), Accuracy::Offline { .. }),
            failure: None,
        };
        let mut y = dynamics.state();
        let mut t = 0.0;
        match context.accuracy() {
            Accuracy::RealTime => {
                let steps = (dt / real_time_step(&settings, dynamics.plant)).ceil();
                let h = dt / steps;
                let mut taken = 0.0;
                let thrust = |plant: &Plant| {
                    plant
                        .thrust(&settings, &fluids.equilibrium, &fluids.frozen)
                        .context("thrust")
                };
                let mut before = thrust(dynamics.plant)?;
                while taken < steps {
                    rk4_step(&mut dynamics, &mut t, &mut y, h).context("thruster RK4 step")?;
                    if let Some(e) = dynamics.failure.take() {
                        return Err(e).context("thruster RK4 step");
                    }
                    dynamics.unpack(&y);
                    dynamics.plant.settle();
                    dynamics.plant.validate()?;
                    let after = thrust(dynamics.plant)?;
                    dynamics.plant.impulse += 0.5 * h * (before + after);
                    before = after;
                    y = dynamics.state();
                    taken += 1.0;
                }
            }
            Accuracy::Offline { rtol } => {
                // States are scaled to order one: the absolute tolerance is the
                // relative one.
                let options = AdaptiveOptions::new(rtol, rtol).context("thruster accuracy")?;
                let result = DormandPrince45::new(options).integrate(
                    &mut dynamics,
                    &mut NoEvents,
                    &mut t,
                    &mut y,
                    dt,
                );
                // A trial stage can land on an impossible state; the
                // integrator rejects it and shrinks the step. Only when it
                // gives up is the derivative's failure the cause.
                if let Err(e) = result {
                    let cause = dynamics.failure.take();
                    return Err(match cause {
                        Some(cause) => cause.context(format!("{e:#}")),
                        None => e,
                    })
                    .context("thruster adaptive integration");
                }
                dynamics.unpack(&y);
                dynamics.plant.settle();
                dynamics.plant.validate()?;
            }
        }
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
            .context("cold-gas sample")?;
        let settings = self.settings()?;
        // Before the first step the tables may not exist yet: sample exactly.
        let fluids = match self.tables {
            Some(tables) if !self.exact() => tables,
            _ => Fluids::exact(self.fluid),
        };
        let p = &self.plant;
        let flows = p
            .flows(&settings, &fluids.equilibrium, &fluids.frozen, true)
            .context("cold-gas sample")?;
        let nozzle = flows.nozzle;
        let (p_exit, t_exit) = (
            nozzle.exit_pressure.get::<pascal>(),
            nozzle.exit_temperature,
        );
        let mass_flow = nozzle.mass_flow.get::<kilogram_per_second>();
        // p / p_sat at the exit: how far the frozen vapour is supersaturated.
        let supersaturation = if mass_flow > 0.0 {
            self.fluid
                .vapour_pressure(t_exit)
                .map_or(0.0, |p_sat| p_exit / p_sat.get::<pascal>())
        } else {
            0.0
        };
        out[channel::P_TANK] = p.tank.pressure().get::<pascal>();
        out[channel::T_TANK] = p.tank.temperature().get::<kelvin>();
        out[channel::T_WALL] = p.wall;
        out[channel::M_TANK] = p.tank.mass().get::<kilogram>();
        out[channel::P_PLENUM] = p.plenum.pressure().get::<pascal>();
        out[channel::T_PLENUM] = p.plenum.temperature().get::<kelvin>();
        out[channel::P_OUTLET] = p.outlet.pressure().get::<pascal>();
        out[channel::LIFT] = settings.effective_lift(p.lift) / regulator::FULL_LIFT;
        out[channel::VALVE] = p.valve;
        out[channel::MASS_FLOW] = mass_flow;
        out[channel::THRUST] = nozzle.thrust.get::<newton>();
        out[channel::ISP] = nozzle.specific_impulse().map_or(0.0, |t| t.get::<second>());
        out[channel::IMPULSE] = p.impulse;
        out[channel::EXIT_MACH] = nozzle.exit_mach;
        out[channel::EXIT_PRESSURE] = p_exit;
        out[channel::EXIT_TEMPERATURE] = t_exit.get::<kelvin>();
        out[channel::SUPERSATURATION] = supersaturation;
        out[channel::SEPARATED] = f64::from(u8::from(nozzle.flow_separation));
        out[channel::REGIME] = regime_index(nozzle.regime);
        out[channel::TANK_VOLUME] = self.tank.volume;
        out[channel::THROAT] = self.controls[THROAT_DIAMETER];
        out[channel::EXIT] = self.controls[EXIT_DIAMETER].max(self.controls[THROAT_DIAMETER]);
        Ok(())
    }

    fn reset(&mut self) -> anyhow::Result<()> {
        let fluid = propellant(&self.controls);
        (self.plant, self.tank) = initial_plant(&self.controls, fluid)?;
        self.scale = self.plant.scales();
        if fluid != self.fluid {
            self.tables = None;
        }
        self.fluid = fluid;
        self.time = 0.0;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
