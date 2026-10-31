//! The thruster's state, its flows and the ODE over them.

use anyhow::Context;
use engine_bench_physics::{
    flow::{Nozzle, NozzleOperatingPoint, NozzleRegime, Orifice, real::RealFluid},
    gas::{ControlVolume, GasPort, RealPort, Transfer, eos::Reference, real_orifice_transfer_in},
    numerics::ode::OdeSystem,
    thermo::ColdGas,
    units::{
        Area, Force, MassRate, Pressure, ThermodynamicTemperature, Velocity, area::square_meter,
        available_energy::joule_per_kilogram, energy::joule, force::newton, mass::kilogram,
        mass_rate::kilogram_per_second, pressure::pascal,
        specific_heat_capacity::joule_per_kilogram_kelvin, thermodynamic_temperature::kelvin,
        velocity::meter_per_second,
    },
};

use super::{
    hardware::{OUTLET_VOLUME, PLENUM_VOLUME, Tank},
    regulator, valve,
    valve::FireMode,
};

/// A propellant volume: one pure fluid of the cold-gas set.
pub(super) type Volume = ControlVolume<ColdGas, 4, Reference<ColdGas, 4>>;

/// ODE states: tank, plenum and valve outlet (4 species masses and internal
/// energy each), wall temperature, poppet lift and velocity, valve position,
/// total impulse.
pub(super) const STATES: usize = 20;
const TANK: usize = 0;
const WALL: usize = 5;
const PLENUM: usize = 6;
const OUTLET: usize = 11;
const LIFT: usize = 16;
const LIFT_RATE: usize = 17;
const VALVE: usize = 18;
const IMPULSE: usize = 19;

/// The lowest back pressure the nozzle flows against, Pa. In vacuum the valve
/// outlet drains to it after the valve shuts instead of towards an empty
/// volume, whose draining tail would expand to ever lower exit pressures;
/// what stays behind (~1e-8 kg) is negligible. Thrust is still taken against
/// the true ambient pressure.
pub(super) const VACUUM_FLOOR: f64 = 1e3;
/// Heat transfer between the plenum or valve outlet gas and their metal
/// bodies (at ambient temperature), W/(m²·K): forced convection.
const BODY_HEAT_TRANSFER: f64 = 200.0;
/// Shortest thermal relaxation time of that gas, s. Nearly empty (a valve
/// outlet drained to vacuum), the gas would follow the wall arbitrarily fast;
/// it then does so in this time, which changes nothing of consequence and
/// keeps the equations from turning stiff.
const FASTEST_RELAXATION: f64 = 1e-3;

/// Everything that changes in time.
#[derive(Debug, Clone)]
pub(super) struct Plant {
    pub(super) tank: Volume,
    /// Wall temperature, K.
    pub(super) wall: f64,
    pub(super) plenum: Volume,
    pub(super) outlet: Volume,
    /// Poppet lift, m, and its rate, m/s.
    pub(super) lift: f64,
    pub(super) lift_rate: f64,
    /// Valve position, 0 shut to 1 open.
    pub(super) valve: f64,
    /// Total impulse, N·s.
    pub(super) impulse: f64,
}

/// Settings the flows depend on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Settings {
    pub(super) tank: Tank,
    /// Pa; 0 holds the regulator open.
    pub(super) setpoint: f64,
    /// Pa.
    pub(super) ambient_pressure: f64,
    /// K.
    pub(super) ambient_temperature: f64,
    /// W/(m²·K).
    pub(super) wall_inside: f64,
    pub(super) wall_outside: f64,
    pub(super) fire: FireMode,
    /// s.
    pub(super) pulse_width: f64,
    pub(super) pulse_period: f64,
    pub(super) nozzle: Nozzle,
    /// m².
    pub(super) throat_area: f64,
    /// m².
    pub(super) exit_area: f64,
}

/// The flows at one instant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Flows {
    /// Tank to plenum.
    pub(super) regulator: Transfer<4>,
    /// Plenum to valve outlet.
    pub(super) valve: Transfer<4>,
    /// Out of the valve outlet through the nozzle.
    pub(super) nozzle: NozzleOperatingPoint,
    /// Enthalpy leaving through the nozzle, W.
    pub(super) nozzle_energy: f64,
    /// Heat from the wall into the gas, W, and from outside into the wall.
    pub(super) wall_to_gas: f64,
    pub(super) outside_to_wall: f64,
    /// Heat from their bodies into the plenum and valve outlet gas, W.
    pub(super) plenum_heat: f64,
    pub(super) outlet_heat: f64,
}

impl Settings {
    /// Lift the regulator acts with: fully open without a setpoint.
    pub(super) fn effective_lift(&self, lift: f64) -> f64 {
        if self.setpoint > 0.0 {
            lift.clamp(0.0, regulator::FULL_LIFT)
        } else {
            regulator::FULL_LIFT
        }
    }
}

impl Plant {
    /// The flows in `fluid` (equilibrium) and `frozen` (nozzle) under
    /// `settings`; the nozzle's full operating point with `thrust`, else its
    /// mass flow alone (no shock solved, thrust and exit state left at rest).
    pub(super) fn flows<E: RealFluid, F: RealFluid>(
        &self,
        settings: &Settings,
        fluid: &E,
        frozen: &F,
        thrust: bool,
    ) -> anyhow::Result<Flows> {
        let regulator_area = regulator::flow_area(settings.effective_lift(self.lift));
        let regulator = transfer(
            fluid,
            regulator_area,
            regulator::DISCHARGE_COEFFICIENT,
            &self.tank,
            &self.plenum,
        )
        .context("regulator")?;
        let valve_area = self.valve.clamp(0.0, 1.0) * valve::AREA_PER_THROAT * settings.throat_area;
        let valve = transfer(
            fluid,
            valve_area,
            valve::DISCHARGE_COEFFICIENT,
            &self.plenum,
            &self.outlet,
        )
        .context("valve")?;
        let back = settings.ambient_pressure.max(VACUUM_FLOOR);
        let p_outlet = self.outlet.pressure().get::<pascal>();
        let (nozzle, nozzle_energy) = if p_outlet > back {
            let chamber = self
                .outlet
                .port_stagnation_in(fluid)
                .context("valve outlet")?;
            let back_pressure = Pressure::new::<pascal>(back);
            let mut point = if thrust {
                settings
                    .nozzle
                    .real_operating_point(frozen, &chamber, back_pressure)
                    .context("nozzle")?
            } else {
                NozzleOperatingPoint {
                    mass_flow: settings
                        .nozzle
                        .real_mass_flow(frozen, &chamber, back_pressure)
                        .context("nozzle")?,
                    ..no_flow(p_outlet, self.outlet.temperature())
                }
            };
            // An attached supersonic exit does not depend on the back pressure:
            // only the pressure term of the thrust does.
            let attached = matches!(
                point.regime,
                NozzleRegime::UnderExpanded | NozzleRegime::Ideal | NozzleRegime::OverExpanded
            );
            if thrust && attached {
                point.thrust +=
                    Force::new::<newton>((back - settings.ambient_pressure) * settings.exit_area);
            }
            let h = self.outlet.specific_enthalpy().get::<joule_per_kilogram>();
            (point, point.mass_flow.get::<kilogram_per_second>() * h)
        } else {
            (no_flow(p_outlet, self.outlet.temperature()), 0.0)
        };
        let t_gas = self.tank.temperature().get::<kelvin>();
        let area = settings.tank.area;
        Ok(Flows {
            regulator,
            valve,
            nozzle,
            nozzle_energy,
            wall_to_gas: settings.wall_inside * area * (self.wall - t_gas),
            outside_to_wall: settings.wall_outside
                * area
                * (settings.ambient_temperature - self.wall),
            plenum_heat: body_heat(&self.plenum, PLENUM_VOLUME, settings.ambient_temperature),
            outlet_heat: body_heat(&self.outlet, OUTLET_VOLUME, settings.ambient_temperature),
        })
    }

    /// The ODE state.
    pub(super) fn state(&self) -> [f64; STATES] {
        let mut y = [0.0; STATES];
        for (offset, volume) in [
            (TANK, &self.tank),
            (PLENUM, &self.plenum),
            (OUTLET, &self.outlet),
        ] {
            y[offset..offset + 4].copy_from_slice(volume.species_masses());
            y[offset + 4] = volume.internal_energy().get::<joule>();
        }
        y[WALL] = self.wall;
        y[LIFT] = self.lift;
        y[LIFT_RATE] = self.lift_rate;
        y[VALVE] = self.valve;
        y[IMPULSE] = self.impulse;
        y
    }

    /// Thrust now, N.
    pub(super) fn thrust<E: RealFluid, F: RealFluid>(
        &self,
        settings: &Settings,
        fluid: &E,
        frozen: &F,
    ) -> anyhow::Result<f64> {
        Ok(self
            .flows(settings, fluid, frozen, true)?
            .nozzle
            .thrust
            .get::<newton>())
    }

    /// Takes the ODE state back.
    pub(super) fn unpack(&mut self, y: &[f64; STATES], tank_volume: f64) {
        for (offset, volume, size) in [
            (TANK, &mut self.tank, tank_volume),
            (PLENUM, &mut self.plenum, PLENUM_VOLUME),
            (OUTLET, &mut self.outlet, OUTLET_VOLUME),
        ] {
            let masses = [y[offset], y[offset + 1], y[offset + 2], y[offset + 3]];
            volume.set_conserved(&masses, y[offset + 4], size);
        }
        self.wall = y[WALL];
        self.lift = y[LIFT];
        self.lift_rate = y[LIFT_RATE];
        self.valve = y[VALVE];
        self.impulse = y[IMPULSE];
    }

    /// Typical magnitude of every state, fixed at reset: the integrators work
    /// on states divided by these, so one tolerance fits masses of 1e-10 kg and
    /// energies of 1e5 J alike, and states starting at zero (the valve, the
    /// lift, the impulse, absent species) get an absolute tolerance.
    pub(super) fn scales(&self) -> [f64; STATES] {
        // Specific internal energy scale, J/kg.
        const ENERGY_PER_MASS: f64 = 1e5;
        let tank = self.tank.mass().get::<kilogram>();
        let plenum = self.plenum.mass().get::<kilogram>();
        // The outlet fills to about plenum density.
        let outlet = plenum * OUTLET_VOLUME / PLENUM_VOLUME;
        let mut scale = [1.0; STATES];
        for (offset, mass, volume) in [
            (TANK, tank, &self.tank),
            (PLENUM, plenum, &self.plenum),
            (OUTLET, outlet, &self.outlet),
        ] {
            scale[offset..offset + 4].fill(mass);
            scale[offset + 4] =
                (mass * ENERGY_PER_MASS).max(volume.internal_energy().get::<joule>().abs());
        }
        scale[WALL] = self.wall;
        scale[LIFT] = regulator::FULL_LIFT;
        // Full lift in a millisecond.
        scale[LIFT_RATE] = regulator::FULL_LIFT * 1e3;
        scale
    }

    /// Puts the poppet and the valve back within their stops after a step.
    pub(super) fn settle(&mut self) {
        (self.lift, self.lift_rate) = regulator::at_stops(self.lift, self.lift_rate);
        self.valve = self.valve.clamp(0.0, 1.0);
    }

    /// Checks every volume.
    pub(super) fn validate(&self) -> anyhow::Result<()> {
        self.tank.validate().context("tank")?;
        self.plenum.validate().context("plenum")?;
        self.outlet.validate().context("valve outlet")
    }
}

/// Heat from a body at `t_body` K into the gas of `volume` (`size` m³), W.
fn body_heat(volume: &Volume, size: f64, t_body: f64) -> f64 {
    let area = 6.0 * size.powf(2.0 / 3.0);
    // c_v is at least 1.5·R (monatomic): a lower bound on the heat capacity.
    let r = volume
        .mixture()
        .gas_constant()
        .get::<joule_per_kilogram_kelvin>();
    let capacity = 1.5 * r * volume.mass().get::<kilogram>();
    let conductance = (BODY_HEAT_TRANSFER * area).min(capacity / FASTEST_RELAXATION);
    conductance * (t_body - volume.temperature().get::<kelvin>())
}

/// Flow through an orifice of `area` m² between two volumes, none when shut.
fn transfer<F: RealFluid>(
    fluid: &F,
    area: f64,
    discharge_coefficient: f64,
    a: &Volume,
    b: &Volume,
) -> anyhow::Result<Transfer<4>> {
    if area <= 0.0 {
        return Ok(Transfer {
            species: [0.0; 4],
            energy: 0.0,
        });
    }
    let orifice = Orifice::new(Area::new::<square_meter>(area), discharge_coefficient)?;
    real_orifice_transfer_in(fluid, &orifice, a, b)
}

/// The nozzle without flow, the outlet at rest.
fn no_flow(p: f64, t: ThermodynamicTemperature) -> NozzleOperatingPoint {
    NozzleOperatingPoint {
        regime: NozzleRegime::NoFlow,
        mass_flow: MassRate::new::<kilogram_per_second>(0.0),
        exit_mach: 0.0,
        exit_pressure: Pressure::new::<pascal>(p),
        exit_temperature: t,
        exit_velocity: Velocity::new::<meter_per_second>(0.0),
        thrust: Force::new::<newton>(0.0),
        flow_separation: false,
    }
}

/// The ODE over a [`Plant`], from absolute time `start` s, in states divided
/// by `scale` (see [`Plant::scales`]).
pub(super) struct Dynamics<'a, E, F> {
    pub(super) plant: &'a mut Plant,
    pub(super) scale: &'a [f64; STATES],
    pub(super) settings: &'a Settings,
    pub(super) fluid: &'a E,
    pub(super) frozen: &'a F,
    pub(super) start: f64,
    /// Whether the impulse is integrated in the ODE (which then needs the
    /// thrust at every stage) or by the caller once per step.
    pub(super) impulse_in_ode: bool,
    /// The first error met in a derivative (derivatives cannot fail; the step
    /// checks this afterwards).
    pub(super) failure: Option<anyhow::Error>,
}

impl<E, F> Dynamics<'_, E, F> {
    /// The plant's state, scaled.
    pub(super) fn state(&self) -> [f64; STATES] {
        let mut y = self.plant.state();
        for (y, scale) in y.iter_mut().zip(self.scale) {
            *y /= scale;
        }
        y
    }

    /// Sets the plant from a scaled state.
    pub(super) fn unpack(&mut self, scaled: &[f64; STATES]) {
        let y = self.unscaled(scaled);
        self.plant.unpack(&y, self.settings.tank.volume);
    }

    fn unscaled(&self, scaled: &[f64; STATES]) -> [f64; STATES] {
        let mut y = *scaled;
        for (y, scale) in y.iter_mut().zip(self.scale) {
            *y *= scale;
        }
        y
    }
}

impl<E: RealFluid, F: RealFluid> OdeSystem<STATES> for Dynamics<'_, E, F> {
    fn derivatives(&mut self, t: f64, scaled: &[f64; STATES], dydt: &mut [f64; STATES]) {
        let s = self.settings;
        let y = &self.unscaled(scaled);
        self.plant.unpack(y, s.tank.volume);
        let flows = match self
            .plant
            .flows(s, self.fluid, self.frozen, self.impulse_in_ode)
        {
            Ok(flows) => flows,
            Err(e) => {
                self.failure.get_or_insert(e);
                *dydt = [f64::NAN; STATES];
                return;
            }
        };
        *dydt = [0.0; STATES];
        let (reg, val) = (&flows.regulator, &flows.valve);
        let nozzle_flow = flows.nozzle.mass_flow.get::<kilogram_per_second>();
        let outlet_fractions = *self.plant.outlet.mass_fractions();
        for i in 0..4 {
            dydt[TANK + i] = -reg.species[i];
            dydt[PLENUM + i] = reg.species[i] - val.species[i];
            dydt[OUTLET + i] = val.species[i] - nozzle_flow * outlet_fractions[i];
        }
        dydt[TANK + 4] = -reg.energy + flows.wall_to_gas;
        dydt[PLENUM + 4] = reg.energy - val.energy + flows.plenum_heat;
        dydt[OUTLET + 4] = val.energy - flows.nozzle_energy + flows.outlet_heat;
        dydt[WALL] = (flows.outside_to_wall - flows.wall_to_gas) / s.tank.heat_capacity;
        if s.setpoint > 0.0 {
            let (p_in, p_out) = (
                self.plant.tank.pressure().get::<pascal>(),
                self.plant.plenum.pressure().get::<pascal>(),
            );
            (dydt[LIFT], dydt[LIFT_RATE]) =
                regulator::motion(y[LIFT], y[LIFT_RATE], p_in, p_out, s.setpoint);
        }
        let open = s.fire.open(self.start + t, s.pulse_width, s.pulse_period);
        dydt[VALVE] = valve::travel(y[VALVE], open);
        if self.impulse_in_ode {
            dydt[IMPULSE] = flows.nozzle.thrust.get::<newton>();
        }
        for (rate, scale) in dydt.iter_mut().zip(self.scale) {
            *rate /= scale;
        }
    }
}
