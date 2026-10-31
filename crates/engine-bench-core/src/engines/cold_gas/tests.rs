//! Thruster behaviour against physics: perfect-gas nozzle limits, impulse,
//! conservation, regulation, pulsing and heat transfer.

use engine_bench_physics::{
    flow::{Nozzle, PerfectGas, Reservoir},
    units::{
        Area, Pressure, SpecificHeatCapacity, area::square_meter, force::newton,
        mass_rate::kilogram_per_second, pressure::pascal,
        specific_heat_capacity::joule_per_kilogram_kelvin,
    },
};

use super::*;

fn seconds(value: f64) -> Time {
    Time::new::<second>(value)
}

fn sample(engine: &ColdGasThruster) -> anyhow::Result<[f64; channel::COUNT]> {
    let mut out = [0.0; channel::COUNT];
    engine.sample(&mut out)?;
    Ok(out)
}

fn set(engine: &mut ColdGasThruster, control: usize, value: f64) -> anyhow::Result<()> {
    engine.set_control(ControlId::new(control), value)
}

/// The thruster firing continuously at default settings, after `time` s of
/// real-time stepping.
fn firing(time: f64) -> anyhow::Result<ColdGasThruster> {
    let mut engine = ColdGasThruster::new()?;
    set(&mut engine, FIRE, 1.0)?;
    if time > 0.0 {
        engine.advance(seconds(time), &StepContext::real_time())?;
    }
    Ok(engine)
}

/// Propellant upstream of the valve (tank and plenum), kg.
fn upstream_mass(engine: &ColdGasThruster) -> f64 {
    let p = &engine.plant;
    p.tank.mass().get::<kilogram>() + p.plenum.mass().get::<kilogram>()
}

#[test]
fn test_thruster_starts_locked_up_and_quiet() -> anyhow::Result<()> {
    let engine = ColdGasThruster::new()?;
    let s = sample(&engine)?;
    anyhow::ensure!((s[channel::P_TANK] / 2e7 - 1.0).abs() < 1e-9, "tank {s:?}");
    // The plenum sits at the regulator's lock-up, just under the setpoint.
    let lock_up = regulator::lock_up(2e7, 1e6);
    anyhow::ensure!(
        (s[channel::P_PLENUM] / lock_up - 1.0).abs() < 1e-9,
        "plenum {s:?}"
    );
    anyhow::ensure!(
        s[channel::THRUST] == 0.0 && s[channel::MASS_FLOW] == 0.0 && s[channel::VALVE] == 0.0,
        "{s:?}"
    );
    Ok(())
}

#[test]
fn test_thruster_matches_the_perfect_gas_nozzle_from_its_valve_outlet() -> anyhow::Result<()> {
    // At 10 bar and room temperature N2 is within ~0.2 % of a perfect gas
    // with γ = 1.4; the frozen expansion to the exit stays dilute.
    let engine = firing(0.3)?;
    let s = sample(&engine)?;
    let gas = PerfectGas::new(
        1.4,
        SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(8.314_462_618 / 0.028_013_4),
    )?;
    let outlet = &engine.plant.outlet;
    let chamber = Reservoir::new(outlet.pressure(), outlet.temperature(), gas)?;
    let nozzle = Nozzle::new(
        Area::new::<square_meter>(hardware::circle_area(1e-3)),
        Area::new::<square_meter>(hardware::circle_area(7e-3)),
    )?;
    let ideal = nozzle.operating_point(&chamber, Pressure::new::<pascal>(0.0))?;
    let close = |a: f64, b: f64, what: &str| -> anyhow::Result<()> {
        anyhow::ensure!((a / b - 1.0).abs() < 1e-2, "{what}: {a} vs {b}");
        Ok(())
    };
    close(s[channel::THRUST], ideal.thrust.get::<newton>(), "thrust")?;
    close(
        s[channel::MASS_FLOW],
        ideal.mass_flow.get::<kilogram_per_second>(),
        "mass flow",
    )?;
    let ideal_isp = ideal.specific_impulse().context("no flow")?.get::<second>();
    close(s[channel::ISP], ideal_isp, "Isp")?;
    // N2 in vacuum through area ratio 49: around 72 s at room temperature.
    anyhow::ensure!(
        (60.0..90.0).contains(&s[channel::ISP]),
        "Isp {}",
        s[channel::ISP]
    );
    // The frozen exit is far supersaturated.
    anyhow::ensure!(
        s[channel::SUPERSATURATION] > 1.0,
        "{}",
        s[channel::SUPERSATURATION]
    );
    Ok(())
}

#[test]
fn test_thruster_regulates_with_droop() -> anyhow::Result<()> {
    let engine = firing(0.3)?;
    let s = sample(&engine)?;
    let p = s[channel::P_PLENUM];
    // Below the setpoint while flowing (the spring relaxes as the poppet
    // lifts), within a few per cent of it.
    anyhow::ensure!(p < 1e6 && p > 0.94e6, "plenum {p}");
    let lift = s[channel::LIFT];
    anyhow::ensure!(lift > 0.0 && lift < 1.0, "lift {lift}");
    Ok(())
}

#[test]
fn test_thruster_impulse_integrates_the_thrust() -> anyhow::Result<()> {
    let mut engine = ColdGasThruster::new()?;
    set(&mut engine, FIRE, 1.0)?;
    let (step, mut integral, mut last) = (1e-3, 0.0, 0.0);
    for _ in 0..150 {
        engine.advance(seconds(step), &StepContext::real_time())?;
        let thrust = sample(&engine)?[channel::THRUST];
        integral += 0.5 * step * (last + thrust);
        last = thrust;
    }
    let impulse = sample(&engine)?[channel::IMPULSE];
    anyhow::ensure!(
        (impulse / integral - 1.0).abs() < 2e-3,
        "{impulse} vs {integral}"
    );
    Ok(())
}

#[test]
fn test_thruster_conserves_propellant_with_the_valve_shut() -> anyhow::Result<()> {
    // Raising the setpoint opens the regulator: the plenum fills from the
    // tank and nothing passes the valve. (The valve outlet's residue at the
    // vacuum floor may seep out through the nozzle.)
    let mut engine = ColdGasThruster::new()?;
    let before = upstream_mass(&engine);
    set(&mut engine, SETPOINT, 2e6)?;
    engine.advance(seconds(0.2), &StepContext::real_time())?;
    let after = upstream_mass(&engine);
    anyhow::ensure!((after / before - 1.0).abs() < 1e-12, "{before} → {after}");
    let s = sample(&engine)?;
    anyhow::ensure!(
        s[channel::P_PLENUM] > 1.5e6,
        "plenum {}",
        s[channel::P_PLENUM]
    );
    anyhow::ensure!(
        s[channel::IMPULSE] < 1e-9,
        "impulse {}",
        s[channel::IMPULSE]
    );
    Ok(())
}

#[test]
fn test_thruster_pulses_deliver_impulse_bits() -> anyhow::Result<()> {
    // 20 ms pulses every 100 ms: the valve opens and shuts, and the impulse
    // grows in equal steps.
    let mut engine = ColdGasThruster::new()?;
    set(&mut engine, FIRE, 2.0)?;
    set(&mut engine, PULSE_WIDTH, 0.02)?;
    set(&mut engine, PULSE_PERIOD, 0.1)?;
    let rt = StepContext::real_time();
    engine.advance(seconds(0.06), &rt)?;
    let first = sample(&engine)?;
    anyhow::ensure!(first[channel::VALVE] == 0.0, "valve still open {first:?}");
    engine.advance(seconds(0.1), &rt)?;
    let later = sample(&engine)?;
    // Thrust × width, plus the tail-off of the valve outlet draining after
    // the valve shuts (2 cm³ at ~9 bar, ~2e-5 kg at ~700 m/s).
    let bit = first[channel::IMPULSE];
    let nominal = 1.28 * 0.02;
    anyhow::ensure!(bit > nominal && bit < nominal + 0.02, "impulse bit {bit}");
    let ratio = later[channel::IMPULSE] / bit;
    anyhow::ensure!((ratio - 2.0).abs() < 0.05, "two pulses give {ratio} bits");
    Ok(())
}

#[test]
fn test_thruster_tank_wall_warms_a_draining_tank() -> anyhow::Result<()> {
    // Unregulated blowdown through a large throat: with heat from the wall the
    // tank stays warmer than adiabatic, and the wall cools.
    let run = |h: f64| -> anyhow::Result<[f64; channel::COUNT]> {
        let mut engine = ColdGasThruster::new()?;
        set(&mut engine, SETPOINT, 0.0)?;
        set(&mut engine, TANK_VOLUME, 2e-4)?;
        set(&mut engine, TANK_PRESSURE, 5e6)?;
        set(&mut engine, THROAT_DIAMETER, 3e-3)?;
        set(&mut engine, EXIT_DIAMETER, 6e-3)?;
        set(&mut engine, WALL_INSIDE, h)?;
        engine.reset()?;
        set(&mut engine, FIRE, 1.0)?;
        engine.advance(seconds(0.3), &StepContext::real_time())?;
        sample(&engine)
    };
    let (adiabatic, heated) = (run(0.0)?, run(200.0)?);
    let start = 293.15;
    anyhow::ensure!(adiabatic[channel::T_TANK] < start - 10.0, "{adiabatic:?}");
    anyhow::ensure!(
        heated[channel::T_TANK] > adiabatic[channel::T_TANK],
        "{} vs {}",
        heated[channel::T_TANK],
        adiabatic[channel::T_TANK]
    );
    anyhow::ensure!(
        (adiabatic[channel::T_WALL] - start).abs() < 1e-12 && heated[channel::T_WALL] < start,
        "wall {} / {}",
        adiabatic[channel::T_WALL],
        heated[channel::T_WALL]
    );
    Ok(())
}

#[test]
fn test_thruster_real_time_agrees_with_offline() -> anyhow::Result<()> {
    let mut fast = firing(0.0)?;
    let mut accurate = firing(0.0)?;
    fast.advance(seconds(0.1), &StepContext::real_time())?;
    accurate.advance(seconds(0.1), &StepContext::offline(1e-8)?)?;
    let (f, a) = (sample(&fast)?, sample(&accurate)?);
    for i in [
        channel::P_PLENUM,
        channel::P_OUTLET,
        channel::THRUST,
        channel::IMPULSE,
    ] {
        anyhow::ensure!(
            (f[i] / a[i] - 1.0).abs() < 1e-3,
            "channel {i}: {} vs {}",
            f[i],
            a[i]
        );
    }
    Ok(())
}

#[test]
fn test_thruster_helium_has_the_higher_specific_impulse() -> anyhow::Result<()> {
    let mut engine = ColdGasThruster::new()?;
    set(&mut engine, PROPELLANT, 2.0)?;
    engine.reset()?;
    set(&mut engine, FIRE, 1.0)?;
    engine.advance(seconds(0.15), &StepContext::real_time())?;
    let isp = sample(&engine)?[channel::ISP];
    // He in vacuum through area ratio 49: about 175 s.
    anyhow::ensure!(isp > 160.0 && isp < 190.0, "Isp {isp}");
    Ok(())
}

#[test]
fn test_thruster_exact_and_tabulated_properties_agree() -> anyhow::Result<()> {
    let mut tables = firing(0.05)?;
    let mut exact = tables.clone();
    set(&mut exact, PROPERTIES, 1.0)?;
    let (t, e) = (sample(&tables)?, sample(&exact)?);
    for i in [
        channel::THRUST,
        channel::MASS_FLOW,
        channel::EXIT_TEMPERATURE,
    ] {
        anyhow::ensure!(
            (t[i] / e[i] - 1.0).abs() < 1e-4,
            "channel {i}: {} vs {}",
            t[i],
            e[i]
        );
    }
    // Switching back keeps the tables.
    set(&mut exact, PROPERTIES, 0.0)?;
    tables.advance(seconds(1e-3), &StepContext::real_time())?;
    Ok(())
}

#[test]
fn test_thruster_rejects_invalid_controls_and_steps() -> anyhow::Result<()> {
    let mut engine = ColdGasThruster::new()?;
    anyhow::ensure!(set(&mut engine, PROPELLANT, 3.0).is_err(), "propellant 3");
    anyhow::ensure!(
        set(&mut engine, THROAT_DIAMETER, 1.0).is_err(),
        "1 m throat"
    );
    anyhow::ensure!(set(&mut engine, 99, 0.0).is_err(), "control 99");
    for dt in [0.0, -1.0, f64::NAN] {
        anyhow::ensure!(
            engine
                .advance(seconds(dt), &StepContext::real_time())
                .is_err(),
            "dt {dt}"
        );
    }
    anyhow::ensure!(engine.sample(&mut [0.0; 3]).is_err(), "short sample buffer");
    Ok(())
}

#[test]
fn test_thruster_reset_refills_and_restarts_time() -> anyhow::Result<()> {
    let mut engine = firing(0.05)?;
    set(&mut engine, TANK_PRESSURE, 1e7)?;
    set(&mut engine, RESET, 1.0)?;
    let s = sample(&engine)?;
    anyhow::ensure!(
        (s[channel::P_TANK] / 1e7 - 1.0).abs() < 1e-9,
        "tank {}",
        s[channel::P_TANK]
    );
    anyhow::ensure!(engine.time().get::<second>() == 0.0, "time");
    anyhow::ensure!(s[channel::IMPULSE] == 0.0, "impulse");
    Ok(())
}
