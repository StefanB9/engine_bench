//! Planar rigid-body flight integrated in time: ballistics, terminal velocity,
//! static stability, pitch damping and weathervaning, through the public API
//! only.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use anyhow::Context;
use engine_bench_physics::{
    flight::{
        aero::{AeroModel, aerodynamic_loads},
        atmosphere::STANDARD_GRAVITY,
        barrowman::{Fins, NoseShape, RocketGeometry, barrowman},
        rigid_body::{Loads, MassProperties, PlanarState, derivatives},
    },
    numerics::ode::{
        AdaptiveOptions, Direction, DormandPrince45, EventSet, NoEvents, OdeSystem, StepOutcome,
    },
    units::{Length, length::meter},
};

const DENSITY: f64 = 1.225;

/// A 1 m, 0.1 m diameter finned rocket (centre of pressure ≈ 0.8 m aft of the
/// nose).
fn rocket() -> anyhow::Result<AeroModel> {
    let m = Length::new::<meter>;
    let geometry = RocketGeometry {
        nose_shape: NoseShape::Ogive,
        nose_length: m(0.2),
        body_diameter: m(0.1),
        fins: Some(Fins {
            count: 4,
            root_chord: m(0.12),
            tip_chord: m(0.06),
            semi_span: m(0.08),
            sweep: m(0.06),
            position: m(0.88),
        }),
    };
    AeroModel::new(barrowman(&geometry)?, 0.4)
}

/// A body of constant mass under gravity `gravity`, with aerodynamics if `aero`
/// is set.
struct Flight {
    aero: Option<AeroModel>,
    mass: MassProperties,
    gravity: f64,
    wind: f64,
}

impl Flight {
    fn new(aero: Option<AeroModel>, cg_from_nose: f64, gravity: f64) -> Self {
        Self {
            aero,
            // A 0.5 kg uniform 1 m rod: I = m·L²/12.
            mass: MassProperties {
                mass: 0.5,
                inertia: 0.5 / 12.0,
                cg_from_nose,
            },
            gravity,
            wind: 0.0,
        }
    }
}

impl OdeSystem<6> for Flight {
    fn derivatives(&mut self, _t: f64, y: &[f64; 6], dydt: &mut [f64; 6]) {
        let state = PlanarState::from_array(y);
        let loads = self.aero.as_ref().map_or_else(Loads::default, |aero| {
            aerodynamic_loads(&state, aero, self.mass.cg_from_nose, DENSITY, self.wind)
        });
        *dydt = derivatives(&state, &self.mass, &loads, self.gravity);
    }
}

/// Fires when the body comes back down to y = 0.
struct Ground;

impl EventSet<6, 1> for Ground {
    fn directions(&self) -> [Direction; 1] {
        [Direction::Falling]
    }

    fn values(&mut self, _t: f64, y: &[f64; 6], g: &mut [f64; 1]) {
        g[0] = y[1];
    }
}

fn integrator() -> anyhow::Result<DormandPrince45<6>> {
    Ok(DormandPrince45::new(
        AdaptiveOptions::new(1e-10, 1e-12)?.with_max_step(0.01)?,
    ))
}

/// Angle of attack: body axis relative to the air-relative velocity, wrapped to
/// (−π, π].
fn angle_of_attack(state: &PlanarState, wind: f64) -> f64 {
    let flight_path = state.vy.atan2(state.vx - wind);
    let alpha = state.theta - flight_path;
    (alpha + PI).rem_euclid(2.0 * PI) - PI
}

#[test]
fn test_flight_without_air_lands_at_projectile_range() -> anyhow::Result<()> {
    let speed = 20.0;
    let start = PlanarState {
        vx: speed * FRAC_PI_4.cos(),
        vy: speed * FRAC_PI_4.sin(),
        theta: FRAC_PI_4,
        y: 1e-9,
        ..Default::default()
    };
    let mut system = Flight::new(None, 0.5, STANDARD_GRAVITY);
    let (mut t, mut y) = (0.0, start.to_array());
    let outcome = integrator()?.integrate(&mut system, &mut Ground, &mut t, &mut y, 60.0)?;
    anyhow::ensure!(outcome == StepOutcome::Event { index: 0 }, "{outcome:?}");
    let range = speed * speed / STANDARD_GRAVITY;
    anyhow::ensure!(
        (y[0] / range - 1.0).abs() < 1e-6,
        "range {} vs {range}",
        y[0]
    );
    Ok(())
}

#[test]
fn test_flight_falls_at_terminal_velocity() -> anyhow::Result<()> {
    // Nose down from rest; a stable rocket stays aligned and reaches √(2mg /
    // (ρ·A·C_A)).
    let start = PlanarState {
        y: 10_000.0,
        theta: -FRAC_PI_2,
        ..Default::default()
    };
    let mut system = Flight::new(Some(rocket()?), 0.4, STANDARD_GRAVITY);
    let (mut t, mut y) = (0.0, start.to_array());
    integrator()?.integrate(&mut system, &mut NoEvents, &mut t, &mut y, 120.0)?;
    let area = FRAC_PI_4 * 0.01;
    let terminal = (2.0 * 0.5 * STANDARD_GRAVITY / (DENSITY * area * 0.4)).sqrt();
    anyhow::ensure!(
        (-y[3] / terminal - 1.0).abs() < 1e-3,
        "v {} vs {terminal}",
        -y[3]
    );
    anyhow::ensure!(y[2].abs() < 1e-6, "drifted sideways at {} m/s", y[2]);
    Ok(())
}

/// Largest |angle of attack| over the last second of a 3 s flight started 10°
/// off the velocity, without gravity, plus the largest over the whole flight.
fn pitch_response(cg_from_nose: f64) -> anyhow::Result<(f64, f64)> {
    let start = PlanarState {
        vx: 50.0,
        theta: 10f64.to_radians(),
        ..Default::default()
    };
    let mut system = Flight::new(Some(rocket()?), cg_from_nose, 0.0);
    let mut dp = integrator()?;
    let (mut t, mut y) = (0.0, start.to_array());
    let (mut late, mut overall) = (0.0f64, 0.0f64);
    while t < 3.0 {
        dp.step(&mut system, &mut NoEvents, &mut t, &mut y, 3.0)
            .context("pitch response step")?;
        let alpha = angle_of_attack(&PlanarState::from_array(&y), 0.0).abs();
        overall = overall.max(alpha);
        if t > 2.0 {
            late = late.max(alpha);
        }
    }
    Ok((late, overall))
}

#[test]
fn test_flight_stable_rocket_damps_its_pitch_oscillation() -> anyhow::Result<()> {
    let (late, overall) = pitch_response(0.4)?;
    anyhow::ensure!(
        late < 0.25 * 10f64.to_radians(),
        "still {}° off after 2 s",
        late.to_degrees()
    );
    anyhow::ensure!(
        overall < 11f64.to_radians(),
        "overshot to {}°",
        overall.to_degrees()
    );
    Ok(())
}

#[test]
fn test_flight_unstable_rocket_tumbles() -> anyhow::Result<()> {
    // Centre of gravity behind the centre of pressure.
    let (_, overall) = pitch_response(0.95)?;
    anyhow::ensure!(
        overall > FRAC_PI_2,
        "never passed 90°: {}°",
        overall.to_degrees()
    );
    Ok(())
}

#[test]
fn test_flight_stable_rocket_turns_into_the_wind() -> anyhow::Result<()> {
    // Climbing straight up at 30 m/s through a 5 m/s wind blowing towards +x,
    // gravity off.
    let start = PlanarState {
        vy: 30.0,
        theta: FRAC_PI_2,
        ..Default::default()
    };
    let mut system = Flight::new(Some(rocket()?), 0.4, 0.0);
    system.wind = 5.0;
    let (mut t, mut y) = (0.0, start.to_array());
    integrator()?.integrate(&mut system, &mut NoEvents, &mut t, &mut y, 3.0)?;
    let state = PlanarState::from_array(&y);
    // Nose tilted towards −x, i.e. upwind, and lined up with the relative wind.
    anyhow::ensure!(
        state.theta > FRAC_PI_2 + 1f64.to_radians(),
        "θ = {}°",
        state.theta.to_degrees()
    );
    let alpha = angle_of_attack(&state, system.wind);
    anyhow::ensure!(
        alpha.abs() < 1f64.to_radians(),
        "α = {}°",
        alpha.to_degrees()
    );
    Ok(())
}
