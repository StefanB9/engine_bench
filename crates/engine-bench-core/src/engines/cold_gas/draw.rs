//! Scene of the cold-gas thruster, built only from a channel sample.
//!
//! The tank is drawn to scale, lying along x with its outlet at x = 0; the
//! regulator, line and valve schematically; the nozzle contour from its
//! diameters at twice scale, and a plume whose length grows with thrust.

use anyhow::Context;

use super::channel;
use crate::{
    CoreError,
    scene::{Fill, Label, Point, Scene, Stroke},
    spec::QuantityKind,
};

/// The tank is this many diameters long (as built).
const TANK_ASPECT: f64 = 3.0;
/// Drawn wall thickness, m.
const WALL: f64 = 0.004;
/// Pipe half-height, m.
const PIPE: f64 = 0.003;
/// x positions, m: regulator box, valve, nozzle inlet, throat and exit.
const REGULATOR: (f64, f64) = (0.03, 0.07);
const VALVE: f64 = 0.15;
const NOZZLE_INLET: f64 = 0.17;
const THROAT_X: f64 = 0.19;
const EXIT_X: f64 = 0.24;
/// Nozzle drawing scale (half-width = this × diameter / 2).
const NOZZLE_SCALE: f64 = 2.0;
/// Plume length per newton of thrust, m, and its cap.
const PLUME_PER_NEWTON: f64 = 0.08;
const PLUME_MAX: f64 = 0.4;
/// Under-expanded plumes spread, over-expanded ones contract.
const UNDER_EXPANDED: f64 = 6.0;

/// Draws the thruster for `sample` (laid out as in `channel`) into `scene`.
pub(crate) fn draw(sample: &[f64], scene: &mut Scene) -> anyhow::Result<()> {
    let s: &[f64; channel::COUNT] = sample
        .try_into()
        .map_err(|_| CoreError::SampleLength {
            expected: channel::COUNT,
            actual: sample.len(),
        })
        .context("cold-gas scene")?;
    scene.clear();
    let (length, radius) = tank(s, scene);
    plumbing(s, scene);
    let exit = nozzle(s, scene);
    plume(s, exit, scene);

    let value = |value, quantity| Label::Value { value, quantity };
    scene.label(
        Point::new(-0.5 * length, radius + 0.03),
        value(s[channel::P_TANK], QuantityKind::Pressure),
    );
    scene.label(
        Point::new(-0.5 * length, radius + 0.06),
        value(s[channel::T_TANK], QuantityKind::Temperature),
    );
    scene.label(
        Point::new(0.11, 0.03),
        value(s[channel::P_PLENUM], QuantityKind::Pressure),
    );
    scene.label(
        Point::new(f64::midpoint(REGULATOR.0, REGULATOR.1), -0.04),
        Label::Text("Regulator"),
    );
    scene.label(Point::new(VALVE, -0.04), Label::Text("Valve"));

    let reach = radius.max(0.05) + 0.1;
    scene.set_view(
        Point::new(-length - 0.04, -reach),
        Point::new(EXIT_X + PLUME_MAX + 0.04, reach),
    );
    Ok(())
}

/// The tank to scale, its gas coloured by temperature; returns its length
/// and radius, m.
fn tank(s: &[f64; channel::COUNT], scene: &mut Scene) -> (f64, f64) {
    let diameter = (4.0 * s[channel::TANK_VOLUME] / (std::f64::consts::PI * TANK_ASPECT)).cbrt();
    let (length, radius) = (TANK_ASPECT * diameter, 0.5 * diameter);
    scene.rectangle(
        Point::new(-length - WALL, -radius - WALL),
        Point::new(0.0, radius + WALL),
        Fill::Metal,
        Stroke::Outline,
    );
    scene.rectangle(
        Point::new(-length, -radius),
        Point::new(-WALL, radius),
        Fill::Gas {
            temperature: s[channel::T_TANK],
        },
        Stroke::None,
    );
    (length, radius)
}

/// The line from the tank through the regulator (its lift a gap) to the
/// valve (two wedges that part as it opens).
fn plumbing(s: &[f64; channel::COUNT], scene: &mut Scene) {
    let plenum_gas = Fill::Gas {
        temperature: s[channel::T_PLENUM],
    };
    scene.rectangle(
        Point::new(0.0, -PIPE),
        Point::new(REGULATOR.0, PIPE),
        Fill::None,
        Stroke::Outline,
    );
    scene.rectangle(
        Point::new(REGULATOR.0, -0.02),
        Point::new(REGULATOR.1, 0.02),
        Fill::Metal,
        Stroke::Outline,
    );
    let gap = 0.015 * s[channel::LIFT].clamp(0.0, 1.0);
    let middle = f64::midpoint(REGULATOR.0, REGULATOR.1);
    if gap > 0.0 {
        scene.rectangle(
            Point::new(middle - 0.5 * gap, -PIPE),
            Point::new(middle + 0.5 * gap, PIPE),
            plenum_gas,
            Stroke::None,
        );
    }
    scene.rectangle(
        Point::new(REGULATOR.1, -PIPE),
        Point::new(NOZZLE_INLET, PIPE),
        plenum_gas,
        Stroke::Outline,
    );
    let opening = s[channel::VALVE].clamp(0.0, 1.0);
    let (w, lift) = (0.012, 0.008 * opening);
    for side in [-1.0, 1.0] {
        scene.polygon(
            &[
                Point::new(VALVE - w, side * (w + lift)),
                Point::new(VALVE + w, side * (w + lift)),
                Point::new(VALVE, side * lift),
            ],
            Fill::Metal,
            Stroke::Outline,
        );
    }
}

/// The nozzle's upper and lower walls; returns the drawn exit half-width, m.
fn nozzle(s: &[f64; channel::COUNT], scene: &mut Scene) -> f64 {
    let half = |d: f64| 0.5 * NOZZLE_SCALE * d;
    let exit = half(s[channel::EXIT]);
    let contour = [
        Point::new(NOZZLE_INLET, 0.012_f64.max(exit)),
        Point::new(THROAT_X, half(s[channel::THROAT])),
        Point::new(EXIT_X, exit),
    ];
    for side in [-1.0, 1.0] {
        let wall: Vec<Point> = contour
            .iter()
            .map(|p| Point::new(p.x, side * p.y))
            .chain(
                contour
                    .iter()
                    .rev()
                    .map(|p| Point::new(p.x, side * (p.y + 0.004))),
            )
            .collect();
        scene.polygon(&wall, Fill::Metal, Stroke::Outline);
    }
    exit
}

/// The plume from an exit of half-width `exit`: longer with thrust, spread
/// when under-expanded, coloured by the exit temperature.
fn plume(s: &[f64; channel::COUNT], exit: f64, scene: &mut Scene) {
    let thrust = s[channel::THRUST];
    if thrust <= 0.0 {
        return;
    }
    let plume = (PLUME_PER_NEWTON * thrust).clamp(0.02, PLUME_MAX);
    let spread = if s[channel::REGIME] >= UNDER_EXPANDED {
        1.6
    } else {
        0.7
    };
    scene.polygon(
        &[
            Point::new(EXIT_X, exit),
            Point::new(EXIT_X + plume, spread * exit),
            Point::new(EXIT_X + plume, -spread * exit),
            Point::new(EXIT_X, -exit),
        ],
        Fill::Gas {
            temperature: s[channel::EXIT_TEMPERATURE],
        },
        Stroke::None,
    );
    scene.label(
        Point::new(EXIT_X + 0.5 * plume, -spread.max(1.0) * exit - 0.03),
        Label::Value {
            value: thrust,
            quantity: QuantityKind::Force,
        },
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene::Shape;

    fn sample(thrust: f64) -> [f64; channel::COUNT] {
        let mut s = [0.0; channel::COUNT];
        s[channel::P_TANK] = 2e7;
        s[channel::T_TANK] = 280.0;
        s[channel::P_PLENUM] = 1e6;
        s[channel::T_PLENUM] = 270.0;
        s[channel::LIFT] = 0.3;
        s[channel::VALVE] = 1.0;
        s[channel::THRUST] = thrust;
        s[channel::EXIT_TEMPERATURE] = 40.0;
        s[channel::REGIME] = 6.0;
        s[channel::TANK_VOLUME] = 2e-3;
        s[channel::THROAT] = 1e-3;
        s[channel::EXIT] = 7e-3;
        s
    }

    fn gas_polygons(scene: &Scene) -> Vec<(f64, Vec<Point>)> {
        scene
            .shapes()
            .iter()
            .filter_map(|shape| match shape {
                Shape::Polygon {
                    points,
                    fill: Fill::Gas { temperature },
                    ..
                } => Some((*temperature, scene.points(*points).to_vec())),
                _ => None,
            })
            .collect()
    }

    fn extent(points: &[Point]) -> f64 {
        let xs = points.iter().map(|p| p.x);
        xs.clone().fold(f64::MIN, f64::max) - xs.fold(f64::MAX, f64::min)
    }

    #[test]
    fn test_draw_colours_the_tank_and_plume_by_temperature() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        draw(&sample(2.0), &mut scene)?;
        let temperatures: Vec<f64> = gas_polygons(&scene).iter().map(|(t, _)| *t).collect();
        anyhow::ensure!(
            temperatures.contains(&280.0) && temperatures.contains(&40.0),
            "{temperatures:?}"
        );
        Ok(())
    }

    #[test]
    fn test_draw_plume_grows_with_thrust_and_vanishes_without() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        let plume = |scene: &Scene| {
            gas_polygons(scene)
                .into_iter()
                .find(|(t, _)| (*t - 40.0).abs() < f64::EPSILON)
                .map(|(_, p)| extent(&p))
        };
        draw(&sample(1.0), &mut scene)?;
        let small = plume(&scene).context("no plume at 1 N")?;
        draw(&sample(3.0), &mut scene)?;
        let large = plume(&scene).context("no plume at 3 N")?;
        anyhow::ensure!(large > small, "{small} vs {large}");
        draw(&sample(0.0), &mut scene)?;
        anyhow::ensure!(plume(&scene).is_none(), "plume without thrust");
        Ok(())
    }

    #[test]
    fn test_draw_tank_length_follows_its_volume() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        let tank_length = |scene: &Scene| {
            gas_polygons(scene)
                .into_iter()
                .find(|(t, _)| (*t - 280.0).abs() < f64::EPSILON)
                .map(|(_, p)| extent(&p) + WALL)
        };
        let mut s = sample(0.0);
        draw(&s, &mut scene)?;
        let two_litres = tank_length(&scene).context("tank")?;
        s[channel::TANK_VOLUME] = 16e-3;
        draw(&s, &mut scene)?;
        let sixteen = tank_length(&scene).context("tank")?;
        // Length ∝ V^(1/3).
        anyhow::ensure!(
            (sixteen / two_litres - 2.0).abs() < 1e-12,
            "{two_litres} vs {sixteen}"
        );
        anyhow::ensure!(scene.view().is_some(), "view");
        Ok(())
    }

    #[test]
    fn test_draw_rejects_short_samples() {
        let mut scene = Scene::new();
        assert!(draw(&[0.0; 3], &mut scene).is_err());
    }
}
