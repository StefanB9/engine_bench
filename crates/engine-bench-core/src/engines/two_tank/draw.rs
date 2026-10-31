//! Scene of the two-tank demo, built only from a channel sample.

use anyhow::Context;

use super::channel;
use crate::{
    CoreError,
    scene::{Fill, Label, Point, Scene, Stroke},
    spec::QuantityKind,
};

/// Tanks are drawn as squares whose area is proportional to their volume.
const DEPTH: f64 = 0.1;
/// Flows below this fraction of [`ARROW_REFERENCE_FLOW`] are not drawn.
const NEGLIGIBLE_FLOW: f64 = 1e-3;
/// Mass flow drawn as a full-length arrow, kg/s.
const ARROW_REFERENCE_FLOW: f64 = 0.02;

/// Draws the demo for `sample` (laid out as in `channel`) into `scene`.
pub(crate) fn draw(sample: &[f64], scene: &mut Scene) -> anyhow::Result<()> {
    // Same order as `channel`.
    let &[
        p_a,
        p_b,
        t_a,
        t_b,
        _m_a,
        _m_b,
        mass_flow,
        _q_a,
        _q_b,
        v_a,
        v_b,
        valve,
    ] = sample
    else {
        return Err(CoreError::SampleLength {
            expected: channel::COUNT,
            actual: sample.len(),
        })
        .context("two-tank scene");
    };
    scene.clear();

    let (side_a, side_b) = ((v_a / DEPTH).sqrt(), (v_b / DEPTH).sqrt());
    let gap = 0.4;
    let (left, right) = (-0.5 * gap, 0.5 * gap);
    scene.rectangle(
        Point::new(left - side_a, 0.0),
        Point::new(left, side_a),
        Fill::Gas { temperature: t_a },
        Stroke::Outline,
    );
    scene.rectangle(
        Point::new(right, 0.0),
        Point::new(right + side_b, side_b),
        Fill::Gas { temperature: t_b },
        Stroke::Outline,
    );

    let pipe = 0.25 * side_a.min(side_b);
    scene.rectangle(
        Point::new(left, pipe - 0.015),
        Point::new(right, pipe + 0.015),
        Fill::None,
        Stroke::Outline,
    );
    let w = 0.04;
    scene.polygon(
        &[
            Point::new(-w, pipe - w),
            Point::new(-w, pipe + w),
            Point::new(0.0, pipe),
        ],
        Fill::Metal,
        Stroke::Outline,
    );
    scene.polygon(
        &[
            Point::new(w, pipe - w),
            Point::new(w, pipe + w),
            Point::new(0.0, pipe),
        ],
        Fill::Metal,
        Stroke::Outline,
    );
    let valve_label = Label::Value {
        value: valve,
        quantity: QuantityKind::Ratio,
    };
    scene.label(Point::new(0.0, pipe + 2.0 * w), valve_label);

    // Near equilibrium the flow is numerically zero with a random sign; drawing
    // it would show a stray arrowhead.
    if mass_flow.abs() > NEGLIGIBLE_FLOW * ARROW_REFERENCE_FLOW {
        let half = 0.5 * 0.8 * gap * (mass_flow.abs() / ARROW_REFERENCE_FLOW).min(1.0);
        let y = pipe + 4.0 * w;
        let sign = mass_flow.signum();
        scene.arrow(Point::new(-sign * half, y), Point::new(sign * half, y));
    }

    for (center, top, p, t, name) in [
        (left - 0.5 * side_a, side_a, p_a, t_a, "A"),
        (right + 0.5 * side_b, side_b, p_b, t_b, "B"),
    ] {
        scene.label(Point::new(center, -0.06), Label::Text(name));
        let pressure = Label::Value {
            value: p,
            quantity: QuantityKind::Pressure,
        };
        scene.label(Point::new(center, top + 0.06), pressure);
        let temperature = Label::Value {
            value: t,
            quantity: QuantityKind::Temperature,
        };
        scene.label(Point::new(center, top + 0.12), temperature);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene::Shape;

    fn demo_sample(mass_flow: f64) -> [f64; channel::COUNT] {
        let mut s = [0.0; channel::COUNT];
        s[channel::P_A] = 6e5;
        s[channel::P_B] = 2e5;
        s[channel::T_A] = 260.0;
        s[channel::T_B] = 340.0;
        s[channel::M_A] = 0.16;
        s[channel::M_B] = 0.1;
        s[channel::MASS_FLOW] = mass_flow;
        s[channel::V_A] = 0.02;
        s[channel::V_B] = 0.05;
        s[channel::VALVE] = 0.4;
        s
    }

    fn gas_temperatures(scene: &Scene) -> Vec<f64> {
        scene
            .shapes()
            .iter()
            .filter_map(|shape| match shape {
                Shape::Polygon {
                    fill: Fill::Gas { temperature },
                    ..
                } => Some(*temperature),
                _ => None,
            })
            .collect()
    }

    fn arrows(scene: &Scene) -> Vec<(Point, Point)> {
        scene
            .shapes()
            .iter()
            .filter_map(|shape| match shape {
                Shape::Arrow { from, to } => Some((*from, *to)),
                _ => None,
            })
            .collect()
    }

    #[test]
    fn test_draw_colours_tanks_by_temperature() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        draw(&demo_sample(0.01), &mut scene)?;
        anyhow::ensure!(
            gas_temperatures(&scene) == [260.0, 340.0],
            "{:?}",
            gas_temperatures(&scene)
        );
        Ok(())
    }

    #[test]
    fn test_draw_tank_size_follows_volume() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        draw(&demo_sample(0.0), &mut scene)?;
        let areas: Vec<f64> = scene
            .shapes()
            .iter()
            .filter_map(|shape| match shape {
                Shape::Polygon {
                    points,
                    fill: Fill::Gas { .. },
                    ..
                } => {
                    let p = scene.points(*points);
                    let (min_x, max_x) = p.iter().fold((f64::MAX, f64::MIN), |(lo, hi), q| {
                        (lo.min(q.x), hi.max(q.x))
                    });
                    let (min_y, max_y) = p.iter().fold((f64::MAX, f64::MIN), |(lo, hi), q| {
                        (lo.min(q.y), hi.max(q.y))
                    });
                    Some((max_x - min_x) * (max_y - min_y))
                }
                _ => None,
            })
            .collect();
        anyhow::ensure!(areas.len() == 2, "{areas:?}");
        anyhow::ensure!(
            (areas[1] / areas[0] - 2.5).abs() < 1e-12,
            "area ratio {areas:?}"
        );
        Ok(())
    }

    #[test]
    fn test_draw_flow_arrow_follows_mass_flow_sign() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        draw(&demo_sample(0.01), &mut scene)?;
        let forward = arrows(&scene);
        anyhow::ensure!(
            forward.len() == 1 && forward[0].1.x > forward[0].0.x,
            "{forward:?}"
        );
        draw(&demo_sample(-0.01), &mut scene)?;
        let backward = arrows(&scene);
        anyhow::ensure!(
            backward.len() == 1 && backward[0].1.x < backward[0].0.x,
            "{backward:?}"
        );
        draw(&demo_sample(0.0), &mut scene)?;
        anyhow::ensure!(arrows(&scene).is_empty(), "arrow without flow");
        Ok(())
    }

    #[test]
    fn test_draw_hides_arrow_for_negligible_flow() -> anyhow::Result<()> {
        // At equilibrium the flow is numerically zero with either sign.
        let mut scene = Scene::new();
        for flow in [2e-11, -2e-11] {
            draw(&demo_sample(flow), &mut scene)?;
            anyhow::ensure!(arrows(&scene).is_empty(), "arrow for {flow} kg/s");
        }
        Ok(())
    }

    #[test]
    fn test_draw_labels_pressures_as_values() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        draw(&demo_sample(0.01), &mut scene)?;
        let pressures: Vec<f64> = scene
            .shapes()
            .iter()
            .filter_map(|shape| match shape {
                Shape::Label {
                    label:
                        Label::Value {
                            value,
                            quantity: QuantityKind::Pressure,
                        },
                    ..
                } => Some(*value),
                _ => None,
            })
            .collect();
        anyhow::ensure!(pressures == [6e5, 2e5], "{pressures:?}");
        Ok(())
    }

    #[test]
    fn test_draw_names_tanks_below_them() -> anyhow::Result<()> {
        // Inside the tanks, light text on a light gas fill is illegible.
        let mut scene = Scene::new();
        draw(&demo_sample(0.01), &mut scene)?;
        let names: Vec<f64> = scene
            .shapes()
            .iter()
            .filter_map(|shape| match shape {
                Shape::Label {
                    at,
                    label: Label::Text("A" | "B"),
                } => Some(at.y),
                _ => None,
            })
            .collect();
        anyhow::ensure!(
            names.len() == 2 && names.iter().all(|y| *y < 0.0),
            "{names:?}"
        );
        Ok(())
    }

    #[test]
    fn test_draw_replaces_previous_contents() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        draw(&demo_sample(0.01), &mut scene)?;
        let first = scene.clone();
        draw(&demo_sample(0.01), &mut scene)?;
        anyhow::ensure!(scene == first, "drawing twice differs");
        Ok(())
    }
}
