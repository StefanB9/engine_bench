//! Turns a renderer-agnostic [`Scene`] into egui shapes.
//!
//! World coordinates (metres, y up) are fitted into the view rectangle with a
//! uniform scale; semantic fills become colours, gas by temperature through a
//! colormap; value labels are formatted in the selected display units.
//!
//! Polygons are filled as convex shapes (egui's `convex_polygon`); every
//! current scene shape is convex. Concave parts need tessellation once an
//! engine draws them.

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Shape, Stroke as EguiStroke, Vec2};
use engine_bench_core::scene::{Fill, Label, Point, Scene, Shape as SceneShape, Stroke};

use crate::units::{UnitSelection, format_value};

/// Maps world coordinates into a screen rectangle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ViewTransform {
    /// Screen points per metre.
    scale: f64,
    /// World point shown at `centre_screen`.
    centre_world: Point,
    centre_screen: Pos2,
}

impl ViewTransform {
    /// Fits world `bounds` (min, max) centred into `rect`, leaving `margin`
    /// points free on every side; degenerate bounds still get a finite scale.
    pub(crate) fn fit(bounds: (Point, Point), rect: Rect, margin: f32) -> Self {
        let (min, max) = bounds;
        // A point or a line still needs a finite scale.
        let width = (max.x - min.x).max(1e-6);
        let height = (max.y - min.y).max(1e-6);
        let available_width = f64::from((rect.width() - 2.0 * margin).max(1.0));
        let available_height = f64::from((rect.height() - 2.0 * margin).max(1.0));
        Self {
            scale: (available_width / width).min(available_height / height),
            centre_world: Point::new(f64::midpoint(min.x, max.x), f64::midpoint(min.y, max.y)),
            centre_screen: rect.center(),
        }
    }

    /// Screen position of world point `p` (y flipped: up on screen is up in the
    /// world).
    pub(crate) fn to_screen(self, p: Point) -> Pos2 {
        Pos2::new(
            self.centre_screen.x + screen((p.x - self.centre_world.x) * self.scale),
            self.centre_screen.y - screen((p.y - self.centre_world.y) * self.scale),
        )
    }

    /// Screen length of `metres`.
    pub(crate) fn length(self, metres: f64) -> f32 {
        screen(metres * self.scale)
    }
}

/// Screen distances are small enough for `f32`.
#[expect(
    clippy::cast_possible_truncation,
    reason = "screen coordinates are within a few thousand points"
)]
fn screen(value: f64) -> f32 {
    value as f32
}

/// (temperature K, colour) stops, interpolated in ln T.
const GAS_STOPS: [(f64, [u8; 3]); 5] = [
    (150.0, [40, 90, 210]),
    (300.0, [120, 160, 200]),
    (600.0, [235, 150, 70]),
    (1200.0, [240, 70, 30]),
    (2500.0, [255, 245, 210]),
];

/// Colour of gas at `temperature` K: blue when cold, grey-blue near ambient,
/// orange to near-white as it gets hot.
pub(crate) fn gas_color(temperature: f64) -> Color32 {
    let (cold, cold_rgb) = GAS_STOPS[0];
    let (hot, hot_rgb) = GAS_STOPS[GAS_STOPS.len() - 1];
    // NaN shows as the coldest colour.
    if temperature.is_nan() || temperature <= cold {
        return Color32::from_rgb(cold_rgb[0], cold_rgb[1], cold_rgb[2]);
    }
    if temperature >= hot {
        return Color32::from_rgb(hot_rgb[0], hot_rgb[1], hot_rgb[2]);
    }
    let upper = GAS_STOPS
        .iter()
        .position(|(t, _)| *t > temperature)
        .unwrap_or(1);
    let ((t0, c0), (t1, c1)) = (GAS_STOPS[upper - 1], GAS_STOPS[upper]);
    let u = (temperature / t0).ln() / (t1 / t0).ln();
    let mix = |a: u8, b: u8| channel(f64::from(a) + u * (f64::from(b) - f64::from(a)));
    Color32::from_rgb(mix(c0[0], c1[0]), mix(c0[1], c1[1]), mix(c0[2], c1[2]))
}

#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "clamped to 0..=255 first"
)]
fn channel(value: f64) -> u8 {
    value.round().clamp(0.0, 255.0) as u8
}

fn fill_color(fill: Fill) -> Color32 {
    match fill {
        Fill::None => Color32::TRANSPARENT,
        Fill::Metal => Color32::from_rgb(150, 150, 160),
        Fill::Liquid => Color32::from_rgb(60, 120, 220),
        Fill::Gas { temperature } => gas_color(temperature),
        Fill::Flame => Color32::from_rgb(255, 140, 30),
    }
}

/// Arrowhead length and half-angle.
const HEAD_LENGTH: f32 = 10.0;
const HEAD_ANGLE: f32 = 0.45;

/// The egui shapes for `scene`; `painter` lays out text in its font and
/// colours.
pub(crate) fn to_shapes(
    scene: &Scene,
    transform: &ViewTransform,
    units: &UnitSelection,
    painter: &egui::Painter,
) -> Vec<Shape> {
    let visuals = painter.ctx().global_style().visuals.clone();
    let stroke = |style: Stroke| match style {
        Stroke::None => EguiStroke::NONE,
        Stroke::Outline => EguiStroke::new(1.5, visuals.text_color()),
        Stroke::Thin => EguiStroke::new(1.0, visuals.weak_text_color()),
    };
    let mut shapes = Vec::with_capacity(scene.shapes().len() + 2);
    for shape in scene.shapes() {
        match *shape {
            SceneShape::Polygon {
                points,
                fill,
                stroke: style,
            } => {
                let points = scene
                    .points(points)
                    .iter()
                    .map(|p| transform.to_screen(*p))
                    .collect();
                shapes.push(Shape::convex_polygon(
                    points,
                    fill_color(fill),
                    stroke(style),
                ));
            }
            SceneShape::Polyline {
                points,
                stroke: style,
            } => {
                let points = scene
                    .points(points)
                    .iter()
                    .map(|p| transform.to_screen(*p))
                    .collect();
                shapes.push(Shape::line(points, stroke(style)));
            }
            SceneShape::Circle {
                center,
                radius,
                fill,
                stroke: style,
            } => {
                let (centre, radius) = (transform.to_screen(center), transform.length(radius));
                shapes.push(Shape::circle_filled(centre, radius, fill_color(fill)));
                shapes.push(Shape::circle_stroke(centre, radius, stroke(style)));
            }
            SceneShape::Arrow { from, to } => {
                let (tail, head) = (transform.to_screen(from), transform.to_screen(to));
                let line = EguiStroke::new(2.0, visuals.text_color());
                shapes.push(Shape::line_segment([tail, head], line));
                // A head longer than half the arrow would dwarf short vectors.
                let head_length = HEAD_LENGTH.min(0.5 * (head - tail).length());
                let back = (tail - head).normalized();
                for angle in [HEAD_ANGLE, -HEAD_ANGLE] {
                    let barb = Vec2::angled(back.angle() + angle) * head_length;
                    shapes.push(Shape::line_segment([head, head + barb], line));
                }
            }
            SceneShape::Label { at, label } => {
                let text = match label {
                    Label::Text(text) => text.to_owned(),
                    Label::Value { value, quantity } => format_value(value, units.unit(quantity)),
                };
                let galley =
                    painter.layout_no_wrap(text, FontId::proportional(13.0), visuals.text_color());
                let rect =
                    Align2::CENTER_CENTER.anchor_size(transform.to_screen(at), galley.size());
                shapes.push(Shape::galley(rect.min, galley, visuals.text_color()));
            }
        }
    }
    shapes
}

#[cfg(test)]
mod tests {
    use eframe::egui::{Context, LayerId, Vec2};
    use proptest::prelude::*;

    use super::*;
    use crate::units::Preset;

    fn rect(w: f32, h: f32) -> Rect {
        Rect::from_min_size(Pos2::new(10.0, 20.0), Vec2::new(w, h))
    }

    #[test]
    fn test_view_transform_fits_and_flips() -> anyhow::Result<()> {
        // A 2 m × 1 m world into a 400 × 400 view: limited by width, centred
        // vertically.
        let bounds = (Point::new(-1.0, 0.0), Point::new(1.0, 1.0));
        let view = ViewTransform::fit(bounds, rect(400.0, 400.0), 0.0);
        let left_bottom = view.to_screen(Point::new(-1.0, 0.0));
        let right_top = view.to_screen(Point::new(1.0, 1.0));
        anyhow::ensure!(
            (left_bottom.x - 10.0).abs() < 1e-3 && (right_top.x - 410.0).abs() < 1e-3,
            "x {left_bottom:?} {right_top:?}"
        );
        anyhow::ensure!(right_top.y < left_bottom.y, "y not flipped");
        anyhow::ensure!(
            (f32::midpoint(left_bottom.y, right_top.y) - 220.0).abs() < 1e-3,
            "not centred"
        );
        anyhow::ensure!(
            (view.length(1.0) - 200.0).abs() < 1e-3,
            "scale {}",
            view.length(1.0)
        );
        Ok(())
    }

    #[test]
    fn test_view_transform_survives_degenerate_bounds() -> anyhow::Result<()> {
        let point = Point::new(3.0, 3.0);
        let view = ViewTransform::fit((point, point), rect(100.0, 100.0), 10.0);
        let centre = view.to_screen(point);
        anyhow::ensure!(
            centre.x.is_finite() && (centre.x - 60.0).abs() < 1e-3,
            "{centre:?}"
        );
        Ok(())
    }

    #[test]
    fn test_gas_color_runs_cold_blue_to_hot_bright() -> anyhow::Result<()> {
        let (cold, hot) = (gas_color(150.0), gas_color(2500.0));
        anyhow::ensure!(cold.b() > cold.r(), "cold {cold:?}");
        anyhow::ensure!(hot.r() > 240 && hot.g() > 200, "hot {hot:?}");
        anyhow::ensure!(
            gas_color(f64::NAN) == gas_color(150.0),
            "NaN not shown as cold"
        );
        anyhow::ensure!(gas_color(1e6) == gas_color(2500.0), "clamped above");
        Ok(())
    }

    #[test]
    fn test_to_shapes_converts_every_scene_shape() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        scene.rectangle(
            Point::new(0.0, 0.0),
            Point::new(1.0, 1.0),
            Fill::Gas { temperature: 900.0 },
            Stroke::Outline,
        );
        scene.polyline(&[Point::new(0.0, 0.0), Point::new(1.0, 1.0)], Stroke::Thin);
        scene.circle(Point::new(0.5, 0.5), 0.1, Fill::Metal, Stroke::None);
        scene.arrow(Point::new(0.0, 0.5), Point::new(1.0, 0.5));
        scene.label(
            Point::new(0.5, 1.2),
            Label::Value {
                value: 8e5,
                quantity: engine_bench_core::spec::QuantityKind::Pressure,
            },
        );
        let units = UnitSelection::new(Preset::Engineering);
        let ctx = Context::default();
        let mut shapes = Vec::new();
        let mut output = ctx.run_pass(egui::RawInput::default(), |ctx| {
            let painter =
                egui::Painter::new(ctx.clone(), LayerId::background(), rect(300.0, 300.0));
            let bounds = scene.bounds().unwrap_or_default();
            let view = ViewTransform::fit(bounds, rect(300.0, 300.0), 10.0);
            shapes = to_shapes(&scene, &view, &units, &painter);
        });
        // No renderer uploads the font atlas here; egui asserts deltas are
        // handled.
        output.textures_delta.clear();
        anyhow::ensure!(shapes.len() >= 5, "{} shapes", shapes.len());
        let has_hot_fill = shapes
            .iter()
            .any(|s| matches!(s, Shape::Path(p) if p.fill == gas_color(900.0)));
        anyhow::ensure!(has_hot_fill, "gas fill missing");
        let label = shapes.iter().find_map(|s| match s {
            Shape::Text(text) => Some(text.galley.text().to_owned()),
            _ => None,
        });
        anyhow::ensure!(label.as_deref() == Some("8.000 bar"), "label {label:?}");
        Ok(())
    }

    #[test]
    fn test_short_arrow_keeps_its_head_inside() -> anyhow::Result<()> {
        // A 4 px arrow must not sprout fixed-size 10 px barbs.
        let mut scene = Scene::new();
        scene.arrow(Point::new(0.0, 0.0), Point::new(0.004, 0.0));
        let view = ViewTransform::fit(
            (Point::new(-1.0, -1.0), Point::new(1.0, 1.0)),
            rect(2000.0, 2000.0),
            0.0,
        );
        let units = UnitSelection::new(Preset::Si);
        let ctx = Context::default();
        let mut shapes = Vec::new();
        let mut output = ctx.run_pass(egui::RawInput::default(), |ctx| {
            let painter =
                egui::Painter::new(ctx.clone(), LayerId::background(), rect(2000.0, 2000.0));
            shapes = to_shapes(&scene, &view, &units, &painter);
        });
        output.textures_delta.clear();
        let lengths: Vec<f32> = shapes
            .iter()
            .filter_map(|s| match s {
                Shape::LineSegment { points, .. } => Some(points[0].distance(points[1])),
                _ => None,
            })
            .collect();
        let shaft = lengths.first().copied().unwrap_or_default();
        anyhow::ensure!((shaft - 4.0).abs() < 0.01, "shaft {shaft}");
        anyhow::ensure!(
            lengths[1..].iter().all(|barb| *barb <= 0.5 * shaft + 1e-3),
            "{lengths:?}"
        );
        Ok(())
    }

    /// egui panics (in debug builds) on text its bundled fonts cannot draw;
    /// every label an engine or the unit table can put on screen must be
    /// drawable.
    #[test]
    fn test_every_engine_and_unit_label_has_glyphs() -> anyhow::Result<()> {
        use engine_bench_core::{engine::EngineId, spec::QuantityKind};

        let mut texts: Vec<&str> = Vec::new();
        for id in EngineId::ALL {
            let descriptor = id.descriptor();
            texts.extend([descriptor.name, descriptor.summary]);
            texts.extend(id.controls().iter().map(|c| c.label));
            for control in id.controls() {
                if let engine_bench_core::spec::ControlKind::Choice { options, .. } = control.kind {
                    texts.extend(options.iter().map(|o| o.label));
                }
            }
            texts.extend(id.channels().iter().map(|c| c.label));
        }
        for kind in QuantityKind::ALL {
            texts.extend(crate::units::units(kind).iter().map(|u| u.symbol));
        }
        texts.extend(["×", "·", "–"]);
        let ctx = Context::default();
        let mut widths = Vec::with_capacity(texts.len());
        let mut output = ctx.run_pass(egui::RawInput::default(), |ctx| {
            let painter =
                egui::Painter::new(ctx.clone(), LayerId::background(), rect(100.0, 100.0));
            for text in &texts {
                let galley = painter.layout_no_wrap(
                    (*text).to_owned(),
                    FontId::proportional(13.0),
                    Color32::WHITE,
                );
                widths.push(galley.size().x);
            }
        });
        output.textures_delta.clear();
        for (text, width) in texts.iter().zip(&widths) {
            anyhow::ensure!(text.is_empty() || *width > 0.0, "{text:?} laid out empty");
        }
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_view_transform_keeps_bounds_inside_the_view(
            x0 in -100.0..100.0f64, y0 in -100.0..100.0f64,
            w in 0.001..100.0f64, h in 0.001..100.0f64,
            vw in 50.0..2000.0f32, vh in 50.0..2000.0f32,
        ) {
            let bounds = (Point::new(x0, y0), Point::new(x0 + w, y0 + h));
            let view_rect = rect(vw, vh);
            let view = ViewTransform::fit(bounds, view_rect, 5.0);
            for corner in [bounds.0, bounds.1, Point::new(x0, y0 + h), Point::new(x0 + w, y0)] {
                let p = view.to_screen(corner);
                prop_assert!(view_rect.expand(0.01).contains(p), "{p:?} outside {view_rect:?}");
            }
        }
    }
}
