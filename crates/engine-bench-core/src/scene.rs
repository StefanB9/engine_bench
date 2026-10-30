//! Renderer-agnostic drawing of an engine: a display list in world coordinates.
//!
//! World units are metres, y pointing up. Fills are semantic (metal, gas at a
//! temperature) and labels can carry a value with its quantity kind, so the
//! renderer picks colours and display units. A `Scene` is reused between
//! frames: [`Scene::clear`] keeps its capacity.

use crate::spec::QuantityKind;

/// A point in world coordinates, metres.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point {
    /// Horizontal, m.
    pub x: f64,
    /// Vertical (up), m.
    pub y: f64,
}

impl Point {
    /// The point (`x`, `y`).
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

/// What fills a closed shape; the renderer maps it to colours.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Fill {
    /// Not filled.
    None,
    /// Solid structure (walls, pistons, casings).
    Metal,
    /// A liquid (water, propellant).
    Liquid,
    /// Gas, coloured by its temperature in K.
    Gas {
        /// Temperature, K.
        temperature: f64,
    },
    /// Combustion zone.
    Flame,
}

/// How a shape is outlined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stroke {
    /// No outline.
    None,
    /// Normal outline.
    Outline,
    /// Thin, subdued line (guides, dimensions).
    Thin,
}

/// Text placed in the scene.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Label {
    /// Fixed text.
    Text(&'static str),
    /// A value in SI units, formatted by the renderer in the chosen display
    /// unit.
    Value {
        /// Value, SI.
        value: f64,
        /// What it measures.
        quantity: QuantityKind,
    },
}

/// A run of points in the scene's shared point buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PointRange {
    start: usize,
    len: usize,
}

/// One drawable item.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Shape {
    /// Closed polygon.
    Polygon {
        /// Vertices.
        points: PointRange,
        /// Interior.
        fill: Fill,
        /// Outline.
        stroke: Stroke,
    },
    /// Open line through points.
    Polyline {
        /// Vertices.
        points: PointRange,
        /// Line style.
        stroke: Stroke,
    },
    /// Circle.
    Circle {
        /// Centre.
        center: Point,
        /// Radius, m.
        radius: f64,
        /// Interior.
        fill: Fill,
        /// Outline.
        stroke: Stroke,
    },
    /// Arrow for a vector quantity (flow, force), drawn from tail to head.
    Arrow {
        /// Tail.
        from: Point,
        /// Head.
        to: Point,
    },
    /// Text anchored at a point.
    Label {
        /// Anchor (centre of the text).
        at: Point,
        /// Content.
        label: Label,
    },
}

/// A reusable display list.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Scene {
    points: Vec<Point>,
    shapes: Vec<Shape>,
}

impl Scene {
    /// An empty scene.
    pub fn new() -> Self {
        Self::default()
    }

    /// Removes all shapes, keeping allocated capacity for the next frame.
    pub fn clear(&mut self) {
        self.points.clear();
        self.shapes.clear();
    }

    fn push_points(&mut self, points: &[Point]) -> PointRange {
        let start = self.points.len();
        self.points.extend_from_slice(points);
        PointRange {
            start,
            len: points.len(),
        }
    }

    /// Adds a closed polygon.
    pub fn polygon(&mut self, points: &[Point], fill: Fill, stroke: Stroke) {
        let points = self.push_points(points);
        self.shapes.push(Shape::Polygon {
            points,
            fill,
            stroke,
        });
    }

    /// Adds an axis-aligned rectangle between corners `min` and `max`.
    pub fn rectangle(&mut self, min: Point, max: Point, fill: Fill, stroke: Stroke) {
        let corners = [min, Point::new(max.x, min.y), max, Point::new(min.x, max.y)];
        self.polygon(&corners, fill, stroke);
    }

    /// Adds an open polyline.
    pub fn polyline(&mut self, points: &[Point], stroke: Stroke) {
        let points = self.push_points(points);
        self.shapes.push(Shape::Polyline { points, stroke });
    }

    /// Adds a circle.
    pub fn circle(&mut self, center: Point, radius: f64, fill: Fill, stroke: Stroke) {
        self.shapes.push(Shape::Circle {
            center,
            radius,
            fill,
            stroke,
        });
    }

    /// Adds an arrow from `from` to `to`.
    pub fn arrow(&mut self, from: Point, to: Point) {
        self.shapes.push(Shape::Arrow { from, to });
    }

    /// Adds a label centred at `at`.
    pub fn label(&mut self, at: Point, label: Label) {
        self.shapes.push(Shape::Label { at, label });
    }

    /// Shapes in drawing order.
    pub fn shapes(&self) -> &[Shape] {
        &self.shapes
    }

    /// The points of a polygon or polyline.
    pub fn points(&self, range: PointRange) -> &[Point] {
        self.points
            .get(range.start..range.start + range.len)
            .unwrap_or_default()
    }

    /// Smallest axis-aligned box containing every shape, as (min, max); `None`
    /// when empty. Labels count by their anchor only (their size depends on the
    /// renderer).
    pub fn bounds(&self) -> Option<(Point, Point)> {
        let mut bounds: Option<(Point, Point)> = None;
        let mut include = |p: Point, r: f64| {
            let (lo, hi) = (Point::new(p.x - r, p.y - r), Point::new(p.x + r, p.y + r));
            bounds = Some(match bounds {
                None => (lo, hi),
                Some((min, max)) => (
                    Point::new(min.x.min(lo.x), min.y.min(lo.y)),
                    Point::new(max.x.max(hi.x), max.y.max(hi.y)),
                ),
            });
        };
        for shape in &self.shapes {
            match *shape {
                Shape::Polygon { points, .. } | Shape::Polyline { points, .. } => {
                    for p in self.points(points) {
                        include(*p, 0.0);
                    }
                }
                Shape::Circle { center, radius, .. } => include(center, radius),
                Shape::Arrow { from, to } => {
                    include(from, 0.0);
                    include(to, 0.0);
                }
                Shape::Label { at, .. } => include(at, 0.0),
            }
        }
        bounds
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    #[test]
    fn test_scene_records_shapes_in_order() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        let triangle = [
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(0.0, 1.0),
        ];
        scene.polygon(&triangle, Fill::Metal, Stroke::Outline);
        scene.polyline(&[Point::new(2.0, 2.0), Point::new(3.0, 2.0)], Stroke::Thin);
        scene.circle(
            Point::new(5.0, 5.0),
            0.5,
            Fill::Gas { temperature: 300.0 },
            Stroke::None,
        );
        scene.arrow(Point::new(0.0, 0.0), Point::new(1.0, 1.0));
        let value = Label::Value {
            value: 1e5,
            quantity: QuantityKind::Pressure,
        };
        scene.label(Point::new(0.5, 0.5), value);

        let shapes = scene.shapes();
        anyhow::ensure!(shapes.len() == 5, "{} shapes", shapes.len());
        let Shape::Polygon {
            points,
            fill: Fill::Metal,
            stroke: Stroke::Outline,
        } = shapes[0]
        else {
            anyhow::bail!("first shape {:?}", shapes[0]);
        };
        anyhow::ensure!(scene.points(points) == triangle, "polygon points");
        let Shape::Polyline { points, .. } = shapes[1] else {
            anyhow::bail!("second shape {:?}", shapes[1]);
        };
        anyhow::ensure!(scene.points(points).len() == 2, "polyline points");
        anyhow::ensure!(
            matches!(shapes[4], Shape::Label { label, .. } if label == value),
            "label"
        );
        Ok(())
    }

    #[test]
    fn test_scene_rectangle_has_four_corners() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        scene.rectangle(
            Point::new(-1.0, 0.0),
            Point::new(1.0, 2.0),
            Fill::Liquid,
            Stroke::None,
        );
        let Some(Shape::Polygon {
            points,
            fill: Fill::Liquid,
            ..
        }) = scene.shapes().first()
        else {
            anyhow::bail!("{:?}", scene.shapes());
        };
        let corners = scene.points(*points);
        anyhow::ensure!(corners.len() == 4, "{corners:?}");
        anyhow::ensure!(
            corners.contains(&Point::new(-1.0, 0.0)) && corners.contains(&Point::new(1.0, 2.0)),
            "{corners:?}"
        );
        Ok(())
    }

    #[test]
    fn test_scene_clear_keeps_capacity() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        for i in 0..100 {
            scene.circle(Point::new(f64::from(i), 0.0), 1.0, Fill::None, Stroke::Thin);
            scene.polyline(&[Point::default(); 8], Stroke::Thin);
        }
        let capacity = (scene.points.capacity(), scene.shapes.capacity());
        scene.clear();
        anyhow::ensure!(
            scene.shapes().is_empty() && scene.bounds().is_none(),
            "not cleared"
        );
        anyhow::ensure!(
            (scene.points.capacity(), scene.shapes.capacity()) == capacity,
            "capacity dropped"
        );
        Ok(())
    }

    #[test]
    fn test_scene_bounds_cover_every_shape() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        anyhow::ensure!(scene.bounds().is_none(), "empty scene has bounds");
        scene.circle(Point::new(0.0, 0.0), 2.0, Fill::None, Stroke::Outline);
        scene.arrow(Point::new(1.0, 1.0), Point::new(5.0, 1.0));
        scene.label(Point::new(0.0, -4.0), Label::Text("A"));
        scene.polyline(&[Point::new(-3.0, 0.5)], Stroke::Thin);
        let bounds = scene.bounds();
        anyhow::ensure!(
            bounds == Some((Point::new(-3.0, -4.0), Point::new(5.0, 2.0))),
            "{bounds:?}"
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_scene_bounds_contain_all_points(
            raw in proptest::collection::vec((-1e3..1e3f64, -1e3..1e3f64), 1..40),
        ) {
            let points: Vec<Point> = raw.iter().map(|&(x, y)| Point::new(x, y)).collect();
            let mut scene = Scene::new();
            scene.polygon(&points, Fill::None, Stroke::Outline);
            let Some((min, max)) = scene.bounds() else {
                return Err(TestCaseError::fail("no bounds"));
            };
            for p in &points {
                prop_assert!(min.x <= p.x && p.x <= max.x && min.y <= p.y && p.y <= max.y);
            }
        }
    }
}
