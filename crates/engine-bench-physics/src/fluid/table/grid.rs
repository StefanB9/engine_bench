//! Uniform grids with cubic Lagrange interpolation: four nodes per axis, so
//! the interpolant is exact for cubics and its error is O(h⁴). A stencil that
//! touches an invalid (NaN) node or leaves the grid gives `None`.

/// A uniform axis: `len` nodes from `start` in steps of `step`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Axis {
    start: f64,
    step: f64,
    len: usize,
}

/// Where a coordinate falls: the first of its four stencil nodes, the
/// weights of those nodes, and the weights of the derivative.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Stencil {
    pub(super) first: usize,
    pub(super) weights: [f64; 4],
    pub(super) slopes: [f64; 4],
}

#[expect(clippy::cast_precision_loss, reason = "node counts are far below 2^52")]
impl Axis {
    /// `len` (≥ 4) nodes evenly spaced from `start` to `end`.
    pub(super) fn new(start: f64, end: f64, len: usize) -> Self {
        let len = len.max(4);
        Self {
            start,
            step: (end - start) / (len - 1) as f64,
            len,
        }
    }

    pub(super) fn len(self) -> usize {
        self.len
    }

    pub(super) fn step(self) -> f64 {
        self.step
    }

    /// The node at or below `x`, if on the axis.
    pub(super) fn below(self, x: f64) -> Option<usize> {
        let u = ((x - self.start) / self.step).floor();
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "u is checked to lie on the axis"
        )]
        (u >= 0.0 && u < (self.len - 1) as f64).then_some(u as usize)
    }

    /// Coordinate of node `i`.
    pub(super) fn node(self, i: usize) -> f64 {
        (i as f64).mul_add(self.step, self.start)
    }

    /// The stencil for `x`, if all four nodes lie on the axis.
    #[expect(
        clippy::many_single_char_names,
        reason = "the stencil's local coordinates"
    )]
    pub(super) fn stencil(self, x: f64) -> Option<Stencil> {
        let u = (x - self.start) / self.step;
        // Cell i spans nodes i and i + 1; its stencil is i − 1 ..= i + 2.
        let cell = u.floor();
        let last_cell = (self.len - 3) as f64;
        if !(cell >= 1.0 && cell <= last_cell) {
            return None;
        }
        let t = u - cell;
        let (a, b, c, d) = (t + 1.0, t, t - 1.0, t - 2.0);
        let weights = [
            -b * c * d / 6.0,
            a * c * d / 2.0,
            -a * b * d / 2.0,
            a * b * c / 6.0,
        ];
        // d/dt of each product, divided by the step for d/dx.
        let slopes = [
            -(c * d + b * d + b * c) / 6.0,
            (c * d + a * d + a * c) / 2.0,
            -(b * d + a * d + a * b) / 2.0,
            (b * c + a * c + a * b) / 6.0,
        ]
        .map(|s| s / self.step);
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "cell is checked to lie in 1 ..= len − 3"
        )]
        let first = cell as usize - 1;
        Some(Stencil {
            first,
            weights,
            slopes,
        })
    }
}

/// `K` quantities on a two-dimensional grid, row-major (x outer, y inner).
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Grid2<const K: usize> {
    pub(super) x: Axis,
    pub(super) y: Axis,
    values: Vec<[f64; K]>,
}

/// Interpolated values and their derivatives along y.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Sample<const K: usize> {
    pub(super) value: [f64; K],
    pub(super) dy: [f64; K],
}

impl<const K: usize> Grid2<K> {
    /// A grid filled row by row from `rows`, each `y.len()` long.
    pub(super) fn from_rows(x: Axis, y: Axis, rows: Vec<Vec<[f64; K]>>) -> Self {
        let mut values = Vec::with_capacity(x.len() * y.len());
        for row in rows {
            values.extend(row);
        }
        Self { x, y, values }
    }

    /// Values at (x, y), or `None` if the stencil is off the grid or touches
    /// an invalid node.
    pub(super) fn sample(&self, x: f64, y: f64) -> Option<Sample<K>> {
        let (sx, sy) = (self.x.stencil(x)?, self.y.stencil(y)?);
        let mut value = [0.0; K];
        let mut dy = [0.0; K];
        for (i, wx) in sx.weights.iter().enumerate() {
            let row = (sx.first + i) * self.y.len() + sy.first;
            let nodes = self.values.get(row..row + 4)?;
            for ((node, wy), sl) in nodes.iter().zip(sy.weights).zip(sy.slopes) {
                for k in 0..K {
                    if node[k].is_nan() {
                        return None;
                    }
                    value[k] += wx * wy * node[k];
                    dy[k] += wx * sl * node[k];
                }
            }
        }
        Some(Sample { value, dy })
    }
}

/// `K` quantities on a one-dimensional grid.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Grid1<const K: usize> {
    pub(super) x: Axis,
    values: Vec<[f64; K]>,
}

impl<const K: usize> Grid1<K> {
    pub(super) fn new(x: Axis, values: Vec<[f64; K]>) -> Self {
        Self { x, values }
    }

    /// Values at `x` and their derivatives, or `None` off the grid.
    pub(super) fn sample(&self, x: f64) -> Option<Sample<K>> {
        let s = self.x.stencil(x)?;
        let nodes = self.values.get(s.first..s.first + 4)?;
        let mut value = [0.0; K];
        let mut dy = [0.0; K];
        for ((node, w), sl) in nodes.iter().zip(s.weights).zip(s.slopes) {
            for k in 0..K {
                value[k] += w * node[k];
                dy[k] += sl * node[k];
            }
        }
        value
            .iter()
            .all(|v| v.is_finite())
            .then_some(Sample { value, dy })
    }

    /// Node values.
    pub(super) fn nodes(&self) -> &[[f64; K]] {
        &self.values
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    fn cubic(x: f64, y: f64) -> f64 {
        1.5 - 2.0 * x + 0.5 * y + 0.3 * x * x * y - 0.2 * y * y * y + 0.1 * x * x * x
    }

    fn grid() -> Grid2<1> {
        let (x, y) = (Axis::new(-1.0, 2.0, 13), Axis::new(0.5, 4.0, 9));
        let rows = (0..x.len())
            .map(|i| {
                (0..y.len())
                    .map(|j| [cubic(x.node(i), y.node(j))])
                    .collect()
            })
            .collect();
        Grid2::from_rows(x, y, rows)
    }

    #[test]
    fn test_axis_has_evenly_spaced_nodes() {
        let axis = Axis::new(1.0, 3.0, 5);
        assert_eq!(axis.len(), 5);
        assert!((axis.node(4) - 3.0).abs() < 1e-15);
        assert!((axis.node(2) - 2.0).abs() < 1e-15);
    }

    #[test]
    fn test_axis_stencil_needs_four_nodes() {
        let axis = Axis::new(0.0, 10.0, 11);
        // The first and last cells lack a node on one side.
        assert!(axis.stencil(0.5).is_none());
        assert!(axis.stencil(9.5).is_none());
        assert!(axis.stencil(1.0).is_some());
        assert!(axis.stencil(8.99).is_some());
        assert!(axis.stencil(f64::NAN).is_none());
    }

    #[test]
    fn test_grid2_rejects_stencils_touching_invalid_nodes() -> anyhow::Result<()> {
        let (x, y) = (Axis::new(0.0, 9.0, 10), Axis::new(0.0, 9.0, 10));
        let mut rows: Vec<Vec<[f64; 1]>> = vec![vec![[1.0]; 10]; 10];
        if let Some(cell) = rows.get_mut(3).and_then(|r| r.get_mut(3)) {
            *cell = [f64::NAN];
        }
        let g = Grid2::from_rows(x, y, rows);
        anyhow::ensure!(g.sample(2.5, 2.5).is_none(), "touches (3, 3)");
        anyhow::ensure!(g.sample(6.5, 2.5).is_some(), "clear of (3, 3)");
        Ok(())
    }

    #[test]
    fn test_grid1_interpolates_a_cubic_and_its_slope() -> anyhow::Result<()> {
        let x = Axis::new(0.0, 3.0, 7);
        let f = |x: f64| 2.0 - x + 0.5 * x * x * x;
        let g = Grid1::new(x, (0..7).map(|i| [f(x.node(i))]).collect());
        let s = g.sample(1.3).ok_or_else(|| anyhow::anyhow!("off grid"))?;
        anyhow::ensure!((s.value[0] - f(1.3)).abs() < 1e-12, "{s:?}");
        anyhow::ensure!((s.dy[0] - (-1.0 + 1.5 * 1.3 * 1.3)).abs() < 1e-12, "{s:?}");
        anyhow::ensure!(g.nodes().len() == 7, "nodes");
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        /// Tensor-product cubic Lagrange reproduces a bicubic polynomial and
        /// its y-derivative exactly.
        #[test]
        fn test_grid2_is_exact_for_cubics(x in -0.75f64..1.75, y in 0.94f64..3.56) {
            let g = grid();
            let s = g.sample(x, y);
            prop_assert!(s.is_some(), "off grid at ({}, {})", x, y);
            if let Some(s) = s {
                prop_assert!((s.value[0] - cubic(x, y)).abs() < 1e-12, "{} vs {}", s.value[0], cubic(x, y));
                let dy = 0.5 + 0.3 * x * x - 0.6 * y * y;
                prop_assert!((s.dy[0] - dy).abs() < 1e-11, "{} vs {}", s.dy[0], dy);
            }
        }
    }
}
