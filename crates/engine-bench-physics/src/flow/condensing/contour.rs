//! Nozzle contours along the axis: the area a condensing flow expands through,
//! from the throat (x = 0) to the exit.
//!
//! Condensation is a rate process, so where it starts depends on how fast the
//! area grows, not only on the exit area. Both shapes have a circular-arc
//! throat (radius of curvature `R`) joined tangentially to a straight
//! diverging wall at half-angle θ: the profile y(x) is y* + R − √(R² − x²) up
//! to x = R·sin θ, then a straight line.

use anyhow::Context;

use crate::{
    PhysicsError,
    units::{Angle, Area, Length, angle::radian, area::square_meter, length::meter},
};

/// How the profile y(x) makes an area.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Shape {
    /// A body of revolution: y is the radius.
    Axisymmetric,
    /// A planar channel of constant `depth` m: y is the half-height.
    Planar { depth: f64 },
}

/// A nozzle's diverging contour: circular-arc throat into a straight wall.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Contour {
    shape: Shape,
    /// Throat half-size (radius or half-height), m.
    throat: f64,
    /// Exit half-size, m.
    exit: f64,
    /// Throat radius of curvature, m.
    curvature: f64,
    /// Wall half-angle, rad.
    angle: f64,
}

impl Contour {
    /// A conical nozzle: throat and exit diameters, throat radius of
    /// curvature, diverging half-angle (0° to 45°, exclusive).
    pub fn conical(
        throat_diameter: Length,
        exit_diameter: Length,
        throat_curvature: Length,
        half_angle: Angle,
    ) -> anyhow::Result<Self> {
        Self::new(
            Shape::Axisymmetric,
            0.5 * throat_diameter.get::<meter>(),
            0.5 * exit_diameter.get::<meter>(),
            throat_curvature.get::<meter>(),
            half_angle.get::<radian>(),
        )
        .context("conical nozzle")
    }

    /// A planar nozzle of constant `depth`: throat and exit heights, throat
    /// radius of curvature, wall angle.
    pub fn planar(
        depth: Length,
        throat_height: Length,
        exit_height: Length,
        throat_curvature: Length,
        wall_angle: Angle,
    ) -> anyhow::Result<Self> {
        let depth = depth.get::<meter>();
        if !(depth > 0.0 && depth.is_finite()) {
            return Err(PhysicsError::InvalidGeometry("depth must be positive"))
                .with_context(|| format!("depth = {depth} m"));
        }
        Self::new(
            Shape::Planar { depth },
            0.5 * throat_height.get::<meter>(),
            0.5 * exit_height.get::<meter>(),
            throat_curvature.get::<meter>(),
            wall_angle.get::<radian>(),
        )
        .context("planar nozzle")
    }

    fn new(
        shape: Shape,
        throat: f64,
        exit: f64,
        curvature: f64,
        angle: f64,
    ) -> anyhow::Result<Self> {
        let positive = |v: f64| v > 0.0 && v.is_finite();
        if !(positive(throat) && positive(curvature)) {
            return Err(PhysicsError::InvalidGeometry(
                "throat and its radius of curvature must be positive",
            )
            .into());
        }
        if !(exit >= throat && exit.is_finite()) {
            return Err(PhysicsError::InvalidGeometry("exit narrower than the throat").into());
        }
        if !(angle > 0.0 && angle < std::f64::consts::FRAC_PI_4) {
            return Err(PhysicsError::InvalidGeometry(
                "half-angle must lie in (0°, 45°)",
            ))
            .with_context(|| format!("{angle} rad"));
        }
        Ok(Self {
            shape,
            throat,
            exit,
            curvature,
            angle,
        })
    }

    /// Distance from the throat to the exit.
    pub fn length(&self) -> Length {
        Length::new::<meter>(self.length_si())
    }

    /// Area at distance `x` downstream of the throat (clamped to the nozzle).
    pub fn area(&self, x: Length) -> Area {
        Area::new::<square_meter>(self.area_si(x.get::<meter>()).0)
    }

    /// Throat area.
    pub fn throat_area(&self) -> Area {
        Area::new::<square_meter>(self.area_from(self.throat))
    }

    /// Exit area.
    pub fn exit_area(&self) -> Area {
        Area::new::<square_meter>(self.area_from(self.exit))
    }

    /// Where the throat arc meets the straight wall, m (or where the arc
    /// already reaches the exit).
    fn tangent(&self) -> (f64, f64) {
        let (r, theta) = (self.curvature, self.angle);
        let x = r * theta.sin();
        let y = self.throat + r * (1.0 - theta.cos());
        if y >= self.exit {
            // The arc reaches the exit first.
            let dy = self.exit - self.throat;
            ((r * r - (r - dy) * (r - dy)).sqrt(), self.exit)
        } else {
            (x, y)
        }
    }

    pub(crate) fn length_si(&self) -> f64 {
        let (x_t, y_t) = self.tangent();
        x_t + (self.exit - y_t) / self.angle.tan()
    }

    /// Profile y(x) and dy/dx, x clamped to [0, length].
    fn profile(&self, x: f64) -> (f64, f64) {
        let x = x.clamp(0.0, self.length_si());
        let (x_t, y_t) = self.tangent();
        if x <= x_t {
            let root = (self.curvature * self.curvature - x * x).sqrt();
            (self.throat + self.curvature - root, x / root)
        } else {
            let slope = self.angle.tan();
            (y_t + (x - x_t) * slope, slope)
        }
    }

    fn area_from(&self, y: f64) -> f64 {
        match self.shape {
            Shape::Axisymmetric => std::f64::consts::PI * y * y,
            Shape::Planar { depth } => 2.0 * y * depth,
        }
    }

    /// Area, m², and dA/dx, m, at `x` m downstream of the throat.
    pub(crate) fn area_si(&self, x: f64) -> (f64, f64) {
        let (y, dy) = self.profile(x);
        match self.shape {
            Shape::Axisymmetric => (
                std::f64::consts::PI * y * y,
                2.0 * std::f64::consts::PI * y * dy,
            ),
            Shape::Planar { depth } => (2.0 * y * depth, 2.0 * depth * dy),
        }
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::units::{angle::degree, length::millimeter};

    fn mm(v: f64) -> Length {
        Length::new::<millimeter>(v)
    }

    fn cone() -> anyhow::Result<Contour> {
        Contour::conical(mm(1.0), mm(7.0), mm(1.5), Angle::new::<degree>(15.0))
    }

    #[test]
    fn test_contour_starts_at_the_throat_and_ends_at_the_exit() -> anyhow::Result<()> {
        let c = cone()?;
        let at = |x: f64| c.area(Length::new::<meter>(x)).get::<square_meter>();
        let throat = std::f64::consts::PI * 0.25e-6;
        let exit = std::f64::consts::PI * 12.25e-6;
        anyhow::ensure!((at(0.0) / throat - 1.0).abs() < 1e-12, "throat {}", at(0.0));
        let l = c.length().get::<meter>();
        anyhow::ensure!((at(l) / exit - 1.0).abs() < 1e-12, "exit {}", at(l));
        anyhow::ensure!(
            (c.exit_area().get::<square_meter>() / exit - 1.0).abs() < 1e-12,
            "exit"
        );
        anyhow::ensure!(
            (c.throat_area().get::<square_meter>() / throat - 1.0).abs() < 1e-12,
            "throat"
        );
        // Beyond the ends the area is clamped.
        anyhow::ensure!(
            (at(-1.0) - at(0.0)).abs() < 1e-18 && (at(1.0) - at(l)).abs() < 1e-18,
            "clamp"
        );
        Ok(())
    }

    #[test]
    fn test_contour_length_follows_the_cone() -> anyhow::Result<()> {
        // Arc to x = R sin θ, then the cone: (r_e − y_t) / tan θ.
        let c = cone()?;
        let theta = 15f64.to_radians();
        let (r, rt, re) = (1.5e-3, 0.5e-3, 3.5e-3);
        let (xt, yt) = (r * theta.sin(), rt + r * (1.0 - theta.cos()));
        let want = xt + (re - yt) / theta.tan();
        let got = c.length().get::<meter>();
        anyhow::ensure!((got / want - 1.0).abs() < 1e-12, "{got} vs {want}");
        Ok(())
    }

    #[test]
    fn test_planar_contour_reproduces_the_yale_nozzle() -> anyhow::Result<()> {
        // Yale nozzle B: h* = 4.0 mm, b = 7.0 mm, R* = 20 mm, divergence
        // α = 6.8° between the walls (3.4° each), h = 13.3 mm: A* = 0.28 cm²,
        // A/A* = 3.33 and l = 8.0 cm (the length fixes how α is meant).
        let c = Contour::planar(
            mm(7.0),
            mm(4.0),
            mm(13.3),
            mm(20.0),
            Angle::new::<degree>(3.4),
        )?;
        let throat = c.throat_area().get::<square_meter>();
        anyhow::ensure!((throat / 0.28e-4 - 1.0).abs() < 1e-12, "A* {throat}");
        let ratio = c.exit_area().get::<square_meter>() / throat;
        anyhow::ensure!((ratio - 3.325).abs() < 1e-12, "A/A* {ratio}");
        let l = c.length().get::<meter>();
        anyhow::ensure!((l - 0.08).abs() < 0.002, "length {l} m");
        Ok(())
    }

    #[test]
    fn test_contour_rejects_impossible_geometry() {
        let deg = Angle::new::<degree>;
        assert!(Contour::conical(mm(2.0), mm(1.0), mm(1.0), deg(15.0)).is_err());
        assert!(Contour::conical(mm(1.0), mm(2.0), mm(0.0), deg(15.0)).is_err());
        assert!(Contour::conical(mm(1.0), mm(2.0), mm(1.0), deg(50.0)).is_err());
        assert!(Contour::conical(mm(1.0), mm(2.0), mm(1.0), deg(0.0)).is_err());
        assert!(Contour::planar(mm(0.0), mm(1.0), mm(2.0), mm(1.0), deg(5.0)).is_err());
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        /// The area grows monotonically and smoothly (the wall's slope is
        /// continuous where the arc meets the cone), and dA/dx matches the area.
        #[test]
        fn test_contour_area_is_smooth_and_growing(
            exit in 1.05f64..20.0,
            curvature in 0.3f64..5.0,
            angle in 2.0f64..40.0,
            f in 0.0f64..1.0,
        ) {
            let c = Contour::conical(mm(1.0), mm(exit), mm(curvature), Angle::new::<degree>(angle))
                .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let l = c.length_si();
            let x = f * l;
            let (a, da) = c.area_si(x);
            prop_assert!(da >= 0.0 && a.is_finite());
            let h = 1e-7 * l;
            if x > h && x < l - h {
                let numeric = (c.area_si(x + h).0 - c.area_si(x - h).0) / (2.0 * h);
                prop_assert!((numeric - da).abs() <= 1e-4 * da.max(1e-9), "{} vs {}", numeric, da);
            }
            if x < l - h {
                prop_assert!(c.area_si(x + h).0 >= a, "not growing at {}", x);
            }
        }
    }
}
