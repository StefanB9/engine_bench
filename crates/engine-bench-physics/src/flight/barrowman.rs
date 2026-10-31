//! Barrowman's method: normal-force slope and centre of pressure of a slender
//! finned rocket at small angles of attack (J. S. Barrowman, "The Practical
//! Calculation of the Aerodynamic Characteristics of Slender Finned Vehicles",
//! 1967).
//!
//! Slopes are per radian, referenced to the body cross-section π·d²/4;
//! positions are measured aft from the nose tip. Components are kept separate
//! because pitch damping sums `slope · (position − CG)²` over them.

use anyhow::Context;

use crate::{
    PhysicsError,
    units::{Length, length::meter},
};

/// Nose-cone profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoseShape {
    /// Straight cone.
    Conical,
    /// Tangent ogive.
    Ogive,
    /// Parabolic series, K = 1.
    Parabolic,
}

/// A set of identical trapezoidal fins equally spaced around the body.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fins {
    /// Number of fins (3, 4 or 6 for the method's interference factor to hold).
    pub count: u32,
    /// Chord at the body.
    pub root_chord: Length,
    /// Chord at the tip.
    pub tip_chord: Length,
    /// Span from body surface to tip.
    pub semi_span: Length,
    /// Axial distance from the root leading edge to the tip leading edge.
    pub sweep: Length,
    /// Distance from the nose tip to the root leading edge.
    pub position: Length,
}

/// Rocket outline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RocketGeometry {
    /// Nose profile.
    pub nose_shape: NoseShape,
    /// Nose length.
    pub nose_length: Length,
    /// Body diameter (the nose base matches it).
    pub body_diameter: Length,
    /// Fins, if any.
    pub fins: Option<Fins>,
}

/// Normal-force slope and centre of pressure of one component, SI (per radian,
/// metres from the nose tip). Kernel data for aerodynamic loads.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Component {
    /// Normal-force slope `dC_N/dα`, per radian.
    pub normal_force_slope: f64,
    /// Centre of pressure, m aft of the nose tip.
    pub centre_of_pressure: f64,
}

/// Barrowman result: components (nose, then fins) and their combination.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stability {
    /// Nose, then fins (if any). Unused slots have zero slope.
    pub components: [Component; 2],
    /// Body diameter, m (sets the reference area).
    pub reference_diameter: f64,
}

impl Stability {
    /// Total normal-force slope, per radian.
    pub fn normal_force_slope(&self) -> f64 {
        self.components.iter().map(|c| c.normal_force_slope).sum()
    }

    /// Centre of pressure of the whole rocket.
    pub fn centre_of_pressure(&self) -> Length {
        let moment: f64 = self
            .components
            .iter()
            .map(|c| c.normal_force_slope * c.centre_of_pressure)
            .sum();
        Length::new::<meter>(moment / self.normal_force_slope())
    }

    /// Reference area π·d²/4, m².
    pub fn reference_area(&self) -> f64 {
        std::f64::consts::FRAC_PI_4 * self.reference_diameter * self.reference_diameter
    }
}

/// Applies Barrowman's equations to `geometry`.
pub fn barrowman(geometry: &RocketGeometry) -> anyhow::Result<Stability> {
    let diameter = positive(geometry.body_diameter, "body diameter must be positive")?;
    let nose_length = positive(geometry.nose_length, "nose length must be positive")?;
    // Every nose shape has C_Nα = 2; only the centre of pressure moves.
    let fraction = match geometry.nose_shape {
        NoseShape::Conical => 2.0 / 3.0,
        NoseShape::Ogive => 0.466,
        NoseShape::Parabolic => 0.5,
    };
    let nose = Component {
        normal_force_slope: 2.0,
        centre_of_pressure: fraction * nose_length,
    };
    let fins = match &geometry.fins {
        Some(fins) => fin_component(fins, diameter, nose_length).context("fins")?,
        None => Component {
            normal_force_slope: 0.0,
            centre_of_pressure: 0.0,
        },
    };
    Ok(Stability {
        components: [nose, fins],
        reference_diameter: diameter,
    })
}

fn fin_component(fins: &Fins, diameter: f64, nose_length: f64) -> anyhow::Result<Component> {
    if !matches!(fins.count, 3 | 4 | 6) {
        return Err(PhysicsError::InvalidGeometry(
            "Barrowman fins need 3, 4 or 6 fins",
        ))
        .with_context(|| format!("{} fins", fins.count));
    }
    let root = positive(fins.root_chord, "root chord must be positive")?;
    let span = positive(fins.semi_span, "fin span must be positive")?;
    let tip = non_negative(fins.tip_chord, "tip chord must not be negative")?;
    let sweep = non_negative(fins.sweep, "fin sweep must not be negative")?;
    let position = fins.position.get::<meter>();
    if !(position >= nose_length && position.is_finite()) {
        return Err(PhysicsError::InvalidGeometry(
            "fins must sit behind the nose",
        ))
        .with_context(|| format!("fins at {position} m, nose {nose_length} m long"));
    }
    let radius = 0.5 * diameter;
    let mid_chord = (span * span + (sweep + 0.5 * tip - 0.5 * root).powi(2)).sqrt();
    // Fin-body interference: the body boosts the fins' lift.
    let interference = 1.0 + radius / (span + radius);
    let fins_alone = 4.0 * f64::from(fins.count) * (span / diameter).powi(2)
        / (1.0 + (1.0 + (2.0 * mid_chord / (root + tip)).powi(2)).sqrt());
    let chord_sum = root + tip;
    let offset =
        sweep * (root + 2.0 * tip) / (3.0 * chord_sum) + (chord_sum - root * tip / chord_sum) / 6.0;
    Ok(Component {
        normal_force_slope: interference * fins_alone,
        centre_of_pressure: position + offset,
    })
}

fn positive(length: Length, message: &'static str) -> anyhow::Result<f64> {
    let value = length.get::<meter>();
    if value > 0.0 && value.is_finite() {
        Ok(value)
    } else {
        Err(PhysicsError::InvalidGeometry(message)).with_context(|| format!("{value} m"))
    }
}

fn non_negative(length: Length, message: &'static str) -> anyhow::Result<f64> {
    let value = length.get::<meter>();
    if value >= 0.0 && value.is_finite() {
        Ok(value)
    } else {
        Err(PhysicsError::InvalidGeometry(message)).with_context(|| format!("{value} m"))
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    fn m(value: f64) -> Length {
        Length::new::<meter>(value)
    }

    fn fins(position: f64) -> Fins {
        Fins {
            count: 4,
            root_chord: m(0.1),
            tip_chord: m(0.05),
            semi_span: m(0.1),
            sweep: m(0.05),
            position: m(position),
        }
    }

    fn geometry(nose_shape: NoseShape, fins: Option<Fins>) -> RocketGeometry {
        RocketGeometry {
            nose_shape,
            nose_length: m(0.3),
            body_diameter: m(0.1),
            fins,
        }
    }

    #[test]
    fn test_barrowman_nose_alone() -> anyhow::Result<()> {
        for (shape, fraction) in [
            (NoseShape::Conical, 2.0 / 3.0),
            (NoseShape::Ogive, 0.466),
            (NoseShape::Parabolic, 0.5),
        ] {
            let s = barrowman(&geometry(shape, None))?;
            anyhow::ensure!(
                (s.normal_force_slope() - 2.0).abs() < 1e-12,
                "{shape:?} slope"
            );
            let cp = s.centre_of_pressure().get::<meter>();
            anyhow::ensure!((cp - fraction * 0.3).abs() < 1e-12, "{shape:?} CP {cp}");
        }
        Ok(())
    }

    /// Worked by hand from the published formulas (see the module
    /// documentation).
    #[test]
    fn test_barrowman_fins_match_hand_calculation() -> anyhow::Result<()> {
        let s = barrowman(&geometry(NoseShape::Conical, Some(fins(1.0))))?;
        let fin = s.components[1];
        anyhow::ensure!((fin.normal_force_slope - 7.902_190).abs() < 1e-5, "{fin:?}");
        anyhow::ensure!((fin.centre_of_pressure - 1.041_667).abs() < 1e-5, "{fin:?}");
        let total = 2.0 + fin.normal_force_slope;
        let cp = (2.0 * 0.2 + fin.normal_force_slope * fin.centre_of_pressure) / total;
        anyhow::ensure!(
            (s.normal_force_slope() - total).abs() < 1e-12,
            "total slope"
        );
        anyhow::ensure!(
            (s.centre_of_pressure().get::<meter>() - cp).abs() < 1e-12,
            "total CP"
        );
        anyhow::ensure!(
            (s.reference_area() - std::f64::consts::PI * 0.0025).abs() < 1e-15,
            "area"
        );
        Ok(())
    }

    #[test]
    fn test_barrowman_rejects_bad_geometry() -> anyhow::Result<()> {
        let mut no_body = geometry(NoseShape::Ogive, None);
        no_body.body_diameter = m(0.0);
        anyhow::ensure!(barrowman(&no_body).is_err(), "zero diameter");
        let mut fins_in_nose = fins(0.1);
        anyhow::ensure!(
            barrowman(&geometry(NoseShape::Ogive, Some(fins_in_nose))).is_err(),
            "fins on the nose"
        );
        fins_in_nose = fins(1.0);
        fins_in_nose.count = 2;
        anyhow::ensure!(
            barrowman(&geometry(NoseShape::Ogive, Some(fins_in_nose))).is_err(),
            "2 fins"
        );
        let mut negative = fins(1.0);
        negative.tip_chord = m(-0.01);
        anyhow::ensure!(
            barrowman(&geometry(NoseShape::Ogive, Some(negative))).is_err(),
            "negative chord"
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_barrowman_fins_move_the_centre_of_pressure_aft(
            position in 0.4..3.0f64, span in 0.01..0.3f64, root in 0.02..0.3f64,
            tip_fraction in 0.0..1.0f64, sweep in 0.0..0.3f64, count in prop::sample::select(vec![3u32, 4, 6]),
        ) {
            let fin = Fins {
                count,
                root_chord: m(root),
                tip_chord: m(root * tip_fraction),
                semi_span: m(span),
                sweep: m(sweep),
                position: m(position),
            };
            let bare = barrowman(&geometry(NoseShape::Ogive, None))
                .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let finned = barrowman(&geometry(NoseShape::Ogive, Some(fin)))
                .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let (nose_cp, total_cp) = (bare.centre_of_pressure().get::<meter>(), finned.centre_of_pressure().get::<meter>());
            let fin_cp = finned.components[1].centre_of_pressure;
            prop_assert!(finned.components[1].normal_force_slope > 0.0);
            prop_assert!(nose_cp < total_cp && total_cp < fin_cp, "{nose_cp} {total_cp} {fin_cp}");
            prop_assert!(fin_cp >= position && fin_cp <= position + root.max(sweep + root * tip_fraction));
        }
    }
}
