//! Tolerances and step limits for the adaptive integrator.

use anyhow::Context;

use crate::PhysicsError;

/// Tolerances and step limits for [`DormandPrince45`](super::DormandPrince45).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveOptions {
    pub(super) rtol: f64,
    pub(super) atol: f64,
    pub(super) max_step: f64,
    pub(super) initial_step: Option<f64>,
}

impl AdaptiveOptions {
    /// Per-component error target `atol + rtol·|y|`. `rtol` must be positive,
    /// `atol` non-negative.
    pub fn new(rtol: f64, atol: f64) -> anyhow::Result<Self> {
        if !(rtol > 0.0 && rtol.is_finite()) {
            return Err(PhysicsError::InvalidSetting(
                "rtol must be positive and finite",
            ))
            .with_context(|| format!("rtol = {rtol}"));
        }
        if !(atol >= 0.0 && atol.is_finite()) {
            return Err(PhysicsError::InvalidSetting(
                "atol must be non-negative and finite",
            ))
            .with_context(|| format!("atol = {atol}"));
        }
        Ok(Self {
            rtol,
            atol,
            max_step: f64::INFINITY,
            initial_step: None,
        })
    }

    /// Caps every step at `max_step` (default: unlimited). Use it so events
    /// cannot be stepped over when `g` may cross zero twice within one step.
    pub fn with_max_step(self, max_step: f64) -> anyhow::Result<Self> {
        if max_step.is_nan() || max_step <= 0.0 {
            return Err(PhysicsError::InvalidSetting("max_step must be positive"))
                .with_context(|| format!("max_step = {max_step}"));
        }
        Ok(Self { max_step, ..self })
    }

    /// Uses `initial_step` for the first step instead of estimating one.
    pub fn with_initial_step(self, initial_step: f64) -> anyhow::Result<Self> {
        if !(initial_step > 0.0 && initial_step.is_finite()) {
            return Err(PhysicsError::InvalidSetting(
                "initial_step must be positive and finite",
            ))
            .with_context(|| format!("initial_step = {initial_step}"));
        }
        Ok(Self {
            initial_step: Some(initial_step),
            ..self
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_adaptive_options_validates_settings() -> anyhow::Result<()> {
        for (rtol, atol) in [(0.0, 1e-9), (-1e-6, 1e-9), (f64::NAN, 0.0), (1e-6, -1.0)] {
            anyhow::ensure!(
                AdaptiveOptions::new(rtol, atol).is_err(),
                "accepted rtol {rtol}, atol {atol}"
            );
        }
        let ok = AdaptiveOptions::new(1e-6, 0.0)?;
        anyhow::ensure!(ok.with_max_step(0.0).is_err(), "accepted max_step 0");
        anyhow::ensure!(ok.with_initial_step(-1.0).is_err(), "accepted negative h0");
        anyhow::ensure!(
            ok.with_max_step(1.0)?.with_initial_step(0.1).is_ok(),
            "rejected valid"
        );
        Ok(())
    }
}
