//! How accurately an engine should advance.

use anyhow::Context;

use crate::CoreError;

/// Accuracy target for one advance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Accuracy {
    /// Interactive: fixed-step integration sized for speed and stability.
    RealTime,
    /// Batch: adaptive integration to relative tolerance `rtol`.
    Offline {
        /// Relative tolerance for adaptive integrators.
        rtol: f64,
    },
}

/// Context passed to every [`advance`](crate::engine::Engine::advance).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepContext {
    accuracy: Accuracy,
}

impl StepContext {
    /// Interactive stepping.
    pub fn real_time() -> Self {
        Self {
            accuracy: Accuracy::RealTime,
        }
    }

    /// Batch stepping to relative tolerance `rtol` in (0, 1).
    pub fn offline(rtol: f64) -> anyhow::Result<Self> {
        if !(rtol > 0.0 && rtol < 1.0) {
            return Err(CoreError::InvalidSetting("rtol must lie in (0, 1)"))
                .with_context(|| format!("rtol = {rtol}"));
        }
        Ok(Self {
            accuracy: Accuracy::Offline { rtol },
        })
    }

    /// The accuracy target.
    pub fn accuracy(&self) -> Accuracy {
        self.accuracy
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_step_context_modes() -> anyhow::Result<()> {
        anyhow::ensure!(
            StepContext::real_time().accuracy() == Accuracy::RealTime,
            "real time"
        );
        let offline = StepContext::offline(1e-8)?;
        anyhow::ensure!(
            offline.accuracy() == Accuracy::Offline { rtol: 1e-8 },
            "offline"
        );
        for bad in [0.0, -1e-6, 1.0, f64::NAN] {
            anyhow::ensure!(StepContext::offline(bad).is_err(), "accepted rtol {bad}");
        }
        Ok(())
    }
}
