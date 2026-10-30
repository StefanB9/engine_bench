//! Recoverable, matchable physics failures.
//!
//! Raised by constructors and solver steps (never by hot-path kernels) and
//! carried inside `anyhow::Error`, so callers can
//! `downcast_ref::<PhysicsError>()` to react to a specific case.

/// A physics or numerics operation could not produce a valid result.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum PhysicsError {
    /// An iterative solver did not reach its tolerance.
    #[error("no convergence after {iterations} iterations (residual {residual:e})")]
    NoConvergence {
        /// Iterations performed.
        iterations: u32,
        /// Last residual magnitude.
        residual: f64,
    },
    /// A root finder was given an interval whose ends have the same sign.
    #[error("root not bracketed: f({a}) = {fa}, f({b}) = {fb}")]
    NotBracketed {
        /// Lower end of the interval.
        a: f64,
        /// Function value at `a`.
        fa: f64,
        /// Upper end of the interval.
        b: f64,
        /// Function value at `b`.
        fb: f64,
    },
    /// The state or its derivative contains NaN or infinity.
    #[error("non-finite state at t = {t}")]
    NonFiniteState {
        /// Independent variable where it was detected.
        t: f64,
    },
    /// An adaptive integrator needed a step smaller than floating point can
    /// resolve.
    #[error("step size {h:e} underflowed at t = {t}")]
    StepSizeUnderflow {
        /// Independent variable where it happened.
        t: f64,
        /// Rejected step size.
        h: f64,
    },
    /// A gas composition is negative, non-finite or does not sum to one.
    #[error("invalid composition: {0}")]
    InvalidComposition(&'static str),
    /// A temperature left the range covered by the thermodynamic data.
    #[error("temperature {temperature} K outside the data range {min}..{max} K")]
    TemperatureOutOfRange {
        /// Offending temperature, K.
        temperature: f64,
        /// Lowest valid temperature, K.
        min: f64,
        /// Highest valid temperature, K.
        max: f64,
    },
    /// A geometric parameter (area, length) is non-physical.
    #[error("invalid geometry: {0}")]
    InvalidGeometry(&'static str),
    /// A physical state or gas property (pressure, temperature, γ) is out of
    /// its valid range.
    #[error("invalid state: {0}")]
    InvalidState(&'static str),
    /// A solver option is out of its valid range.
    #[error("invalid solver setting: {0}")]
    InvalidSetting(&'static str),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_physics_error_display_for_every_variant() -> anyhow::Result<()> {
        let cases = [
            (
                PhysicsError::NoConvergence {
                    iterations: 50,
                    residual: 1.5e-3,
                },
                "no convergence after 50 iterations (residual 1.5e-3)",
            ),
            (
                PhysicsError::NotBracketed {
                    a: 0.0,
                    fa: 1.0,
                    b: 2.0,
                    fb: 3.0,
                },
                "root not bracketed: f(0) = 1, f(2) = 3",
            ),
            (
                PhysicsError::NonFiniteState { t: 0.5 },
                "non-finite state at t = 0.5",
            ),
            (
                PhysicsError::StepSizeUnderflow { t: 1.0, h: 1e-15 },
                "step size 1e-15 underflowed at t = 1",
            ),
            (
                PhysicsError::InvalidComposition("fractions sum to 0.5"),
                "invalid composition: fractions sum to 0.5",
            ),
            (
                PhysicsError::TemperatureOutOfRange {
                    temperature: 150.5,
                    min: 200.0,
                    max: 6000.0,
                },
                "temperature 150.5 K outside the data range 200..6000 K",
            ),
            (
                PhysicsError::InvalidGeometry("exit area below throat area"),
                "invalid geometry: exit area below throat area",
            ),
            (
                PhysicsError::InvalidState("gamma must exceed 1"),
                "invalid state: gamma must exceed 1",
            ),
            (
                PhysicsError::InvalidSetting("rtol must be positive"),
                "invalid solver setting: rtol must be positive",
            ),
        ];
        for (error, expected) in cases {
            let shown = error.to_string();
            anyhow::ensure!(shown == expected, "got {shown:?}, expected {expected:?}");
        }
        Ok(())
    }
}
