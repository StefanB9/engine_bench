//! ODE integrators on fixed-size state `[f64; N]` (no heap, fully unrollable).
//!
//! - [`rk4_step`]: classic fixed-step RK4 for the real-time loop's fixed
//!   physics timestep.
//! - [`DormandPrince45`]: adaptive 5(4) pair with PI step control, dense output
//!   and event location, for offline runs and for problems with discrete
//!   transitions (a tank running dry, a rocket reaching apogee).
//!
//! Layout: [`OdeSystem`] and shared stage arithmetic here; [`rk4_step`] in
//! `rk4`; the adaptive integrator in `dormand_prince` with its options and
//! event types in `options` / `events`.

mod dormand_prince;
mod events;
mod options;
mod rk4;
#[cfg(test)]
mod test_systems;

pub use dormand_prince::{DormandPrince45, Statistics, StepOutcome};
pub use events::{Direction, EventSet, NoEvents};
pub use options::AdaptiveOptions;
pub use rk4::rk4_step;

/// A first-order system `dy/dt = f(t, y)`.
///
/// `derivatives` is a hot-path kernel: it must not fail and must not allocate.
/// Validation of the state happens in the integrator's `step`, which rejects
/// non-finite results.
pub trait OdeSystem<const N: usize> {
    /// Writes `f(t, y)` into `dydt`.
    fn derivatives(&mut self, t: f64, y: &[f64; N], dydt: &mut [f64; N]);
}

/// `out = base + Σ cᵢ·kᵢ`.
#[inline]
fn combine<const N: usize>(out: &mut [f64; N], base: &[f64; N], terms: &[(f64, &[f64; N])]) {
    *out = *base;
    for &(c, k) in terms {
        for (o, ki) in out.iter_mut().zip(k) {
            *o += c * ki;
        }
    }
}

fn all_finite<const N: usize>(y: &[f64; N]) -> bool {
    y.iter().all(|v| v.is_finite())
}
