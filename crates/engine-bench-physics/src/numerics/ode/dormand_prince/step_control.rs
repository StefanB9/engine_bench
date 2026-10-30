//! Step-size selection: PI controller constants, error scale and initial-step
//! estimate.

use super::DormandPrince45;
use crate::numerics::ode::{OdeSystem, combine};

// PI step-size controller (Hairer, Nørsett & Wanner, Solving ODEs I, §II.4).
pub(super) const SAFETY: f64 = 0.9;
pub(super) const BETA: f64 = 0.04;
pub(super) const EXPO1: f64 = 0.2 - BETA * 0.75;
pub(super) const MIN_SHRINK: f64 = 0.2;
pub(super) const MAX_GROWTH: f64 = 10.0;

impl<const N: usize> DormandPrince45<N> {
    #[inline]
    pub(super) fn evaluate<S: OdeSystem<N>>(
        &mut self,
        system: &mut S,
        t: f64,
        y: &[f64; N],
        out: &mut [f64; N],
    ) {
        self.statistics.evaluations += 1;
        system.derivatives(t, y, out);
    }

    pub(super) fn scale(&self, y0: f64, y1: f64) -> f64 {
        self.options.atol + self.options.rtol * y0.abs().max(y1.abs())
    }

    /// Hairer's `hinit`: makes the explicit Euler step's local error roughly
    /// match the tolerance, using one extra derivative evaluation.
    pub(super) fn initial_step<S: OdeSystem<N>>(
        &mut self,
        system: &mut S,
        t: f64,
        y: &[f64; N],
        k1: &[f64; N],
    ) -> f64 {
        let max_step = self.options.max_step;
        if let Some(h) = self.options.initial_step {
            return h.min(max_step);
        }
        // A zero scale (atol = 0 on a component that is exactly 0) says nothing
        // about magnitudes, so such components are left out of the estimate.
        let (mut norm_f, mut norm_y) = (0.0, 0.0);
        for i in 0..N {
            let sk = self.scale(y[i], y[i]);
            if sk > 0.0 {
                norm_f += (k1[i] / sk).powi(2);
                norm_y += (y[i] / sk).powi(2);
            }
        }
        let mut h = if norm_f <= 1e-10 || norm_y <= 1e-10 {
            1e-6
        } else {
            (norm_y / norm_f).sqrt() * 0.01
        };
        h = h.min(max_step);

        let mut euler = [0.0; N];
        combine(&mut euler, y, &[(h, k1)]);
        let mut f1 = [0.0; N];
        self.evaluate(system, t + h, &euler, &mut f1);
        let mut second = 0.0;
        for i in 0..N {
            let sk = self.scale(y[i], y[i]);
            if sk > 0.0 {
                second += ((f1[i] - k1[i]) / sk).powi(2);
            }
        }
        let second = second.sqrt() / h;
        let curvature = second.max(norm_f.sqrt());
        let h1 = if curvature <= 1e-15 {
            (h * 1e-3).max(1e-6)
        } else {
            (0.01 / curvature).powf(0.2)
        };
        (100.0 * h).min(h1).min(max_step)
    }
}
