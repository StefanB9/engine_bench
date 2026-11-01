//! Adaptive Dormand–Prince 5(4) integrator.

mod event_location;
mod step_control;
mod tableau;

use anyhow::Context;

use self::{
    step_control::{BETA, EXPO1, MAX_GROWTH, MIN_SHRINK, SAFETY},
    tableau::Dense,
};
use super::{AdaptiveOptions, EventSet, OdeSystem, all_finite};
use crate::PhysicsError;

/// Result of one [`DormandPrince45::step`] call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepOutcome {
    /// One step was accepted; `t` is still before `t_end`.
    Stepped,
    /// `t` reached `t_end` exactly.
    ReachedEnd,
    /// Event `index` fired; `t` and `y` are at the located crossing, on its far
    /// side.
    Event {
        /// Index into the event set.
        index: usize,
    },
}

/// Work counters for an integrator instance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Statistics {
    /// Calls to [`OdeSystem::derivatives`].
    pub evaluations: u64,
    /// Accepted steps.
    pub accepted: u64,
    /// Rejected trial steps.
    pub rejected: u64,
}

/// Adaptive Dormand–Prince 5(4) integrator (Hairer's DOPRI5 step control and
/// dense output).
///
/// The last stage of an accepted step is the first stage of the next (FSAL), so
/// a step costs six evaluations. That cache assumes `y` and the system only
/// change through this integrator: call [`reset`](Self::reset) after any
/// outside change.
#[derive(Debug, Clone)]
pub struct DormandPrince45<const N: usize> {
    options: AdaptiveOptions,
    statistics: Statistics,
    /// `f(t, y)` at the current point, if still valid.
    first_stage: Option<[f64; N]>,
    /// Proposed next step; `None` until estimated.
    next_step: Option<f64>,
    previous_error: f64,
    dense: Option<Dense<N>>,
}

impl<const N: usize> DormandPrince45<N> {
    /// Creates an integrator; the first step size is estimated on the first
    /// call.
    pub fn new(options: AdaptiveOptions) -> Self {
        Self {
            options,
            statistics: Statistics::default(),
            first_stage: None,
            next_step: None,
            previous_error: 1e-4,
            dense: None,
        }
    }

    /// Forgets the cached derivative and step history. Call it whenever the
    /// state or the system changes outside the integrator (new control input,
    /// discrete jump in `y`).
    pub fn reset(&mut self) {
        self.first_stage = None;
        self.next_step = None;
        self.previous_error = 1e-4;
    }

    /// Forgets only the cached derivative, keeping the step size and its
    /// history: for an outside change that rearranges the state without
    /// changing the solution it represents.
    pub fn forget_derivative(&mut self) {
        self.first_stage = None;
    }

    /// Work counters since construction.
    pub fn statistics(&self) -> Statistics {
        self.statistics
    }

    /// Takes one accepted step towards `t_end` (never past it), retrying
    /// internally with smaller steps as needed, and stops early at the first
    /// event crossing.
    pub fn step<S, E, const M: usize>(
        &mut self,
        system: &mut S,
        events: &mut E,
        t: &mut f64,
        y: &mut [f64; N],
        t_end: f64,
    ) -> anyhow::Result<StepOutcome>
    where
        S: OdeSystem<N>,
        E: EventSet<N, M>,
    {
        let t0 = *t;
        if t_end.is_nan() || t_end < t0 {
            return Err(PhysicsError::InvalidSetting("t_end must not be before t"))
                .with_context(|| format!("t = {t0}, t_end = {t_end}"));
        }
        if t_end <= t0 {
            return Ok(StepOutcome::ReachedEnd);
        }

        let k1 = if let Some(k1) = self.first_stage {
            k1
        } else {
            let mut k1 = [0.0; N];
            self.evaluate(system, t0, y, &mut k1);
            if !all_finite(&k1) {
                return Err(PhysicsError::NonFiniteState { t: t0 })
                    .context("derivative at step start");
            }
            k1
        };
        let proposed = match self.next_step {
            Some(h) => h,
            None => self.initial_step(system, t0, y, &k1),
        };

        let mut h = proposed.min(self.options.max_step);
        let mut rejected_before = false;
        loop {
            let remaining = t_end - t0;
            let last = h * 1.01 >= remaining;
            if last {
                h = remaining;
            }
            if 0.1 * h <= t0.abs() * f64::EPSILON || h <= 0.0 {
                return Err(PhysicsError::StepSizeUnderflow { t: t0, h }).context("dormand-prince");
            }

            let trial = self.trial_step(system, t0, y, &k1, h);
            let error = trial.error;
            // NaN fails both comparisons and is treated as a large error.
            let accept = error <= 1.0;
            let shrink = error.powf(EXPO1);
            if !accept {
                self.statistics.rejected += 1;
                rejected_before = true;
                let factor = if shrink.is_finite() {
                    (shrink / SAFETY).min(1.0 / MIN_SHRINK)
                } else {
                    1.0 / MIN_SHRINK
                };
                h /= factor;
                continue;
            }

            self.statistics.accepted += 1;
            let factor = (shrink / self.previous_error.powf(BETA) / SAFETY)
                .clamp(1.0 / MAX_GROWTH, 1.0 / MIN_SHRINK);
            self.previous_error = error.max(1e-4);
            let mut h_new = (h / factor).min(self.options.max_step);
            if rejected_before {
                h_new = h_new.min(h);
            }
            // A step clamped to hit t_end says little about the natural step
            // size.
            self.next_step = Some(if last { h_new.max(proposed) } else { h_new });
            self.dense = Some(trial.dense);

            let t1 = if last { t_end } else { t0 + h };
            if let Some((index, t_event)) = self.locate_event(events, t0, y, t1, &trial.y) {
                self.dense_at(t_event, y);
                *t = t_event;
                self.first_stage = None;
                return Ok(StepOutcome::Event { index });
            }
            *t = t1;
            *y = trial.y;
            self.first_stage = Some(trial.k7);
            return Ok(if last {
                StepOutcome::ReachedEnd
            } else {
                StepOutcome::Stepped
            });
        }
    }

    /// Steps until `t_end` or the first event.
    pub fn integrate<S, E, const M: usize>(
        &mut self,
        system: &mut S,
        events: &mut E,
        t: &mut f64,
        y: &mut [f64; N],
        t_end: f64,
    ) -> anyhow::Result<StepOutcome>
    where
        S: OdeSystem<N>,
        E: EventSet<N, M>,
    {
        loop {
            match self.step(system, events, t, y, t_end)? {
                StepOutcome::Stepped => {}
                done => return Ok(done),
            }
        }
    }

    /// Evaluates the dense-output polynomial of the last accepted step at `t`,
    /// which must lie within that step.
    pub fn interpolate(&self, t: f64, out: &mut [f64; N]) -> anyhow::Result<()> {
        let Some(dense) = &self.dense else {
            return Err(PhysicsError::InvalidSetting("no step taken yet")).context("interpolate");
        };
        let slack = 4.0 * f64::EPSILON * (dense.t0.abs() + dense.h.abs());
        if !(t >= dense.t0 - slack && t <= dense.t0 + dense.h + slack) {
            return Err(PhysicsError::InvalidSetting("t outside the last step"))
                .with_context(|| format!("t = {t}, step [{}, {}]", dense.t0, dense.t0 + dense.h));
        }
        dense.eval(t, out);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::TAU;

    use proptest::prelude::*;

    use super::*;
    use crate::numerics::ode::{
        Direction, NoEvents,
        test_systems::{
            AlwaysNan, ApogeeAndGround, Ball, BlowsUp, Decay, G, Ground, Oscillator, physics_error,
        },
    };

    #[test]
    fn test_dormand_prince_hits_end_exactly_with_requested_accuracy() -> anyhow::Result<()> {
        let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-10, 1e-14)?);
        let (mut t, mut y) = (0.0, [1.0]);
        let outcome = dp.integrate(&mut Decay { rate: 1.0 }, &mut NoEvents, &mut t, &mut y, 5.0)?;
        anyhow::ensure!(outcome == StepOutcome::ReachedEnd, "{outcome:?}");
        anyhow::ensure!(t.to_bits() == 5.0f64.to_bits(), "t = {t}");
        let rel = (y[0] / (-5.0f64).exp() - 1.0).abs();
        anyhow::ensure!(rel < 1e-8, "relative error {rel}");
        Ok(())
    }

    #[test]
    fn test_dormand_prince_error_shrinks_with_tolerance() -> anyhow::Result<()> {
        let mut errors = [0.0; 2];
        for (error, rtol) in errors.iter_mut().zip([1e-5, 1e-9]) {
            let mut dp = DormandPrince45::new(AdaptiveOptions::new(rtol, 0.0)?);
            let (mut t, mut y) = (0.0, [1.0, 0.0]);
            dp.integrate(&mut Oscillator, &mut NoEvents, &mut t, &mut y, 10.0 * TAU)?;
            *error = (y[0] - 1.0).abs().max(y[1].abs());
        }
        anyhow::ensure!(errors[1] < 1e-7, "tight error {}", errors[1]);
        anyhow::ensure!(errors[0] / errors[1] > 1e3, "errors {errors:?}");
        Ok(())
    }

    #[test]
    fn test_dormand_prince_is_fifth_order_at_fixed_step() -> anyhow::Result<()> {
        // A huge tolerance accepts every step, so max_step fixes the step size.
        let error_at = |h: f64| -> anyhow::Result<f64> {
            let options = AdaptiveOptions::new(1e6, 0.0)?
                .with_max_step(h)?
                .with_initial_step(h)?;
            let mut dp = DormandPrince45::new(options);
            let (mut t, mut y) = (0.0, [1.0]);
            dp.integrate(&mut Decay { rate: 1.0 }, &mut NoEvents, &mut t, &mut y, 2.0)?;
            Ok((y[0] - (-2.0f64).exp()).abs())
        };
        let ratio = error_at(0.2)? / error_at(0.1)?;
        anyhow::ensure!(
            (26.0..38.0).contains(&ratio),
            "error ratio {ratio}, expected ≈ 32"
        );
        Ok(())
    }

    #[test]
    fn test_dormand_prince_reuses_last_stage_between_steps() -> anyhow::Result<()> {
        let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-8, 1e-12)?);
        let (mut t, mut y) = (0.0, [1.0, 0.0]);
        dp.integrate(&mut Oscillator, &mut NoEvents, &mut t, &mut y, TAU)?;
        let s = dp.statistics();
        anyhow::ensure!(s.accepted > 10, "{s:?}");
        // FSAL: 6 new evaluations per trial step, plus initial derivative and
        // step estimate.
        let expected = 6 * (s.accepted + s.rejected) + 2;
        anyhow::ensure!(
            s.evaluations == expected,
            "{s:?}, expected {expected} evaluations"
        );
        Ok(())
    }

    #[test]
    fn test_dormand_prince_reset_discards_cached_derivative() -> anyhow::Result<()> {
        let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-8, 1e-12)?);
        let (mut t, mut y) = (0.0, [1.0]);
        let mut slow = Decay { rate: 1.0 };
        dp.integrate(&mut slow, &mut NoEvents, &mut t, &mut y, 1.0)?;
        // The system changes externally: rate doubles for t > 1.
        dp.reset();
        let mut fast = Decay { rate: 2.0 };
        dp.integrate(&mut fast, &mut NoEvents, &mut t, &mut y, 2.0)?;
        let exact = (-1.0f64 - 2.0).exp();
        anyhow::ensure!(
            (y[0] / exact - 1.0).abs() < 1e-6,
            "y = {}, exact {exact}",
            y[0]
        );
        Ok(())
    }

    #[test]
    fn test_dormand_prince_forget_derivative_keeps_the_step_size() -> anyhow::Result<()> {
        // One step, then a rearrangement that keeps the solution: the next
        // step re-evaluates the derivative (one extra evaluation) but needs no
        // new step estimate (which reset would cost: two).
        let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-8, 1e-12)?);
        let (mut t, mut y) = (0.0, [1.0]);
        let mut decay = Decay { rate: 1.0 };
        dp.step(&mut decay, &mut NoEvents, &mut t, &mut y, 2.0)?;
        let before = dp.statistics();
        dp.forget_derivative();
        dp.step(&mut decay, &mut NoEvents, &mut t, &mut y, 2.0)?;
        let after = dp.statistics();
        let trials = after.accepted + after.rejected - before.accepted - before.rejected;
        anyhow::ensure!(
            after.evaluations - before.evaluations == 6 * trials + 1,
            "{before:?} then {after:?}"
        );
        anyhow::ensure!((y[0] / (-t).exp() - 1.0).abs() < 1e-7, "y = {}", y[0]);
        Ok(())
    }

    #[test]
    fn test_dormand_prince_locates_ground_contact() -> anyhow::Result<()> {
        let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-10, 1e-12)?);
        let (mut t, mut y) = (0.0, [10.0, 0.0]);
        let outcome = dp.integrate(
            &mut Ball,
            &mut Ground(Direction::Falling),
            &mut t,
            &mut y,
            100.0,
        )?;
        anyhow::ensure!(outcome == StepOutcome::Event { index: 0 }, "{outcome:?}");
        let exact = (2.0 * 10.0 / G).sqrt();
        anyhow::ensure!((t - exact).abs() < 1e-9, "t = {t}, exact {exact}");
        anyhow::ensure!(
            y[0] <= 0.0 && y[0] > -1e-8,
            "height {} not just past zero",
            y[0]
        );
        Ok(())
    }

    #[test]
    fn test_dormand_prince_does_not_refire_event_after_resuming() -> anyhow::Result<()> {
        let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-10, 1e-12)?);
        let (mut t, mut y) = (0.0, [10.0, 0.0]);
        let mut ground = Ground(Direction::Falling);
        dp.integrate(&mut Ball, &mut ground, &mut t, &mut y, 100.0)?;
        let outcome = dp.integrate(&mut Ball, &mut ground, &mut t, &mut y, 5.0)?;
        anyhow::ensure!(outcome == StepOutcome::ReachedEnd, "refired: {outcome:?}");
        Ok(())
    }

    #[test]
    fn test_dormand_prince_ignores_crossing_in_other_direction() -> anyhow::Result<()> {
        let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-10, 1e-12)?);
        let (mut t, mut y) = (0.0, [10.0, 0.0]);
        let outcome = dp.integrate(
            &mut Ball,
            &mut Ground(Direction::Rising),
            &mut t,
            &mut y,
            3.0,
        )?;
        anyhow::ensure!(outcome == StepOutcome::ReachedEnd, "{outcome:?}");
        Ok(())
    }

    #[test]
    fn test_dormand_prince_reports_earliest_of_several_events() -> anyhow::Result<()> {
        let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-10, 1e-12)?);
        // Thrown upward at 20 m/s from 1 m: apogee at 20/g, ground later.
        let (mut t, mut y) = (0.0, [1.0, 20.0]);
        let mut events = ApogeeAndGround;
        let first = dp.integrate(&mut Ball, &mut events, &mut t, &mut y, 100.0)?;
        anyhow::ensure!(first == StepOutcome::Event { index: 0 }, "{first:?}");
        anyhow::ensure!((t - 20.0 / G).abs() < 1e-9, "apogee at t = {t}");
        let second = dp.integrate(&mut Ball, &mut events, &mut t, &mut y, 100.0)?;
        anyhow::ensure!(second == StepOutcome::Event { index: 1 }, "{second:?}");
        Ok(())
    }

    #[test]
    fn test_dormand_prince_interpolates_within_last_step() -> anyhow::Result<()> {
        let options = AdaptiveOptions::new(1e-9, 1e-12)?.with_max_step(0.5)?;
        let mut dp = DormandPrince45::new(options);
        let (mut t0, mut y) = (0.0, [1.0, 0.0]);
        dp.step(&mut Oscillator, &mut NoEvents, &mut t0, &mut y, 10.0)?;
        let mut t1 = t0;
        dp.step(&mut Oscillator, &mut NoEvents, &mut t1, &mut y, 10.0)?;
        let mid = f64::midpoint(t0, t1);
        let mut out = [0.0; 2];
        dp.interpolate(mid, &mut out)?;
        anyhow::ensure!((out[0] - mid.cos()).abs() < 1e-8, "x({mid}) = {}", out[0]);
        anyhow::ensure!((out[1] + mid.sin()).abs() < 1e-8, "v({mid}) = {}", out[1]);
        anyhow::ensure!(
            dp.interpolate(t1 + 0.1, &mut out).is_err(),
            "extrapolated past step"
        );
        anyhow::ensure!(
            dp.interpolate(t0 - 0.1, &mut out).is_err(),
            "extrapolated before"
        );
        Ok(())
    }

    #[test]
    fn test_dormand_prince_interpolate_requires_a_step() -> anyhow::Result<()> {
        let dp = DormandPrince45::<1>::new(AdaptiveOptions::new(1e-6, 0.0)?);
        anyhow::ensure!(
            dp.interpolate(0.0, &mut [0.0]).is_err(),
            "no step taken yet"
        );
        Ok(())
    }

    #[test]
    fn test_dormand_prince_rejects_backward_end() -> anyhow::Result<()> {
        let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-6, 0.0)?);
        let (mut t, mut y) = (1.0, [1.0]);
        let result = dp.step(&mut Decay { rate: 1.0 }, &mut NoEvents, &mut t, &mut y, 0.5);
        anyhow::ensure!(result.is_err(), "integrated backwards: {result:?}");
        let same = dp.step(&mut Decay { rate: 1.0 }, &mut NoEvents, &mut t, &mut y, 1.0)?;
        anyhow::ensure!(same == StepOutcome::ReachedEnd, "{same:?}");
        Ok(())
    }

    #[test]
    fn test_dormand_prince_rejects_non_finite_start() -> anyhow::Result<()> {
        let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-6, 0.0)?);
        let (mut t, mut y) = (0.0, [1.0]);
        let Err(e) = dp.step(&mut AlwaysNan, &mut NoEvents, &mut t, &mut y, 1.0) else {
            anyhow::bail!("accepted NaN derivative");
        };
        anyhow::ensure!(
            matches!(physics_error(&e), Some(PhysicsError::NonFiniteState { .. })),
            "wrong error: {e:#}"
        );
        Ok(())
    }

    #[test]
    fn test_dormand_prince_fails_cleanly_at_singularity() -> anyhow::Result<()> {
        let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-8, 1e-12)?);
        let (mut t, mut y) = (0.0, [0.0]);
        let Err(e) = dp.integrate(&mut BlowsUp, &mut NoEvents, &mut t, &mut y, 2.0) else {
            anyhow::bail!("integrated through a singularity");
        };
        anyhow::ensure!(
            matches!(
                physics_error(&e),
                Some(PhysicsError::StepSizeUnderflow { .. } | PhysicsError::NonFiniteState { .. })
            ),
            "wrong error: {e:#}"
        );
        anyhow::ensure!(
            t < 1.0 && y[0].is_finite(),
            "state corrupted: t = {t}, y = {y:?}"
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_dormand_prince_decay_matches_exponential(
            rate in 0.1..10.0f64, y0 in -100.0..100.0f64, t_end in 0.1..3.0f64,
        ) {
            let options = AdaptiveOptions::new(1e-9, 1e-12)
                .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let mut dp = DormandPrince45::new(options);
            let (mut t, mut y) = (0.0, [y0]);
            dp.integrate(&mut Decay { rate }, &mut NoEvents, &mut t, &mut y, t_end)
                .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let exact = y0 * (-rate * t_end).exp();
            prop_assert!((y[0] - exact).abs() <= 1e-7 * exact.abs() + 1e-10,
                "y = {}, exact {exact}", y[0]);
        }

    }
}
