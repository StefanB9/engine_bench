//! Real-time pacing: how many fixed physics steps to run for the wall time that
//! passed.
//!
//! Pure logic driven by timestamps from a [`Clock`](crate::clock::Clock); the
//! simulation thread asks for a [`FramePlan`] each tick and runs that many
//! steps.

use std::time::Duration;

use anyhow::Context;

use crate::CoreError;

/// Smoothing time constant for the achieved real-time factor, s.
const FACTOR_SMOOTHING: f64 = 1.0;

/// Fastest allowed time scale (simulated seconds per wall second).
const MAX_TIME_SCALE: f64 = 1000.0;

/// Fixed physics step and catch-up limit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PacerSettings {
    step: f64,
    max_steps_per_frame: u32,
}

impl PacerSettings {
    /// Physics step `step` in simulated seconds (> 0) and at most
    /// `max_steps_per_frame` (≥ 1) steps per tick.
    pub fn new(step: f64, max_steps_per_frame: u32) -> anyhow::Result<Self> {
        if !(step > 0.0 && step.is_finite()) {
            return Err(CoreError::InvalidSetting("physics step must be positive"))
                .with_context(|| format!("step = {step} s"));
        }
        if max_steps_per_frame == 0 {
            return Err(CoreError::InvalidSetting("at least one step per frame").into());
        }
        Ok(Self {
            step,
            max_steps_per_frame,
        })
    }

    /// Physics step, simulated seconds.
    pub fn step(&self) -> f64 {
        self.step
    }
}

/// What to do this tick.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FramePlan {
    /// Physics steps to run now.
    pub steps: u32,
    /// The catch-up limit cut this tick short: the simulation falls behind real
    /// time and the backlog is dropped rather than carried (which would make
    /// every later tick slower).
    pub capped: bool,
}

/// Turns elapsed wall time into fixed physics steps.
#[derive(Debug, Clone)]
pub struct RealTimePacer {
    settings: PacerSettings,
    time_scale: f64,
    paused: bool,
    pending_single_steps: u32,
    last_wall: Option<Duration>,
    /// Simulated time owed but not yet stepped, s.
    owed: f64,
    /// Exponentially smoothed simulated and wall time, for the achieved factor.
    smoothed_simulated: f64,
    smoothed_wall: f64,
}

impl RealTimePacer {
    /// A running pacer at time scale 1.
    pub fn new(settings: PacerSettings) -> Self {
        Self {
            settings,
            time_scale: 1.0,
            paused: false,
            pending_single_steps: 0,
            last_wall: None,
            owed: 0.0,
            smoothed_simulated: 0.0,
            smoothed_wall: 0.0,
        }
    }

    /// Simulated seconds per wall second, in (0, 1000]; below 1 is slow motion.
    pub fn set_time_scale(&mut self, scale: f64) -> anyhow::Result<()> {
        if !(scale > 0.0 && scale <= MAX_TIME_SCALE) {
            return Err(CoreError::InvalidSetting(
                "time scale must lie in (0, 1000]",
            ))
            .with_context(|| format!("scale = {scale}"));
        }
        self.time_scale = scale;
        Ok(())
    }

    /// Current time scale.
    pub fn time_scale(&self) -> f64 {
        self.time_scale
    }

    /// Stops advancing (single steps still run).
    pub fn pause(&mut self) {
        self.paused = true;
    }

    /// Continues advancing from now, without catching up on the paused time.
    pub fn resume(&mut self) {
        self.paused = false;
        self.owed = 0.0;
    }

    /// Whether paused.
    pub fn is_paused(&self) -> bool {
        self.paused
    }

    /// Runs exactly one physics step at the next tick (meant for use while
    /// paused).
    pub fn request_step(&mut self) {
        self.pending_single_steps = self.pending_single_steps.saturating_add(1);
    }

    /// Steps to run for a tick at wall time `now`.
    pub fn plan(&mut self, now: Duration) -> FramePlan {
        let Some(last) = self.last_wall.replace(now) else {
            return FramePlan {
                steps: 0,
                capped: false,
            };
        };
        let elapsed = now.saturating_sub(last).as_secs_f64();
        let step = self.settings.step;
        let single_steps = std::mem::take(&mut self.pending_single_steps);
        let plan = if self.paused {
            self.owed = 0.0;
            FramePlan {
                steps: single_steps,
                capped: false,
            }
        } else {
            self.owed += elapsed * self.time_scale;
            // The epsilon keeps exact multiples of the step from rounding down.
            let wanted = (self.owed / step + 1e-9).floor();
            let max = self.settings.max_steps_per_frame;
            if wanted > f64::from(max) {
                self.owed = 0.0;
                FramePlan {
                    steps: max,
                    capped: true,
                }
            } else {
                self.owed -= wanted * step;
                FramePlan {
                    steps: whole_steps(wanted),
                    capped: false,
                }
            }
        };
        let decay = (-elapsed / FACTOR_SMOOTHING).exp();
        self.smoothed_simulated = self.smoothed_simulated * decay + f64::from(plan.steps) * step;
        self.smoothed_wall = self.smoothed_wall * decay + elapsed;
        plan
    }

    /// Simulated seconds per wall second actually achieved recently (smoothed
    /// over about a second); `None` before any wall time has passed.
    pub fn achieved_time_scale(&self) -> Option<f64> {
        (self.smoothed_wall > 0.0).then(|| self.smoothed_simulated / self.smoothed_wall)
    }
}

/// Converts a step count already checked to lie in `0..=u32::MAX` and be whole.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "callers pass whole, non-negative counts no larger than a u32 frame cap"
)]
fn whole_steps(count: f64) -> u32 {
    count as u32
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    const STEP: f64 = 1e-3;

    fn pacer(max_steps: u32) -> anyhow::Result<RealTimePacer> {
        Ok(RealTimePacer::new(PacerSettings::new(STEP, max_steps)?))
    }

    fn ms(value: f64) -> Duration {
        Duration::from_secs_f64(value * 1e-3)
    }

    #[test]
    fn test_pacer_settings_validate() -> anyhow::Result<()> {
        for (step, max) in [(0.0, 10), (-1e-3, 10), (f64::NAN, 10), (1e-3, 0)] {
            anyhow::ensure!(
                PacerSettings::new(step, max).is_err(),
                "accepted {step}, {max}"
            );
        }
        anyhow::ensure!(
            (PacerSettings::new(2e-3, 5)?.step() - 2e-3).abs() < 1e-18,
            "step"
        );
        Ok(())
    }

    #[test]
    fn test_pacer_first_tick_only_starts_the_clock() -> anyhow::Result<()> {
        let mut p = pacer(1000)?;
        anyhow::ensure!(p.plan(ms(500.0)).steps == 0, "stepped on first tick");
        anyhow::ensure!(p.plan(ms(516.0)).steps == 16, "second tick");
        Ok(())
    }

    #[test]
    fn test_pacer_carries_fractional_steps() -> anyhow::Result<()> {
        let mut p = pacer(1000)?;
        p.plan(ms(0.0));
        let steps: Vec<u32> = [16.5, 33.0, 49.5, 66.0]
            .iter()
            .map(|t| p.plan(ms(*t)).steps)
            .collect();
        anyhow::ensure!(steps == [16, 17, 16, 17], "{steps:?}");
        Ok(())
    }

    #[test]
    fn test_pacer_time_scale_slows_and_speeds() -> anyhow::Result<()> {
        let mut p = pacer(1000)?;
        p.set_time_scale(0.25)?;
        p.plan(ms(0.0));
        anyhow::ensure!(p.plan(ms(40.0)).steps == 10, "slow motion");
        p.set_time_scale(4.0)?;
        anyhow::ensure!(p.plan(ms(50.0)).steps == 40, "fast forward");
        for bad in [0.0, -1.0, 1001.0, f64::NAN] {
            anyhow::ensure!(p.set_time_scale(bad).is_err(), "accepted scale {bad}");
        }
        anyhow::ensure!(
            (p.time_scale() - 4.0).abs() < 1e-15,
            "scale changed by a rejected value"
        );
        Ok(())
    }

    #[test]
    fn test_pacer_pause_step_and_resume() -> anyhow::Result<()> {
        let mut p = pacer(1000)?;
        p.plan(ms(0.0));
        p.pause();
        anyhow::ensure!(p.is_paused() && p.plan(ms(100.0)).steps == 0, "paused");
        p.request_step();
        p.request_step();
        anyhow::ensure!(p.plan(ms(110.0)).steps == 2, "single steps");
        anyhow::ensure!(p.plan(ms(120.0)).steps == 0, "single steps repeated");
        p.resume();
        // Paused wall time is not owed after resuming.
        anyhow::ensure!(
            !p.is_paused() && p.plan(ms(130.0)).steps == 10,
            "resume caught up"
        );
        Ok(())
    }

    #[test]
    fn test_pacer_caps_and_drops_backlog() -> anyhow::Result<()> {
        let mut p = pacer(50)?;
        p.plan(ms(0.0));
        let plan = p.plan(ms(500.0));
        anyhow::ensure!(
            plan == FramePlan {
                steps: 50,
                capped: true
            },
            "{plan:?}"
        );
        // The 450 ms that did not fit are dropped, not carried into the next
        // tick.
        let next = p.plan(ms(510.0));
        anyhow::ensure!(
            next == FramePlan {
                steps: 10,
                capped: false
            },
            "{next:?}"
        );
        Ok(())
    }

    #[test]
    fn test_pacer_reports_achieved_time_scale() -> anyhow::Result<()> {
        let mut keeping_up = pacer(1000)?;
        anyhow::ensure!(
            keeping_up.achieved_time_scale().is_none(),
            "factor before any time"
        );
        let mut behind = pacer(8)?;
        for frame in 0..200_u32 {
            let now = ms(f64::from(frame) * 16.0);
            keeping_up.plan(now);
            behind.plan(now);
        }
        let ok = keeping_up.achieved_time_scale().unwrap_or(0.0);
        anyhow::ensure!((ok - 1.0).abs() < 0.05, "keeping up: {ok}");
        let slow = behind.achieved_time_scale().unwrap_or(0.0);
        anyhow::ensure!((slow - 0.5).abs() < 0.05, "8 of 16 steps per frame: {slow}");
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_pacer_accounts_for_every_wall_second(
            frames in proptest::collection::vec(0.0..40.0f64, 1..200),
            scale in 0.1..10.0f64,
        ) {
            let mut p = pacer(u32::MAX).map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            p.set_time_scale(scale).map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let mut now = 0.0;
            p.plan(ms(now));
            let mut steps = 0.0;
            for frame in frames {
                now += frame;
                let plan = p.plan(ms(now));
                prop_assert!(!plan.capped);
                steps += f64::from(plan.steps);
            }
            // Owed time stays below one step; nothing is lost or invented.
            let owed = scale * ms(now).as_secs_f64() - steps * STEP;
            prop_assert!((-1e-9..STEP + 1e-9).contains(&owed), "owed {owed}");
        }

        #[test]
        fn test_pacer_never_exceeds_cap(
            frames in proptest::collection::vec(0.0..1000.0f64, 1..100),
            cap in 1..64u32,
        ) {
            let mut p = pacer(cap).map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let mut now = 0.0;
            for frame in frames {
                now += frame;
                let plan = p.plan(ms(now));
                prop_assert!(plan.steps <= cap, "{} steps > cap {cap}", plan.steps);
            }
        }
    }
}
