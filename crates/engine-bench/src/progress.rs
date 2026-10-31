//! Progress bars for offline runs and sweeps on the terminal.
//!
//! A row reports as one atomic store; simulated time, speed and the time left
//! are computed only when the bar redraws (at most 20 times a second), and a
//! steady tick keeps the elapsed time moving through slow rows. On anything
//! but a terminal (pipes, logs, tests) the bars draw nothing.

use std::{
    fmt::Write as _,
    io::Write as _,
    sync::atomic::{AtomicUsize, Ordering},
    time::Duration,
};

use anyhow::Context;
use engine_bench_core::offline::{Progress, RunSettings};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressState, ProgressStyle};

/// How often the elapsed time redraws without new rows.
const TICK: Duration = Duration::from_millis(100);
const BAR_CHARACTERS: &str = "█▉▊▋▌▍▎▏ ";

/// The bar of one run: simulated time out of the total, speed against real
/// time, time left.
pub(crate) struct RunBar {
    bar: ProgressBar,
}

impl RunBar {
    /// A bar for a run with `settings`, drawn to `target`.
    pub(crate) fn new(settings: RunSettings, target: ProgressDrawTarget) -> anyhow::Result<Self> {
        let rows = settings.rows();
        let total = settings.row_time(rows.saturating_sub(1));
        let style = ProgressStyle::with_template(
            "{spinner:.cyan} {elapsed_precise} {bar:30.cyan/blue} {percent:>3}% · t = {simulated} \
             of {total}{speed} · {left}",
        )
        .context("progress template")?
        .progress_chars(BAR_CHARACTERS)
        .with_key(
            "simulated",
            move |state: &ProgressState, w: &mut dyn std::fmt::Write| {
                let _ = w.write_str(&seconds(simulated(state.pos(), settings)));
            },
        )
        .with_key(
            "total",
            move |_: &ProgressState, w: &mut dyn std::fmt::Write| {
                let _ = w.write_str(&seconds(total));
            },
        )
        .with_key(
            "speed",
            move |state: &ProgressState, w: &mut dyn std::fmt::Write| {
                let _ = w.write_str(&speed(
                    simulated(state.pos(), settings),
                    state.elapsed().as_secs_f64(),
                ));
            },
        )
        .with_key("left", time_left_key);
        let bar = ProgressBar::with_draw_target(Some(u64::from(rows)), target).with_style(style);
        bar.enable_steady_tick(TICK);
        Ok(Self { bar })
    }

    /// Leaves the full bar on screen and ends its line, so output that
    /// follows starts on a line of its own.
    pub(crate) fn finish(self) {
        end_line(&self.bar);
    }

    /// Leaves the bar where the run failed.
    pub(crate) fn abandon(&self) {
        self.bar.abandon();
    }

    #[cfg(test)]
    fn position(&self) -> u64 {
        self.bar.position()
    }
}

impl Progress for RunBar {
    fn recorded(&self, _run: usize, done: u32, _total: u32) {
        self.bar.set_position(u64::from(done));
    }
}

/// The bar of a sweep: every row of every run, and how many runs are done.
pub(crate) struct SweepBar {
    bar: ProgressBar,
    runs: usize,
    finished: AtomicUsize,
}

impl SweepBar {
    /// A bar for `runs` runs with `settings`, drawn to `target`.
    pub(crate) fn new(
        runs: usize,
        settings: RunSettings,
        target: ProgressDrawTarget,
    ) -> anyhow::Result<Self> {
        let style = ProgressStyle::with_template(
            "{spinner:.cyan} {elapsed_precise} {bar:30.cyan/blue} {percent:>3}% · {msg} · {left}",
        )
        .context("progress template")?
        .progress_chars(BAR_CHARACTERS)
        .with_key("left", time_left_key);
        let rows = u64::try_from(runs)?.saturating_mul(u64::from(settings.rows()));
        let bar = ProgressBar::with_draw_target(Some(rows), target).with_style(style);
        bar.set_message(runs_done(0, runs));
        bar.enable_steady_tick(TICK);
        Ok(Self {
            bar,
            runs,
            finished: AtomicUsize::new(0),
        })
    }

    /// Leaves the full bar on screen and ends its line, so output that
    /// follows starts on a line of its own.
    pub(crate) fn finish(self) {
        end_line(&self.bar);
    }

    /// Leaves the bar where the sweep failed.
    pub(crate) fn abandon(&self) {
        self.bar.abandon();
    }

    #[cfg(test)]
    fn position(&self) -> u64 {
        self.bar.position()
    }
}

impl Progress for SweepBar {
    fn recorded(&self, _run: usize, done: u32, total: u32) {
        // Each report is one more row.
        self.bar.inc(1);
        if done == total {
            let finished = self.finished.fetch_add(1, Ordering::Relaxed) + 1;
            self.bar.set_message(runs_done(finished, self.runs));
        }
    }
}

/// Finishes `bar` and, where it is drawn, ends its line.
fn end_line(bar: &ProgressBar) {
    bar.finish();
    if !bar.is_hidden() {
        let _ = writeln!(std::io::stderr());
    }
}

/// Simulated time when `rows` rows are recorded, s.
fn simulated(rows: u64, settings: RunSettings) -> f64 {
    let last = u32::try_from(rows.saturating_sub(1)).unwrap_or(u32::MAX);
    settings.row_time(last.min(settings.rows().saturating_sub(1)))
}

/// A time in seconds with three significant figures, e.g. `0.120 s`, `12.0 s`.
fn seconds(t: f64) -> String {
    let mut out = String::with_capacity(10);
    let _ = if t < 10.0 {
        write!(out, "{t:.3} s")
    } else if t < 100.0 {
        write!(out, "{t:.2} s")
    } else {
        write!(out, "{t:.1} s")
    };
    out
}

/// Simulated over wall time with its separator, e.g. ` · 2.4× real time` or,
/// slower than real time, ` · 1/640 real time`; blank before either has
/// passed.
fn speed(simulated: f64, wall: f64) -> String {
    if simulated <= 0.0 || wall <= 0.0 {
        return String::new();
    }
    let ratio = simulated / wall;
    if ratio >= 10.0 {
        format!(" · {ratio:.0}× real time")
    } else if ratio >= 0.995 {
        format!(" · {ratio:.1}× real time")
    } else {
        format!(" · 1/{:.0} real time", 1.0 / ratio)
    }
}

/// Writes the time left from the bar's elapsed time and fraction done.
fn time_left_key(state: &ProgressState, w: &mut dyn std::fmt::Write) {
    let fraction = state.len().map_or(0.0, |len| {
        #[expect(clippy::cast_precision_loss, reason = "row counts are far below 2^52")]
        let fraction = state.pos() as f64 / len.max(1) as f64;
        fraction
    });
    let _ = w.write_str(&time_left(state.elapsed().as_secs_f64(), fraction));
}

/// Time left after `elapsed` s with `fraction` done, assuming the pace so far:
/// `12s left`, `3m 05s left`, `1h 02m left`; `… left` until something is done.
fn time_left(elapsed: f64, fraction: f64) -> String {
    if fraction.is_nan() || fraction <= 0.0 || elapsed <= 0.0 {
        return "… left".to_owned();
    }
    let left = elapsed * (1.0 - fraction.min(1.0)) / fraction;
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a non-negative, finite number of seconds"
    )]
    let secs = left.round().min(f64::from(u32::MAX)) as u64;
    let (h, m, s) = (secs / 3600, secs / 60 % 60, secs % 60);
    if h > 0 {
        format!("{h}h {m:02}m left")
    } else if m > 0 {
        format!("{m}m {s:02}s left")
    } else {
        format!("{s}s left")
    }
}

/// `3 of 10 runs done`.
fn runs_done(finished: usize, runs: usize) -> String {
    format!("{finished} of {runs} runs done")
}

#[cfg(test)]
mod tests {
    use engine_bench_core::spec::StepContext;
    use engine_bench_physics::units::{Frequency, Time, frequency::hertz, time::second};

    use super::*;

    fn settings() -> anyhow::Result<RunSettings> {
        RunSettings::new(
            Time::new::<second>(0.5),
            Frequency::new::<hertz>(10.0),
            StepContext::offline(1e-6)?,
        )
    }

    #[test]
    fn test_run_bar_follows_the_recorded_rows() -> anyhow::Result<()> {
        let bar = RunBar::new(settings()?, ProgressDrawTarget::hidden())?;
        bar.recorded(0, 3, 6);
        anyhow::ensure!(bar.position() == 3, "{}", bar.position());
        bar.finish();
        let failed = RunBar::new(settings()?, ProgressDrawTarget::hidden())?;
        failed.abandon();
        Ok(())
    }

    #[test]
    fn test_sweep_bar_counts_rows_and_finished_runs() -> anyhow::Result<()> {
        let bar = SweepBar::new(2, settings()?, ProgressDrawTarget::hidden())?;
        for run in 0..2 {
            for done in 1..=6 {
                bar.recorded(run, done, 6);
            }
        }
        anyhow::ensure!(bar.position() == 12, "{}", bar.position());
        anyhow::ensure!(bar.finished.load(Ordering::Relaxed) == 2, "runs");
        bar.finish();
        SweepBar::new(1, settings()?, ProgressDrawTarget::hidden())?.abandon();
        Ok(())
    }

    #[test]
    fn test_simulated_time_is_the_last_recorded_row() -> anyhow::Result<()> {
        let s = settings()?;
        anyhow::ensure!(simulated(0, s).abs() < f64::EPSILON, "nothing yet");
        anyhow::ensure!(simulated(1, s).abs() < f64::EPSILON, "t = 0 row");
        anyhow::ensure!((simulated(4, s) - 0.3).abs() < 1e-12, "fourth row");
        anyhow::ensure!((simulated(99, s) - 0.5).abs() < 1e-12, "capped");
        Ok(())
    }

    #[test]
    fn test_progress_text() {
        assert_eq!(seconds(0.12), "0.120 s");
        assert_eq!(seconds(12.0), "12.00 s");
        assert_eq!(seconds(120.0), "120.0 s");
        assert_eq!(speed(2.4, 1.0), " · 2.4× real time");
        assert_eq!(speed(50.0, 2.0), " · 25× real time");
        assert_eq!(speed(0.05, 32.0), " · 1/640 real time");
        assert_eq!(speed(0.0, 1.0), "");
        assert_eq!(runs_done(3, 10), "3 of 10 runs done");
        assert_eq!(time_left(10.0, 0.5), "10s left");
        assert_eq!(time_left(60.0, 0.25), "3m 00s left");
        assert_eq!(time_left(3700.0, 0.5), "1h 01m left");
        assert_eq!(time_left(5.0, 0.0), "… left");
        assert_eq!(time_left(5.0, 1.0), "0s left");
    }
}
