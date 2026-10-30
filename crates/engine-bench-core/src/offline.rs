//! Offline (batch) runs: advance an engine at offline accuracy and record its
//! channels at a fixed sample rate; sweeps run several variants of one control
//! in parallel.

use anyhow::Context;
use engine_bench_physics::units::{Frequency, Time, frequency::hertz, time::second};

use crate::{
    CoreError,
    engine::{AnyEngine, Engine, EngineId},
    recording::Recording,
    spec::{ControlId, StepContext},
};

/// Refuses runs whose recording would exceed this many rows (about 1 GB per
/// channel).
pub const MAX_ROWS: u32 = 100_000_000;

/// Length, sampling and accuracy of a batch run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RunSettings {
    /// Hz.
    sample_rate: f64,
    rows: u32,
    context: StepContext,
}

impl RunSettings {
    /// Runs for `duration` (> 0), sampling at `sample_rate` (> 0) from t = 0,
    /// stepping at `context` accuracy. The last sample is the last whole
    /// interval within `duration`.
    pub fn new(
        duration: Time,
        sample_rate: Frequency,
        context: StepContext,
    ) -> anyhow::Result<Self> {
        let (duration, sample_rate) = (duration.get::<second>(), sample_rate.get::<hertz>());
        if !(duration > 0.0 && duration.is_finite()) {
            return Err(CoreError::InvalidSetting("duration must be positive"))
                .with_context(|| format!("duration = {duration} s"));
        }
        if !(sample_rate > 0.0 && sample_rate.is_finite()) {
            return Err(CoreError::InvalidSetting("sample rate must be positive"))
                .with_context(|| format!("rate = {sample_rate} Hz"));
        }
        // The epsilon keeps an exact final sample (1 s at 3 Hz → t = 1 s) from
        // rounding away.
        let intervals = (duration * sample_rate + 1e-9).floor();
        if intervals >= f64::from(MAX_ROWS) {
            return Err(CoreError::InvalidSetting("run would record too many rows"))
                .with_context(|| format!("{duration} s at {sample_rate} Hz, limit {MAX_ROWS}"));
        }
        Ok(Self {
            sample_rate,
            rows: whole_rows(intervals) + 1,
            context,
        })
    }

    /// Recorded rows, including t = 0.
    pub fn rows(&self) -> u32 {
        self.rows
    }

    /// Time of row `index`, s: exactly `index / sample_rate`, so long runs do
    /// not drift.
    pub fn row_time(&self, index: u32) -> f64 {
        f64::from(index) / self.sample_rate
    }
}

/// Runs engine `id` with `controls` applied (then reset, so initial conditions
/// take effect) and records every channel.
///
/// Each advance spans the difference of consecutive exact row times rather than
/// a fixed interval, so the engine's own clock stays locked to the row times
/// over long runs.
pub fn run(
    id: EngineId,
    controls: &[(ControlId, f64)],
    settings: &RunSettings,
) -> anyhow::Result<Recording> {
    let mut engine = AnyEngine::new(id)?;
    for &(control, value) in controls {
        engine.set_control(control, value)?;
    }
    engine.reset()?;
    let mut recording = Recording::new(id, usize::try_from(settings.rows)?);
    let mut sample = vec![0.0; id.channels().len()];
    engine.sample(&mut sample)?;
    recording.push(settings.row_time(0), &sample)?;
    for index in 1..settings.rows {
        let dt = settings.row_time(index) - settings.row_time(index - 1);
        engine
            .advance(Time::new::<second>(dt), &settings.context)
            .with_context(|| format!("advancing to t = {} s", settings.row_time(index)))?;
        engine.sample(&mut sample)?;
        recording.push(settings.row_time(index), &sample)?;
    }
    Ok(recording)
}

/// Runs engine `id` once per value of control `varied`, on top of `controls`,
/// in parallel. Recordings come back in the order of `values`.
pub fn sweep(
    id: EngineId,
    controls: &[(ControlId, f64)],
    varied: ControlId,
    values: &[f64],
    settings: &RunSettings,
) -> anyhow::Result<Vec<Recording>> {
    if values.is_empty() {
        return Err(CoreError::InvalidSetting("sweep needs at least one value").into());
    }
    let workers = std::thread::available_parallelism()
        .map_or(1, std::num::NonZero::get)
        .min(values.len());
    let one = |value: f64| -> anyhow::Result<Recording> {
        let mut variant = Vec::with_capacity(controls.len() + 1);
        variant.extend_from_slice(controls);
        variant.push((varied, value));
        run(id, &variant, settings).with_context(|| format!("sweep value {value}"))
    };
    // Static striding: worker w runs values w, w + workers, …; no shared
    // mutable state.
    let mut results: Vec<Option<anyhow::Result<Recording>>> =
        (0..values.len()).map(|_| None).collect();
    std::thread::scope(|scope| -> anyhow::Result<()> {
        let handles: Vec<_> = (0..workers)
            .map(|worker| {
                scope.spawn(move || {
                    (worker..values.len())
                        .step_by(workers)
                        .map(|index| (index, one(values[index])))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        for handle in handles {
            let finished = handle
                .join()
                .map_err(|_| CoreError::InvalidSetting("sweep worker panicked"))?;
            for (index, result) in finished {
                results[index] = Some(result);
            }
        }
        Ok(())
    })?;
    results
        .into_iter()
        .map(|result| {
            result.unwrap_or_else(|| Err(CoreError::InvalidSetting("sweep value not run").into()))
        })
        .collect()
}

/// Converts a row count already checked to be whole, non-negative and below
/// [`MAX_ROWS`].
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "callers pass whole, non-negative counts below MAX_ROWS"
)]
fn whole_rows(count: f64) -> u32 {
    count as u32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::find_control;

    const DEMO: EngineId = EngineId::TwoTankDemo;

    fn settings(duration: f64, rate: f64) -> anyhow::Result<RunSettings> {
        RunSettings::new(
            Time::new::<second>(duration),
            Frequency::new::<hertz>(rate),
            StepContext::offline(1e-8)?,
        )
    }

    fn control(key: &str) -> anyhow::Result<ControlId> {
        find_control(DEMO.controls(), key).with_context(|| format!("no control {key}"))
    }

    fn column(recording: &Recording, key: &str) -> anyhow::Result<Vec<f64>> {
        let index = crate::spec::find_channel(recording.channels(), key)
            .with_context(|| format!("no channel {key}"))?;
        Ok(recording.column(index).context("column")?.to_vec())
    }

    #[test]
    fn test_run_settings_validate_and_count_rows() -> anyhow::Result<()> {
        anyhow::ensure!(settings(1.0, 100.0)?.rows() == 101, "1 s at 100 Hz");
        anyhow::ensure!(settings(1.0, 3.0)?.rows() == 4, "exact end included");
        anyhow::ensure!(settings(0.5, 3.0)?.rows() == 2, "partial interval dropped");
        for (duration, rate) in [
            (0.0, 10.0),
            (1.0, 0.0),
            (f64::NAN, 10.0),
            (1.0, -5.0),
            (1e9, 1.0),
        ] {
            anyhow::ensure!(
                settings(duration, rate).is_err(),
                "accepted {duration} s at {rate} Hz"
            );
        }
        Ok(())
    }

    #[test]
    fn test_run_samples_at_exact_times() -> anyhow::Result<()> {
        let s = settings(0.3, 10.0)?;
        let recording = run(DEMO, &[], &s)?;
        anyhow::ensure!(
            recording.engine() == DEMO && recording.len() == 4,
            "{}",
            recording.len()
        );
        let expected: Vec<f64> = (0..4).map(|i| s.row_time(i)).collect();
        anyhow::ensure!(
            recording.time() == expected.as_slice(),
            "{:?}",
            recording.time()
        );
        anyhow::ensure!(
            s.row_time(3).to_bits() == (3.0f64 / 10.0).to_bits(),
            "row time"
        );
        Ok(())
    }

    #[test]
    fn test_run_matches_manual_stepping() -> anyhow::Result<()> {
        let s = settings(0.5, 20.0)?;
        let recording = run(DEMO, &[], &s)?;
        let mut engine = AnyEngine::new(DEMO)?;
        let mut sample = vec![0.0; DEMO.channels().len()];
        let mut row = sample.clone();
        for i in 0..s.rows() {
            if i > 0 {
                let dt = s.row_time(i) - s.row_time(i - 1);
                engine.advance(Time::new::<second>(dt), &StepContext::offline(1e-8)?)?;
            }
            engine.sample(&mut sample)?;
            recording.row(usize::try_from(i)?, &mut row)?;
            anyhow::ensure!(row == sample, "row {i} differs");
        }
        Ok(())
    }

    #[test]
    fn test_run_applies_controls_including_initial_conditions() -> anyhow::Result<()> {
        let controls = [(control("valve")?, 0.0), (control("p_a0")?, 1.2e6)];
        let recording = run(DEMO, &controls, &settings(1.0, 10.0)?)?;
        let p_a = column(&recording, "p_a")?;
        anyhow::ensure!(
            (p_a[0] / 1.2e6 - 1.0).abs() < 1e-12,
            "initial p_a {}",
            p_a[0]
        );
        anyhow::ensure!(
            p_a.iter().all(|p| p.to_bits() == p_a[0].to_bits()),
            "closed valve: {p_a:?}"
        );
        Ok(())
    }

    #[test]
    fn test_run_rejects_invalid_control() -> anyhow::Result<()> {
        let result = run(DEMO, &[(control("valve")?, 2.0)], &settings(1.0, 10.0)?);
        anyhow::ensure!(result.is_err(), "valve 2.0 accepted");
        Ok(())
    }

    #[test]
    fn test_sweep_matches_individual_runs_in_order() -> anyhow::Result<()> {
        let valve = control("valve")?;
        let values = [0.2, 0.6, 1.0, 0.4];
        let s = settings(2.0, 5.0)?;
        let recordings = sweep(DEMO, &[], valve, &values, &s)?;
        anyhow::ensure!(
            recordings.len() == values.len(),
            "{} recordings",
            recordings.len()
        );
        for (recording, value) in recordings.iter().zip(values) {
            let single = run(DEMO, &[(valve, value)], &s)?;
            anyhow::ensure!(
                recording == &single,
                "valve {value} differs from a single run"
            );
        }
        let final_p_a: Vec<f64> = recordings
            .iter()
            .map(|r| column(r, "p_a").map(|c| c[c.len() - 1]))
            .collect::<anyhow::Result<_>>()?;
        // Wider valve, faster discharge: order follows the valve values 0.2,
        // 0.6, 1.0, 0.4.
        anyhow::ensure!(
            final_p_a[2] < final_p_a[1]
                && final_p_a[1] < final_p_a[3]
                && final_p_a[3] < final_p_a[0],
            "{final_p_a:?}"
        );
        Ok(())
    }

    #[test]
    fn test_sweep_rejects_empty_or_invalid_values() -> anyhow::Result<()> {
        let valve = control("valve")?;
        let s = settings(0.2, 10.0)?;
        anyhow::ensure!(sweep(DEMO, &[], valve, &[], &s).is_err(), "empty sweep");
        let Err(e) = sweep(DEMO, &[], valve, &[0.5, 3.0], &s) else {
            anyhow::bail!("valve 3.0 accepted");
        };
        anyhow::ensure!(
            format!("{e:#}").contains('3'),
            "error does not name the value: {e:#}"
        );
        Ok(())
    }
}
