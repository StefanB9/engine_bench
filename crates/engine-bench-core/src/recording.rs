//! In-memory time series of an engine's channels, in SI units.
//!
//! File formats live in the binary; this is the shape they read into and write
//! from. Rows can be drawn with
//! [`EngineId::draw`](crate::engine::EngineId::draw) for animated replay.

use anyhow::Context;

use crate::{CoreError, engine::EngineId, spec::ChannelSpec};

/// Samples of every channel of one engine over time.
#[derive(Debug, Clone, PartialEq)]
pub struct Recording {
    engine: EngineId,
    /// s, strictly increasing.
    time: Vec<f64>,
    /// One column per channel, in channel order.
    columns: Vec<Vec<f64>>,
}

impl Recording {
    /// An empty recording of `engine`, with room for `capacity` rows.
    pub fn new(engine: EngineId, capacity: usize) -> Self {
        Self {
            engine,
            time: Vec::with_capacity(capacity),
            columns: (0..engine.channels().len())
                .map(|_| Vec::with_capacity(capacity))
                .collect(),
        }
    }

    /// A recording from existing columns (as read from a file); every column
    /// must have as many rows as `time`, one column per channel of `engine`,
    /// and time must increase.
    pub fn from_columns(
        engine: EngineId,
        time: Vec<f64>,
        columns: Vec<Vec<f64>>,
    ) -> anyhow::Result<Self> {
        let expected = engine.channels().len();
        if columns.len() != expected {
            return Err(CoreError::SampleLength {
                expected,
                actual: columns.len(),
            })
            .context("recording columns");
        }
        if let Some(bad) = columns.iter().position(|c| c.len() != time.len()) {
            return Err(CoreError::InvalidSetting("column length differs from time"))
                .with_context(|| format!("column {bad} of {}", engine.key()));
        }
        if !strictly_increasing(&time) {
            return Err(CoreError::InvalidSetting("recording time must increase"))
                .context("recording time");
        }
        Ok(Self {
            engine,
            time,
            columns,
        })
    }

    /// Appends one row: `sample` must have one value per channel and `time`
    /// must exceed the previous row's time.
    pub fn push(&mut self, time: f64, sample: &[f64]) -> anyhow::Result<()> {
        if sample.len() != self.columns.len() {
            return Err(CoreError::SampleLength {
                expected: self.columns.len(),
                actual: sample.len(),
            })
            .context("recording row");
        }
        let after_last = self.time.last().is_none_or(|last| time > *last);
        if !(time.is_finite() && after_last) {
            return Err(CoreError::InvalidSetting("recording time must increase"))
                .with_context(|| format!("t = {time} s after {:?}", self.time.last()));
        }
        self.time.push(time);
        for (column, value) in self.columns.iter_mut().zip(sample) {
            column.push(*value);
        }
        Ok(())
    }

    /// The recorded engine.
    pub fn engine(&self) -> EngineId {
        self.engine
    }

    /// The recorded channels, in column order.
    pub fn channels(&self) -> &'static [ChannelSpec] {
        self.engine.channels()
    }

    /// Number of rows.
    pub fn len(&self) -> usize {
        self.time.len()
    }

    /// Whether no rows were recorded.
    pub fn is_empty(&self) -> bool {
        self.time.is_empty()
    }

    /// Row times, s.
    pub fn time(&self) -> &[f64] {
        &self.time
    }

    /// Values of channel `index`, if it exists.
    pub fn column(&self, index: usize) -> Option<&[f64]> {
        self.columns.get(index).map(Vec::as_slice)
    }

    /// Copies row `index` into `out` (one value per channel).
    pub fn row(&self, index: usize, out: &mut [f64]) -> anyhow::Result<()> {
        if out.len() != self.columns.len() {
            return Err(CoreError::SampleLength {
                expected: self.columns.len(),
                actual: out.len(),
            })
            .context("recording row");
        }
        if index >= self.len() {
            return Err(CoreError::InvalidSetting("row index past the end"))
                .with_context(|| format!("row {index} of {}", self.len()));
        }
        for (value, column) in out.iter_mut().zip(&self.columns) {
            *value = column[index];
        }
        Ok(())
    }
}

fn strictly_increasing(time: &[f64]) -> bool {
    time.iter().all(|t| t.is_finite()) && time.windows(2).all(|w| w[0] < w[1])
}

#[cfg(test)]
mod tests {
    use engine_bench_physics::units::{Time, time::second};

    use super::*;
    use crate::{
        engine::{AnyEngine, Engine},
        scene::Scene,
        spec::StepContext,
    };

    const ENGINE: EngineId = EngineId::TwoTankDemo;

    fn width() -> usize {
        ENGINE.channels().len()
    }

    /// Row `i` holds 100·i, 100·i + 1, … across the channels.
    fn numbered_row(i: u32) -> Vec<f64> {
        let mut row = vec![0.0; width()];
        let mut value = f64::from(i * 100);
        for slot in &mut row {
            *slot = value;
            value += 1.0;
        }
        row
    }

    #[test]
    fn test_recording_push_and_read_back() -> anyhow::Result<()> {
        let mut recording = Recording::new(ENGINE, 4);
        anyhow::ensure!(
            recording.is_empty() && recording.channels() == ENGINE.channels(),
            "new"
        );
        for i in 0..3 {
            recording.push(f64::from(i) * 0.1, &numbered_row(i))?;
        }
        anyhow::ensure!(
            recording.len() == 3 && recording.time() == [0.0, 0.1, 0.2],
            "{:?}",
            recording.time()
        );
        anyhow::ensure!(
            recording.column(2) == Some(&[2.0, 102.0, 202.0][..]),
            "{:?}",
            recording.column(2)
        );
        anyhow::ensure!(recording.column(width()).is_none(), "column past the end");
        let mut row = vec![0.0; width()];
        recording.row(1, &mut row)?;
        anyhow::ensure!(row == numbered_row(1), "{row:?}");
        anyhow::ensure!(recording.row(3, &mut row).is_err(), "row past the end");
        Ok(())
    }

    #[test]
    fn test_recording_rejects_bad_rows() -> anyhow::Result<()> {
        let mut recording = Recording::new(ENGINE, 2);
        let row = vec![0.0; width()];
        recording.push(1.0, &row)?;
        anyhow::ensure!(recording.push(1.0, &row).is_err(), "repeated time");
        anyhow::ensure!(recording.push(0.5, &row).is_err(), "time went back");
        anyhow::ensure!(recording.push(f64::NAN, &row).is_err(), "NaN time");
        let Err(e) = recording.push(2.0, &row[1..]) else {
            anyhow::bail!("short row accepted");
        };
        anyhow::ensure!(
            matches!(
                e.downcast_ref::<CoreError>(),
                Some(CoreError::SampleLength { .. })
            ),
            "{e:#}"
        );
        anyhow::ensure!(recording.len() == 1, "bad rows were stored");
        Ok(())
    }

    #[test]
    fn test_recording_from_columns_validates_shape() -> anyhow::Result<()> {
        let columns = vec![vec![1.0, 2.0]; width()];
        let ok = Recording::from_columns(ENGINE, vec![0.0, 1.0], columns.clone())?;
        anyhow::ensure!(ok.len() == 2, "len");
        anyhow::ensure!(
            Recording::from_columns(ENGINE, vec![0.0], columns.clone()).is_err(),
            "row mismatch"
        );
        anyhow::ensure!(
            Recording::from_columns(ENGINE, vec![0.0, 1.0], columns[1..].to_vec()).is_err(),
            "column count"
        );
        anyhow::ensure!(
            Recording::from_columns(ENGINE, vec![1.0, 0.0], columns).is_err(),
            "time order"
        );
        Ok(())
    }

    /// The replay guarantee: drawing a recorded row reproduces the live scene
    /// exactly.
    #[test]
    fn test_recording_row_redraws_the_live_scene() -> anyhow::Result<()> {
        let mut engine = AnyEngine::new(ENGINE)?;
        let mut recording = Recording::new(ENGINE, 8);
        let mut sample = vec![0.0; width()];
        let mut live = Vec::new();
        for _ in 0..5 {
            engine.advance(Time::new::<second>(0.1), &StepContext::real_time())?;
            engine.sample(&mut sample)?;
            recording.push(engine.time().get::<second>(), &sample)?;
            let mut scene = Scene::new();
            ENGINE.draw(&sample, &mut scene)?;
            live.push(scene);
        }
        let mut replayed = Scene::new();
        for (i, live_scene) in live.iter().enumerate() {
            recording.row(i, &mut sample)?;
            ENGINE.draw(&sample, &mut replayed)?;
            anyhow::ensure!(&replayed == live_scene, "row {i} redraws differently");
        }
        Ok(())
    }
}
