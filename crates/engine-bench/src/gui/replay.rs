//! Playing back a recording: a time cursor over one run, drawn with the live
//! code path.

use anyhow::Context;
use engine_bench_core::{engine::EngineId, recording::Recording};

use crate::recording_io::RecordingFile;

/// Playback state of a loaded recording.
#[derive(Debug, Clone)]
pub(crate) struct Replay {
    file: RecordingFile,
    run: usize,
    cursor: f64,
    playing: bool,
    speed: f64,
    row: Vec<f64>,
}

impl Replay {
    /// Starts at the beginning of run 0, paused, at normal speed.
    pub(crate) fn new(file: RecordingFile) -> anyhow::Result<Self> {
        let first = file.runs.first().context("recording has no runs")?;
        anyhow::ensure!(!first.is_empty(), "recording has no rows");
        let row = vec![0.0; first.channels().len()];
        let mut replay = Self {
            file,
            run: 0,
            cursor: 0.0,
            playing: false,
            speed: 1.0,
            row,
        };
        replay.seek(0.0);
        Ok(replay)
    }

    fn recording(&self) -> &Recording {
        &self.file.runs[self.run]
    }

    /// The recorded engine.
    pub(crate) fn engine(&self) -> EngineId {
        self.recording().engine()
    }

    /// The loaded file.
    pub(crate) fn file(&self) -> &RecordingFile {
        &self.file
    }

    /// Index of the run shown.
    pub(crate) fn run(&self) -> usize {
        self.run
    }

    /// Shows run `index` from its start.
    pub(crate) fn select_run(&mut self, index: usize) -> anyhow::Result<()> {
        let runs = self.file.runs.len();
        let recording = self
            .file
            .runs
            .get(index)
            .with_context(|| format!("run {index} of {runs}"))?;
        anyhow::ensure!(!recording.is_empty(), "run {index} has no rows");
        self.run = index;
        self.seek(f64::NEG_INFINITY);
        Ok(())
    }

    /// Time of the first row of the current run, s.
    fn start(&self) -> f64 {
        self.recording().time().first().copied().unwrap_or_default()
    }

    /// Time of the last row of the current run, s.
    pub(crate) fn duration(&self) -> f64 {
        self.recording().time().last().copied().unwrap_or_default()
    }

    /// Current time, s.
    pub(crate) fn cursor(&self) -> f64 {
        self.cursor
    }

    /// Moves the cursor to `time`, clamped to the run.
    pub(crate) fn seek(&mut self, time: f64) {
        self.cursor = time.clamp(self.start(), self.duration());
    }

    /// Starts or stops playback; starting at the end rewinds first.
    pub(crate) fn set_playing(&mut self, playing: bool) {
        if playing && self.cursor >= self.duration() {
            self.seek(f64::NEG_INFINITY);
        }
        self.playing = playing;
    }

    /// Whether playing.
    pub(crate) fn is_playing(&self) -> bool {
        self.playing
    }

    /// Recorded seconds per wall second (> 0).
    pub(crate) fn set_speed(&mut self, speed: f64) {
        if speed > 0.0 && speed.is_finite() {
            self.speed = speed;
        }
    }

    /// Playback speed.
    pub(crate) fn speed(&self) -> f64 {
        self.speed
    }

    /// Moves the cursor by `wall_seconds` of playback; stops at the end.
    pub(crate) fn advance(&mut self, wall_seconds: f64) {
        if !self.playing {
            return;
        }
        self.seek(self.cursor + wall_seconds * self.speed);
        if self.cursor >= self.duration() {
            self.playing = false;
        }
    }

    /// Index of the last row at or before the cursor.
    pub(crate) fn row_index(&self) -> usize {
        self.recording()
            .time()
            .partition_point(|t| *t <= self.cursor)
            .saturating_sub(1)
    }

    /// Channel values at the cursor (the row at or before it).
    pub(crate) fn sample(&mut self) -> &[f64] {
        let index = self.row_index();
        // `index` is always a valid row of the current run, so this cannot
        // fail.
        let _ = self.file.runs[self.run].row(index, &mut self.row);
        &self.row
    }
}

#[cfg(test)]
mod tests {
    use engine_bench_core::{
        offline::{RunSettings, run},
        scene::Scene,
        spec::{StepContext, find_control},
    };
    use engine_bench_physics::units::{Frequency, Time, frequency::hertz, time::second};

    use super::*;
    use crate::recording_io::Sweep;

    const DEMO: EngineId = EngineId::TwoTankDemo;

    fn recording(valve: f64) -> anyhow::Result<Recording> {
        let settings = RunSettings::new(
            Time::new::<second>(1.0),
            Frequency::new::<hertz>(10.0),
            StepContext::offline(1e-8)?,
        )?;
        run(
            DEMO,
            &[(
                find_control(DEMO.controls(), "valve").context("valve")?,
                valve,
            )],
            &settings,
        )
    }

    fn replay() -> anyhow::Result<Replay> {
        Replay::new(RecordingFile {
            runs: vec![recording(0.3)?, recording(0.9)?],
            sweep: Some(Sweep {
                control: "valve".to_owned(),
                values: vec![0.3, 0.9],
            }),
        })
    }

    #[test]
    fn test_replay_starts_paused_at_zero() -> anyhow::Result<()> {
        let r = replay()?;
        anyhow::ensure!(
            r.engine() == DEMO && r.run() == 0 && r.cursor() == 0.0 && !r.is_playing(),
            "start"
        );
        anyhow::ensure!(
            (r.duration() - 1.0).abs() < 1e-12 && (r.speed() - 1.0).abs() < 1e-15,
            "{}",
            r.duration()
        );
        anyhow::ensure!(
            Replay::new(RecordingFile {
                runs: vec![],
                sweep: None
            })
            .is_err(),
            "empty file"
        );
        Ok(())
    }

    #[test]
    fn test_replay_seek_clamps_and_finds_rows() -> anyhow::Result<()> {
        let mut r = replay()?;
        r.seek(0.55);
        anyhow::ensure!(r.row_index() == 5, "row {}", r.row_index());
        r.seek(0.6);
        anyhow::ensure!(r.row_index() == 6, "exact row time");
        r.seek(-3.0);
        anyhow::ensure!(r.cursor() == 0.0 && r.row_index() == 0, "below");
        r.seek(99.0);
        anyhow::ensure!(
            (r.cursor() - 1.0).abs() < 1e-12 && r.row_index() == 10,
            "above"
        );
        Ok(())
    }

    #[test]
    fn test_replay_sample_matches_the_recorded_row_and_draws() -> anyhow::Result<()> {
        let mut r = replay()?;
        r.seek(0.3);
        let mut expected = vec![0.0; DEMO.channels().len()];
        r.file().runs[0].row(3, &mut expected)?;
        anyhow::ensure!(r.sample() == expected.as_slice(), "sample differs");
        let mut scene = Scene::new();
        DEMO.draw(r.sample(), &mut scene)?;
        anyhow::ensure!(!scene.shapes().is_empty(), "nothing drawn");
        Ok(())
    }

    #[test]
    fn test_replay_plays_at_speed_and_stops_at_the_end() -> anyhow::Result<()> {
        let mut r = replay()?;
        r.set_speed(2.0);
        r.set_playing(true);
        r.advance(0.25);
        anyhow::ensure!((r.cursor() - 0.5).abs() < 1e-12, "cursor {}", r.cursor());
        r.advance(10.0);
        anyhow::ensure!((r.cursor() - 1.0).abs() < 1e-12 && !r.is_playing(), "end");
        r.set_playing(true);
        anyhow::ensure!(r.cursor() == 0.0 && r.is_playing(), "play at end rewinds");
        r.set_playing(false);
        r.advance(0.3);
        anyhow::ensure!(r.cursor() == 0.0, "moved while paused");
        r.set_speed(-1.0);
        anyhow::ensure!((r.speed() - 2.0).abs() < 1e-15, "bad speed accepted");
        Ok(())
    }

    #[test]
    fn test_replay_switches_runs() -> anyhow::Result<()> {
        let mut r = replay()?;
        r.seek(0.5);
        r.select_run(1)?;
        anyhow::ensure!(r.run() == 1 && r.cursor() == 0.0, "switch");
        let mut expected = vec![0.0; DEMO.channels().len()];
        r.file().runs[1].row(0, &mut expected)?;
        anyhow::ensure!(r.sample() == expected.as_slice(), "run 1 sample");
        anyhow::ensure!(r.select_run(2).is_err(), "run 2 does not exist");
        Ok(())
    }
}
