//! Recordings on disk: Parquet files (one or many runs) and CSV export in
//! display units.
//!
//! Layout: column `run` (u32), optional `sweep_value` (the swept control's
//! value, SI), `time` (s), then one column per channel key, all SI. File
//! metadata names the engine, its channel-schema version and, for sweeps, the
//! swept control.

use std::{fs::File, ops::Range, path::Path};

use anyhow::Context;
use engine_bench_core::{
    engine::EngineId,
    recording::Recording,
    spec::{ControlKind, QuantityKind},
};
use polars::prelude::{
    Column, CsvWriter, DataFrame, KeyValueMetadata, ParquetReader, ParquetWriter, SerReader,
    SerWriter,
};

use crate::units::UnitSelection;

/// Metadata key: engine key (`EngineId::key`).
const META_ENGINE: &str = "engine_bench.engine";
/// Metadata key: `EngineDescriptor::schema_version` at write time.
const META_SCHEMA: &str = "engine_bench.schema_version";
/// Metadata key: command-line key of the swept control (sweeps only).
const META_SWEEP: &str = "engine_bench.sweep_control";

/// The varied control of a sweep and its value per run.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Sweep {
    /// Control key.
    pub(crate) control: String,
    /// Value (SI) of run `i`.
    pub(crate) values: Vec<f64>,
}

/// Everything in one recording file.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RecordingFile {
    /// One recording per run, in run order.
    pub(crate) runs: Vec<Recording>,
    /// Present for sweeps.
    pub(crate) sweep: Option<Sweep>,
}

/// Writes `file` to `path` as Parquet. All runs must record the same engine; a
/// sweep needs one value per run.
pub(crate) fn write_parquet(file: &RecordingFile, path: &Path) -> anyhow::Result<()> {
    let engine = file.engine()?;
    let mut frame = file.frame(None)?;
    let mut metadata = vec![
        (META_ENGINE.to_owned(), engine.key().to_owned()),
        (
            META_SCHEMA.to_owned(),
            engine.descriptor().schema_version.to_string(),
        ),
    ];
    if let Some(sweep) = &file.sweep {
        metadata.push((META_SWEEP.to_owned(), sweep.control.clone()));
    }
    let out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    ParquetWriter::new(out)
        .with_key_value_metadata(Some(KeyValueMetadata::from_static(metadata)))
        .finish(&mut frame)
        .with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

/// Reads a Parquet file written by [`write_parquet`].
pub(crate) fn read_parquet(path: &Path) -> anyhow::Result<RecordingFile> {
    let input = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut reader = ParquetReader::new(input);
    let metadata: Vec<(String, String)> = reader
        .get_metadata()
        .with_context(|| format!("{} is not a Parquet file", path.display()))?
        .key_value_metadata()
        .iter()
        .flatten()
        .filter_map(|kv| Some((kv.key.clone(), kv.value.clone()?)))
        .collect();
    let lookup = |key: &str| {
        metadata
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    };
    let engine_key = lookup(META_ENGINE).context("not an engine-bench recording")?;
    let engine =
        EngineId::from_key(engine_key).with_context(|| format!("unknown engine {engine_key:?}"))?;
    let schema: u32 = lookup(META_SCHEMA)
        .context("missing schema version")?
        .parse()?;
    let current = engine.descriptor().schema_version;
    anyhow::ensure!(
        schema == current,
        "recording uses {engine_key} channel schema {schema}, this build has {current}"
    );
    let sweep_control = lookup(META_SWEEP).map(str::to_owned);

    let frame = reader
        .finish()
        .with_context(|| format!("reading {}", path.display()))?;
    let run = frame
        .column("run")?
        .u32()?
        .into_no_null_iter()
        .collect::<Vec<u32>>();
    let time = f64_column(&frame, "time")?;
    let channels = engine
        .channels()
        .iter()
        .map(|c| f64_column(&frame, c.key))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let ranges = run_ranges(&run)?;

    let mut runs = Vec::with_capacity(ranges.len());
    for range in &ranges {
        let columns = channels.iter().map(|c| c[range.clone()].to_vec()).collect();
        runs.push(Recording::from_columns(
            engine,
            time[range.clone()].to_vec(),
            columns,
        )?);
    }
    let sweep = match sweep_control {
        Some(control) => {
            let values = f64_column(&frame, "sweep_value")?;
            Some(Sweep {
                control,
                values: ranges.iter().map(|r| values[r.start]).collect(),
            })
        }
        None => None,
    };
    Ok(RecordingFile { runs, sweep })
}

/// Writes `file` as CSV with values in `units`; headers carry the unit, e.g.
/// `p_a [bar]`.
pub(crate) fn export_csv(
    file: &RecordingFile,
    units: &UnitSelection,
    path: &Path,
) -> anyhow::Result<()> {
    let mut frame = file.frame(Some(units))?;
    let out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    CsvWriter::new(out)
        .include_header(true)
        .finish(&mut frame)
        .with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

impl RecordingFile {
    /// The engine all runs recorded.
    fn engine(&self) -> anyhow::Result<EngineId> {
        let first = self
            .runs
            .first()
            .context("recording file has no runs")?
            .engine();
        anyhow::ensure!(
            self.runs.iter().all(|r| r.engine() == first),
            "runs record different engines"
        );
        if let Some(sweep) = &self.sweep {
            anyhow::ensure!(
                sweep.values.len() == self.runs.len(),
                "sweep has {} values for {} runs",
                sweep.values.len(),
                self.runs.len()
            );
        }
        Ok(first)
    }

    /// All runs stacked into one table: SI with plain keys, or converted to
    /// `units` with unit-labelled headers.
    fn frame(&self, units: Option<&UnitSelection>) -> anyhow::Result<DataFrame> {
        let engine = self.engine()?;
        let rows: usize = self.runs.iter().map(Recording::len).sum();
        let header = |key: &str, kind: QuantityKind| match units {
            None => key.to_owned(),
            Some(units) => match units.unit(kind).symbol {
                "" => key.to_owned(),
                symbol => format!("{key} [{symbol}]"),
            },
        };
        let convert = |value: f64, kind: QuantityKind| match units {
            None => value,
            Some(units) => units.unit(kind).to_display(value),
        };

        let mut run = Vec::with_capacity(rows);
        for (index, recording) in self.runs.iter().enumerate() {
            run.extend(std::iter::repeat_n(u32::try_from(index)?, recording.len()));
        }
        let mut columns = vec![Column::new("run".into(), run)];
        if let Some(sweep) = &self.sweep {
            let kind = control_kind(engine, &sweep.control);
            let mut values = Vec::with_capacity(rows);
            for (recording, value) in self.runs.iter().zip(&sweep.values) {
                values.extend(std::iter::repeat_n(convert(*value, kind), recording.len()));
            }
            columns.push(Column::new(header("sweep_value", kind).into(), values));
        }
        let stacked = |column: &dyn Fn(&Recording) -> &[f64], kind: QuantityKind| {
            let mut values = Vec::with_capacity(rows);
            for recording in &self.runs {
                values.extend(column(recording).iter().map(|v| convert(*v, kind)));
            }
            values
        };
        columns.push(Column::new(
            header("time", QuantityKind::Time).into(),
            stacked(&|r: &Recording| r.time(), QuantityKind::Time),
        ));
        for (index, channel) in engine.channels().iter().enumerate() {
            let values = stacked(
                &|r: &Recording| r.column(index).unwrap_or_default(),
                channel.quantity,
            );
            columns.push(Column::new(
                header(channel.key, channel.quantity).into(),
                values,
            ));
        }
        Ok(DataFrame::new(rows, columns)?)
    }
}

/// The quantity a control measures (ratio for toggles, triggers and unknown
/// keys).
fn control_kind(engine: EngineId, key: &str) -> QuantityKind {
    match engine
        .controls()
        .iter()
        .find(|c| c.key == key)
        .map(|c| c.kind)
    {
        Some(ControlKind::Continuous { quantity, .. }) => quantity,
        _ => QuantityKind::Ratio,
    }
}

fn f64_column(frame: &DataFrame, name: &str) -> anyhow::Result<Vec<f64>> {
    Ok(frame
        .column(name)
        .with_context(|| format!("missing column {name:?}"))?
        .f64()
        .with_context(|| format!("column {name:?} is not f64"))?
        .into_no_null_iter()
        .collect())
}

/// Row ranges of runs 0, 1, 2, … which must be stored contiguously and in
/// order.
fn run_ranges(run: &[u32]) -> anyhow::Result<Vec<Range<usize>>> {
    let mut ranges: Vec<Range<usize>> = Vec::new();
    for (row, &index) in run.iter().enumerate() {
        let index = usize::try_from(index)?;
        if index + 1 == ranges.len() {
            if let Some(last) = ranges.last_mut() {
                last.end = row + 1;
            }
        } else {
            anyhow::ensure!(index == ranges.len(), "runs out of order at row {row}");
            ranges.push(row..row + 1);
        }
    }
    anyhow::ensure!(!ranges.is_empty(), "recording file has no rows");
    Ok(ranges)
}

#[cfg(test)]
mod tests {
    use engine_bench_core::{
        offline::{RunSettings, run},
        spec::{QuantityKind, find_channel},
    };
    use engine_bench_physics::units::{Frequency, Time, frequency::hertz, time::second};

    use super::*;
    use crate::units::Preset;

    const DEMO: EngineId = EngineId::TwoTankDemo;

    /// A unique path in the system temp directory, removed on drop.
    struct TempPath(std::path::PathBuf);

    impl TempPath {
        fn new(name: &str) -> Self {
            Self(std::env::temp_dir().join(format!("engine-bench-{}-{name}", std::process::id())))
        }
    }

    impl Drop for TempPath {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    fn demo_run(valve: f64) -> anyhow::Result<Recording> {
        let settings = RunSettings::new(
            Time::new::<second>(0.5),
            Frequency::new::<hertz>(20.0),
            engine_bench_core::spec::StepContext::offline(1e-8)?,
        )?;
        let valve_id =
            engine_bench_core::spec::find_control(DEMO.controls(), "valve").context("valve")?;
        run(DEMO, &[(valve_id, valve)], &settings)
    }

    #[test]
    fn test_parquet_round_trips_a_single_run() -> anyhow::Result<()> {
        let file = RecordingFile {
            runs: vec![demo_run(0.3)?],
            sweep: None,
        };
        let path = TempPath::new("single.parquet");
        write_parquet(&file, &path.0)?;
        let back = read_parquet(&path.0)?;
        anyhow::ensure!(back == file, "round trip differs");
        Ok(())
    }

    #[test]
    fn test_parquet_round_trips_a_sweep() -> anyhow::Result<()> {
        let values = vec![0.2, 0.7, 1.0];
        let runs = values
            .iter()
            .map(|v| demo_run(*v))
            .collect::<anyhow::Result<Vec<_>>>()?;
        let file = RecordingFile {
            runs,
            sweep: Some(Sweep {
                control: "valve".to_owned(),
                values,
            }),
        };
        let path = TempPath::new("sweep.parquet");
        write_parquet(&file, &path.0)?;
        anyhow::ensure!(read_parquet(&path.0)? == file, "sweep round trip differs");
        Ok(())
    }

    #[test]
    fn test_write_rejects_inconsistent_files() -> anyhow::Result<()> {
        let path = TempPath::new("bad.parquet");
        let empty = RecordingFile {
            runs: vec![],
            sweep: None,
        };
        anyhow::ensure!(
            write_parquet(&empty, &path.0).is_err(),
            "empty file written"
        );
        let short_sweep = RecordingFile {
            runs: vec![demo_run(0.5)?, demo_run(0.6)?],
            sweep: Some(Sweep {
                control: "valve".to_owned(),
                values: vec![0.5],
            }),
        };
        anyhow::ensure!(
            write_parquet(&short_sweep, &path.0).is_err(),
            "sweep without values"
        );
        Ok(())
    }

    #[test]
    fn test_read_rejects_missing_or_foreign_files() -> anyhow::Result<()> {
        let path = TempPath::new("missing.parquet");
        anyhow::ensure!(read_parquet(&path.0).is_err(), "read a missing file");
        std::fs::write(&path.0, b"not parquet")?;
        anyhow::ensure!(read_parquet(&path.0).is_err(), "read garbage");
        Ok(())
    }

    #[test]
    fn test_csv_export_converts_units_and_labels_headers() -> anyhow::Result<()> {
        let recording = demo_run(0.4)?;
        let file = RecordingFile {
            runs: vec![recording.clone()],
            sweep: None,
        };
        let path = TempPath::new("export.csv");
        export_csv(&file, &UnitSelection::new(Preset::Engineering), &path.0)?;
        let text = std::fs::read_to_string(&path.0)?;
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().context("no header")?.split(',').collect();
        anyhow::ensure!(
            header[..3] == ["run", "time [s]", "p_a [bar]"],
            "{header:?}"
        );
        anyhow::ensure!(
            header.contains(&"t_b [°C]") && header.contains(&"valve [%]"),
            "{header:?}"
        );
        let first: Vec<f64> = lines
            .next()
            .context("no rows")?
            .split(',')
            .map(str::parse)
            .collect::<Result<_, _>>()?;
        let p_a = recording
            .column(find_channel(recording.channels(), "p_a").context("p_a")?)
            .context("col")?[0];
        anyhow::ensure!(
            (first[2] - p_a * 1e-5).abs() < 1e-9,
            "p_a {} vs {}",
            first[2],
            p_a * 1e-5
        );
        let t_b = recording
            .column(find_channel(recording.channels(), "t_b").context("t_b")?)
            .context("col")?[0];
        let t_index = header
            .iter()
            .position(|h| *h == "t_b [°C]")
            .context("t_b column")?;
        anyhow::ensure!(
            (first[t_index] - (t_b - 273.15)).abs() < 1e-9,
            "t_b {}",
            first[t_index]
        );
        anyhow::ensure!(text.lines().count() == 1 + recording.len(), "row count");
        let _ = QuantityKind::Ratio;
        Ok(())
    }
}
