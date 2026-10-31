//! What each subcommand does.

use std::io::Write;

use anyhow::Context;
use engine_bench_core::{
    engine::EngineId,
    offline::{RunSettings, run_with_progress, sweep_with_progress},
    spec::{ControlId, ControlKind, ControlSpec, QuantityKind, StepContext, find_choice_option},
};
use engine_bench_physics::units::{Frequency, Time, frequency::hertz, time::second};
use indicatif::ProgressDrawTarget;
use tracing::info;

use crate::{
    cli::{Command, ExportArgs, RunArgs, SimulationArgs, SweepArgs},
    progress::{RunBar, SweepBar},
    recording_io::{RecordingFile, Sweep, export_csv, read_parquet, write_parquet},
    units::{UnitSelection, format_value, parse_value, quantity_name, units},
};

/// Runs a batch command, writing human-readable output to `out`.
pub(crate) fn execute(command: Command, out: &mut impl Write) -> anyhow::Result<()> {
    match command {
        Command::List { engine: None } => list_engines(out),
        Command::List { engine: Some(key) } => describe_engine(resolve_engine(&key)?, out),
        Command::Run(args) => run_command(&args, out),
        Command::Sweep(args) => sweep_command(&args, out),
        Command::Export(args) => export_command(&args, out),
    }
}

fn list_engines(out: &mut impl Write) -> anyhow::Result<()> {
    for id in EngineId::ALL {
        let descriptor = id.descriptor();
        writeln!(
            out,
            "{:<20} {}: {}",
            id.key(),
            descriptor.name,
            descriptor.summary
        )?;
    }
    Ok(())
}

fn describe_engine(engine: EngineId, out: &mut impl Write) -> anyhow::Result<()> {
    let descriptor = engine.descriptor();
    writeln!(
        out,
        "{} ({})\n{}\n\ncontrols (SI):",
        descriptor.name,
        engine.key(),
        descriptor.summary
    )?;
    for control in engine.controls() {
        let detail = match control.kind {
            ControlKind::Continuous {
                quantity,
                min,
                max,
                default,
            } => {
                let si = &units(quantity)[0];
                format!(
                    "{}: {} to {}, default {}",
                    quantity_name(quantity),
                    format_value(min, si),
                    format_value(max, si),
                    format_value(default, si)
                )
            }
            ControlKind::Toggle { default } => {
                format!("on/off, default {}", if default { "on" } else { "off" })
            }
            ControlKind::Trigger => "trigger".to_owned(),
            ControlKind::Choice { options, default } => {
                let keys: Vec<&str> = options.iter().map(|o| o.key).collect();
                let default = options.get(default).map_or("", |o| o.key);
                format!("one of {}, default {default}", keys.join(", "))
            }
        };
        writeln!(out, "  {:<15} {:<34} {detail}", control.key, control.label)?;
    }
    writeln!(out, "\nchannels (SI):")?;
    for channel in engine.channels() {
        writeln!(
            out,
            "  {:<15} {:<34} {} {:?}",
            channel.key,
            channel.label,
            quantity_name(channel.quantity),
            channel.role
        )?;
    }
    Ok(())
}

/// What a `run` or `sweep` simulates.
struct Setup {
    engine: EngineId,
    controls: Vec<(ControlId, f64)>,
    settings: RunSettings,
}

/// Engine, controls and run settings from the shared simulation options.
fn simulation_setup(args: &SimulationArgs) -> anyhow::Result<Setup> {
    let engine = resolve_engine(&args.engine)?;
    let controls = args
        .settings
        .iter()
        .map(|text| parse_setting(engine, text))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let duration = parse_value(QuantityKind::Time, &args.duration).context("--duration")?;
    let settings = RunSettings::new(
        Time::new::<second>(duration),
        Frequency::new::<hertz>(args.rate),
        StepContext::offline(args.rtol).context("--rtol")?,
    )?;
    Ok(Setup {
        engine,
        controls,
        settings,
    })
}

fn run_command(args: &RunArgs, out: &mut impl Write) -> anyhow::Result<()> {
    let Setup {
        engine,
        controls,
        settings,
    } = simulation_setup(&args.simulation)?;
    info!(engine = engine.key(), rows = settings.rows(), "offline run");
    let bar = RunBar::new(settings, ProgressDrawTarget::stderr())?;
    let recording =
        run_with_progress(engine, &controls, &settings, &bar).inspect_err(|_| bar.abandon())?;
    bar.finish();
    let file = RecordingFile {
        runs: vec![recording],
        sweep: None,
    };
    write_parquet(&file, &args.simulation.output)?;
    writeln!(
        out,
        "wrote {} rows of {} channels to {}",
        settings.rows(),
        engine.channels().len(),
        args.simulation.output.display()
    )?;
    Ok(())
}

fn sweep_command(args: &SweepArgs, out: &mut impl Write) -> anyhow::Result<()> {
    let Setup {
        engine,
        controls,
        settings,
    } = simulation_setup(&args.simulation)?;
    let (varied, values) = parse_vary(engine, &args.vary).context("--vary")?;
    info!(engine = engine.key(), runs = values.len(), "offline sweep");
    let bar = SweepBar::new(values.len(), settings, ProgressDrawTarget::stderr())?;
    let runs = sweep_with_progress(engine, &controls, varied, &values, &settings, &bar)
        .inspect_err(|_| bar.abandon())?;
    bar.finish();
    let control = engine.controls()[varied.index()].key.to_owned();
    let file = RecordingFile {
        runs,
        sweep: Some(Sweep { control, values }),
    };
    write_parquet(&file, &args.simulation.output)?;
    writeln!(
        out,
        "wrote {} runs of {} rows to {}",
        file.runs.len(),
        settings.rows(),
        args.simulation.output.display()
    )?;
    Ok(())
}

fn export_command(args: &ExportArgs, out: &mut impl Write) -> anyhow::Result<()> {
    let file = read_parquet(&args.input)?;
    let mut units = UnitSelection::new(args.units);
    for text in &args.unit_overrides {
        parse_unit_override(&mut units, text).context("--unit")?;
    }
    export_csv(&file, &units, &args.csv)?;
    writeln!(out, "wrote {}", args.csv.display())?;
    Ok(())
}

/// The engine with command-line key `key`.
pub(crate) fn resolve_engine(key: &str) -> anyhow::Result<EngineId> {
    EngineId::from_key(key).with_context(|| {
        let known: Vec<&str> = EngineId::ALL.iter().map(|id| id.key()).collect();
        format!("unknown engine {key:?}; known: {}", known.join(", "))
    })
}

/// Parses `key=value` for one of `engine`'s controls; the value may carry a
/// unit.
pub(crate) fn parse_setting(engine: EngineId, text: &str) -> anyhow::Result<(ControlId, f64)> {
    let (key, value) = text
        .split_once('=')
        .with_context(|| format!("expected KEY=VALUE, got {text:?}"))?;
    let (id, spec) = control(engine, key.trim())?;
    let value = parse_control_value(spec, value.trim())
        .with_context(|| format!("value of {}", spec.key))?;
    Ok((id, value))
}

/// The control `key` of `engine`.
fn control(engine: EngineId, key: &str) -> anyhow::Result<(ControlId, &'static ControlSpec)> {
    let controls = engine.controls();
    let index = controls
        .iter()
        .position(|c| c.key == key)
        .with_context(|| {
            let known: Vec<&str> = controls.iter().map(|c| c.key).collect();
            format!(
                "{} has no control {key:?}; controls: {}",
                engine.key(),
                known.join(", ")
            )
        })?;
    Ok((ControlId::new(index), &controls[index]))
}

/// Parses a control value in the control's quantity: `on`/`off` for toggles,
/// an option's key or label for choices.
pub(crate) fn parse_control_value(spec: &ControlSpec, text: &str) -> anyhow::Result<f64> {
    match spec.kind {
        ControlKind::Continuous { quantity, .. } => parse_value(quantity, text),
        ControlKind::Toggle { .. } => match text.to_ascii_lowercase().as_str() {
            "on" | "true" | "1" => Ok(1.0),
            "off" | "false" | "0" => Ok(0.0),
            _ => anyhow::bail!("expected on or off, got {text:?}"),
        },
        ControlKind::Trigger => Ok(1.0),
        ControlKind::Choice { options, .. } => {
            let index = find_choice_option(options, text).with_context(|| {
                let keys: Vec<&str> = options.iter().map(|o| o.key).collect();
                format!("expected one of {}, got {text:?}", keys.join(", "))
            })?;
            Ok(f64::from(u32::try_from(index)?))
        }
    }
}

/// Parses `key=start:stop:count` into the control and `count` evenly spaced
/// values.
pub(crate) fn parse_vary(engine: EngineId, text: &str) -> anyhow::Result<(ControlId, Vec<f64>)> {
    let (key, range) = text
        .split_once('=')
        .with_context(|| format!("expected KEY=START:STOP:COUNT, got {text:?}"))?;
    let (id, spec) = control(engine, key.trim())?;
    let parts: Vec<&str> = range.split(':').collect();
    let [start, stop, count] = parts.as_slice() else {
        anyhow::bail!("expected START:STOP:COUNT, got {range:?}");
    };
    let start = parse_control_value(spec, start.trim()).context("sweep start")?;
    let stop = parse_control_value(spec, stop.trim()).context("sweep stop")?;
    let count: u32 = count
        .trim()
        .parse()
        .with_context(|| format!("sweep count {count:?}"))?;
    anyhow::ensure!(count >= 1, "sweep needs at least one value");
    let values = if count == 1 {
        vec![start]
    } else {
        let last = f64::from(count - 1);
        (0..count)
            .map(|i| start + (stop - start) * f64::from(i) / last)
            .collect()
    };
    Ok((id, values))
}

/// Parses `quantity=symbol` overrides on top of `preset`.
pub(crate) fn parse_unit_override(selection: &mut UnitSelection, text: &str) -> anyhow::Result<()> {
    let (name, symbol) = text
        .split_once('=')
        .with_context(|| format!("expected QUANTITY=SYMBOL, got {text:?}"))?;
    let name = name.trim();
    let kind = QuantityKind::ALL
        .into_iter()
        .find(|k| quantity_name(*k) == name)
        .with_context(|| {
            let known: Vec<&str> = QuantityKind::ALL
                .iter()
                .map(|k| quantity_name(*k))
                .collect();
            format!("unknown quantity {name:?}; known: {}", known.join(", "))
        })?;
    selection.set_override(kind, symbol.trim())
}

#[cfg(test)]
mod tests {
    use clap::Parser;

    use super::*;
    use crate::{cli::Cli, units::Preset};

    const DEMO: EngineId = EngineId::TwoTankDemo;

    struct TempPath(std::path::PathBuf);

    impl TempPath {
        fn new(name: &str) -> Self {
            Self(
                std::env::temp_dir()
                    .join(format!("engine-bench-cmd-{}-{name}", std::process::id())),
            )
        }

        fn arg(&self) -> &str {
            self.0.to_str().unwrap_or_default()
        }
    }

    impl Drop for TempPath {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    fn command(args: &[&str]) -> anyhow::Result<Command> {
        Cli::try_parse_from(std::iter::once("engine-bench").chain(args.iter().copied()))?
            .command
            .context("no subcommand")
    }

    fn run_quietly(args: &[&str]) -> anyhow::Result<String> {
        let mut out = Vec::new();
        execute(command(args)?, &mut out)?;
        Ok(String::from_utf8(out)?)
    }

    #[test]
    fn test_resolve_engine_lists_known_keys_on_error() -> anyhow::Result<()> {
        anyhow::ensure!(resolve_engine("two-tank-demo")? == DEMO, "demo");
        let Err(e) = resolve_engine("warp-drive") else {
            anyhow::bail!("resolved an unknown engine");
        };
        anyhow::ensure!(format!("{e:#}").contains("two-tank-demo"), "{e:#}");
        Ok(())
    }

    #[test]
    fn test_parse_setting_converts_units() -> anyhow::Result<()> {
        let (id, value) = parse_setting(DEMO, "p_a0=8bar")?;
        anyhow::ensure!(
            DEMO.controls()[id.index()].key == "p_a0" && (value - 8e5).abs() < 1e-9,
            "{value}"
        );
        let (_, opening) = parse_setting(DEMO, "valve = 50%")?;
        anyhow::ensure!((opening - 0.5).abs() < 1e-15, "{opening}");
        for bad in ["valve", "nozzle=1", "valve=lots", "p_a0=8 kg"] {
            anyhow::ensure!(parse_setting(DEMO, bad).is_err(), "accepted {bad:?}");
        }
        let Err(e) = parse_setting(DEMO, "nozzle=1") else {
            anyhow::bail!("accepted nozzle");
        };
        anyhow::ensure!(
            format!("{e:#}").contains("valve"),
            "error does not list controls: {e:#}"
        );
        Ok(())
    }

    #[test]
    fn test_parse_control_value_handles_toggles_and_triggers() -> anyhow::Result<()> {
        let toggle = ControlSpec {
            key: "pump",
            label: "Pump",
            kind: ControlKind::Toggle { default: false },
        };
        for (text, expected) in [("on", 1.0), ("off", 0.0), ("true", 1.0), ("0", 0.0)] {
            let value = parse_control_value(&toggle, text)?;
            anyhow::ensure!(value.to_bits() == f64::to_bits(expected), "{text}: {value}");
        }
        anyhow::ensure!(parse_control_value(&toggle, "maybe").is_err(), "maybe");
        let trigger = ControlSpec {
            key: "reset",
            label: "Reset",
            kind: ControlKind::Trigger,
        };
        anyhow::ensure!(
            parse_control_value(&trigger, "now")?.to_bits() == 1.0f64.to_bits(),
            "trigger fires"
        );
        Ok(())
    }

    #[test]
    fn test_parse_control_value_reads_choice_keys_and_labels() -> anyhow::Result<()> {
        use engine_bench_core::spec::ChoiceOption;

        const OPTIONS: [ChoiceOption; 2] = [
            ChoiceOption {
                key: "n2",
                label: "Nitrogen",
            },
            ChoiceOption {
                key: "he",
                label: "Helium",
            },
        ];
        let choice = ControlSpec {
            key: "propellant",
            label: "Propellant",
            kind: ControlKind::Choice {
                options: &OPTIONS,
                default: 0,
            },
        };
        for (text, expected) in [("he", 1.0), ("N2", 0.0), ("helium", 1.0)] {
            let value = parse_control_value(&choice, text)?;
            anyhow::ensure!(value.to_bits() == f64::to_bits(expected), "{text}: {value}");
        }
        let Err(e) = parse_control_value(&choice, "xenon") else {
            anyhow::bail!("accepted xenon");
        };
        anyhow::ensure!(format!("{e:#}").contains("n2, he"), "{e:#}");
        Ok(())
    }

    #[test]
    fn test_parse_vary_spaces_values_evenly() -> anyhow::Result<()> {
        let (id, values) = parse_vary(DEMO, "valve=10%:100%:4")?;
        anyhow::ensure!(DEMO.controls()[id.index()].key == "valve", "control");
        let expected = [0.1, 0.4, 0.7, 1.0];
        anyhow::ensure!(
            values.len() == 4
                && values
                    .iter()
                    .zip(expected)
                    .all(|(a, b)| (a - b).abs() < 1e-12),
            "{values:?}"
        );
        let (_, single) = parse_vary(DEMO, "p_a0=5bar:5bar:1")?;
        anyhow::ensure!(
            single.len() == 1 && (single[0] - 5e5).abs() < 1e-9,
            "{single:?}"
        );
        for bad in [
            "valve=0:1",
            "valve=0:1:0",
            "valve=0:1:x",
            "valve",
            "nozzle=0:1:2",
        ] {
            anyhow::ensure!(parse_vary(DEMO, bad).is_err(), "accepted {bad:?}");
        }
        Ok(())
    }

    #[test]
    fn test_parse_unit_override() -> anyhow::Result<()> {
        let mut units = UnitSelection::new(Preset::Si);
        parse_unit_override(&mut units, "pressure=psi")?;
        parse_unit_override(&mut units, "angular_velocity=rpm")?;
        anyhow::ensure!(
            units.unit(QuantityKind::Pressure).symbol == "psi",
            "pressure"
        );
        anyhow::ensure!(
            units.unit(QuantityKind::AngularVelocity).symbol == "rpm",
            "rpm"
        );
        for bad in ["pressure", "colour=red", "pressure=kg"] {
            anyhow::ensure!(
                parse_unit_override(&mut units, bad).is_err(),
                "accepted {bad:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_list_shows_engines_and_their_controls() -> anyhow::Result<()> {
        let all = run_quietly(&["list"])?;
        anyhow::ensure!(all.contains("two-tank-demo"), "{all}");
        let one = run_quietly(&["list", "two-tank-demo"])?;
        // Ranges are listed in SI.
        anyhow::ensure!(
            one.contains("valve") && one.contains("p_a") && !one.contains("bar"),
            "{one}"
        );
        Ok(())
    }

    #[test]
    fn test_run_command_writes_the_offline_recording() -> anyhow::Result<()> {
        let path = TempPath::new("run.parquet");
        let report = run_quietly(&[
            "run",
            "two-tank-demo",
            "--duration",
            "300ms",
            "--rate",
            "20",
            "--set",
            "valve=50%",
            "-o",
            path.arg(),
        ])?;
        anyhow::ensure!(report.contains("7 rows"), "{report}");
        let file = read_parquet(&path.0)?;
        let settings = RunSettings::new(
            Time::new::<second>(0.3),
            Frequency::new::<hertz>(20.0),
            StepContext::offline(1e-8)?,
        )?;
        let (valve, _) = parse_setting(DEMO, "valve=0.5")?;
        let expected = engine_bench_core::offline::run(DEMO, &[(valve, 0.5)], &settings)?;
        anyhow::ensure!(
            file.runs == [expected] && file.sweep.is_none(),
            "recording differs"
        );
        Ok(())
    }

    #[test]
    fn test_sweep_and_export_commands() -> anyhow::Result<()> {
        let parquet = TempPath::new("sweep.parquet");
        let csv = TempPath::new("sweep.csv");
        run_quietly(&[
            "sweep",
            "two-tank-demo",
            "--duration",
            "0.2",
            "--rate",
            "10",
            "--vary",
            "valve=20%:80%:3",
            "-o",
            parquet.arg(),
        ])?;
        let file = read_parquet(&parquet.0)?;
        let Some(Sweep { control, values }) = &file.sweep else {
            anyhow::bail!("no sweep metadata");
        };
        anyhow::ensure!(
            control == "valve" && values.len() == 3 && file.runs.len() == 3,
            "{values:?}"
        );
        run_quietly(&[
            "export",
            parquet.arg(),
            "--csv",
            csv.arg(),
            "--units",
            "engineering",
        ])?;
        let text = std::fs::read_to_string(&csv.0)?;
        let header = text.lines().next().context("empty csv")?;
        anyhow::ensure!(
            header.starts_with("run,sweep_value [%],time [s],p_a [bar]"),
            "{header}"
        );
        anyhow::ensure!(
            text.lines().count() == 1 + 3 * 3,
            "rows: {}",
            text.lines().count()
        );
        Ok(())
    }

    #[test]
    fn test_commands_report_bad_input() -> anyhow::Result<()> {
        let path = TempPath::new("bad.parquet");
        for args in [
            vec!["run", "warp-drive", "--duration", "1", "-o", path.arg()],
            vec!["run", "two-tank-demo", "--duration", "-1", "-o", path.arg()],
            vec![
                "run",
                "two-tank-demo",
                "--duration",
                "1",
                "--set",
                "valve=200%",
                "-o",
                path.arg(),
            ],
            vec!["export", path.arg(), "--csv", path.arg()],
        ] {
            anyhow::ensure!(run_quietly(&args).is_err(), "accepted {args:?}");
        }
        Ok(())
    }
}
