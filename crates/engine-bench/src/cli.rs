//! Command-line interface. Without a subcommand the GUI opens.

use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

use crate::units::Preset;

/// Physically simulate and visualise engines.
#[derive(Debug, Parser)]
#[command(name = "engine-bench", version, about)]
pub(crate) struct Cli {
    /// What to do; omit to open the GUI.
    #[command(subcommand)]
    pub(crate) command: Option<Command>,
}

/// Batch commands.
#[derive(Debug, Subcommand)]
pub(crate) enum Command {
    /// List engines, or one engine's controls and channels.
    List {
        /// Engine key, e.g. `two-tank-demo`.
        engine: Option<String>,
    },
    /// Run one engine offline and record every channel to Parquet.
    Run(RunArgs),
    /// Run one engine for several values of a control (in parallel) into one
    /// Parquet file.
    Sweep(SweepArgs),
    /// Convert a Parquet recording to CSV in chosen display units.
    Export(ExportArgs),
}

/// Options shared by `run` and `sweep`.
#[derive(Debug, Args)]
pub(crate) struct SimulationArgs {
    /// Engine key, e.g. `two-tank-demo` (see `list`).
    pub(crate) engine: String,
    /// Simulated duration, with unit (`10 s`, `500ms`); a bare number is
    /// seconds.
    #[arg(long)]
    pub(crate) duration: String,
    /// Sample rate in Hz.
    #[arg(long, default_value_t = 100.0)]
    pub(crate) rate: f64,
    /// Relative tolerance of the adaptive integrator.
    #[arg(long, default_value_t = 1e-8)]
    pub(crate) rtol: f64,
    /// Control value, repeatable: `--set valve=50%`, `--set p_a0=8bar`.
    #[arg(long = "set", value_name = "KEY=VALUE")]
    pub(crate) settings: Vec<String>,
    /// Output Parquet file.
    #[arg(short, long)]
    pub(crate) output: PathBuf,
}

/// `run` options.
#[derive(Debug, Args)]
pub(crate) struct RunArgs {
    /// Engine, duration, rate, controls, output.
    #[command(flatten)]
    pub(crate) simulation: SimulationArgs,
}

/// `sweep` options.
#[derive(Debug, Args)]
pub(crate) struct SweepArgs {
    /// Engine, duration, rate, controls, output.
    #[command(flatten)]
    pub(crate) simulation: SimulationArgs,
    /// Control to vary over evenly spaced values: `valve=10%:100%:10`.
    #[arg(long, value_name = "KEY=START:STOP:COUNT")]
    pub(crate) vary: String,
}

/// `export` options.
#[derive(Debug, Args)]
pub(crate) struct ExportArgs {
    /// Parquet recording to read.
    pub(crate) input: PathBuf,
    /// CSV file to write.
    #[arg(long)]
    pub(crate) csv: PathBuf,
    /// Display-unit preset.
    #[arg(long, value_enum, default_value_t = Preset::Si)]
    pub(crate) units: Preset,
    /// Per-quantity override, repeatable: `--unit pressure=psi`.
    #[arg(long = "unit", value_name = "QUANTITY=SYMBOL")]
    pub(crate) unit_overrides: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> anyhow::Result<Cli> {
        Ok(Cli::try_parse_from(
            std::iter::once("engine-bench").chain(args.iter().copied()),
        )?)
    }

    #[test]
    fn test_cli_without_subcommand_means_gui() -> anyhow::Result<()> {
        anyhow::ensure!(parse(&[])?.command.is_none(), "expected no command");
        Ok(())
    }

    #[test]
    fn test_cli_parses_run_with_repeated_settings() -> anyhow::Result<()> {
        let cli = parse(&[
            "run",
            "two-tank-demo",
            "--duration",
            "10 s",
            "--set",
            "valve=50%",
            "--set",
            "p_a0=8bar",
            "-o",
            "out.parquet",
        ])?;
        let Some(Command::Run(RunArgs { simulation })) = cli.command else {
            anyhow::bail!("{:?}", cli.command);
        };
        anyhow::ensure!(simulation.engine == "two-tank-demo", "{simulation:?}");
        anyhow::ensure!(
            simulation.settings == ["valve=50%", "p_a0=8bar"],
            "{simulation:?}"
        );
        anyhow::ensure!((simulation.rate - 100.0).abs() < 1e-12, "default rate");
        anyhow::ensure!((simulation.rtol - 1e-8).abs() < 1e-20, "default rtol");
        anyhow::ensure!(
            simulation.output.as_path() == std::path::Path::new("out.parquet"),
            "output"
        );
        Ok(())
    }

    #[test]
    fn test_cli_parses_sweep_and_export() -> anyhow::Result<()> {
        let sweep = parse(&[
            "sweep",
            "two-tank-demo",
            "--duration",
            "5",
            "--rate",
            "50",
            "--vary",
            "valve=10%:100%:10",
            "-o",
            "s.parquet",
        ])?;
        let Some(Command::Sweep(args)) = sweep.command else {
            anyhow::bail!("{:?}", sweep.command);
        };
        anyhow::ensure!(
            args.vary == "valve=10%:100%:10" && (args.simulation.rate - 50.0).abs() < 1e-12,
            "{args:?}"
        );

        let export = parse(&[
            "export",
            "s.parquet",
            "--csv",
            "s.csv",
            "--units",
            "imperial",
            "--unit",
            "pressure=bar",
        ])?;
        let Some(Command::Export(args)) = export.command else {
            anyhow::bail!("{:?}", export.command);
        };
        anyhow::ensure!(
            args.units == Preset::Imperial && args.unit_overrides == ["pressure=bar"],
            "{args:?}"
        );
        Ok(())
    }

    #[test]
    fn test_cli_rejects_missing_or_bad_arguments() -> anyhow::Result<()> {
        anyhow::ensure!(
            parse(&["run", "two-tank-demo"]).is_err(),
            "run without duration/output"
        );
        anyhow::ensure!(
            parse(&["export", "x.parquet", "--csv", "x.csv", "--units", "cubits"]).is_err(),
            "bad preset"
        );
        anyhow::ensure!(parse(&["fly"]).is_err(), "unknown subcommand");
        Ok(())
    }
}
