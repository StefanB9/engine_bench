//! Engine Bench entry point: batch commands, or the GUI when no subcommand is
//! given.

mod cli;
mod commands;
mod gui;
mod logging;
mod recording_io;
mod units;

use std::{io::Write, process::ExitCode};

use clap::Parser;
use tracing::info;

use crate::cli::Cli;

fn main() -> ExitCode {
    // Parse first: --help and --version print and exit before logging starts.
    let cli = Cli::parse();
    logging::setup_logging();
    match run(cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            let _ = writeln!(std::io::stderr(), "error: {error:#}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: Cli) -> anyhow::Result<()> {
    info!(
        version = engine_bench_core::VERSION,
        "engine-bench starting"
    );
    match cli.command {
        Some(command) => commands::execute(command, &mut std::io::stdout().lock()),
        None => gui::run(),
    }
}
