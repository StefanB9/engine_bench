//! The real binary, driven as a user would: list, run to Parquet, export to
//! CSV, errors.

use std::{path::PathBuf, process::Command};

use anyhow::Context;

fn engine_bench(args: &[&str]) -> anyhow::Result<std::process::Output> {
    Command::new(env!("CARGO_BIN_EXE_engine-bench"))
        .args(args)
        .env("RUST_LOG", "off")
        .output()
        .context("starting engine-bench")
}

/// A unique path in the system temp directory, removed on drop.
struct TempPath(PathBuf);

impl TempPath {
    fn new(name: &str) -> Self {
        Self(std::env::temp_dir().join(format!("engine-bench-cli-{}-{name}", std::process::id())))
    }

    fn arg(&self) -> anyhow::Result<&str> {
        self.0.to_str().context("non-UTF-8 temp path")
    }
}

impl Drop for TempPath {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

#[test]
fn test_cli_list_names_the_demo() -> anyhow::Result<()> {
    let output = engine_bench(&["list"])?;
    anyhow::ensure!(output.status.success(), "{output:?}");
    anyhow::ensure!(
        String::from_utf8(output.stdout)?.contains("two-tank-demo"),
        "missing demo"
    );
    Ok(())
}

#[test]
fn test_cli_run_then_export_in_engineering_units() -> anyhow::Result<()> {
    let parquet = TempPath::new("run.parquet");
    let csv = TempPath::new("run.csv");
    let run = engine_bench(&[
        "run",
        "two-tank-demo",
        "--duration",
        "0.5 s",
        "--rate",
        "10",
        "--set",
        "valve=50%",
        "--set",
        "p_a0=12 bar",
        "-o",
        parquet.arg()?,
    ])?;
    anyhow::ensure!(
        run.status.success(),
        "run failed: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let export = engine_bench(&[
        "export",
        parquet.arg()?,
        "--csv",
        csv.arg()?,
        "--units",
        "engineering",
        "--unit",
        "temperature=K",
    ])?;
    anyhow::ensure!(
        export.status.success(),
        "export failed: {}",
        String::from_utf8_lossy(&export.stderr)
    );

    let text = std::fs::read_to_string(&csv.0)?;
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().context("empty csv")?.split(',').collect();
    anyhow::ensure!(
        header.contains(&"p_a [bar]") && header.contains(&"t_a [K]"),
        "{header:?}"
    );
    let p_a = header
        .iter()
        .position(|h| *h == "p_a [bar]")
        .context("p_a")?;
    let first: Vec<&str> = lines.next().context("no rows")?.split(',').collect();
    let initial: f64 = first[p_a].parse()?;
    anyhow::ensure!((initial - 12.0).abs() < 1e-9, "initial p_a {initial} bar");
    anyhow::ensure!(
        text.lines().count() == 1 + 6,
        "rows: {}",
        text.lines().count()
    );
    Ok(())
}

#[test]
fn test_cli_reports_errors_with_nonzero_exit() -> anyhow::Result<()> {
    let out = TempPath::new("never.parquet");
    let output = engine_bench(&["run", "warp-drive", "--duration", "1", "-o", out.arg()?])?;
    anyhow::ensure!(!output.status.success(), "unknown engine succeeded");
    let stderr = String::from_utf8(output.stderr)?;
    anyhow::ensure!(
        stderr.contains("warp-drive") && stderr.contains("two-tank-demo"),
        "{stderr}"
    );
    anyhow::ensure!(!out.0.exists(), "wrote output despite the error");
    Ok(())
}

#[test]
fn test_cli_without_subcommand_explains_the_gui_is_not_built_yet() -> anyhow::Result<()> {
    let output = engine_bench(&[])?;
    anyhow::ensure!(!output.status.success(), "no-subcommand run succeeded");
    anyhow::ensure!(
        String::from_utf8(output.stderr)?.contains("GUI"),
        "no GUI hint"
    );
    Ok(())
}
