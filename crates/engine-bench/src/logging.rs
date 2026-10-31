//! Log output to the console, filtered by `RUST_LOG` when set.

use tracing_subscriber::{EnvFilter, fmt, prelude::*};

/// Filter used when `RUST_LOG` is unset: our crates log more in debug builds.
pub(crate) fn default_directives() -> &'static str {
    if cfg!(debug_assertions) {
        "engine_bench_core=debug,engine_bench=debug"
    } else {
        "engine_bench_core=info,engine_bench=info"
    }
}

/// Installs the global tracing subscriber.
pub(crate) fn setup_logging() {
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| {
        EnvFilter::builder()
            .with_default_directive(tracing::level_filters::LevelFilter::INFO.into())
            .parse_lossy(default_directives())
    });

    tracing_subscriber::registry()
        // stderr, so command output on stdout stays clean for pipes and files.
        .with(fmt::layer().with_writer(std::io::stderr))
        .with(filter)
        .init();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_directives_parse() -> anyhow::Result<()> {
        let parsed = EnvFilter::builder().parse(default_directives());
        anyhow::ensure!(
            parsed.is_ok(),
            "default directives do not parse: {parsed:?}"
        );
        Ok(())
    }
}
