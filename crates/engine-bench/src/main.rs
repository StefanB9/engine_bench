//! Engine Bench entry point.

mod logging;

use tracing::info;

fn main() {
    logging::setup_logging();
    info!(
        version = engine_bench_core::VERSION,
        "engine-bench starting"
    );
}
