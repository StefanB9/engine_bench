//! Engine Bench core library: engine models and the renderer-agnostic runtime
//! shared by the binaries.
//!
//! Engines describe their controls and observable channels ([`spec`]), advance
//! themselves in time, and draw a [`scene`] from a sample of
//! their channels, so the live view and replay of a recording run the same
//! code.

pub mod engine;
pub mod engines;
pub mod error;
pub mod recording;
pub mod scene;
pub mod spec;

pub use error::CoreError;

/// Version of the engine-bench crates, from `Cargo.toml`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_version_matches_manifest() -> anyhow::Result<()> {
        anyhow::ensure!(VERSION == "0.1.0", "unexpected version {VERSION}");
        Ok(())
    }
}
