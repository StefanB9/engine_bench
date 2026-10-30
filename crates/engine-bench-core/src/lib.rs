//! Engine Bench core library: domain logic shared by the binaries.

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
