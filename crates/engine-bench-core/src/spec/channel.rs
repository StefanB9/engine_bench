//! Values an engine reports each sample.

use super::QuantityKind;

/// Why a channel exists.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChannelRole {
    /// A result worth plotting (pressure, thrust, efficiency).
    Observable,
    /// State the scene needs to be drawn from a recording (geometry, valve
    /// positions); hidden from plot pickers by default.
    SceneState,
}

/// One reported value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelSpec {
    /// Stable column name in recordings, e.g. `p_a`.
    pub key: &'static str,
    /// Human-readable label.
    pub label: &'static str,
    /// What it measures (values are SI).
    pub quantity: QuantityKind,
    /// Why it exists.
    pub role: ChannelRole,
}

/// Index of the channel with recording key `key`, if any.
pub fn find_channel(channels: &[ChannelSpec], key: &str) -> Option<usize> {
    channels.iter().position(|c| c.key == key)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_channel_by_key() -> anyhow::Result<()> {
        let channels = [
            ChannelSpec {
                key: "p_a",
                label: "Pressure A",
                quantity: QuantityKind::Pressure,
                role: ChannelRole::Observable,
            },
            ChannelSpec {
                key: "v_a",
                label: "Volume A",
                quantity: QuantityKind::Volume,
                role: ChannelRole::SceneState,
            },
        ];
        anyhow::ensure!(find_channel(&channels, "v_a") == Some(1), "v_a");
        anyhow::ensure!(find_channel(&channels, "x").is_none(), "x");
        Ok(())
    }
}
