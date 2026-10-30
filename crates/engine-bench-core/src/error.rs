//! Recoverable, matchable runtime failures, carried inside `anyhow::Error`.

/// An engine or runtime request could not be honoured.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum CoreError {
    /// A control id or key does not exist on this engine.
    #[error("unknown control {0}")]
    UnknownControl(usize),
    /// A control value lies outside its allowed range.
    #[error("control {key} = {value} outside {min}..={max}")]
    ControlOutOfRange {
        /// Control key.
        key: &'static str,
        /// Rejected value.
        value: f64,
        /// Lowest allowed value.
        min: f64,
        /// Highest allowed value.
        max: f64,
    },
    /// A sample buffer has the wrong number of channels.
    #[error("sample has {actual} values, engine has {expected} channels")]
    SampleLength {
        /// Channels of the engine.
        expected: usize,
        /// Values supplied.
        actual: usize,
    },
    /// A runtime setting is out of its valid range.
    #[error("invalid setting: {0}")]
    InvalidSetting(&'static str),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_core_error_display_for_every_variant() -> anyhow::Result<()> {
        let cases = [
            (CoreError::UnknownControl(7), "unknown control 7"),
            (
                CoreError::ControlOutOfRange {
                    key: "valve",
                    value: 1.5,
                    min: 0.0,
                    max: 1.0,
                },
                "control valve = 1.5 outside 0..=1",
            ),
            (
                CoreError::SampleLength {
                    expected: 9,
                    actual: 3,
                },
                "sample has 3 values, engine has 9 channels",
            ),
            (
                CoreError::InvalidSetting("rtol must be positive"),
                "invalid setting: rtol must be positive",
            ),
        ];
        for (error, expected) in cases {
            let shown = error.to_string();
            anyhow::ensure!(shown == expected, "got {shown:?}, expected {expected:?}");
        }
        Ok(())
    }
}
