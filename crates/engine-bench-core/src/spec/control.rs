//! Controls an engine accepts from the UI or the command line.

use anyhow::Context;

use super::QuantityKind;
use crate::CoreError;

/// Index of a control in its engine's [`ControlSpec`] list.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ControlId(usize);

impl ControlId {
    /// The control at position `index`.
    pub const fn new(index: usize) -> Self {
        Self(index)
    }

    /// Position in the engine's control list.
    pub const fn index(self) -> usize {
        self.0
    }
}

/// How a control is operated. Every value travels as `f64`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ControlKind {
    /// A continuous value in SI units.
    Continuous {
        /// What it measures.
        quantity: QuantityKind,
        /// Lowest value.
        min: f64,
        /// Highest value.
        max: f64,
        /// Value at start and after reset.
        default: f64,
    },
    /// On (1.0) or off (0.0).
    Toggle {
        /// State at start and after reset.
        default: bool,
    },
    /// A one-shot action such as reset; any value fires it.
    Trigger,
}

/// A control: stable key for the command line, label for the UI, and its kind.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ControlSpec {
    /// Command-line key, e.g. `valve` in `--set valve=0.5`.
    pub key: &'static str,
    /// Human-readable label.
    pub label: &'static str,
    /// How it is operated.
    pub kind: ControlKind,
}

impl ControlSpec {
    /// Value at start and after reset (0 for triggers).
    pub fn default_value(&self) -> f64 {
        match self.kind {
            ControlKind::Continuous { default, .. } => default,
            ControlKind::Toggle { default } => f64::from(u8::from(default)),
            ControlKind::Trigger => 0.0,
        }
    }

    /// Accepts `value` if it is valid for this control: inside the range for
    /// continuous controls, exactly 0 or 1 for toggles, any finite value for
    /// triggers.
    pub fn validate(&self, value: f64) -> anyhow::Result<f64> {
        let (min, max) = match self.kind {
            ControlKind::Continuous { min, max, .. } => (min, max),
            ControlKind::Toggle { .. } => (0.0, 1.0),
            ControlKind::Trigger => (f64::MIN, f64::MAX),
        };
        let toggle_ok = !matches!(self.kind, ControlKind::Toggle { .. })
            || value.to_bits() == 0.0f64.to_bits()
            || value.to_bits() == 1.0f64.to_bits();
        if (min..=max).contains(&value) && toggle_ok {
            Ok(value)
        } else {
            Err(CoreError::ControlOutOfRange {
                key: self.key,
                value,
                min,
                max,
            })
            .context(self.label)
        }
    }
}

/// The control with command-line `key`, if any.
pub fn find_control(controls: &[ControlSpec], key: &str) -> Option<ControlId> {
    controls
        .iter()
        .position(|c| c.key == key)
        .map(ControlId::new)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `validate` passes accepted values through unchanged, so compare bit for
    /// bit.
    fn same(a: f64, b: f64) -> bool {
        a.to_bits() == b.to_bits()
    }

    const CONTROLS: [ControlSpec; 3] = [
        ControlSpec {
            key: "valve",
            label: "Valve opening",
            kind: ControlKind::Continuous {
                quantity: QuantityKind::Ratio,
                min: 0.0,
                max: 1.0,
                default: 0.25,
            },
        },
        ControlSpec {
            key: "heater",
            label: "Heater",
            kind: ControlKind::Toggle { default: true },
        },
        ControlSpec {
            key: "reset",
            label: "Reset",
            kind: ControlKind::Trigger,
        },
    ];

    #[test]
    fn test_control_id_round_trips_index() -> anyhow::Result<()> {
        anyhow::ensure!(ControlId::new(4).index() == 4, "index");
        Ok(())
    }

    #[test]
    fn test_control_spec_default_values() -> anyhow::Result<()> {
        let defaults = CONTROLS.map(|c| c.default_value());
        anyhow::ensure!(
            defaults
                .iter()
                .zip([0.25, 1.0, 0.0])
                .all(|(a, b)| same(*a, b)),
            "{defaults:?}"
        );
        Ok(())
    }

    #[test]
    fn test_control_spec_validates_continuous_range() -> anyhow::Result<()> {
        let valve = CONTROLS[0];
        anyhow::ensure!(
            same(valve.validate(0.0)?, 0.0) && same(valve.validate(1.0)?, 1.0),
            "ends"
        );
        for bad in [-0.1, 1.1, f64::NAN] {
            let Err(e) = valve.validate(bad) else {
                anyhow::bail!("accepted {bad}");
            };
            anyhow::ensure!(
                matches!(
                    e.downcast_ref::<CoreError>(),
                    Some(CoreError::ControlOutOfRange { key: "valve", .. })
                ),
                "wrong error for {bad}: {e:#}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_control_spec_validates_toggle_and_trigger() -> anyhow::Result<()> {
        let (toggle, trigger) = (CONTROLS[1], CONTROLS[2]);
        anyhow::ensure!(
            same(toggle.validate(1.0)?, 1.0) && same(toggle.validate(0.0)?, 0.0),
            "toggle"
        );
        anyhow::ensure!(toggle.validate(0.5).is_err(), "toggle accepted 0.5");
        anyhow::ensure!(same(trigger.validate(42.0)?, 42.0), "trigger");
        anyhow::ensure!(
            trigger.validate(f64::INFINITY).is_err(),
            "trigger accepted inf"
        );
        Ok(())
    }

    #[test]
    fn test_find_control_by_key() -> anyhow::Result<()> {
        anyhow::ensure!(
            find_control(&CONTROLS, "heater") == Some(ControlId::new(1)),
            "heater"
        );
        anyhow::ensure!(find_control(&CONTROLS, "missing").is_none(), "missing");
        Ok(())
    }
}
