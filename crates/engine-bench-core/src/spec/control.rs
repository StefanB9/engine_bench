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
    /// One of a fixed list of options; the value is its index (0.0, 1.0, …).
    Choice {
        /// The options, in index order.
        options: &'static [ChoiceOption],
        /// Index at start and after reset.
        default: usize,
    },
}

/// One option of a [`ControlKind::Choice`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChoiceOption {
    /// Command-line key, e.g. `n2` in `--set propellant=n2`.
    pub key: &'static str,
    /// Human-readable label.
    pub label: &'static str,
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
            ControlKind::Choice { default, .. } => index_value(default),
        }
    }

    /// Accepts `value` if it is valid for this control: inside the range for
    /// continuous controls, exactly 0 or 1 for toggles, an option's index for
    /// choices, any finite value for triggers.
    pub fn validate(&self, value: f64) -> anyhow::Result<f64> {
        let (min, max) = match self.kind {
            ControlKind::Continuous { min, max, .. } => (min, max),
            ControlKind::Toggle { .. } => (0.0, 1.0),
            ControlKind::Trigger => (f64::MIN, f64::MAX),
            ControlKind::Choice { options, .. } => {
                (0.0, index_value(options.len().saturating_sub(1)))
            }
        };
        let discrete_ok = match self.kind {
            ControlKind::Toggle { .. } => {
                value.to_bits() == 0.0f64.to_bits() || value.to_bits() == 1.0f64.to_bits()
            }
            ControlKind::Choice { options, .. } => choice_index(value, options.len()).is_some(),
            ControlKind::Continuous { .. } | ControlKind::Trigger => true,
        };
        if (min..=max).contains(&value) && discrete_ok {
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

/// The option index a validated choice `value` holds, among `len` options.
pub fn choice_index(value: f64, len: usize) -> Option<usize> {
    (0..len).find(|&i| index_value(i).to_bits() == value.to_bits())
}

/// The option whose key or label is `text` (ignoring case), if any.
pub fn find_choice_option(options: &[ChoiceOption], text: &str) -> Option<usize> {
    options
        .iter()
        .position(|o| o.key.eq_ignore_ascii_case(text) || o.label.eq_ignore_ascii_case(text))
}

/// An option index as a control value.
fn index_value(index: usize) -> f64 {
    // Option lists are short: the index is exact in f64.
    f64::from(u32::try_from(index).unwrap_or(u32::MAX))
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

    const GASES: [ChoiceOption; 3] = [
        ChoiceOption {
            key: "n2",
            label: "Nitrogen",
        },
        ChoiceOption {
            key: "ar",
            label: "Argon",
        },
        ChoiceOption {
            key: "he",
            label: "Helium",
        },
    ];

    const PROPELLANT: ControlSpec = ControlSpec {
        key: "propellant",
        label: "Propellant",
        kind: ControlKind::Choice {
            options: &GASES,
            default: 1,
        },
    };

    #[test]
    fn test_control_spec_choice_defaults_to_its_option_index() -> anyhow::Result<()> {
        anyhow::ensure!(same(PROPELLANT.default_value(), 1.0), "default");
        Ok(())
    }

    #[test]
    fn test_control_spec_validates_choice_indices() -> anyhow::Result<()> {
        for good in [0.0, 1.0, 2.0] {
            anyhow::ensure!(same(PROPELLANT.validate(good)?, good), "{good}");
        }
        for bad in [-1.0, 0.5, 3.0, f64::NAN] {
            anyhow::ensure!(PROPELLANT.validate(bad).is_err(), "accepted {bad}");
        }
        Ok(())
    }

    #[test]
    fn test_choice_index_reads_a_validated_value() -> anyhow::Result<()> {
        anyhow::ensure!(choice_index(2.0, GASES.len()) == Some(2), "2");
        anyhow::ensure!(choice_index(3.0, GASES.len()).is_none(), "3");
        anyhow::ensure!(choice_index(1.5, GASES.len()).is_none(), "1.5");
        Ok(())
    }

    #[test]
    fn test_find_choice_option_by_key_or_label() -> anyhow::Result<()> {
        anyhow::ensure!(find_choice_option(&GASES, "he") == Some(2), "key");
        anyhow::ensure!(find_choice_option(&GASES, "ARGON") == Some(1), "label");
        anyhow::ensure!(find_choice_option(&GASES, "xenon").is_none(), "unknown");
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
