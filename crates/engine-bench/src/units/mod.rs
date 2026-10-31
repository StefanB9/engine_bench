//! Display units: what the user sees and types. Values are SI everywhere else.
//!
//! `table` lists the units of every quantity kind and the presets; this module
//! converts, selects (preset plus per-quantity overrides), parses `"8 bar"` and
//! formats values.

mod table;

use anyhow::Context;
use engine_bench_core::spec::QuantityKind;

pub(crate) use self::table::units;

/// A unit for display: `display = si · scale + offset`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct DisplayUnit {
    /// Printed after the number; empty for plain ratios.
    pub(crate) symbol: &'static str,
    scale: f64,
    offset: f64,
}

impl DisplayUnit {
    /// Converts an SI value to this unit.
    pub(crate) fn to_display(self, value: f64) -> f64 {
        value * self.scale + self.offset
    }

    /// Converts a value in this unit to SI.
    pub(crate) fn to_si(self, value: f64) -> f64 {
        (value - self.offset) / self.scale
    }
}

/// A coherent set of display units.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, clap::ValueEnum)]
pub(crate) enum Preset {
    /// Pa, K, N, rad/s, …
    #[default]
    Si,
    /// bar, °C, kN, rpm, kW, L, …
    Engineering,
    /// psi, °F, lbf, rpm, hp, ft³, …
    Imperial,
}

/// The unit in use for every quantity kind: a preset plus per-kind overrides.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct UnitSelection {
    preset: Preset,
    overrides: [Option<&'static DisplayUnit>; QuantityKind::ALL.len()],
}

impl UnitSelection {
    /// The preset's units, no overrides.
    pub(crate) fn new(preset: Preset) -> Self {
        Self {
            preset,
            overrides: [None; QuantityKind::ALL.len()],
        }
    }

    /// The preset in use.
    pub(crate) fn preset(&self) -> Preset {
        self.preset
    }

    /// Switches to `preset`, dropping per-quantity overrides.
    pub(crate) fn set_preset(&mut self, preset: Preset) {
        *self = Self::new(preset);
    }

    /// Uses the unit spelled `symbol` (or an ASCII alias) for `kind`.
    pub(crate) fn set_override(&mut self, kind: QuantityKind, symbol: &str) -> anyhow::Result<()> {
        let unit = find_unit(kind, symbol)?;
        self.overrides[kind_index(kind)] = Some(unit);
        Ok(())
    }

    /// The unit in use for `kind`.
    pub(crate) fn unit(&self, kind: QuantityKind) -> &'static DisplayUnit {
        if let Some(unit) = self.overrides[kind_index(kind)] {
            return unit;
        }
        let list = units(kind);
        table::preset_symbol(self.preset, kind)
            .and_then(|symbol| list.iter().find(|u| u.symbol == symbol))
            .unwrap_or(&list[0])
    }
}

fn kind_index(kind: QuantityKind) -> usize {
    QuantityKind::ALL
        .iter()
        .position(|k| *k == kind)
        .unwrap_or_default()
}

/// The unit of `kind` spelled `symbol` or one of its ASCII aliases.
fn find_unit(kind: QuantityKind, symbol: &str) -> anyhow::Result<&'static DisplayUnit> {
    let symbol = table::ALIASES
        .iter()
        .find(|(alias, _)| *alias == symbol)
        .map_or(symbol, |(_, canonical)| canonical);
    units(kind)
        .iter()
        .find(|u| u.symbol == symbol)
        .with_context(|| {
            let known: Vec<&str> = units(kind).iter().map(|u| u.symbol).collect();
            format!("unknown {kind:?} unit {symbol:?}; known: {known:?}")
        })
}

/// Command-line name of a quantity kind, e.g. `angular_velocity`.
pub(crate) fn quantity_name(kind: QuantityKind) -> &'static str {
    match kind {
        QuantityKind::Ratio => "ratio",
        QuantityKind::Time => "time",
        QuantityKind::Length => "length",
        QuantityKind::Area => "area",
        QuantityKind::Volume => "volume",
        QuantityKind::Mass => "mass",
        QuantityKind::MassFlow => "mass_flow",
        QuantityKind::Pressure => "pressure",
        QuantityKind::Temperature => "temperature",
        QuantityKind::Energy => "energy",
        QuantityKind::Power => "power",
        QuantityKind::Force => "force",
        QuantityKind::Velocity => "velocity",
        QuantityKind::AngularVelocity => "angular_velocity",
        QuantityKind::HeatTransferCoefficient => "heat_transfer_coefficient",
        QuantityKind::Impulse => "impulse",
    }
}

/// Parses `"8 bar"`, `"8bar"`, `"20 degC"` or a bare SI number into an SI value
/// of `kind`.
pub(crate) fn parse_value(kind: QuantityKind, text: &str) -> anyhow::Result<f64> {
    let text = text.trim();
    // Longest prefix that is a number; the rest must be a unit of `kind` (or
    // nothing: SI).
    let split = text
        .char_indices()
        .map(|(i, _)| i)
        .chain(std::iter::once(text.len()))
        .rev()
        .find(|&i| i > 0 && text[..i].trim().parse::<f64>().is_ok())
        .with_context(|| format!("no number in {text:?}"))?;
    let number: f64 = text[..split].trim().parse()?;
    if !number.is_finite() {
        anyhow::bail!("{text:?} is not a finite number");
    }
    let rest = text[split..].trim();
    let unit = if rest.is_empty() {
        &units(kind)[0]
    } else {
        find_unit(kind, rest)?
    };
    Ok(unit.to_si(number))
}

/// `value` (SI) in `unit` to four significant digits, with its symbol: `"8.000
/// bar"`.
pub(crate) fn format_value(value: f64, unit: &DisplayUnit) -> String {
    let shown = unit.to_display(value);
    let magnitude = shown.abs();
    let number = if magnitude == 0.0 {
        format!("{shown:.3}")
    } else if !(1e-3..1e6).contains(&magnitude) {
        format!("{shown:.3e}")
    } else {
        let decimals = [1e3, 1e2, 1e1, 1.0, 1e-1, 1e-2, 1e-3]
            .iter()
            .position(|limit| magnitude >= *limit)
            .unwrap_or(6);
        format!("{shown:.decimals$}")
    };
    if unit.symbol.is_empty() {
        number
    } else {
        let mut out = String::with_capacity(number.len() + 1 + unit.symbol.len());
        out.push_str(&number);
        out.push(' ');
        out.push_str(unit.symbol);
        out
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    fn si(kind: QuantityKind, symbol: &str, value: f64) -> anyhow::Result<f64> {
        let unit = units(kind)
            .iter()
            .find(|u| u.symbol == symbol)
            .with_context(|| format!("{kind:?} has no unit {symbol}"))?;
        Ok(unit.to_si(value))
    }

    fn close(got: f64, want: f64, rel: f64) -> anyhow::Result<()> {
        anyhow::ensure!(
            (got - want).abs() <= rel * want.abs().max(1e-300),
            "{got} vs {want}"
        );
        Ok(())
    }

    #[test]
    fn test_units_every_kind_lists_si_first_without_duplicates() -> anyhow::Result<()> {
        for kind in QuantityKind::ALL {
            let list = units(kind);
            let first = list
                .first()
                .with_context(|| format!("{kind:?} has no units"))?;
            anyhow::ensure!(
                first.symbol == kind.si_symbol()
                    && first.to_display(2.5).to_bits() == 2.5f64.to_bits(),
                "{kind:?} does not start with SI"
            );
            for (i, a) in list.iter().enumerate() {
                anyhow::ensure!(
                    list[i + 1..].iter().all(|b| b.symbol != a.symbol),
                    "{kind:?} repeats {}",
                    a.symbol
                );
            }
        }
        Ok(())
    }

    #[test]
    fn test_units_golden_conversions() -> anyhow::Result<()> {
        use QuantityKind as K;
        close(si(K::Pressure, "bar", 1.0)?, 1e5, 1e-15)?;
        close(
            si(K::Pressure, "psi", 14.695_948_775_513_4)?,
            101_325.0,
            1e-12,
        )?;
        close(si(K::Pressure, "atm", 1.0)?, 101_325.0, 1e-15)?;
        close(si(K::Temperature, "°C", 0.0)?, 273.15, 1e-15)?;
        close(si(K::Temperature, "°F", 32.0)?, 273.15, 1e-12)?;
        close(si(K::Temperature, "°F", 212.0)?, 373.15, 1e-12)?;
        close(si(K::Power, "hp", 1.0)?, 745.699_871_582_270_2, 1e-12)?;
        close(si(K::Force, "lbf", 1.0)?, 4.448_221_615_260_5, 1e-12)?;
        close(si(K::Impulse, "lbf·s", 1.0)?, 4.448_221_615_260_5, 1e-12)?;
        close(
            si(K::AngularVelocity, "rpm", 3000.0)?,
            314.159_265_358_979_3,
            1e-12,
        )?;
        close(si(K::Volume, "ft³", 1.0)?, 0.028_316_846_592, 1e-12)?;
        close(si(K::Volume, "L", 1.0)?, 1e-3, 1e-15)?;
        close(si(K::Energy, "BTU", 1.0)?, 1_055.055_852_62, 1e-12)?;
        close(si(K::Velocity, "mph", 60.0)?, 26.8224, 1e-12)?;
        close(
            si(K::HeatTransferCoefficient, "BTU/(h·ft²·°F)", 1.0)?,
            5.678_263_340_8,
            1e-9,
        )?;
        close(si(K::Ratio, "%", 50.0)?, 0.5, 1e-15)
    }

    #[test]
    fn test_unit_selection_presets_and_overrides() -> anyhow::Result<()> {
        let mut units = UnitSelection::new(Preset::Engineering);
        anyhow::ensure!(
            units.unit(QuantityKind::Pressure).symbol == "bar",
            "engineering pressure"
        );
        anyhow::ensure!(
            units.unit(QuantityKind::Temperature).symbol == "°C",
            "engineering T"
        );
        anyhow::ensure!(
            UnitSelection::new(Preset::Imperial)
                .unit(QuantityKind::Force)
                .symbol
                == "lbf",
            "imperial force"
        );
        anyhow::ensure!(
            UnitSelection::new(Preset::Si)
                .unit(QuantityKind::Power)
                .symbol
                == "W",
            "SI power"
        );
        units.set_override(QuantityKind::Pressure, "psi")?;
        units.set_override(QuantityKind::Temperature, "degF")?;
        anyhow::ensure!(
            units.unit(QuantityKind::Pressure).symbol == "psi",
            "override"
        );
        anyhow::ensure!(
            units.unit(QuantityKind::Temperature).symbol == "°F",
            "alias override"
        );
        anyhow::ensure!(
            units.unit(QuantityKind::Force).symbol == "kN",
            "other kinds kept"
        );
        anyhow::ensure!(
            units
                .set_override(QuantityKind::Pressure, "furlong")
                .is_err(),
            "bogus"
        );
        anyhow::ensure!(
            units.set_override(QuantityKind::Pressure, "kg").is_err(),
            "wrong kind"
        );
        units.set_preset(Preset::Imperial);
        anyhow::ensure!(units.preset() == Preset::Imperial, "preset");
        anyhow::ensure!(
            units.unit(QuantityKind::Pressure).symbol == "psi",
            "overrides kept"
        );
        Ok(())
    }

    #[test]
    fn test_quantity_names_are_unique_snake_case() -> anyhow::Result<()> {
        let names: Vec<&str> = QuantityKind::ALL
            .iter()
            .map(|k| quantity_name(*k))
            .collect();
        for (i, name) in names.iter().enumerate() {
            anyhow::ensure!(
                name.chars().all(|c| c.is_ascii_lowercase() || c == '_'),
                "{name}"
            );
            anyhow::ensure!(!names[i + 1..].contains(name), "duplicate {name}");
        }
        anyhow::ensure!(
            quantity_name(QuantityKind::MassFlow) == "mass_flow",
            "mass_flow"
        );
        anyhow::ensure!(
            quantity_name(QuantityKind::HeatTransferCoefficient) == "heat_transfer_coefficient",
            "heat_transfer_coefficient"
        );
        Ok(())
    }

    #[test]
    fn test_presets_cover_every_kind() -> anyhow::Result<()> {
        for preset in [Preset::Si, Preset::Engineering, Preset::Imperial] {
            let selection = UnitSelection::new(preset);
            for kind in QuantityKind::ALL {
                let unit = selection.unit(kind);
                anyhow::ensure!(units(kind).contains(unit), "{preset:?} {kind:?}");
            }
        }
        Ok(())
    }

    #[test]
    fn test_parse_value_accepts_units_aliases_and_bare_si() -> anyhow::Result<()> {
        use QuantityKind as K;
        close(parse_value(K::Pressure, "8 bar")?, 8e5, 1e-15)?;
        close(parse_value(K::Pressure, "8bar")?, 8e5, 1e-15)?;
        close(parse_value(K::Pressure, "250000")?, 2.5e5, 1e-15)?;
        close(parse_value(K::Pressure, "1.5e5 Pa")?, 1.5e5, 1e-15)?;
        close(parse_value(K::Temperature, "20 degC")?, 293.15, 1e-15)?;
        close(parse_value(K::Temperature, "-40 °F")?, 233.15, 1e-12)?;
        close(parse_value(K::Time, "500 ms")?, 0.5, 1e-15)?;
        close(parse_value(K::Volume, " 2 l ")?, 2e-3, 1e-15)?;
        close(parse_value(K::Ratio, "40%")?, 0.4, 1e-15)?;
        for bad in ["", "bar", "8 furlong", "8 kg", "nan", "1e400"] {
            anyhow::ensure!(parse_value(K::Pressure, bad).is_err(), "accepted {bad:?}");
        }
        Ok(())
    }

    #[test]
    fn test_format_value_uses_four_significant_digits() -> anyhow::Result<()> {
        let units = UnitSelection::new(Preset::Engineering);
        let bar = units.unit(QuantityKind::Pressure);
        for (value, expected) in [
            (8e5, "8.000 bar"),
            (123_456.0, "1.235 bar"),
            (12_345_678.0, "123.5 bar"),
            (0.0, "0.000 bar"),
            (-2e5, "-2.000 bar"),
            (1.0, "1.000e-5 bar"),
        ] {
            let shown = format_value(value, bar);
            anyhow::ensure!(
                shown == expected,
                "{value}: {shown:?}, expected {expected:?}"
            );
        }
        let ratio = units.unit(QuantityKind::Ratio);
        anyhow::ensure!(
            format_value(0.25, ratio) == "25.00 %",
            "{}",
            format_value(0.25, ratio)
        );
        let plain = UnitSelection::new(Preset::Si).unit(QuantityKind::Ratio);
        anyhow::ensure!(
            format_value(0.25, plain) == "0.2500",
            "{}",
            format_value(0.25, plain)
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_units_round_trip_through_si(kind_index in 0..14usize, value in -1e9..1e9f64) {
            let kind = QuantityKind::ALL[kind_index];
            for unit in units(kind) {
                let back = unit.to_si(unit.to_display(value));
                prop_assert!((back - value).abs() <= 1e-12 * value.abs().max(1.0),
                    "{kind:?} {}: {value} → {back}", unit.symbol);
            }
        }
    }
}
