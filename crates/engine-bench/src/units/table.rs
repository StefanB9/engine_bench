//! Display units per quantity kind and the presets that pick among them.

use engine_bench_core::spec::QuantityKind;

use super::DisplayUnit;

const fn unit(symbol: &'static str, scale: f64) -> DisplayUnit {
    DisplayUnit {
        symbol,
        scale,
        offset: 0.0,
    }
}

const fn affine(symbol: &'static str, scale: f64, offset: f64) -> DisplayUnit {
    DisplayUnit {
        symbol,
        scale,
        offset,
    }
}

// Exact definitions where they exist (inch, pound, foot, gallon-free volumes,
// psi from pound-force and inch, international BTU, mechanical horsepower).
const POUND: f64 = 0.453_592_37;
const INCH: f64 = 0.0254;
const FOOT: f64 = 0.3048;
const POUND_FORCE: f64 = POUND * 9.806_65;

const RATIO: [DisplayUnit; 2] = [unit("", 1.0), unit("%", 100.0)];
const TIME: [DisplayUnit; 4] = [
    unit("s", 1.0),
    unit("ms", 1e3),
    unit("min", 1.0 / 60.0),
    unit("h", 1.0 / 3600.0),
];
const LENGTH: [DisplayUnit; 5] = [
    unit("m", 1.0),
    unit("mm", 1e3),
    unit("cm", 1e2),
    unit("in", 1.0 / INCH),
    unit("ft", 1.0 / FOOT),
];
const AREA: [DisplayUnit; 5] = [
    unit("m²", 1.0),
    unit("cm²", 1e4),
    unit("mm²", 1e6),
    unit("in²", 1.0 / (INCH * INCH)),
    unit("ft²", 1.0 / (FOOT * FOOT)),
];
const VOLUME: [DisplayUnit; 5] = [
    unit("m³", 1.0),
    unit("L", 1e3),
    unit("cm³", 1e6),
    unit("in³", 1.0 / (INCH * INCH * INCH)),
    unit("ft³", 1.0 / (FOOT * FOOT * FOOT)),
];
const MASS: [DisplayUnit; 3] = [unit("kg", 1.0), unit("g", 1e3), unit("lb", 1.0 / POUND)];
const MASS_FLOW: [DisplayUnit; 4] = [
    unit("kg/s", 1.0),
    unit("g/s", 1e3),
    unit("kg/h", 3600.0),
    unit("lb/s", 1.0 / POUND),
];
const PRESSURE: [DisplayUnit; 6] = [
    unit("Pa", 1.0),
    unit("kPa", 1e-3),
    unit("bar", 1e-5),
    unit("MPa", 1e-6),
    unit("atm", 1.0 / 101_325.0),
    unit("psi", INCH * INCH / POUND_FORCE),
];
const TEMPERATURE: [DisplayUnit; 3] = [
    unit("K", 1.0),
    affine("°C", 1.0, -273.15),
    affine("°F", 1.8, -459.67),
];
const ENERGY: [DisplayUnit; 5] = [
    unit("J", 1.0),
    unit("kJ", 1e-3),
    unit("MJ", 1e-6),
    unit("kWh", 1.0 / 3.6e6),
    unit("BTU", 1.0 / BTU),
];
const POWER: [DisplayUnit; 4] = [
    unit("W", 1.0),
    unit("kW", 1e-3),
    unit("MW", 1e-6),
    unit("hp", 1.0 / (550.0 * FOOT * POUND_FORCE)),
];
const FORCE: [DisplayUnit; 3] = [
    unit("N", 1.0),
    unit("kN", 1e-3),
    unit("lbf", 1.0 / POUND_FORCE),
];
const VELOCITY: [DisplayUnit; 4] = [
    unit("m/s", 1.0),
    unit("km/h", 3.6),
    unit("ft/s", 1.0 / FOOT),
    unit("mph", 3600.0 / (5280.0 * FOOT)),
];
const ANGULAR_VELOCITY: [DisplayUnit; 3] = [
    unit("rad/s", 1.0),
    unit("rpm", 60.0 / std::f64::consts::TAU),
    unit("deg/s", 180.0 / std::f64::consts::PI),
];

/// International BTU, J.
const BTU: f64 = 1_055.055_852_62;
const HEAT_TRANSFER_COEFFICIENT: [DisplayUnit; 2] = [
    unit("W/(m²·K)", 1.0),
    // 1 BTU/(h·ft²·°F) = BTU / (3600 s · ft² · 5/9 K).
    unit("BTU/(h·ft²·°F)", 3600.0 * FOOT * FOOT * (5.0 / 9.0) / BTU),
];

/// Units for `kind`, SI first.
pub(crate) fn units(kind: QuantityKind) -> &'static [DisplayUnit] {
    match kind {
        QuantityKind::Ratio => &RATIO,
        QuantityKind::Time => &TIME,
        QuantityKind::Length => &LENGTH,
        QuantityKind::Area => &AREA,
        QuantityKind::Volume => &VOLUME,
        QuantityKind::Mass => &MASS,
        QuantityKind::MassFlow => &MASS_FLOW,
        QuantityKind::Pressure => &PRESSURE,
        QuantityKind::Temperature => &TEMPERATURE,
        QuantityKind::Energy => &ENERGY,
        QuantityKind::Power => &POWER,
        QuantityKind::Force => &FORCE,
        QuantityKind::Velocity => &VELOCITY,
        QuantityKind::AngularVelocity => &ANGULAR_VELOCITY,
        QuantityKind::HeatTransferCoefficient => &HEAT_TRANSFER_COEFFICIENT,
    }
}

/// Symbols chosen by the engineering preset, in `QuantityKind::ALL` order.
const ENGINEERING: [&str; 15] = [
    "%",
    "s",
    "mm",
    "cm²",
    "L",
    "kg",
    "kg/s",
    "bar",
    "°C",
    "kJ",
    "kW",
    "kN",
    "m/s",
    "rpm",
    "W/(m²·K)",
];

/// Symbols chosen by the imperial preset, in `QuantityKind::ALL` order.
const IMPERIAL: [&str; 15] = [
    "%",
    "s",
    "in",
    "in²",
    "ft³",
    "lb",
    "lb/s",
    "psi",
    "°F",
    "BTU",
    "hp",
    "lbf",
    "ft/s",
    "rpm",
    "BTU/(h·ft²·°F)",
];

/// Preset symbol for `kind`: `None` means SI.
pub(super) fn preset_symbol(preset: super::Preset, kind: QuantityKind) -> Option<&'static str> {
    let index = QuantityKind::ALL.iter().position(|k| *k == kind)?;
    match preset {
        super::Preset::Si => None,
        super::Preset::Engineering => ENGINEERING.get(index).copied(),
        super::Preset::Imperial => IMPERIAL.get(index).copied(),
    }
}

/// ASCII spellings accepted when parsing, mapped to their symbols.
pub(super) const ALIASES: [(&str, &str); 14] = [
    ("m2", "m²"),
    ("cm2", "cm²"),
    ("mm2", "mm²"),
    ("in2", "in²"),
    ("ft2", "ft²"),
    ("m3", "m³"),
    ("cm3", "cm³"),
    ("in3", "in³"),
    ("ft3", "ft³"),
    ("l", "L"),
    ("degC", "°C"),
    ("C", "°C"),
    ("degF", "°F"),
    ("F", "°F"),
];
