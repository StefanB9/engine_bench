//! Controls and channels of the cold-gas thruster.

use crate::spec::{ChannelRole, ChannelSpec, ChoiceOption, ControlKind, ControlSpec, QuantityKind};

pub(super) const PROPELLANT: usize = 0;
pub(super) const FIRE: usize = 1;
pub(super) const PULSE_WIDTH: usize = 2;
pub(super) const PULSE_PERIOD: usize = 3;
pub(super) const SETPOINT: usize = 4;
pub(super) const AMBIENT_PRESSURE: usize = 5;
pub(super) const THROAT_DIAMETER: usize = 6;
pub(super) const EXIT_DIAMETER: usize = 7;
pub(super) const TANK_VOLUME: usize = 8;
pub(super) const TANK_PRESSURE: usize = 9;
pub(super) const TEMPERATURE: usize = 10;
pub(super) const WALL_INSIDE: usize = 11;
pub(super) const WALL_OUTSIDE: usize = 12;
pub(super) const PROPERTIES: usize = 13;
pub(super) const RESET: usize = 14;
pub(super) const CONTROL_COUNT: usize = 15;

/// Propellants, in option order (CO₂ waits for a condensation model).
pub(super) const PROPELLANTS: [ChoiceOption; 3] = [
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

/// How fluid properties are evaluated, in option order.
const PROPERTY_MODELS: [ChoiceOption; 2] = [
    ChoiceOption {
        key: "tables",
        label: "Tables",
    },
    ChoiceOption {
        key: "exact",
        label: "Exact",
    },
];

const FIRE_MODES: [ChoiceOption; 3] = [
    ChoiceOption {
        key: "off",
        label: "Off",
    },
    ChoiceOption {
        key: "on",
        label: "Continuous",
    },
    ChoiceOption {
        key: "pulse",
        label: "Pulsed",
    },
];

const fn continuous(
    key: &'static str,
    label: &'static str,
    quantity: QuantityKind,
    (min, max, default): (f64, f64, f64),
) -> ControlSpec {
    ControlSpec {
        key,
        label,
        kind: ControlKind::Continuous {
            quantity,
            min,
            max,
            default,
        },
    }
}

/// Settings marked "on reset" shape the hardware or its initial state.
pub(crate) const CONTROLS: [ControlSpec; CONTROL_COUNT] = [
    ControlSpec {
        key: "propellant",
        label: "Propellant (on reset)",
        kind: ControlKind::Choice {
            options: &PROPELLANTS,
            default: 0,
        },
    },
    ControlSpec {
        key: "fire",
        label: "Fire",
        kind: ControlKind::Choice {
            options: &FIRE_MODES,
            default: 0,
        },
    },
    continuous(
        "pulse_width",
        "Pulse width",
        QuantityKind::Time,
        (0.005, 1.0, 0.05),
    ),
    continuous(
        "pulse_period",
        "Pulse period",
        QuantityKind::Time,
        (0.01, 5.0, 0.5),
    ),
    // 0 holds the regulator open: pure blowdown.
    continuous(
        "setpoint",
        "Regulator setpoint (0: blowdown)",
        QuantityKind::Pressure,
        (0.0, 5e6, 1e6),
    ),
    // 0 to one standard atmosphere.
    continuous(
        "p_ambient",
        "Ambient pressure",
        QuantityKind::Pressure,
        (0.0, 101_325.0, 0.0),
    ),
    continuous(
        "throat",
        "Throat diameter",
        QuantityKind::Length,
        (2e-4, 3e-3, 1e-3),
    ),
    // At most 10 mm: area ratios to a few hundred keep the exit above 10 Pa.
    continuous(
        "exit",
        "Exit diameter",
        QuantityKind::Length,
        (2e-4, 1e-2, 7e-3),
    ),
    continuous(
        "tank_volume",
        "Tank volume (on reset)",
        QuantityKind::Volume,
        (1e-4, 2e-2, 2e-3),
    ),
    continuous(
        "tank_pressure",
        "Tank pressure (on reset)",
        QuantityKind::Pressure,
        (1e6, 3e7, 2e7),
    ),
    continuous(
        "temperature",
        "Ambient and initial temperature",
        QuantityKind::Temperature,
        (200.0, 350.0, 293.15),
    ),
    // Natural convection inside a draining tank is some 10 to 100 W/(m²·K).
    continuous(
        "h_inside",
        "Tank wall heat transfer, inside",
        QuantityKind::HeatTransferCoefficient,
        (0.0, 200.0, 30.0),
    ),
    // 0 in space (radiation is not modelled); 5 to 25 in still air.
    continuous(
        "h_outside",
        "Tank wall heat transfer, outside",
        QuantityKind::HeatTransferCoefficient,
        (0.0, 50.0, 0.0),
    ),
    // Tables agree with the exact equations to ~5e-6 and run in real time;
    // exact costs ~500 times more.
    ControlSpec {
        key: "properties",
        label: "Fluid properties",
        kind: ControlKind::Choice {
            options: &PROPERTY_MODELS,
            default: 0,
        },
    },
    ControlSpec {
        key: "reset",
        label: "Reset",
        kind: ControlKind::Trigger,
    },
];

/// Sample layout; `draw` reads the same indices.
pub(crate) mod channel {
    pub(crate) const P_TANK: usize = 0;
    pub(crate) const T_TANK: usize = 1;
    pub(crate) const T_WALL: usize = 2;
    pub(crate) const M_TANK: usize = 3;
    pub(crate) const P_PLENUM: usize = 4;
    pub(crate) const T_PLENUM: usize = 5;
    pub(crate) const P_OUTLET: usize = 6;
    pub(crate) const LIFT: usize = 7;
    pub(crate) const VALVE: usize = 8;
    pub(crate) const MASS_FLOW: usize = 9;
    pub(crate) const THRUST: usize = 10;
    pub(crate) const ISP: usize = 11;
    pub(crate) const IMPULSE: usize = 12;
    pub(crate) const EXIT_MACH: usize = 13;
    pub(crate) const EXIT_PRESSURE: usize = 14;
    pub(crate) const EXIT_TEMPERATURE: usize = 15;
    pub(crate) const SUPERSATURATION: usize = 16;
    pub(crate) const SEPARATED: usize = 17;
    pub(crate) const REGIME: usize = 18;
    pub(crate) const TANK_VOLUME: usize = 19;
    pub(crate) const THROAT: usize = 20;
    pub(crate) const EXIT: usize = 21;
    pub(crate) const COUNT: usize = 22;
}

const fn observable(key: &'static str, label: &'static str, quantity: QuantityKind) -> ChannelSpec {
    ChannelSpec {
        key,
        label,
        quantity,
        role: ChannelRole::Observable,
    }
}

const fn scene_state(
    key: &'static str,
    label: &'static str,
    quantity: QuantityKind,
) -> ChannelSpec {
    ChannelSpec {
        key,
        label,
        quantity,
        role: ChannelRole::SceneState,
    }
}

pub(crate) const CHANNELS: [ChannelSpec; channel::COUNT] = [
    observable("p_tank", "Tank pressure", QuantityKind::Pressure),
    observable("t_tank", "Tank temperature", QuantityKind::Temperature),
    observable("t_wall", "Tank wall temperature", QuantityKind::Temperature),
    observable("m_tank", "Propellant in the tank", QuantityKind::Mass),
    observable("p_plenum", "Plenum pressure", QuantityKind::Pressure),
    observable("t_plenum", "Plenum temperature", QuantityKind::Temperature),
    observable("p_outlet", "Valve outlet pressure", QuantityKind::Pressure),
    observable("lift", "Regulator lift", QuantityKind::Ratio),
    observable("valve", "Valve opening", QuantityKind::Ratio),
    observable("mdot", "Mass flow", QuantityKind::MassFlow),
    observable("thrust", "Thrust", QuantityKind::Force),
    observable("isp", "Specific impulse", QuantityKind::Time),
    observable("impulse", "Total impulse", QuantityKind::Impulse),
    observable("exit_mach", "Exit Mach number", QuantityKind::Ratio),
    observable("p_exit", "Exit pressure", QuantityKind::Pressure),
    observable("t_exit", "Exit temperature", QuantityKind::Temperature),
    // p / p_sat at the exit: above 1 the frozen vapour is supersaturated.
    observable(
        "supersaturation",
        "Exit supersaturation",
        QuantityKind::Ratio,
    ),
    observable("separated", "Flow separation", QuantityKind::Ratio),
    scene_state("regime", "Nozzle regime", QuantityKind::Ratio),
    scene_state("v_tank", "Tank volume", QuantityKind::Volume),
    scene_state("d_throat", "Throat diameter", QuantityKind::Length),
    scene_state("d_exit", "Exit diameter", QuantityKind::Length),
];
