//! Controls and channels of the two-tank demo.

use crate::spec::{ChannelRole, ChannelSpec, ControlKind, ControlSpec, QuantityKind};

pub(super) const VALVE: usize = 0;
pub(super) const PRESSURE_A: usize = 1;
pub(super) const PRESSURE_B: usize = 2;
pub(super) const HEATER: usize = 3;
pub(super) const RESET: usize = 4;

/// Initial-pressure limits keep tank A, which cools as it expands, above the
/// 200 K floor of the gas data for every setting (worst case ≈ 213 K).
pub(crate) const CONTROLS: [ControlSpec; 5] = [
    ControlSpec {
        key: "valve",
        label: "Valve opening",
        kind: ControlKind::Continuous {
            quantity: QuantityKind::Ratio,
            min: 0.0,
            max: 1.0,
            default: 0.2,
        },
    },
    ControlSpec {
        key: "p_a0",
        label: "Initial pressure A (on reset)",
        kind: ControlKind::Continuous {
            quantity: QuantityKind::Pressure,
            min: 1e5,
            max: 2e6,
            default: 8e5,
        },
    },
    ControlSpec {
        key: "p_b0",
        label: "Initial pressure B (on reset)",
        kind: ControlKind::Continuous {
            quantity: QuantityKind::Pressure,
            min: 5e4,
            max: 2e6,
            default: 1e5,
        },
    },
    ControlSpec {
        key: "heater",
        label: "Heater power, tank B",
        kind: ControlKind::Continuous {
            quantity: QuantityKind::Power,
            min: 0.0,
            max: 2000.0,
            default: 0.0,
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
    pub(crate) const P_A: usize = 0;
    pub(crate) const P_B: usize = 1;
    pub(crate) const T_A: usize = 2;
    pub(crate) const T_B: usize = 3;
    pub(crate) const M_A: usize = 4;
    pub(crate) const M_B: usize = 5;
    pub(crate) const MASS_FLOW: usize = 6;
    pub(crate) const V_A: usize = 7;
    pub(crate) const V_B: usize = 8;
    pub(crate) const VALVE: usize = 9;
    pub(crate) const COUNT: usize = 10;
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
    observable("p_a", "Pressure A", QuantityKind::Pressure),
    observable("p_b", "Pressure B", QuantityKind::Pressure),
    observable("t_a", "Temperature A", QuantityKind::Temperature),
    observable("t_b", "Temperature B", QuantityKind::Temperature),
    observable("m_a", "Mass A", QuantityKind::Mass),
    observable("m_b", "Mass B", QuantityKind::Mass),
    observable("mdot", "Mass flow A to B", QuantityKind::MassFlow),
    scene_state("v_a", "Volume A", QuantityKind::Volume),
    scene_state("v_b", "Volume B", QuantityKind::Volume),
    scene_state("valve", "Valve opening", QuantityKind::Ratio),
];
