//! Physical quantity tags for controls and channels.

/// What a control or channel measures. Values are always SI; the tag lets
/// displays convert units and plots group axes without per-engine code.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum QuantityKind {
    /// Dimensionless (fractions, ratios, Mach number).
    Ratio,
    /// s.
    Time,
    /// m.
    Length,
    /// m².
    Area,
    /// m³.
    Volume,
    /// kg.
    Mass,
    /// kg/s.
    MassFlow,
    /// Pa.
    Pressure,
    /// K.
    Temperature,
    /// J.
    Energy,
    /// W.
    Power,
    /// N.
    Force,
    /// m/s.
    Velocity,
    /// rad/s.
    AngularVelocity,
    /// W/(m²·K).
    HeatTransferCoefficient,
}

impl QuantityKind {
    /// Every kind, for unit pickers and exhaustive tests.
    pub const ALL: [Self; 15] = [
        Self::Ratio,
        Self::Time,
        Self::Length,
        Self::Area,
        Self::Volume,
        Self::Mass,
        Self::MassFlow,
        Self::Pressure,
        Self::Temperature,
        Self::Energy,
        Self::Power,
        Self::Force,
        Self::Velocity,
        Self::AngularVelocity,
        Self::HeatTransferCoefficient,
    ];

    /// Symbol of the SI unit values are stored in (empty for dimensionless).
    pub fn si_symbol(self) -> &'static str {
        match self {
            Self::Ratio => "",
            Self::Time => "s",
            Self::Length => "m",
            Self::Area => "m²",
            Self::Volume => "m³",
            Self::Mass => "kg",
            Self::MassFlow => "kg/s",
            Self::Pressure => "Pa",
            Self::Temperature => "K",
            Self::Energy => "J",
            Self::Power => "W",
            Self::Force => "N",
            Self::Velocity => "m/s",
            Self::AngularVelocity => "rad/s",
            Self::HeatTransferCoefficient => "W/(m²·K)",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quantity_kind_si_symbols() -> anyhow::Result<()> {
        let expected = [
            "",
            "s",
            "m",
            "m²",
            "m³",
            "kg",
            "kg/s",
            "Pa",
            "K",
            "J",
            "W",
            "N",
            "m/s",
            "rad/s",
            "W/(m²·K)",
        ];
        anyhow::ensure!(expected.len() == QuantityKind::ALL.len(), "kind count");
        for (kind, symbol) in QuantityKind::ALL.into_iter().zip(expected) {
            anyhow::ensure!(
                kind.si_symbol() == symbol,
                "{kind:?}: {:?}",
                kind.si_symbol()
            );
        }
        Ok(())
    }
}
