//! Critical constants and acentric factors for the Peng–Robinson EOS, from
//! `CoolProp` 8.0.0's reference equations of state
//! (`tools/reference/pr_reference.py`).

/// Critical point and acentric factor of a species.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Critical {
    /// Critical temperature, K.
    pub temperature: f64,
    /// Critical pressure, Pa.
    pub pressure: f64,
    /// Pitzer acentric factor ω.
    pub acentric_factor: f64,
}

#[allow(
    clippy::unreadable_literal,
    clippy::excessive_precision,
    reason = "values exactly as CoolProp prints them"
)]
const TABLE: [(&str, Critical); 8] = [
    (
        "N2",
        Critical {
            temperature: 126.19199999958556,
            pressure: 3395800.444647145,
            acentric_factor: 0.0372,
        },
    ),
    (
        "O2",
        Critical {
            temperature: 154.59938983528554,
            pressure: 5046410.521187217,
            acentric_factor: 0.0222,
        },
    ),
    (
        "Ar",
        Critical {
            temperature: 150.68700000001058,
            pressure: 4863000.544876669,
            acentric_factor: -0.00219,
        },
    ),
    (
        "He",
        Critical {
            temperature: 5.195300013635951,
            pressure: 228322.7892147868,
            acentric_factor: -0.38353992538323034,
        },
    ),
    (
        "CO2",
        Critical {
            temperature: 304.1282000029807,
            pressure: 7377298.373446752,
            acentric_factor: 0.22394,
        },
    ),
    (
        "H2O",
        Critical {
            temperature: 647.0959999999873,
            pressure: 22063999.999997754,
            acentric_factor: 0.3442920843,
        },
    ),
    (
        "H2",
        Critical {
            temperature: 33.1443326883113,
            pressure: 1296357.6060553084,
            acentric_factor: -0.219,
        },
    ),
    (
        "CO",
        Critical {
            temperature: 132.85989463386605,
            pressure: 3498194.6661988837,
            acentric_factor: 0.0497,
        },
    ),
];

/// Critical constants of the species named `name` (NASA spelling), if known.
pub fn critical_constants(name: &str) -> Option<Critical> {
    TABLE.iter().find(|(n, _)| *n == name).map(|(_, c)| *c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_critical_constants_cover_the_molecular_species() -> anyhow::Result<()> {
        for name in ["N2", "O2", "Ar", "He", "CO2", "H2O", "H2", "CO"] {
            let c = critical_constants(name).ok_or_else(|| anyhow::anyhow!("{name} missing"))?;
            anyhow::ensure!(c.temperature > 0.0 && c.pressure > 0.0, "{name}: {c:?}");
        }
        anyhow::ensure!(
            critical_constants("OH").is_none(),
            "radicals have no critical point"
        );
        let co2 = critical_constants("CO2").ok_or_else(|| anyhow::anyhow!("CO2"))?;
        anyhow::ensure!((co2.temperature - 304.1282).abs() < 1e-6, "{co2:?}");
        Ok(())
    }
}
