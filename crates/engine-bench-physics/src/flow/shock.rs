//! Normal-shock relations of a calorically perfect gas.

/// State change across a normal shock with upstream Mach number M₁.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NormalShock {
    /// M₂.
    pub downstream_mach: f64,
    /// p₂/p₁.
    pub pressure_ratio: f64,
    /// ρ₂/ρ₁.
    pub density_ratio: f64,
    /// T₂/T₁.
    pub temperature_ratio: f64,
    /// p₀₂/p₀₁ (stagnation pressure loss).
    pub stagnation_pressure_ratio: f64,
}

/// Normal shock at upstream Mach `mach`. M₁ < 1 cannot sustain a shock (it
/// would decrease entropy) and is treated as M₁ = 1, the identity.
#[inline]
pub fn normal_shock(gamma: f64, mach: f64) -> NormalShock {
    // Not `mach.max(1.0)`: that would turn a NaN into a silent identity shock.
    let m2 = if mach < 1.0 { 1.0 } else { mach * mach };
    let g = gamma;
    let pressure_ratio = 1.0 + 2.0 * g / (g + 1.0) * (m2 - 1.0);
    let density_ratio = (g + 1.0) * m2 / ((g - 1.0) * m2 + 2.0);
    NormalShock {
        downstream_mach: ((1.0 + 0.5 * (g - 1.0) * m2) / (g * m2 - 0.5 * (g - 1.0))).sqrt(),
        pressure_ratio,
        density_ratio,
        temperature_ratio: pressure_ratio / density_ratio,
        stagnation_pressure_ratio: density_ratio.powf(g / (g - 1.0))
            * pressure_ratio.powf(-1.0 / (g - 1.0)),
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    /// Anderson, Modern Compressible Flow, Table A.2 (γ = 1.4, M₁ = 2).
    #[test]
    fn test_normal_shock_matches_table_at_mach_2() -> anyhow::Result<()> {
        let s = normal_shock(1.4, 2.0);
        for (what, got, want) in [
            ("M2", s.downstream_mach, 0.5774),
            ("p2/p1", s.pressure_ratio, 4.5),
            ("rho2/rho1", s.density_ratio, 2.6667),
            ("T2/T1", s.temperature_ratio, 1.6875),
            ("p02/p01", s.stagnation_pressure_ratio, 0.7209),
        ] {
            anyhow::ensure!((got - want).abs() < 1e-4, "{what}: {got} vs {want}");
        }
        Ok(())
    }

    #[test]
    fn test_normal_shock_is_identity_at_and_below_sonic() -> anyhow::Result<()> {
        for m in [1.0, 0.5] {
            let s = normal_shock(1.4, m);
            for (what, got) in [
                ("M2", s.downstream_mach),
                ("p2/p1", s.pressure_ratio),
                ("rho2/rho1", s.density_ratio),
                ("T2/T1", s.temperature_ratio),
                ("p02/p01", s.stagnation_pressure_ratio),
            ] {
                anyhow::ensure!((got - 1.0).abs() < 1e-12, "M1 = {m}: {what} = {got}");
            }
        }
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_normal_shock_conserves_and_loses_stagnation_pressure(
            gamma in 1.05..1.67f64, mach in 1.0001..10.0f64,
        ) {
            let s = normal_shock(gamma, mach);
            // Mass: ρ₁u₁ = ρ₂u₂ → (ρ₂/ρ₁)·M₂·√(T₂/T₁) = M₁.
            let mass = s.density_ratio * s.downstream_mach * s.temperature_ratio.sqrt();
            prop_assert!((mass / mach - 1.0).abs() < 1e-12, "mass {mass} vs {mach}");
            // Momentum: p₁(1 + γM₁²) = p₂(1 + γM₂²).
            let lhs = 1.0 + gamma * mach * mach;
            let rhs = s.pressure_ratio * (1.0 + gamma * s.downstream_mach.powi(2));
            prop_assert!((lhs / rhs - 1.0).abs() < 1e-12, "momentum {lhs} vs {rhs}");
            prop_assert!(s.downstream_mach < 1.0, "M2 = {}", s.downstream_mach);
            prop_assert!(s.stagnation_pressure_ratio < 1.0 && s.stagnation_pressure_ratio > 0.0);
        }
    }
}
