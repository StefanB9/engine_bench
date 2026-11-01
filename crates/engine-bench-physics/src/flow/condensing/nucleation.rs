//! Homogeneous nucleation: the classical (Volmer–Frenkel–Zeldovich) rate,
//! with Kantrowitz's non-isothermal correction, as in the Wet Steam Modeling
//! Project's equations (2) to (6).

use std::f64::consts::PI;

/// Boltzmann's constant, J/K.
pub(crate) const BOLTZMANN: f64 = 1.380_649e-23;

/// The vapour and liquid at the vapour's temperature, raw SI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Supersaturated {
    /// Vapour pressure, Pa, and temperature, K.
    pub(crate) p: f64,
    pub(crate) t: f64,
    /// Vapour density, kg/m³, and its cp/cv.
    pub(crate) rho: f64,
    pub(crate) gamma: f64,
    /// Saturation pressure at `t`, Pa.
    pub(crate) p_sat: f64,
    /// Flat-film surface tension, N/m, liquid density, kg/m³, and latent
    /// heat, J/kg, at `t`.
    pub(crate) sigma: f64,
    pub(crate) rho_liquid: f64,
    pub(crate) latent: f64,
    /// Mass of a molecule, kg.
    pub(crate) molecule: f64,
}

/// Nucleation rate, 1/(m³·s), and critical radius, m (infinite at or below
/// saturation, where the rate is zero).
pub(crate) fn nucleation(
    v: &Supersaturated,
    condensation_coefficient: f64,
    kantrowitz: bool,
) -> (f64, f64) {
    let ln_s = (v.p / v.p_sat).ln();
    if ln_s.is_nan() || ln_s <= 0.0 {
        return (0.0, f64::INFINITY);
    }
    let gas_constant = BOLTZMANN / v.molecule;
    let r_crit = 2.0 * v.sigma / (v.rho_liquid * gas_constant * v.t * ln_s);
    let barrier = 4.0 * PI * r_crit * r_crit * v.sigma / (3.0 * BOLTZMANN * v.t);
    let prefactor = condensation_coefficient * v.rho * v.rho / v.rho_liquid
        * (2.0 * v.sigma / (PI * v.molecule.powi(3))).sqrt();
    let correction = if kantrowitz {
        let x = v.latent / (gas_constant * v.t);
        1.0 + 2.0 * (v.gamma - 1.0) / (v.gamma + 1.0) * x * (x - 0.5)
    } else {
        1.0
    };
    (prefactor * (-barrier).exp() / correction, r_crit)
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    /// Steam at 30 kPa, 300 K: about 42 K below saturation (342 K).
    fn steam(p: f64) -> Supersaturated {
        Supersaturated {
            p,
            t: 300.0,
            rho: p / (461.5 * 300.0),
            gamma: 1.32,
            p_sat: 3_537.0,
            sigma: 0.0717,
            rho_liquid: 996.5,
            latent: 2.4375e6,
            molecule: 2.9915e-26,
        }
    }

    #[test]
    fn test_nucleation_needs_supersaturation() {
        let (j, r) = nucleation(&steam(3_537.0), 1.0, true);
        assert!(j.abs() < f64::MIN_POSITIVE && r.is_infinite());
        let (j, _) = nucleation(&steam(2_000.0), 1.0, true);
        assert!(j.abs() < f64::MIN_POSITIVE);
    }

    #[test]
    fn test_nucleation_kantrowitz_divides_by_about_a_hundred_for_steam() -> anyhow::Result<()> {
        // The project: the correction "typically reduces nucleation by two
        // orders of magnitude". Here L/(RT) = 17.6: 1 + φ = 1 +
        // 0.276·17.6·17.1.
        let v = steam(15_000.0);
        let (with, _) = nucleation(&v, 1.0, true);
        let (without, _) = nucleation(&v, 1.0, false);
        let ratio = without / with;
        anyhow::ensure!((ratio / 84.0 - 1.0).abs() < 0.02, "ratio {ratio}");
        // J scales with q_c.
        let (half, _) = nucleation(&v, 0.5, true);
        anyhow::ensure!((half / with - 0.5).abs() < 1e-12, "q_c");
        Ok(())
    }

    #[test]
    fn test_nucleation_critical_radius_is_kelvins() -> anyhow::Result<()> {
        // A cluster of radius r* is in equilibrium with the vapour:
        // p = pₛ·exp(2σ / (ρ_l R T r*)).
        let v = steam(12_000.0);
        let (_, r) = nucleation(&v, 1.0, true);
        let p = v.p_sat * (2.0 * v.sigma / (v.rho_liquid * BOLTZMANN / v.molecule * v.t * r)).exp();
        anyhow::ensure!((p / v.p - 1.0).abs() < 1e-12, "{p}");
        // About a nanometre at these supersaturations (S ≈ 3.4).
        anyhow::ensure!(r > 0.5e-9 && r < 1.5e-9, "r* = {r}");
        Ok(())
    }

    #[test]
    #[expect(clippy::many_single_char_names, reason = "the symbols p, t, v, x, j")]
    fn test_nucleation_rate_of_steam_matches_the_wet_steam_project() -> anyhow::Result<()> {
        use crate::{
            flow::condensing::{Condensable, kernel::VapourKernel},
            steam::{Water, saturation_temperature},
            units::{
                Pressure, ThermodynamicTemperature, available_energy::joule_per_kilogram,
                mass::kilogram, mass_density::kilogram_per_cubic_meter, pressure::pascal,
                surface_tension::newton_per_meter, thermodynamic_temperature::kelvin,
            },
        };
        // Starzmann et al. (2018), Fig. 14, curve "CL, NISO, f_σ = 1.00"
        // (classical rate with Kantrowitz's correction, IF97, IAPWS surface
        // tension) at 50 kPa: (subcooling K, log₁₀ J) read from the plot to
        // about 0.1 K and 0.1 decade.
        const READ: [(f64, f64); 5] = [
            (23.83, 12.38),
            (26.78, 16.01),
            (30.47, 19.23),
            (34.89, 21.83),
            (37.84, 23.14),
        ];
        let p = 5e4;
        let t_sat = saturation_temperature(Pressure::new::<pascal>(p))?.get::<kelvin>();
        for (subcooling, log_j) in READ {
            let t = t_sat - subcooling;
            let v = Water.vapour(p, t);
            let x = Water.condensation(ThermodynamicTemperature::new::<kelvin>(t))?;
            let (j, _) = nucleation(
                &Supersaturated {
                    p,
                    t,
                    rho: v.rho,
                    gamma: v.cp / v.cv,
                    p_sat: x.vapour_pressure.get::<pascal>(),
                    sigma: x.surface_tension.get::<newton_per_meter>(),
                    rho_liquid: x.liquid_density.get::<kilogram_per_cubic_meter>(),
                    latent: x.latent_heat.get::<joule_per_kilogram>(),
                    molecule: Water.molecular_mass().get::<kilogram>(),
                },
                1.0,
                true,
            );
            anyhow::ensure!(
                (j.log10() - log_j).abs() < 0.3,
                "{subcooling} K: log J {} vs {log_j}",
                j.log10()
            );
        }
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        /// The rate rises steeply and monotonically with supersaturation.
        #[test]
        fn test_nucleation_rate_rises_with_supersaturation(s in 1.01f64..20.0, ds in 0.001f64..1.0) {
            let (j1, r1) = nucleation(&steam(3_537.0 * s), 1.0, true);
            let (j2, r2) = nucleation(&steam(3_537.0 * (s + ds)), 1.0, true);
            prop_assert!(j2 >= j1 && j1 >= 0.0, "{} then {}", j1, j2);
            prop_assert!(r2 < r1);
        }
    }
}
