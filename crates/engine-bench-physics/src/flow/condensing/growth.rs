//! Droplet growth and evaporation: Young's (1982) law, bridging the free-
//! molecular and continuum regimes, with the droplet temperature from
//! Gyarmathy's quasi-steady capillary relation (the Wet Steam Modeling
//! Project's equations (8) to (11)).

use super::CondensationModel;

/// Mean free path factor: l = 1.5·μ·√(R·T) / p (Gyarmathy).
const MEAN_FREE_PATH: f64 = 1.5;

/// The vapour around the droplets, raw SI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Surroundings {
    /// Vapour pressure, Pa, and temperature, K.
    pub(crate) p: f64,
    pub(crate) t: f64,
    /// Saturation temperature at `p`, K.
    pub(crate) t_sat: f64,
    /// Surface tension, N/m, liquid density, kg/m³, and latent heat, J/kg, at
    /// `t_sat`.
    pub(crate) sigma: f64,
    pub(crate) rho_liquid: f64,
    pub(crate) latent: f64,
    /// Vapour viscosity, Pa·s, conductivity, W/(m·K), cp, J/(kg·K), cp/cv
    /// and gas constant, J/(kg·K).
    pub(crate) viscosity: f64,
    pub(crate) conductivity: f64,
    pub(crate) cp: f64,
    pub(crate) gamma: f64,
    pub(crate) gas_constant: f64,
}

/// dr/dt of a droplet of radius r, with everything that does not depend on r
/// computed once per state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Growth {
    /// λ / (ρₗ·L), m²/(s·K).
    scale: f64,
    /// Subcooling Tₛ − T, K.
    subcooling: f64,
    /// Capillary term 2σTₛ / (ρₗ·L), K·m: Tₗ − T = ΔT − this / r.
    capillary: f64,
    /// Mean free path, m.
    free_path: f64,
    beta: f64,
    /// 3.78·(1 − ν) / Pr.
    molecular: f64,
}

impl Growth {
    pub(crate) fn new(v: &Surroundings, model: &CondensationModel) -> Self {
        let q = model.condensation_coefficient;
        let latent_ratio = v.gas_constant * v.t_sat / v.latent;
        let nu = latent_ratio
            * (model.alpha
                - 0.5
                - (2.0 - q) / (2.0 * q) * (v.gamma + 1.0) / (2.0 * v.gamma) * v.cp * v.t_sat
                    / v.latent);
        let prandtl = v.viscosity * v.cp / v.conductivity;
        Self {
            scale: model.growth_factor * v.conductivity / (v.rho_liquid * v.latent),
            subcooling: v.t_sat - v.t,
            capillary: 2.0 * v.sigma * v.t_sat / (v.rho_liquid * v.latent),
            free_path: MEAN_FREE_PATH * v.viscosity * (v.gas_constant * v.t).sqrt() / v.p,
            beta: model.beta,
            molecular: 3.78 * (1.0 - nu) / prandtl,
        }
    }

    /// dr/dt, m/s, of a droplet of radius `r` m: negative while it
    /// evaporates.
    pub(crate) fn rate(&self, r: f64) -> f64 {
        let l = self.free_path;
        let resistance = r / (1.0 + self.beta * l / r) + 0.5 * self.molecular * l;
        self.scale * (self.subcooling - self.capillary / r) / resistance
    }

    /// Subcooling Tₛ − T, K.
    #[cfg(test)]
    pub(crate) fn subcooling(&self) -> f64 {
        self.subcooling
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Steam at 10 kPa, 20 K below saturation (319 K).
    fn steam(p: f64) -> Surroundings {
        Surroundings {
            p,
            t: 299.0,
            t_sat: 319.0,
            sigma: 0.0685,
            rho_liquid: 989.8,
            latent: 2.392e6,
            viscosity: 9.8e-6,
            conductivity: 0.019,
            cp: 1_950.0,
            gamma: 1.32,
            gas_constant: 461.5,
        }
    }

    #[test]
    fn test_growth_stops_at_the_critical_radius_and_reverses_below() -> anyhow::Result<()> {
        let g = Growth::new(&steam(1e4), &CondensationModel::default());
        let r_crit = g.capillary / g.subcooling;
        anyhow::ensure!(g.rate(r_crit).abs() < 1e-15, "{}", g.rate(r_crit));
        anyhow::ensure!(
            g.rate(2.0 * r_crit) > 0.0 && g.rate(0.5 * r_crit) < 0.0,
            "sign"
        );
        // Superheated vapour evaporates every droplet.
        let mut hot = steam(1e4);
        hot.t = 330.0;
        let g = Growth::new(&hot, &CondensationModel::default());
        anyhow::ensure!(g.rate(1e-6) < 0.0, "evaporation");
        Ok(())
    }

    #[test]
    fn test_growth_tends_to_the_continuum_law_for_large_droplets() -> anyhow::Result<()> {
        // r ≫ l: dr/dt = λ (T_l − T) / (ρ_l L r), Maxwell's heat conduction.
        let v = steam(1e4);
        let g = Growth::new(&v, &CondensationModel::default());
        let r = 1e-3;
        let maxwell = v.conductivity * (20.0 - g.capillary / r) / (v.rho_liquid * v.latent * r);
        let rate = g.rate(r);
        anyhow::ensure!((rate / maxwell - 1.0).abs() < 1e-3, "{rate} vs {maxwell}");
        Ok(())
    }

    #[test]
    fn test_growth_rate_is_size_independent_in_the_free_molecular_limit() -> anyhow::Result<()> {
        // r ≪ l = 0.55 mm (but ≫ r*): dr/dt tends to a constant, Hill's regime.
        let g = Growth::new(&steam(10.0), &CondensationModel::default());
        let (a, b) = (g.rate(1e-6), g.rate(2e-6));
        anyhow::ensure!((a / b - 1.0).abs() < 0.01, "{a} vs {b}");
        // Larger α means faster growth (ν rises, 1 − ν falls).
        let fast = Growth::new(&steam(10.0), &CondensationModel::young(11.0, 0.0)?);
        anyhow::ensure!(fast.rate(1e-6) > a, "α");
        Ok(())
    }

    #[test]
    fn test_growth_free_path_follows_gyarmathy() -> anyhow::Result<()> {
        let v = steam(1e4);
        let g = Growth::new(&v, &CondensationModel::default());
        let l = 1.5 * 9.8e-6 * (461.5f64 * 299.0).sqrt() / 1e4;
        anyhow::ensure!((g.free_path / l - 1.0).abs() < 1e-12, "{}", g.free_path);
        anyhow::ensure!((g.subcooling() - 20.0).abs() < 1e-12, "ΔT");
        Ok(())
    }
}
