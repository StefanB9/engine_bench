//! Mean-field kinetic nucleation theory (Kalikmanov, J. Chem. Phys. 124,
//! 124505, 2006): the kinetic approach to nucleation with a cluster
//! distribution from the revised Fisher droplet model, valid down to
//! clusters of a few molecules. Large clusters carry the planar surface
//! tension (the classical limit); small ones a microscopic surface tension
//! fixed by the saturated vapour's compressibility. Classical theory
//! underestimates the rates of argon and nitrogen by 10 to 25 orders of
//! magnitude where this theory is within one to three (Figs. 3 and 4).

use std::f64::consts::PI;

use super::nucleation::{BOLTZMANN, Supersaturated};

/// Clusters beyond the critical one are summed until their terms have
/// fallen by e⁻⁴⁶ ≈ 1e-20.
const TAIL: f64 = 46.0;
/// No cluster sum runs past this size: nucleation barriers that far out
/// make the rate nil.
const MAX_CLUSTER: u32 = 200_000;

/// What the theory needs beyond the classical rate's inputs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ClusterData {
    /// Excess compressibility Z − 1 of the saturated vapour (negative).
    pub(crate) z_excess: f64,
    /// Effective molecular (hard-sphere) diameter, m.
    pub(crate) diameter: f64,
}

/// The coordination number of the liquid from its packing fraction η
/// (Cahoon, Can. J. Phys. 82, 291, 2004; the paper's Sec. II G).
fn coordination(packing: f64) -> f64 {
    5.5116 * packing * packing + 6.1383 * packing + 1.275
}

/// Surface molecules of an n cluster with `n1` molecules in its first
/// shell: all of them up to n1, else n − x³ with x the positive root of
/// x³ + 3ωx² + 3ωλx − (n − ωλ²) = 0 (Eq. 28).
fn surface(n: f64, n1: f64, omega: f64, lambda: f64) -> f64 {
    if n <= n1 {
        return n;
    }
    let c = n - omega * lambda * lambda;
    let mut x = c.max(1.0).cbrt();
    for _ in 0..50 {
        let f = x.mul_add(x.mul_add(x + 3.0 * omega, 3.0 * omega * lambda), -c);
        let df = 3.0 * x * x + 6.0 * omega * x + 3.0 * omega * lambda;
        let step = f / df;
        x -= step;
        if step.abs() <= 1e-14 * x {
            break;
        }
    }
    n - x * x * x
}

/// Nucleation rate, 1/(m³·s), and the radius, m, of the cluster that
/// dominates it (the minimum of H(n), Eq. 34); zero rate at or below
/// saturation.
#[expect(
    clippy::many_single_char_names,
    reason = "the paper's symbols: n, h, x, s, a"
)]
pub(crate) fn mean_field_kinetic(v: &Supersaturated, c: &ClusterData) -> (f64, f64) {
    let ln_s = (v.p / v.p_sat).ln();
    if ln_s.is_nan() || ln_s <= 0.0 {
        return (0.0, f64::INFINITY);
    }
    let kt = BOLTZMANN * v.t;
    let liquid = v.rho_liquid / v.molecule;
    let r_l = (3.0 / (4.0 * PI * liquid)).cbrt();
    let s1 = 4.0 * PI * r_l * r_l;
    let theta_planar = v.sigma * s1 / kt;
    let theta_micro = -(-c.z_excess).ln();
    let omega = theta_planar / (3.0 * theta_micro);
    let n1 = coordination(PI / 6.0 * liquid * c.diameter.powi(3));
    let lambda = (n1 / omega - 0.75).sqrt() - 1.5;
    let vapour_sat = v.p_sat / ((1.0 + c.z_excess) * kt);
    let f1 = v.p_sat * s1 / (2.0 * PI * v.molecule * kt).sqrt();
    let a = vapour_sat * f1 * ln_s.exp();
    // Σ e^(−H(n)), scaled by e^(H_min) as it goes so nothing overflows.
    let h_of = |n: f64| {
        2.0 / 3.0 * n.ln() + n * ln_s - theta_micro * (surface(n, n1, omega, lambda) - 1.0)
    };
    let (mut h_min, mut n_min, mut sum) = (f64::INFINITY, 1.0, 0.0);
    for k in 1..=MAX_CLUSTER {
        let n = f64::from(k);
        let h = h_of(n);
        if h < h_min {
            sum = sum * (h - h_min).exp() + 1.0;
            h_min = h;
            n_min = n;
        } else {
            sum += (h_min - h).exp();
            if h - h_min > TAIL {
                break;
            }
        }
    }
    let rate = a * h_min.exp() / sum;
    let radius = (3.0 * n_min * v.molecule / (4.0 * PI * v.rho_liquid)).cbrt();
    (if rate.is_finite() { rate } else { 0.0 }, radius)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Steam at 50 kPa, `subcooling` K below saturation.
    fn steam(subcooling: f64) -> anyhow::Result<(Supersaturated, ClusterData)> {
        use crate::{
            flow::condensing::{Condensable, kernel::VapourKernel},
            steam::Water,
            units::{
                ThermodynamicTemperature, available_energy::joule_per_kilogram, mass::kilogram,
                mass_density::kilogram_per_cubic_meter, pressure::pascal,
                surface_tension::newton_per_meter, thermodynamic_temperature::kelvin,
            },
        };
        let (p, t) = (5e4, 354.47 - subcooling);
        let x = Water.condensation(ThermodynamicTemperature::new::<kelvin>(t))?;
        let v = Water.vapour(p, t);
        let p_sat = x.vapour_pressure.get::<pascal>();
        Ok((
            Supersaturated {
                p,
                t,
                rho: v.rho,
                gamma: v.cp / v.cv,
                p_sat,
                sigma: x.surface_tension.get::<newton_per_meter>(),
                rho_liquid: x.liquid_density.get::<kilogram_per_cubic_meter>(),
                latent: x.latent_heat.get::<joule_per_kilogram>(),
                molecule: Water.molecular_mass().get::<kilogram>(),
            },
            ClusterData {
                z_excess: Water.saturated_excess(t, p_sat),
                diameter: Water.molecular_diameter(),
            },
        ))
    }

    #[test]
    fn test_mknt_tends_to_classical_theory_for_large_critical_clusters() -> anyhow::Result<()> {
        // For large critical clusters the planar surface tension rules (the
        // paper's "macroscopic regime"), up to a Tolman-like correction that
        // fades as n^(-1/3) (its Sec. III C): the barriers' relative
        // difference shrinks as the subcooling, and so S, falls.
        let relative = |subcooling: f64| -> anyhow::Result<(f64, f64, f64)> {
            let (v, c) = steam(subcooling)?;
            let (mknt, radius) = mean_field_kinetic(&v, &c);
            let (classical, critical) = super::super::nucleation::nucleation(&v, 1.0, false);
            Ok(((mknt.ln() / classical.ln() - 1.0).abs(), radius, critical))
        };
        let (small, r, r_classical) = relative(6.0)?;
        let (larger, _, _) = relative(15.0)?;
        anyhow::ensure!(
            small < larger && small < 0.15,
            "{small} at 6 K, {larger} at 15 K"
        );
        anyhow::ensure!(
            (r / r_classical - 1.0).abs() < 0.1,
            "r {r} vs {r_classical}"
        );
        let (v, c) = steam(30.0)?;
        let (none, _) = mean_field_kinetic(&Supersaturated { p: v.p_sat, ..v }, &c);
        anyhow::ensure!(none.abs() < f64::MIN_POSITIVE, "saturated");
        Ok(())
    }

    #[test]
    fn test_mknt_surface_follows_the_cluster_structure() {
        // Up to the first shell every molecule is at the surface; then the
        // core grows and, far out, the surface scales as 3ω n^(2/3).
        let (n1, omega) = (6.0_f64, 0.4_f64);
        let lambda = (n1 / omega - 0.75).sqrt() - 1.5;
        assert!((surface(5.0, n1, omega, lambda) - 5.0).abs() < 1e-12);
        // One core molecule at n = 1 + N1 (x = 1).
        assert!((surface(n1 + 1.0, n1, omega, lambda) - n1).abs() < 1e-9);
        let n = 1e9;
        let ratio = surface(n, n1, omega, lambda) / (3.0 * omega * n.powf(2.0 / 3.0));
        assert!((ratio - 1.0).abs() < 0.01, "{ratio}");
    }
}
