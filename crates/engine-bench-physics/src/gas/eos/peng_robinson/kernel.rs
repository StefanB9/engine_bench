//! Peng–Robinson on a molar basis: pressure and its derivatives, and the
//! departures of u, s and `c_v` from the ideal gas at the same (T, v).
//!
//! van der Waals one-fluid mixing with `k_ij` = 0, so a = (Σ xᵢ √aᵢ)², and
//! √aᵢ(T) = √`a_cᵢ` (1 + κᵢ (1 − √(`T/T_cᵢ`))) is linear in √T. Infallible
//! kernels on raw SI values; callers validate.

use std::f64::consts::{PI, SQRT_2};

use super::critical::Critical;
use crate::thermo::MOLAR_GAS_CONSTANT as R;

/// Per-species constants.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Parameters<const N: usize> {
    /// √`a_c`, (Pa·m⁶/mol²)^½.
    sqrt_ac: [f64; N],
    kappa: [f64; N],
    /// Critical temperature, K.
    tc: [f64; N],
    /// Co-volume, m³/mol.
    b: [f64; N],
}

/// PR state at (T, v) for one composition, molar.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Molar {
    /// Pa.
    pub(crate) p: f64,
    /// (∂p/∂T) at constant v, Pa/K.
    pub(crate) dp_dt: f64,
    /// (∂p/∂v) at constant T, Pa·mol/m³.
    pub(crate) dp_dv: f64,
    /// u − `u_ig` at the same T, J/mol.
    pub(crate) u_dep: f64,
    /// s − `s_ig` at the same (T, v), J/(mol·K).
    pub(crate) s_dep: f64,
    /// `c_v` − `c_v,ig`, J/(mol·K).
    pub(crate) cv_dep: f64,
}

/// `Ω_a` and `Ω_b` of Peng and Robinson (1976): the critical-point conditions
/// of the cubic, with `Z_c` = 0.307 401.
const OMEGA_A: f64 = 0.457_235_528_921_382_7;
const OMEGA_B: f64 = 0.077_796_073_903_888_45;

impl<const N: usize> Parameters<N> {
    pub(crate) fn new(critical: &[Critical; N]) -> Self {
        let mut pr = Self {
            sqrt_ac: [0.0; N],
            kappa: [0.0; N],
            tc: [0.0; N],
            b: [0.0; N],
        };
        for (i, c) in critical.iter().enumerate() {
            let w = c.acentric_factor;
            pr.sqrt_ac[i] = (OMEGA_A / c.pressure).sqrt() * R * c.temperature;
            pr.kappa[i] = (-0.269_92 * w + 1.542_26).mul_add(w, 0.374_64);
            pr.tc[i] = c.temperature;
            pr.b[i] = OMEGA_B * R * c.temperature / c.pressure;
        }
        pr
    }

    /// Mixture co-volume, m³/mol.
    pub(crate) fn covolume(&self, x: &[f64; N]) -> f64 {
        x.iter().zip(&self.b).map(|(x, b)| x * b).sum()
    }

    /// a(T) and its first two temperature derivatives. With `k_ij` = 0,
    /// a = Q², Q = Σ xᵢ √`a_cᵢ` (1 + κᵢ − κᵢ √(`T/T_cᵢ`)).
    #[expect(
        clippy::many_single_char_names,
        reason = "symbols of the Peng-Robinson equations"
    )]
    fn attraction(&self, x: &[f64; N], t: f64) -> (f64, f64, f64) {
        let (mut q, mut dq, mut ddq) = (0.0, 0.0, 0.0);
        for i in 0..N {
            let w = x[i] * self.sqrt_ac[i];
            let k = self.kappa[i];
            let root = (t * self.tc[i]).sqrt();
            q += w * (1.0 + k - k * t / root);
            dq -= w * k / (2.0 * root);
            ddq += w * k / (4.0 * t * root);
        }
        (q * q, 2.0 * q * dq, 2.0 * dq.mul_add(dq, q * ddq))
    }

    /// State at temperature `t` K and molar volume `v` m³/mol (> b).
    #[expect(
        clippy::many_single_char_names,
        reason = "symbols of the Peng-Robinson equations"
    )]
    pub(crate) fn evaluate(&self, x: &[f64; N], t: f64, v: f64) -> Molar {
        let (a, da, dda) = self.attraction(x, t);
        let b = self.covolume(x);
        let d = v.mul_add(v, b * (2.0 * v - b));
        // ∫_∞^v dv / (v² + 2bv − b²).
        let l = ((v + (1.0 - SQRT_2) * b) / (v + (1.0 + SQRT_2) * b)).ln() / (2.0 * SQRT_2 * b);
        let free = v - b;
        Molar {
            p: R * t / free - a / d,
            dp_dt: R / free - da / d,
            dp_dv: -R * t / (free * free) + 2.0 * a * (v + b) / (d * d),
            u_dep: (a - t * da) * l,
            s_dep: R * (free / v).ln() - da * l,
            cv_dep: -t * dda * l,
        }
    }

    /// Molar volume, m³/mol, of the stable phase at `t` K and `p` Pa: of the
    /// cubic's roots, the one with the lowest Gibbs energy.
    #[expect(
        clippy::many_single_char_names,
        reason = "symbols of the Peng-Robinson equations"
    )]
    pub(crate) fn molar_volume(&self, x: &[f64; N], t: f64, p: f64) -> f64 {
        let (a, _, _) = self.attraction(x, t);
        let rt = R * t;
        let big_a = a * p / (rt * rt);
        let big_b = self.covolume(x) * p / rt;
        let coefficients = [
            -(1.0 - big_b),
            big_a - big_b * (3.0 * big_b + 2.0),
            -(big_a * big_b - big_b * big_b * (1.0 + big_b)),
        ];
        let gibbs = |z: f64| {
            z - 1.0
                - (z - big_b).ln()
                - big_a / (2.0 * SQRT_2 * big_b)
                    * ((z + (1.0 + SQRT_2) * big_b) / (z + (1.0 - SQRT_2) * big_b)).ln()
        };
        let z = cubic_roots(coefficients)
            .into_iter()
            .flatten()
            .filter(|z| *z > big_b)
            .min_by(|a, b| gibbs(*a).total_cmp(&gibbs(*b)))
            .unwrap_or(f64::NAN);
        z * rt / p
    }
}

/// Real roots of z³ + c₂z² + c₁z + c₀ (Cardano, polished by Newton).
fn cubic_roots([c2, c1, c0]: [f64; 3]) -> [Option<f64>; 3] {
    let shift = c2 / 3.0;
    let p = c1 - c2 * shift;
    let q = (2.0 * shift * shift - c1) * shift + c0;
    let discriminant = (q / 2.0).powi(2) + (p / 3.0).powi(3);
    let polish = |mut z: f64| {
        for _ in 0..3 {
            let f = ((z + c2) * z + c1) * z + c0;
            let df = (3.0 * z + 2.0 * c2) * z + c1;
            if df != 0.0 {
                z -= f / df;
            }
        }
        z
    };
    if discriminant > 0.0 {
        let root = discriminant.sqrt();
        let y = (-q / 2.0 + root).cbrt() + (-q / 2.0 - root).cbrt();
        [Some(polish(y - shift)), None, None]
    } else {
        let r = (-p / 3.0).sqrt();
        let phi = (q / (-2.0 * r * r * r)).clamp(-1.0, 1.0).acos();
        let y = |k: f64| 2.0 * r * ((phi - 2.0 * PI * k) / 3.0).cos();
        [
            Some(polish(y(0.0) - shift)),
            Some(polish(y(1.0) - shift)),
            Some(polish(y(2.0) - shift)),
        ]
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::gas::eos::peng_robinson::critical::critical_constants;

    fn pure(name: &str) -> anyhow::Result<(Parameters<1>, Critical)> {
        let c = critical_constants(name).ok_or_else(|| anyhow::anyhow!("{name}"))?;
        Ok((Parameters::new(&[c]), c))
    }

    #[test]
    fn test_peng_robinson_critical_compressibility_is_0_3074() -> anyhow::Result<()> {
        for name in ["N2", "CO2", "He", "H2O"] {
            let (pr, c) = pure(name)?;
            let v = pr.molar_volume(&[1.0], c.temperature, c.pressure);
            let z = c.pressure * v / (R * c.temperature);
            anyhow::ensure!((z - 0.307_401).abs() < 1e-5, "{name}: Z_c = {z}");
        }
        Ok(())
    }

    #[test]
    fn test_peng_robinson_picks_the_stable_root() -> anyhow::Result<()> {
        // CO2 at 250 K: saturation near 17.9 bar. Above it the liquid root is
        // stable, below it the vapour root.
        let (pr, _) = pure("CO2")?;
        let liquid = pr.molar_volume(&[1.0], 250.0, 30e5);
        let vapour = pr.molar_volume(&[1.0], 250.0, 10e5);
        let z_l = 30e5 * liquid / (R * 250.0);
        let z_v = 10e5 * vapour / (R * 250.0);
        anyhow::ensure!(z_l < 0.15 && z_v > 0.8, "Z_l = {z_l}, Z_v = {z_v}");
        Ok(())
    }

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_peng_robinson_tends_to_the_ideal_gas_at_low_density(t in 100.0f64..1000.0, log_v in 1.0f64..3.0) {
            let (pr, _) = pure("N2").map_err(|e| fail(&e))?;
            let v = 10f64.powf(log_v); // 10–1000 m³/mol: below 1e-3 bar
            let m = pr.evaluate(&[1.0], t, v);
            prop_assert!((m.p * v / (R * t) - 1.0).abs() < 1e-4, "Z = {}", m.p * v / (R * t));
            prop_assert!(m.u_dep.abs() < 1e-3 * R * t && m.s_dep.abs() < 1e-3 * R, "{:?}", m);
        }

        #[test]
        fn test_peng_robinson_derivatives_are_consistent(
            t in 150.0f64..700.0,
            z_b in 1.3f64..30.0,
            y in 0.0f64..1.0,
        ) {
            // A CO2–N2 mixture at molar volume z_b · b.
            let c = [critical_constants("CO2"), critical_constants("N2")];
            let [Some(co2), Some(n2)] = c else { return Err(TestCaseError::fail("constants")) };
            let pr = Parameters::new(&[co2, n2]);
            let x = [y, 1.0 - y];
            let v = z_b * pr.covolume(&x);
            let m = pr.evaluate(&x, t, v);
            let (dt, dv) = (1e-4 * t, 1e-6 * v);
            let at = |t: f64, v: f64| pr.evaluate(&x, t, v);
            let dp_dt = (at(t + dt, v).p - at(t - dt, v).p) / (2.0 * dt);
            let dp_dv = (at(t, v + dv).p - at(t, v - dv).p) / (2.0 * dv);
            let cv = (at(t + dt, v).u_dep - at(t - dt, v).u_dep) / (2.0 * dt);
            // Maxwell: (∂s/∂v)_T = (∂p/∂T)_v, and s_ig(T, v) contributes R/v.
            let ds_dv = (at(t, v + dv).s_dep - at(t, v - dv).s_dep) / (2.0 * dv) + R / v;
            let close = |a: f64, b: f64, scale: f64| (a - b).abs() <= 1e-5 * scale.abs().max(a.abs());
            prop_assert!(close(m.dp_dt, dp_dt, R / v), "dp/dT {} vs {}", m.dp_dt, dp_dt);
            prop_assert!(close(m.dp_dv, dp_dv, R * t / (v * v)), "dp/dv {} vs {}", m.dp_dv, dp_dv);
            prop_assert!(close(m.cv_dep, cv, R), "cv {} vs {}", m.cv_dep, cv);
            prop_assert!(close(m.dp_dt, ds_dv, R / v), "Maxwell {} vs {}", m.dp_dt, ds_dv);
        }
    }
}
