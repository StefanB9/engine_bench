//! Forward equations of regions 1, 2 and 5: the dimensionless Gibbs free energy
//! γ(π, τ) and its derivatives, and the properties that follow from them.
//!
//! Infallible kernels on raw SI values; callers have validated (p, T) and
//! picked the region.

use super::{
    R,
    coefficients::{
        Term,
        gibbs::{REGION1, REGION2_IDEAL, REGION2_RESIDUAL, REGION5_IDEAL, REGION5_RESIDUAL},
    },
};

/// Per-mass properties at a (p, T) state, raw SI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct RawProperties {
    /// Specific volume, m³/kg.
    pub(crate) v: f64,
    /// Specific enthalpy, J/kg.
    pub(crate) h: f64,
    /// Specific internal energy, J/kg.
    pub(crate) u: f64,
    /// Specific entropy, J/(kg·K).
    pub(crate) s: f64,
    /// Isobaric heat capacity, J/(kg·K).
    pub(crate) cp: f64,
    /// Isochoric heat capacity, J/(kg·K).
    pub(crate) cv: f64,
    /// Speed of sound, m/s.
    pub(crate) w: f64,
}

/// Region 1 (compressed liquid) properties.
pub(crate) fn region1(p: f64, t: f64) -> RawProperties {
    const P_STAR: f64 = 16.53e6;
    const T_STAR: f64 = 1386.0;
    let (pi, tau) = (p / P_STAR, T_STAR / t);
    // γ = Σ n (7.1 − π)^I (τ − 1.222)^J: d/dπ = −d/d(7.1 − π).
    let mut g = residual(&REGION1, 7.1 - pi, tau - 1.222);
    g.g_p = -g.g_p;
    g.g_pt = -g.g_pt;
    g.properties(pi, tau, P_STAR, t)
}

/// Region 2 (vapour) properties.
pub(crate) fn region2(p: f64, t: f64) -> RawProperties {
    const P_STAR: f64 = 1e6;
    const T_STAR: f64 = 540.0;
    let (pi, tau) = (p / P_STAR, T_STAR / t);
    let g = ideal(&REGION2_IDEAL, pi, tau).plus(&residual(&REGION2_RESIDUAL, pi, tau - 0.5));
    g.properties(pi, tau, P_STAR, t)
}

/// Region 5 (high-temperature vapour) properties.
pub(crate) fn region5(p: f64, t: f64) -> RawProperties {
    const P_STAR: f64 = 1e6;
    const T_STAR: f64 = 1000.0;
    let (pi, tau) = (p / P_STAR, T_STAR / t);
    let g = ideal(&REGION5_IDEAL, pi, tau).plus(&residual(&REGION5_RESIDUAL, pi, tau));
    g.properties(pi, tau, P_STAR, t)
}

/// Dimensionless Gibbs free energy γ = g / (RT) and its partial derivatives in
/// reduced pressure π and inverse reduced temperature τ.
#[derive(Debug, Clone, Copy, Default)]
struct Gibbs {
    g: f64,
    g_p: f64,
    g_pp: f64,
    g_t: f64,
    g_tt: f64,
    g_pt: f64,
}

impl Gibbs {
    fn plus(&self, other: &Self) -> Self {
        Self {
            g: self.g + other.g,
            g_p: self.g_p + other.g_p,
            g_pp: self.g_pp + other.g_pp,
            g_t: self.g_t + other.g_t,
            g_tt: self.g_tt + other.g_tt,
            g_pt: self.g_pt + other.g_pt,
        }
    }

    /// Properties from the full γ at π = p / `p_star`, τ = T* / `t` (IAPWS
    /// R7-97(2012) Table 3; the relations hold for regions 2 and 5 too when γ
    /// includes the ideal-gas part).
    fn properties(&self, pi: f64, tau: f64, p_star: f64, t: f64) -> RawProperties {
        let rt = R * t;
        let tau_g_t = tau * self.g_t;
        let tau2_g_tt = tau * tau * self.g_tt;
        let x = self.g_p - tau * self.g_pt;
        RawProperties {
            v: rt * self.g_p / p_star,
            h: rt * tau_g_t,
            u: rt * (tau_g_t - pi * self.g_p),
            s: R * (tau_g_t - self.g),
            cp: -R * tau2_g_tt,
            cv: R * (x * x / self.g_pp - tau2_g_tt),
            w: (rt * self.g_p * self.g_p / (x * x / tau2_g_tt - self.g_pp)).sqrt(),
        }
    }
}

/// Σ n a^I b^J and its derivatives in a and b. Each term is evaluated once and
/// the derivatives scaled by I/a, J/b afterwards, so a and b must be non-zero,
/// which they are everywhere inside IF97's range.
#[expect(clippy::many_single_char_names, reason = "symbols of IAPWS R7-97")]
fn residual(terms: &[Term<i32>], a: f64, b: f64) -> Gibbs {
    let mut sum = Gibbs::default();
    for term in terms {
        let (i, j) = (f64::from(term.i), f64::from(term.j));
        let v = term.n * a.powi(term.i) * b.powi(term.j);
        sum.g += v;
        sum.g_p += i * v;
        sum.g_pp += i * (i - 1.0) * v;
        sum.g_t += j * v;
        sum.g_tt += j * (j - 1.0) * v;
        sum.g_pt += i * j * v;
    }
    Gibbs {
        g: sum.g,
        g_p: sum.g_p / a,
        g_pp: sum.g_pp / (a * a),
        g_t: sum.g_t / b,
        g_tt: sum.g_tt / (b * b),
        g_pt: sum.g_pt / (a * b),
    }
}

/// Ideal-gas part γ° = ln π + Σ n τ^J of regions 2 and 5.
fn ideal(terms: &[(i32, f64)], pi: f64, tau: f64) -> Gibbs {
    let (mut g, mut g_t, mut g_tt) = (pi.ln(), 0.0, 0.0);
    for &(j, n) in terms {
        let jf = f64::from(j);
        let v = n * tau.powi(j);
        g += v;
        g_t += jf * v;
        g_tt += jf * (jf - 1.0) * v;
    }
    Gibbs {
        g,
        g_p: 1.0 / pi,
        g_pp: -1.0 / (pi * pi),
        g_t: g_t / tau,
        g_tt: g_tt / (tau * tau),
        g_pt: 0.0,
    }
}

#[cfg(test)]
#[expect(
    clippy::unreadable_literal,
    reason = "verification values exactly as IAPWS R7-97(2012) prints them"
)]
mod tests {
    use super::*;

    /// IAPWS R7-97(2012) verification rows: T K, p MPa, v m³/kg, h kJ/kg,
    /// u kJ/kg, s kJ/(kg·K), cp kJ/(kg·K), w m/s.
    type Row = [f64; 8];

    /// Table 5 (region 1).
    const TABLE5: [Row; 3] = [
        [
            300.0,
            3.0,
            0.100215168e-2,
            115.331273,
            112.324818,
            0.392294792,
            4.17301218,
            1507.73921,
        ],
        [
            300.0,
            80.0,
            0.971180894e-3,
            184.142828,
            106.448356,
            0.368563852,
            4.01008987,
            1634.69054,
        ],
        [
            500.0,
            3.0,
            0.120241800e-2,
            975.542239,
            971.934985,
            2.58041912,
            4.65580682,
            1240.71337,
        ],
    ];

    /// Table 15 (region 2).
    const TABLE15: [Row; 3] = [
        [
            300.0, 0.0035, 39.4913866, 2549.91145, 2411.69160, 8.52238967, 1.91300162, 427.920172,
        ],
        [
            700.0, 0.0035, 92.3015898, 3335.68375, 3012.62819, 10.1749996, 2.08141274, 644.289068,
        ],
        [
            700.0,
            30.0,
            0.542946619e-2,
            2631.49474,
            2468.61076,
            5.17540298,
            10.3505092,
            480.386523,
        ],
    ];

    /// Table 42 (region 5).
    const TABLE42: [Row; 3] = [
        [
            1500.0, 0.5, 1.38455090, 5219.76855, 4527.49310, 9.65408875, 2.61609445, 917.068690,
        ],
        [
            1500.0,
            30.0,
            0.230761299e-1,
            5167.23514,
            4474.95124,
            7.72970133,
            2.72724317,
            928.548002,
        ],
        [
            2000.0,
            30.0,
            0.311385219e-1,
            6571.22604,
            5637.07038,
            8.53640523,
            2.88569882,
            1067.36948,
        ],
    ];

    /// The tables give 9 significant digits.
    const TOLERANCE: f64 = 5e-8;

    fn check(table: &[Row], kernel: fn(f64, f64) -> RawProperties) -> anyhow::Result<()> {
        for &[t, p_mpa, v, h, u, s, cp, w] in table {
            let got = kernel(p_mpa * 1e6, t);
            let pairs = [
                ("v", got.v, v),
                ("h", got.h, h * 1e3),
                ("u", got.u, u * 1e3),
                ("s", got.s, s * 1e3),
                ("cp", got.cp, cp * 1e3),
                ("w", got.w, w),
            ];
            for (name, got, want) in pairs {
                let error = (got / want - 1.0).abs();
                anyhow::ensure!(
                    error < TOLERANCE,
                    "{name} at {t} K, {p_mpa} MPa: {got} vs {want} (rel {error:e})"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn test_region1_matches_iapws_table_5() -> anyhow::Result<()> {
        check(&TABLE5, region1)
    }

    #[test]
    fn test_region2_matches_iapws_table_15() -> anyhow::Result<()> {
        check(&TABLE15, region2)
    }

    #[test]
    fn test_region5_matches_iapws_table_42() -> anyhow::Result<()> {
        check(&TABLE42, region5)
    }

    #[test]
    fn test_gibbs_regions_satisfy_u_equals_h_minus_pv() -> anyhow::Result<()> {
        for (kernel, p, t) in [
            (region1 as fn(f64, f64) -> RawProperties, 5e6, 400.0),
            (region2, 1e5, 500.0),
            (region5, 1e7, 1800.0),
        ] {
            let x = kernel(p, t);
            let error = (x.u - (x.h - p * x.v)).abs() / x.h.abs();
            anyhow::ensure!(error < 1e-13, "u ≠ h − pv at {p} Pa, {t} K: rel {error:e}");
        }
        Ok(())
    }

    #[test]
    fn test_region2_approaches_an_ideal_gas_at_low_pressure() -> anyhow::Result<()> {
        let (p, t) = (1.0, 500.0);
        let x = region2(p, t);
        let z = p * x.v / (R * t);
        anyhow::ensure!((z - 1.0).abs() < 1e-7, "Z = {z}");
        anyhow::ensure!(
            (x.cp - x.cv - R).abs() < 1e-6 * R,
            "cp − cv = {}",
            x.cp - x.cv
        );
        Ok(())
    }
}
