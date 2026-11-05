//! Cold metastable vapour: below the triple point the reference equations'
//! residual parts are extrapolations, and their temperature curvature
//! diverges a few tens of kelvin down even where the vapour is dilute (cp
//! turns negative: argon below ~35 K, nitrogen ~28 K, CO₂ ~155 K). There the
//! vapour is the equation's ideal gas plus the second virial coefficient of
//! a Lennard-Jones fluid plus a constant, B = b₀·B*(T/ε) + B₀ with b₀ =
//! (2π/3)·N·σ³ (N Avogadro's number), its σ, ε and B₀ chosen so that B, dB/dT
//! and d²B/dT² equal the equation's own at the triple point (cv continuous
//! across it). B*(T*) is the exact series Σⱼ bⱼ·T*^(−(2j+1)/4), bⱼ =
//! −2^(j+½)/(4·j!)· Γ((2j − 1)/4) (Hirschfelder, Curtiss and Bird, Molecular
//! Theory of Gases and Liquids, eq. 3.6-8), finite at any temperature, where
//! corresponding states correlations such as Tsonopoulos's diverge (their
//! (T꜀/T)⁸ term makes B hundreds of times too negative at T ≈ 0.15·T꜀). At
//! these densities (|Z − 1| ≲ 0.1) the second virial term carries the
//! non-ideality.

use super::{
    Fluid,
    helmholtz::{Residual, residual},
};

/// Densities, as δ, up to which the cold vapour is virial: well past any
/// metastable vapour below the triple point, far short of the liquid.
pub(crate) const MAX_DELTA: f64 = 1.0;
/// The equation's higher-order part is held within this fraction of τ
/// beyond the triple point (T ≳ 0.9·Tₜ, where every equation behaves).
const HOLD: f64 = 0.1;
/// Terms of the B* series kept at most (it converges like 2ʲ/j!·T*^(−j/2)).
const SERIES_TERMS: usize = 200;
/// Γ(−1/4) and Γ(1/4): the series' Gamma values follow from Γ(x + 1) = xΓ(x).
const GAMMA_MINUS_QUARTER: f64 = -4.901_666_809_860_711;
const GAMMA_QUARTER: f64 = 3.625_609_908_221_908;

/// Σⱼ bⱼ·x^mⱼ, Σ mⱼ·bⱼ·x^mⱼ and Σ mⱼ(mⱼ − 1)·bⱼ·x^mⱼ with mⱼ = (2j + 1)/4,
/// x = 1/T*: B* and its first two τ-derivatives scaled by τ and τ².
fn series(x: f64) -> (f64, f64, f64) {
    let (mut sum, mut d1, mut d2) = (0.0, 0.0, 0.0);
    // Γ((2j − 1)/4) for j and j + 1, and 2^(j+½)/j!.
    let (mut gamma, mut gamma_next) = (GAMMA_MINUS_QUARTER, GAMMA_QUARTER);
    let mut weight = std::f64::consts::SQRT_2;
    for j in 0..SERIES_TERMS {
        #[expect(clippy::cast_precision_loss, reason = "a small term index")]
        let jf = j as f64;
        let m = (2.0 * jf + 1.0) / 4.0;
        let term = -weight / 4.0 * gamma * x.powf(m);
        sum += term;
        d1 += m * term;
        d2 += m * (m - 1.0) * term;
        if j > 4 && term.abs() <= 1e-16 * sum.abs() {
            break;
        }
        let x_j = (2.0 * jf - 1.0) / 4.0;
        (gamma, gamma_next) = (gamma_next, x_j * gamma);
        weight *= 2.0 / (jf + 1.0);
    }
    (sum, d1, d2)
}

/// The cold vapour's second virial coefficient as βᵥ(τ) = ρ꜀·B: βᵥ =
/// `scale`·B*(T*) + `offset`, T* = `t_star_tau`/τ (the offset, a constant
/// excluded volume, lets B match the equation's in curvature too).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Virial {
    /// τ at the triple point: colder than this the vapour is virial.
    pub(crate) tau_switch: f64,
    scale: f64,
    offset: f64,
    t_star_tau: f64,
}

impl Virial {
    /// Matched to `fluid`'s equation at its triple point in value, slope and
    /// curvature (so cv is continuous there); `None` if no Lennard-Jones B
    /// fits (helium: its B curves the other way).
    pub(crate) fn fitted(fluid: &Fluid) -> Option<Self> {
        const DILUTE: f64 = 1e-10;
        let tau = fluid.data.t_crit / fluid.data.t_triple;
        let r = residual(fluid, DILUTE, tau);
        let (beta, slope, curvature) = (r.d / DILUTE, r.dt / DILUTE, r.tt / DILUTE);
        // τ²β''/(τβ') depends on T* alone and falls monotonically with it:
        // bisect in ln T*.
        let target = curvature / slope;
        let f = |ln_t: f64| {
            let (_, d1, d2) = series((-ln_t).exp());
            d2 / d1 - target
        };
        let (mut lo, mut hi) = (0.05f64.ln(), 3.0f64.ln());
        if !(f(lo) > 0.0 && f(hi) < 0.0) {
            return None;
        }
        for _ in 0..200 {
            let mid = f64::midpoint(lo, hi);
            if f(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let t_star = f64::midpoint(lo, hi).exp();
        let (b_star, d1, _) = series(1.0 / t_star);
        let scale = slope / d1;
        Some(Self {
            tau_switch: tau,
            scale,
            offset: beta - scale * b_star,
            t_star_tau: t_star * tau,
        })
    }

    /// ρ꜀·B and its τ-derivatives scaled by τ and τ² at `tau`.
    fn beta(&self, tau: f64) -> (f64, f64, f64) {
        let (b, d1, d2) = series(tau / self.t_star_tau);
        (
            self.scale * b + self.offset,
            self.scale * d1,
            self.scale * d2,
        )
    }

    /// The scaled residual derivatives of the cold vapour: δ·βᵥ(τ) plus the
    /// equation's own higher-order density part (beyond the second virial),
    /// evaluated at τ* = τₛ + Δ·tanh((τ − τₛ)/Δ): equal to τ to second order
    /// at the triple point, levelling off within [`HOLD`] of it, where the
    /// equation still behaves. With βᵥ matched to the equation's B to second
    /// order, αʳ and its τ-derivatives are continuous there at every density,
    /// and every derivative is exact (chain rule), so the vapour stays
    /// thermodynamically consistent.
    #[expect(
        clippy::many_single_char_names,
        reason = "the symbols of the chain rule: b, m, x, e, l"
    )]
    pub(crate) fn residual(&self, fluid: &Fluid, delta: f64, tau: f64) -> Residual {
        const DILUTE: f64 = 1e-10;
        let (b, d1, d2) = self.beta(tau);
        let ts = self.tau_switch;
        let width = HOLD * ts;
        let x = (tau - ts) / width;
        let (th, sech2) = (x.tanh(), 1.0 / x.cosh().powi(2));
        let m = ts + width * th;
        let (m1, m2) = (sech2, -2.0 * sech2 * th / width);
        // The equation and its second virial part at τ*.
        let e = residual(fluid, delta, m);
        let l = residual(fluid, DILUTE, m);
        let (be, be1, be2) = (l.d / DILUTE, l.dt / DILUTE, l.tt / DILUTE);
        // The higher-order part h and its τ*-scaled derivatives.
        let (ha, hd, hdd) = (e.a - delta * be, e.d - delta * be, e.dd);
        let (ht, htt, hdt) = (e.t - delta * be1, e.tt - delta * be2, e.dt - delta * be1);
        // To τ: h_τ = h_m·m′, h_ττ = h_mm·m′² + h_m·m″.
        let (slope, curvature, cross) = (ht / m, htt / (m * m), hdt / m);
        Residual {
            a: delta * b + ha,
            d: delta * b + hd,
            dd: hdd,
            t: delta * d1 + tau * slope * m1,
            tt: delta * d2 + tau * tau * (curvature * m1 * m1 + slope * m2),
            dt: delta * d1 + tau * cross * m1,
        }
    }

    /// ρ꜀·B at `tau` (colder than the switch).
    pub(crate) fn second_virial(&self, tau: f64) -> f64 {
        self.beta(tau).0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid::{PureFluid, helmholtz::molar};

    /// B*(T*) = −3∫₀^∞ (exp(−4(r⁻¹² − r⁻⁶)/T*) − 1)·r² dr by Simpson's rule
    /// (r in σ, below 0.5 the integrand is −r²).
    fn integrated(t_star: f64) -> f64 {
        let integrand = |r: f64| {
            let r6 = r.powi(-6);
            ((-4.0 * (r6 * r6 - r6) / t_star).exp() - 1.0) * r * r
        };
        let (a, b, n) = (0.5, 40.0, 4_000_000);
        let h = (b - a) / f64::from(n);
        let mut sum = integrand(a) + integrand(b);
        for k in 1..n {
            sum += integrand(a + h * f64::from(k)) * if k % 2 == 1 { 4.0 } else { 2.0 };
        }
        // Inside 0.5σ the integrand is −r²; beyond 40σ, 4r⁻⁴/T*.
        let core = -0.5f64.powi(3) / 3.0;
        let tail = 4.0 / (3.0 * t_star * b.powi(3));
        -3.0 * (core + sum * h / 3.0 + tail)
    }

    #[test]
    fn test_lennard_jones_series_matches_the_integral() {
        for t_star in [0.15, 0.3, 0.6, 1.0, 2.0, 5.0] {
            let (series_b, _, _) = series(1.0 / t_star);
            let direct = integrated(t_star);
            assert!(
                (series_b / direct - 1.0).abs() < 1e-6,
                "B*({t_star}) = {series_b}, integral {direct}"
            );
        }
    }

    /// Molar density, mol/m³, of an ideal gas at `p` Pa and `t` K.
    fn dilute(p: f64, t: f64) -> f64 {
        p / (8.314_462_618 * t)
    }

    #[test]
    fn test_cold_vapour_stays_physical_down_to_a_fifth_of_the_triple_point() {
        // At 2 kPa (CO₂ at 20 kPa) from the triple point down to 0.2 of it,
        // wherever the vapour is trusted (|Z − 1| ≤ 0.1): cp above the ideal
        // gas's by at most twenty times |Z − 1| (B's curvature grows as the
        // vapour cools; the reference equations' own extrapolations turn cp
        // negative).
        for (fluid, p) in [
            (PureFluid::Nitrogen, 2e3),
            (PureFluid::Argon, 2e3),
            (PureFluid::CarbonDioxide, 2e4),
        ] {
            let f = fluid.data();
            for k in 0..=16 {
                let t = f.data.t_triple * (1.0 - 0.05 * f64::from(k));
                let ideal_cp = molar(f, 1e-12, t).cp;
                let x = molar(f, dilute(p, t), t);
                let z = x.p / (dilute(p, t) * f.data.gas_constant * t);
                // Vapour this cold is trusted only while nearly ideal.
                if z < 0.9 {
                    continue;
                }
                let excess = x.cp / ideal_cp - 1.0;
                assert!(
                    excess > 0.0 && excess < 20.0 * (1.0 - z) && z < 1.0,
                    "{fluid:?} at {t} K: cp {} (ideal {ideal_cp}), Z {z}",
                    x.cp
                );
            }
        }
    }

    #[test]
    fn test_cold_vapour_meets_the_equation_at_the_triple_point() {
        // Just below and just above the switch the properties agree.
        for fluid in [
            PureFluid::Nitrogen,
            PureFluid::Argon,
            PureFluid::CarbonDioxide,
        ] {
            let f = fluid.data();
            let t = f.data.t_triple;
            let rho = dilute(1e3, t);
            let (below, above) = (
                molar(f, rho, t * (1.0 - 1e-9)),
                molar(f, rho, t * (1.0 + 1e-9)),
            );
            for (name, a, b) in [
                ("p", below.p, above.p),
                ("h", below.h, above.h),
                ("s", below.s, above.s),
                ("cp", below.cp, above.cp),
                ("w", below.w, above.w),
            ] {
                assert!((a / b - 1.0).abs() < 1e-5, "{fluid:?} {name}: {a} vs {b}");
            }
        }
    }
}
