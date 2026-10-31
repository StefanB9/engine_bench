//! The reduced Helmholtz energy α = α⁰ + αʳ and its derivatives in scaled
//! form (δ^i τ^j ∂^(i+j)α/∂δ^i∂τ^j), and the molar properties at (ρ, T).
//!
//! Infallible kernels on raw SI values; callers have validated the state.

use super::{
    Fluid,
    coefficients::{FluidData, NonAnalytic},
};

/// Scaled residual derivatives: αʳ, `δαʳ_δ`, `δ²αʳ_δδ`, `ταʳ_τ`, `τ²αʳ_ττ`,
/// `δταʳ_δτ`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct Residual {
    pub(crate) a: f64,
    pub(crate) d: f64,
    pub(crate) dd: f64,
    pub(crate) t: f64,
    pub(crate) tt: f64,
    pub(crate) dt: f64,
}

/// Molar properties at (ρ, T).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Molar {
    /// Pa.
    pub(crate) p: f64,
    /// (∂p/∂ρ) at constant T, Pa·m³/mol.
    pub(crate) dp_drho: f64,
    /// (∂p/∂T) at constant ρ, Pa/K.
    pub(crate) dp_dt: f64,
    /// J/mol.
    pub(crate) u: f64,
    /// J/mol.
    pub(crate) h: f64,
    /// J/(mol·K).
    pub(crate) s: f64,
    /// J/(mol·K).
    pub(crate) cv: f64,
    /// J/(mol·K).
    pub(crate) cp: f64,
    /// Speed of sound, m/s.
    pub(crate) w: f64,
}

/// Scaled residual derivatives of `fluid` at reduced density `delta` and
/// inverse reduced temperature `tau`.
pub(crate) fn residual(fluid: &Fluid, delta: f64, tau: f64) -> Residual {
    let data = fluid.data;
    let mut r = Residual::default();
    let ln_tau = tau.ln();
    for term in data.power {
        let (d, l) = (f64::from(term.d), f64::from(term.l));
        let delta_l = if term.l > 0 { delta.powi(term.l) } else { 0.0 };
        let decay = if term.l > 0 { (-delta_l).exp() } else { 1.0 };
        let v = term.n * delta.powi(term.d) * (term.t * ln_tau).exp() * decay;
        // δ·∂ln(term)/∂δ and its own δ-derivative part.
        let gd = d - l * delta_l;
        r.a += v;
        r.d += v * gd;
        r.dd += v * (gd * gd - d - l * (l - 1.0) * delta_l);
        r.t += v * term.t;
        r.tt += v * term.t * (term.t - 1.0);
        r.dt += v * term.t * gd;
    }
    for term in data.gaussian {
        let d = f64::from(term.d);
        let (dd_eps, dt_gam) = (delta - term.epsilon, tau - term.gamma);
        let v = term.n
            * delta.powi(term.d)
            * (term.t * ln_tau - term.eta * dd_eps * dd_eps - term.beta * dt_gam * dt_gam).exp();
        let gd = d - 2.0 * term.eta * delta * dd_eps;
        let gt = term.t - 2.0 * term.beta * tau * dt_gam;
        r.a += v;
        r.d += v * gd;
        r.dd += v * (gd * gd - d - 2.0 * term.eta * delta * delta);
        r.t += v * gt;
        r.tt += v * (gt * gt - term.t - 2.0 * term.beta * tau * tau);
        r.dt += v * gd * gt;
    }
    for term in data.non_analytic {
        non_analytic(term, delta, tau, &mut r);
    }
    r
}

/// Adds Span and Wagner's non-analytic term (IAPWS-95 Table 5 derivatives),
/// scaled. The term is singular exactly at δ = 1, so δ is nudged off it there
/// as `CoolProp` does.
#[expect(clippy::similar_names, reason = "symbols of IAPWS-95 Table 5")]
fn non_analytic(term: &NonAnalytic, delta: f64, tau: f64, r: &mut Residual) {
    let delta = if (delta - 1.0).abs() < 1e-12 {
        1.0 + 1e-12
    } else {
        delta
    };
    let NonAnalytic {
        n,
        a,
        b,
        beta,
        big_a,
        big_b,
        big_c,
        big_d,
    } = *term;
    let dm1 = delta - 1.0;
    let sq = dm1 * dm1;
    let tm1 = tau - 1.0;
    let theta = (1.0 - tau) + big_a * sq.powf(1.0 / (2.0 * beta));
    let big_delta = theta * theta + big_b * sq.powf(a);
    let psi = (-big_c * sq - big_d * tm1 * tm1).exp();
    let psi_d = -2.0 * big_c * dm1 * psi;
    let psi_dd = (2.0 * big_c * sq - 1.0) * 2.0 * big_c * psi;
    let psi_t = -2.0 * big_d * tm1 * psi;
    let psi_tt = (2.0 * big_d * tm1 * tm1 - 1.0) * 2.0 * big_d * psi;
    let psi_dt = 4.0 * big_c * big_d * dm1 * tm1 * psi;
    let e = 1.0 / (2.0 * beta);
    let delta_d = dm1
        * (big_a * theta * (2.0 / beta) * sq.powf(e - 1.0) + 2.0 * big_b * a * sq.powf(a - 1.0));
    let delta_dd = delta_d / dm1
        + sq * (4.0 * big_b * a * (a - 1.0) * sq.powf(a - 2.0)
            + 2.0 * big_a * big_a * (1.0 / beta).powi(2) * sq.powf(e - 1.0).powi(2)
            + big_a * theta * (4.0 / beta) * (e - 1.0) * sq.powf(e - 2.0));
    let db = big_delta.powf(b);
    let db1 = big_delta.powf(b - 1.0);
    let db2 = big_delta.powf(b - 2.0);
    let db_d = b * db1 * delta_d;
    let db_dd = b * (db1 * delta_dd + (b - 1.0) * db2 * delta_d * delta_d);
    let db_t = -2.0 * theta * b * db1;
    let db_tt = 2.0 * b * db1 + 4.0 * theta * theta * b * (b - 1.0) * db2;
    let db_dt = -big_a * b * (2.0 / beta) * db1 * dm1 * sq.powf(e - 1.0)
        - 2.0 * theta * b * (b - 1.0) * db2 * delta_d;
    let phi = n * db * delta * psi;
    let phi_d = n * (db * (psi + delta * psi_d) + db_d * delta * psi);
    let phi_dd = n
        * (db * (2.0 * psi_d + delta * psi_dd)
            + 2.0 * db_d * (psi + delta * psi_d)
            + db_dd * delta * psi);
    let phi_t = n * delta * (db_t * psi + db * psi_t);
    let phi_tt = n * delta * (db_tt * psi + 2.0 * db_t * psi_t + db * psi_tt);
    let phi_dt = n
        * (db * (psi_t + delta * psi_dt)
            + delta * db_d * psi_t
            + db_t * (psi + delta * psi_d)
            + db_dt * delta * psi);
    r.a += phi;
    r.d += delta * phi_d;
    r.dd += delta * delta * phi_dd;
    r.t += tau * phi_t;
    r.tt += tau * tau * phi_tt;
    r.dt += delta * tau * phi_dt;
}

/// Ideal part without its reference constants: α⁰ − ln δ − c₁ − c₂τ, `τα⁰_τ` −
/// c₂τ, `τ²α⁰_ττ`.
pub(crate) fn ideal_shape(data: &FluidData, tau: f64) -> (f64, f64, f64) {
    let (mut a, mut t, mut tt) = (data.log_tau * tau.ln(), data.log_tau, -data.log_tau);
    for &(n, exponent) in data.ideal_power {
        let v = n * tau.powf(exponent);
        a += v;
        t += v * exponent;
        tt += v * exponent * (exponent - 1.0);
    }
    for &(n, theta) in data.planck_einstein {
        let x = theta * tau;
        let em = (-x).exp();
        a += n * (-em).ln_1p();
        // x / (e^x − 1) and x² e^x / (e^x − 1)², written in e^−x to stay
        // finite.
        let q = em / (1.0 - em);
        t += n * x * q;
        tt -= n * x * x * q / (1.0 - em);
    }
    (a, t, tt)
}

/// Molar properties at density `rho` mol/m³ and temperature `t` K.
pub(crate) fn molar(fluid: &Fluid, rho: f64, t: f64) -> Molar {
    let data = fluid.data;
    let (delta, tau) = (rho / data.rho_crit, data.t_crit / t);
    let r = residual(fluid, delta, tau);
    let (a0, t0, tt0) = ideal_shape(data, tau);
    let a0 = a0 + delta.ln() + fluid.c1 + fluid.c2 * tau;
    let t0 = t0 + fluid.c2 * tau;
    let gas = data.gas_constant;
    let rt = gas * t;
    let x = 1.0 + r.d - r.dt;
    let y = 1.0 + 2.0 * r.d + r.dd;
    let cv_r = -(tt0 + r.tt);
    Molar {
        p: rho * rt * (1.0 + r.d),
        dp_drho: rt * y,
        dp_dt: rho * gas * x,
        u: rt * (t0 + r.t),
        h: rt * (1.0 + t0 + r.t + r.d),
        s: gas * (t0 + r.t - a0 - r.a),
        cv: gas * cv_r,
        cp: gas * (cv_r + x * x / y),
        w: (rt / data.molar_mass * (y + x * x / cv_r)).sqrt(),
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::fluid::PureFluid;

    fn fluids() -> impl Strategy<Value = PureFluid> {
        prop_oneof![
            Just(PureFluid::Nitrogen),
            Just(PureFluid::Argon),
            Just(PureFluid::Helium),
            Just(PureFluid::CarbonDioxide),
        ]
    }

    #[test]
    fn test_helmholtz_near_critical_curvature_matches_coolprop() -> anyhow::Result<()> {
        // CO2 at δ = 0.966, τ = 1.000 15, inside Span–Wagner's non-analytic
        // region: τ²αʳ_ττ from CoolProp's HEOS (homogeneous phase).
        let r = residual(
            PureFluid::CarbonDioxide.data(),
            0.965_689_538_175_702_5,
            1.000_150_975_621_335_1,
        );
        anyhow::ensure!(
            (r.tt / -9.503_751_387_434_187 - 1.0).abs() < 1e-11,
            "τ²αʳ_ττ = {}",
            r.tt
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_helmholtz_scaled_derivatives_match_central_differences(
            fluid in fluids(),
            delta in 0.05f64..2.5,
            tau in 0.3f64..2.0,
        ) {
            // Finite differences cannot follow the non-analytic terms' cusp at
            // the critical point; that corner is pinned to CoolProp below.
            prop_assume!(!((delta - 1.0).abs() < 0.05 && (tau - 1.0).abs() < 0.01));
            let f = fluid.data();
            let at = |d: f64, t: f64| residual(f, d, t);
            let r = at(delta, tau);
            // Small steps: CO2's Gaussian terms (β up to 1250) bend a coarser
            // central difference by ~1e-5.
            let (hd, ht) = (1e-6 * delta, 1e-6 * tau);
            let a_d = (at(delta + hd, tau).a - at(delta - hd, tau).a) / (2.0 * hd) * delta;
            let a_t = (at(delta, tau + ht).a - at(delta, tau - ht).a) / (2.0 * ht) * tau;
            let d_d = (at(delta + hd, tau).d - at(delta - hd, tau).d) / (2.0 * hd) * delta;
            let t_t = (at(delta, tau + ht).t - at(delta, tau - ht).t) / (2.0 * ht) * tau;
            let d_t = (at(delta, tau + ht).d - at(delta, tau - ht).d) / (2.0 * ht) * tau;
            let close = |x: f64, y: f64| (x - y).abs() <= 1e-5 * (1.0 + x.abs());
            prop_assert!(close(r.d, a_d), "δα_δ {} vs {}", r.d, a_d);
            prop_assert!(close(r.t, a_t), "τα_τ {} vs {}", r.t, a_t);
            // δ ∂(δα_δ)/∂δ = δα_δ + δ²α_δδ, and likewise in τ.
            // Sums of terms that can nearly cancel: scale by the terms.
            let close2 = |x: f64, y: f64, scale: f64| (x - y).abs() <= 1e-5 * (1.0 + scale);
            prop_assert!(close2(r.d + r.dd, d_d, r.d.abs() + r.dd.abs()), "δ²α_δδ {} vs {}", r.dd, d_d - r.d);
            prop_assert!(close2(r.t + r.tt, t_t, r.t.abs() + r.tt.abs()), "τ²α_ττ {} vs {}", r.tt, t_t - r.t);
            prop_assert!(close(r.dt, d_t), "δτα_δτ {} vs {}", r.dt, d_t);
        }
    }
}
