//! Region 3 forward equation: the dimensionless Helmholtz free energy
//! φ(δ, τ) = f / (RT) and the properties that follow from it at (ρ, T).
//!
//! Infallible kernel on raw SI values; callers have validated the state.

use super::{
    R, RHO_CRIT, T_CRIT,
    coefficients::region3::{HELMHOLTZ, HELMHOLTZ_LN},
    gibbs::RawProperties,
};

/// Region 3 state at density `rho` kg/m³ and temperature `t` K: pressure, Pa,
/// ∂p/∂ρ at constant T, Pa·m³/kg, and the per-mass properties.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Region3State {
    pub(crate) p: f64,
    pub(crate) dp_drho: f64,
    pub(crate) properties: RawProperties,
}

/// Region 3 (near-critical) state at (ρ, T) (IAPWS R7-97(2012) Table 31).
#[expect(
    clippy::many_single_char_names,
    reason = "symbols of IAPWS R7-97 Table 31"
)]
pub(crate) fn region3(rho: f64, t: f64) -> Region3State {
    let (delta, tau) = (rho / RHO_CRIT, T_CRIT / t);
    // Scaled derivatives δφ_δ, δ²φ_δδ, τφ_τ, τ²φ_ττ and δτφ_δτ need no
    // division: each term contributes its exponents times its value.
    let mut phi = HELMHOLTZ_LN * delta.ln();
    let (mut d, mut dd, mut t1, mut tt, mut dt) = (HELMHOLTZ_LN, -HELMHOLTZ_LN, 0.0, 0.0, 0.0);
    for term in &HELMHOLTZ {
        let (i, j) = (f64::from(term.i), f64::from(term.j));
        let v = term.n * delta.powi(term.i) * tau.powi(term.j);
        phi += v;
        d += i * v;
        dd += i * (i - 1.0) * v;
        t1 += j * v;
        tt += j * (j - 1.0) * v;
        dt += i * j * v;
    }
    let rt = R * t;
    let x = d - dt;
    let y = 2.0 * d + dd;
    Region3State {
        p: rho * rt * d,
        dp_drho: rt * y,
        properties: RawProperties {
            v: 1.0 / rho,
            h: rt * (t1 + d),
            u: rt * t1,
            s: R * (t1 - phi),
            cp: R * (x * x / y - tt),
            cv: -R * tt,
            w: (rt * (y - x * x / tt)).sqrt(),
            // ∂ρ/∂p = 1/(∂p/∂ρ) and ∂p/∂T at constant ρ = ρRx.
            dv_dp: -1.0 / (rho * rho * rt * y),
            dv_dt: R * x / (rho * rt * y),
        },
    }
}

#[cfg(test)]
#[expect(
    clippy::unreadable_literal,
    reason = "verification values exactly as IAPWS R7-97(2012) prints them"
)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    #[test]
    fn test_region3_matches_iapws_table_33() -> anyhow::Result<()> {
        // T K, ρ kg/m³, p MPa, h kJ/kg, u kJ/kg, s kJ/(kg·K), cp kJ/(kg·K), w
        // m/s.
        let table: [[f64; 8]; 3] = [
            [
                650.0, 500.0, 25.5837018, 1863.43019, 1812.26279, 4.05427273, 13.8935717,
                502.005554,
            ],
            [
                650.0, 200.0, 22.2930643, 2375.12401, 2263.65868, 4.85438792, 44.6579342,
                383.444594,
            ],
            [
                750.0, 500.0, 78.3095639, 2258.68845, 2102.06932, 4.46971906, 6.34165359,
                760.696041,
            ],
        ];
        for [t, rho, p, h, u, s, cp, w] in table {
            let x = region3(rho, t);
            let pairs = [
                ("p", x.p, p * 1e6),
                ("v", x.properties.v, 1.0 / rho),
                ("h", x.properties.h, h * 1e3),
                ("u", x.properties.u, u * 1e3),
                ("s", x.properties.s, s * 1e3),
                ("cp", x.properties.cp, cp * 1e3),
                ("w", x.properties.w, w),
            ];
            for (name, got, want) in pairs {
                let error = (got / want - 1.0).abs();
                anyhow::ensure!(
                    error < 5e-8,
                    "{name} at {t} K, {rho} kg/m³: {got} vs {want}"
                );
            }
        }
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_region3_dp_drho_matches_a_central_difference(
            rho in 150.0f64..750.0,
            t in 660.0f64..860.0,
        ) {
            let d = 1e-4 * rho;
            let slope = (region3(rho + d, t).p - region3(rho - d, t).p) / (2.0 * d);
            let x = region3(rho, t);
            prop_assert!((x.dp_drho - slope).abs() < 1e-6 * x.dp_drho.abs().max(1e3),
                "{} vs {}", x.dp_drho, slope);
            prop_assert!((x.properties.u - (x.properties.h - x.p / rho)).abs() < 1e-8 * x.properties.h);
            // (∂v/∂p)_T and (∂v/∂T)_p from p(ρ, T): ∂ρ/∂p = 1/(∂p/∂ρ) and
            // ∂ρ/∂T = −(∂p/∂T)/(∂p/∂ρ), each by central differences.
            let dt = 1e-4 * t;
            let dp_dt = (region3(rho, t + dt).p - region3(rho, t - dt).p) / (2.0 * dt);
            let v = 1.0 / rho;
            let dv_dp = -v * v / slope;
            let dv_dt = v * v * dp_dt / slope;
            prop_assert!((x.properties.dv_dp / dv_dp - 1.0).abs() < 1e-5, "dv/dp {} vs {}", x.properties.dv_dp, dv_dp);
            prop_assert!((x.properties.dv_dt / dv_dt - 1.0).abs() < 1e-5, "dv/dT {} vs {}", x.properties.dv_dt, dv_dt);
        }
    }
}
