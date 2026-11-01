//! Water's vapour for condensing flow, with the pressure as natural
//! variable: IF97's region 2 equation where the vapour is superheated, its
//! supplementary metastable-vapour equation (Eq. 18) where it is
//! supersaturated, up to 10 MPa. The two agree on the saturation line to
//! 0.014 % in v and 0.78 % in cp; the basic equation itself diverges a few
//! kelvin into the metastable range (cp triples 30 K below saturation at
//! 12 kPa).

use super::{
    T_CRIT, Water,
    gibbs::{metastable_vapour, region2},
    saturation::saturation_temperature_si,
};
use crate::{
    PhysicsError,
    flow::condensing::{VapourPoint, kernel::VapourKernel},
};

/// Highest pressure of the metastable-vapour equation, Pa.
const METASTABLE_MAX_PRESSURE: f64 = 1e7;

impl VapourKernel for Water {
    /// The pressure itself, Pa.
    fn vapour_variable(&self, p: f64, t: f64) -> anyhow::Result<f64> {
        if p > 0.0 && p.is_finite() && t > 0.0 && t.is_finite() {
            Ok(p)
        } else {
            Err(PhysicsError::InvalidState("vapour needs positive p and T").into())
        }
    }

    fn vapour(&self, a: f64, t: f64) -> VapourPoint {
        let x = if a <= METASTABLE_MAX_PRESSURE && t < saturation_temperature_si(a) {
            metastable_vapour(a, t)
        } else {
            region2(a, t)
        };
        let rho = 1.0 / x.v;
        VapourPoint {
            p: a,
            rho,
            h: x.h,
            s: x.s,
            cp: x.cp,
            cv: x.cv,
            w: x.w,
            dp_da: 1.0,
            dp_dt: 0.0,
            drho_da: -x.dv_dp * rho * rho,
            drho_dt: -x.dv_dt * rho * rho,
            dh_da: x.v - t * x.dv_dt,
            dh_dt: x.cp,
        }
    }

    fn dome_top(&self) -> f64 {
        T_CRIT
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    #[test]
    fn test_water_vapour_is_the_metastable_equation_below_saturation() -> anyhow::Result<()> {
        // 30 K below saturation at 11.56 kPa (321.8 K): cp stays within ~10 %
        // of the ideal gas's 1.87 kJ/(kg·K) (the basic equation gives 5.9), and
        // the isentrope from 40.05 kPa, 374.3 K reaches this pressure at about
        // 277 K (CoolProp IF97 with cp from 330 K).
        let v = Water.vapour(11_560.0, 291.8);
        anyhow::ensure!(v.cp > 1_850.0 && v.cp < 2_200.0, "cp {}", v.cp);
        let s0 = Water.vapour(40_050.0, 374.3).s;
        let (cold, warm) = (
            Water.vapour(11_560.0, 277.0).s,
            Water.vapour(11_560.0, 283.0).s,
        );
        anyhow::ensure!(cold < s0 && s0 < warm, "s {cold}..{warm} vs {s0}");
        // Superheated: the basic equation, as IAPWS Table 15 (300 K, 3.5 kPa).
        let v = Water.vapour(3_500.0, 300.0);
        anyhow::ensure!((v.h / 2_549_911.45 - 1.0).abs() < 1e-8, "h {}", v.h);
        Ok(())
    }

    #[test]
    fn test_water_vapour_rejects_impossible_states() {
        assert!(Water.vapour_variable(0.0, 300.0).is_err());
        assert!(Water.vapour_variable(1e4, f64::NAN).is_err());
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        /// The derivatives match central differences, down into supersaturated
        /// vapour (40 K below saturation).
        #[test]
        fn test_water_vapour_kernel_derivatives_match_differences(
            t in 280.0f64..500.0,
            p in 1e3f64..2e5,
        ) {
            // Supersaturated nozzle vapour is at most ~45 K below saturation.
            prop_assume!(t > crate::steam::saturation::saturation_temperature_si(p) - 45.0);
            let a = Water.vapour_variable(p, t).map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let v = Water.vapour(a, t);
            let (da, dt) = (1e-6 * a, 1e-6 * t);
            let d_a = |f: fn(&VapourPoint) -> f64| {
                (f(&Water.vapour(a + da, t)) - f(&Water.vapour(a - da, t))) / (2.0 * da)
            };
            let d_t = |f: fn(&VapourPoint) -> f64| {
                (f(&Water.vapour(a, t + dt)) - f(&Water.vapour(a, t - dt))) / (2.0 * dt)
            };
            let close = |x: f64, y: f64, scale: f64| (x - y).abs() <= 1e-5 * scale;
            prop_assert!(close(v.drho_da, d_a(|x| x.rho), v.rho / a), "dρ/dp");
            prop_assert!(close(v.drho_dt, d_t(|x| x.rho), v.rho / t), "dρ/dT");
            prop_assert!(close(v.dh_da, d_a(|x| x.h), v.cp * t / a), "dh/dp {} vs {}", v.dh_da, d_a(|x| x.h));
            prop_assert!(close(v.dh_dt, d_t(|x| x.h), v.cp), "dh/dT");
        }
    }
}
