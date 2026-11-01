//! The reference fluids' metastable vapour for condensing flow: the Helmholtz
//! equation at (ρ, T), with ρ the natural variable.

use super::{
    PureFluid,
    frozen::vapour_density,
    helmholtz::{molar, residual},
    saturation::dome_top,
    transport::lennard_jones_diameter,
};
use crate::flow::condensing::{VapourPoint, kernel::VapourKernel};

impl VapourKernel for PureFluid {
    /// The vapour's molar density, mol/m³.
    fn vapour_variable(&self, p: f64, t: f64) -> anyhow::Result<f64> {
        vapour_density(self.data(), p, t)
    }

    fn vapour(&self, a: f64, t: f64) -> VapourPoint {
        let fluid = self.data();
        let m = fluid.data.molar_mass;
        let x = molar(fluid, a, t);
        VapourPoint {
            p: x.p,
            rho: a * m,
            h: x.h / m,
            s: x.s / m,
            cp: x.cp / m,
            cv: x.cv / m,
            w: x.w,
            dp_da: x.dp_drho,
            dp_dt: x.dp_dt,
            drho_da: m,
            drho_dt: 0.0,
            // h = u + p/ρ with (∂u/∂ρ)_T = (p − T (∂p/∂T)_ρ)/ρ².
            dh_da: (x.dp_drho / a - t * x.dp_dt / (a * a)) / m,
            dh_dt: (x.cv + x.dp_dt / a) / m,
        }
    }

    fn dome_top(&self) -> f64 {
        dome_top(self.data())
    }

    /// From the second virial coefficient B(T) = lim (∂αʳ/∂δ)/ρ꜀ (δ → 0): Z − 1
    /// ≈ B·p/(RT), as Kalikmanov evaluates it (his Eq. 26).
    fn saturated_excess(&self, t: f64, p_sat: f64) -> f64 {
        const DILUTE: f64 = 1e-10;
        let fluid = self.data();
        let r = residual(fluid, DILUTE, fluid.data.t_crit / t);
        let b = r.d / DILUTE / fluid.data.rho_crit;
        b * p_sat / (fluid.data.gas_constant * t)
    }

    /// The Lennard-Jones σ of the dilute-gas viscosity correlation.
    fn molecular_diameter(&self) -> f64 {
        lennard_jones_diameter(*self)
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        /// The derivatives match central differences of the kernel's own
        /// values, on vapour from dilute to near saturation.
        #[test]
        fn test_vapour_kernel_derivatives_match_differences(
            fluid in prop_oneof![
                Just(PureFluid::Nitrogen),
                Just(PureFluid::Argon),
                Just(PureFluid::CarbonDioxide),
            ],
            t_fraction in 0.6f64..2.5,
            ln_p in 2.0f64..13.0,
        ) {
            let t = t_fraction * fluid.data().data.t_triple;
            let a = fluid.vapour_variable(ln_p.exp(), t);
            prop_assume!(a.is_ok());
            let a = a.map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let v = fluid.vapour(a, t);
            let (da, dt) = (1e-6 * a, 1e-6 * t);
            let d_a = |f: fn(&VapourPoint) -> f64| {
                (f(&fluid.vapour(a + da, t)) - f(&fluid.vapour(a - da, t))) / (2.0 * da)
            };
            let d_t = |f: fn(&VapourPoint) -> f64| {
                (f(&fluid.vapour(a, t + dt)) - f(&fluid.vapour(a, t - dt))) / (2.0 * dt)
            };
            let close = |x: f64, y: f64, scale: f64| (x - y).abs() <= 1e-5 * scale;
            prop_assert!(close(v.dp_da, d_a(|x| x.p), v.p / a), "dp/da {} vs {}", v.dp_da, d_a(|x| x.p));
            prop_assert!(close(v.dp_dt, d_t(|x| x.p), v.p / t), "dp/dT");
            prop_assert!(close(v.drho_da, d_a(|x| x.rho), v.rho / a), "dρ/da");
            prop_assert!(close(v.dh_da, d_a(|x| x.h), v.cp * t / a), "dh/da {} vs {}", v.dh_da, d_a(|x| x.h));
            prop_assert!(close(v.dh_dt, d_t(|x| x.h), v.cp), "dh/dT {} vs {}", v.dh_dt, d_t(|x| x.h));
            prop_assert!((v.p / ln_p.exp() - 1.0).abs() < 1e-8, "p");
        }
    }
}
