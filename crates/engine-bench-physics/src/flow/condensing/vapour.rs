//! The vapour a condensing flow carries, supersaturated (metastable): its
//! state and first derivatives at the equation of state's own variables, for
//! the march's Newton solves.

use super::Condensable;
use crate::units::{pressure::pascal, thermodynamic_temperature::kelvin};

/// The vapour at one state, per unit mass, with its derivatives in the
/// equation's natural variables (a, T): a is the molar density for a
/// Helmholtz equation, the pressure for IF97's Gibbs equation.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct VapourPoint {
    /// Pa.
    pub(crate) p: f64,
    /// kg/m³.
    pub(crate) rho: f64,
    /// J/kg.
    pub(crate) h: f64,
    /// J/(kg·K).
    pub(crate) s: f64,
    /// J/(kg·K).
    pub(crate) cp: f64,
    /// J/(kg·K).
    pub(crate) cv: f64,
    /// Speed of sound, m/s.
    pub(crate) w: f64,
    /// ∂p/∂a at constant T.
    pub(crate) dp_da: f64,
    /// ∂p/∂T at constant a.
    pub(crate) dp_dt: f64,
    /// ∂ρ/∂a at constant T.
    pub(crate) drho_da: f64,
    /// ∂ρ/∂T at constant a.
    pub(crate) drho_dt: f64,
    /// ∂h/∂a at constant T.
    pub(crate) dh_da: f64,
    /// ∂h/∂T at constant a.
    pub(crate) dh_dt: f64,
}

pub(crate) mod kernel {
    use super::VapourPoint;

    /// The equation-of-state side of a
    /// [`CondensingFluid`](super::CondensingFluid).
    pub trait VapourKernel {
        /// The natural variable a of the vapour at `p` Pa and `t` K (an
        /// iterative solve: used to start a march, not inside it).
        fn vapour_variable(&self, p: f64, t: f64) -> anyhow::Result<f64>;

        /// The vapour at (`a`, `t` K): an infallible kernel on a state the
        /// caller keeps on the vapour branch.
        fn vapour(&self, a: f64, t: f64) -> VapourPoint;

        /// The highest temperature, K, at which liquid and vapour coexist.
        fn dome_top(&self) -> f64;

        /// Excess compressibility Z − 1 (negative) of the vapour saturated
        /// at `t` K, where its pressure is `p_sat` Pa.
        fn saturated_excess(&self, t: f64, p_sat: f64) -> f64;

        /// Effective molecular diameter, m: the position of the first peak
        /// of the liquid's pair correlation, or a hard-sphere diameter.
        fn molecular_diameter(&self) -> f64;
    }
}

/// A fluid whose vapour can be marched through a nozzle while it condenses:
/// its condensation properties and its metastable vapour.
pub trait CondensingFluid: Condensable + kernel::VapourKernel {}

impl<F: Condensable + kernel::VapourKernel> CondensingFluid for F {}

/// Saturation temperature, K, at `p` Pa by Newton's method on ln pₛ(T) with
/// the Clausius–Clapeyron slope L/(R·T²), from `guess` K; NaN if it fails.
pub(crate) fn saturation_temperature<F: Condensable + ?Sized>(
    fluid: &F,
    p: f64,
    guess: f64,
    gas_constant: f64,
) -> f64 {
    let ln_p = p.ln();
    let mut t = guess;
    for _ in 0..40 {
        let Ok(x) = fluid.condensation(crate::units::ThermodynamicTemperature::new::<kelvin>(t))
        else {
            return f64::NAN;
        };
        let p_sat = x.vapour_pressure.get::<pascal>();
        let latent = x
            .latent_heat
            .get::<crate::units::available_energy::joule_per_kilogram>();
        let slope = latent / (gas_constant * t * t);
        let step = (ln_p - p_sat.ln()) / slope;
        // Steps are capped to a tenth of T: far from saturation the slope
        // is only a guide.
        let step = step.clamp(-0.1 * t, 0.1 * t);
        t += step;
        if step.abs() <= 1e-12 * t {
            return t;
        }
    }
    f64::NAN
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{fluid::PureFluid, steam::Water, units::ThermodynamicTemperature};

    #[test]
    fn test_saturation_temperature_inverts_the_vapour_pressure() -> anyhow::Result<()> {
        // Nitrogen at 1 bar boils at 77.3 K; water at 101.325 kPa at 373.12 K
        // (IF97); argon below its triple point (supercooled liquid).
        let n2 = PureFluid::Nitrogen;
        let r_n2 = 8.314_462_618 / 0.028_013_48;
        let t = saturation_temperature(&n2, 1e5, 60.0, r_n2);
        anyhow::ensure!((t - 77.244).abs() < 0.01, "N2 {t}");
        let t = saturation_temperature(&Water, 101_325.0, 330.0, 461.526);
        anyhow::ensure!((t - 373.124).abs() < 0.01, "water {t}");
        let ar = PureFluid::Argon;
        let t_low = 70.0;
        let p = ar
            .condensation(ThermodynamicTemperature::new::<kelvin>(t_low))?
            .vapour_pressure
            .get::<pascal>();
        let t = saturation_temperature(&ar, p, 90.0, 8.314_462_618 / 0.039_948);
        anyhow::ensure!((t - t_low).abs() < 1e-9, "argon {t}");
        // Beyond the dome there is no saturation temperature.
        anyhow::ensure!(
            saturation_temperature(&n2, 1e8, 100.0, r_n2).is_nan(),
            "N2 at 1000 bar"
        );
        Ok(())
    }
}
