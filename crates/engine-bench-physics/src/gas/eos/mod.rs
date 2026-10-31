//! Equations of state for lumped gas volumes: how pressure, temperature and
//! enthalpy follow from what a volume conserves.
//!
//! [`IdealGas`] is the default and costs nothing extra; [`PengRobinson`] adds
//! real-gas departures for dense mixtures.

pub mod peng_robinson;
mod reference;

use std::fmt::Debug;

use anyhow::Context;
pub use peng_robinson::PengRobinson;
pub use reference::Reference;

use crate::{
    PhysicsError,
    flow::PerfectGas,
    thermo::{Mixture, SpeciesSet},
};

/// Per-mass state of a fluid from an equation of state, raw SI.
///
/// `h`, `cp`, `cv` and `gas_constant` describe what leaves through a port at
/// the top of the volume: the bulk for one phase, saturated vapour inside the
/// two-phase dome.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluidState {
    /// K.
    pub(crate) t: f64,
    /// Pa.
    pub(crate) p: f64,
    /// m³/kg.
    pub(crate) v: f64,
    /// J/kg.
    pub(crate) u: f64,
    /// J/kg.
    pub(crate) h: f64,
    /// J/(kg·K).
    pub(crate) cp: f64,
    /// J/(kg·K).
    pub(crate) cv: f64,
    /// p·v / T, J/(kg·K): the mixture gas constant for an ideal gas, Z times it
    /// for a real one.
    pub(crate) gas_constant: f64,
    /// Vapour mass fraction inside the dome; NaN for one phase.
    pub(crate) quality: f64,
    /// Fraction of the volume filled with liquid (1 for a dense single phase).
    pub(crate) liquid_fraction: f64,
}

impl FluidState {
    /// No state: every field NaN.
    pub(crate) const INVALID: Self = Self {
        t: f64::NAN,
        p: f64::NAN,
        v: f64::NAN,
        u: f64::NAN,
        h: f64::NAN,
        cp: f64::NAN,
        cv: f64::NAN,
        gas_constant: f64::NAN,
        quality: f64::NAN,
        liquid_fraction: f64::NAN,
    };

    /// Frozen perfect-gas stand-in for orifice and nozzle flow: γ = `c_p/c_v`
    /// and R = p·v/T. Exact for an ideal gas; for a real gas, used until
    /// real-fluid flow replaces it.
    pub(crate) fn perfect_gas(&self) -> PerfectGas {
        PerfectGas {
            gamma: self.cp / self.cv,
            gas_constant: self.gas_constant,
        }
    }
}

/// An equation of state for a mixture of the species set `S`.
pub trait EquationOfState<S: SpeciesSet<N>, const N: usize>: Copy + Debug + PartialEq {
    /// State at pressure `p` Pa and temperature `t` K. Edge: fails where the
    /// equation has no (stable) state.
    fn state_pt(&self, mixture: &Mixture<S, N>, p: f64, t: f64) -> anyhow::Result<FluidState>;

    /// State at specific internal energy `u` J/kg and specific volume `v`
    /// m³/kg, warm-started from `previous`. Kernel: never fails; returns
    /// [`FluidState`] with NaN fields when there is no stable state.
    fn state_uv(
        &self,
        mixture: &Mixture<S, N>,
        u: f64,
        v: f64,
        previous: &FluidState,
    ) -> FluidState;

    /// Specific entropy, J/(kg·K), of `state`.
    fn entropy(&self, mixture: &Mixture<S, N>, state: &FluidState) -> f64;

    /// Saturated liquid and vapour at `t` K, the liquid filling
    /// `liquid_fraction` of the volume. Only equations with a two-phase
    /// region have one.
    fn state_saturated(
        &self,
        mixture: &Mixture<S, N>,
        t: f64,
        liquid_fraction: f64,
    ) -> anyhow::Result<FluidState> {
        let _ = (mixture, t, liquid_fraction);
        Err(PhysicsError::InvalidState("this equation of state has no two-phase states").into())
    }
}

/// The ideal gas: p = `ρRT`, u and h from the NASA polynomials alone.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct IdealGas;

impl<S: SpeciesSet<N>, const N: usize> EquationOfState<S, N> for IdealGas {
    fn state_pt(&self, mixture: &Mixture<S, N>, p: f64, t: f64) -> anyhow::Result<FluidState> {
        check_pt(p, t)?;
        Ok(ideal_state(mixture, t, mixture.gas_constant_si() * t / p))
    }

    #[inline]
    fn state_uv(
        &self,
        mixture: &Mixture<S, N>,
        u: f64,
        v: f64,
        previous: &FluidState,
    ) -> FluidState {
        let t = mixture.temperature_from_internal_energy_si(u, previous.t);
        ideal_state(mixture, t, v)
    }

    fn entropy(&self, mixture: &Mixture<S, N>, state: &FluidState) -> f64 {
        mixture.entropy_si(state.t, state.p)
    }
}

/// Ideal gas at `t` K and `v` m³/kg.
#[inline]
#[expect(
    clippy::many_single_char_names,
    reason = "thermodynamic symbols p, v, T, u, h, R"
)]
fn ideal_state<S: SpeciesSet<N>, const N: usize>(
    mixture: &Mixture<S, N>,
    t: f64,
    v: f64,
) -> FluidState {
    let r = mixture.gas_constant_si();
    let (h, cp) = mixture.enthalpy_and_cp_si(t);
    let p = r * t / v;
    FluidState {
        t,
        p,
        v,
        u: h - r * t,
        h,
        cp,
        cv: cp - r,
        gas_constant: r,
        quality: f64::NAN,
        liquid_fraction: 0.0,
    }
}

/// Checks that `p` Pa and `t` K are positive and finite.
pub(crate) fn check_pt(p: f64, t: f64) -> anyhow::Result<()> {
    if !(p > 0.0 && p.is_finite()) {
        return Err(PhysicsError::InvalidState("pressure must be positive"))
            .with_context(|| format!("p = {p} Pa"));
    }
    if !(t > 0.0 && t.is_finite()) {
        return Err(PhysicsError::InvalidState("temperature must be positive"))
            .with_context(|| format!("T = {t} K"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermo::Air3;

    #[test]
    fn test_ideal_gas_state_is_the_mixture_arithmetic_of_control_volumes() -> anyhow::Result<()> {
        let air = Mixture::<Air3, 3>::dry_air()?;
        let state = IdealGas.state_pt(&air, 2e5, 400.0)?;
        let r = air.gas_constant_si();
        let (h, cp) = air.enthalpy_and_cp_si(400.0);
        anyhow::ensure!(state.v.to_bits() == (r * 400.0 / 2e5).to_bits(), "v");
        anyhow::ensure!(
            state.h.to_bits() == h.to_bits() && state.cp.to_bits() == cp.to_bits(),
            "h, cp"
        );
        let gas = state.perfect_gas();
        anyhow::ensure!(gas.gamma.to_bits() == (cp / (cp - r)).to_bits(), "γ");
        anyhow::ensure!(gas.gas_constant.to_bits() == r.to_bits(), "R");
        let back = IdealGas.state_uv(&air, state.u, state.v, &state);
        anyhow::ensure!(
            (back.t / 400.0 - 1.0).abs() < 1e-12 && (back.p / 2e5 - 1.0).abs() < 1e-12,
            "{back:?}"
        );
        let s = IdealGas.entropy(&air, &state);
        anyhow::ensure!(s.to_bits() == air.entropy_si(400.0, 2e5).to_bits(), "s");
        Ok(())
    }

    #[test]
    fn test_ideal_gas_rejects_non_positive_pressure_or_temperature() -> anyhow::Result<()> {
        let air = Mixture::<Air3, 3>::dry_air()?;
        for (p, t) in [(0.0, 300.0), (1e5, 0.0), (f64::NAN, 300.0)] {
            anyhow::ensure!(
                IdealGas.state_pt(&air, p, t).is_err(),
                "p {p}, T {t} accepted"
            );
        }
        Ok(())
    }
}
