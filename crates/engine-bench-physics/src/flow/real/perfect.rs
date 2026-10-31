//! The calorically perfect gas as a [`RealFluid`]: h = cp·T and
//! s = cp·ln(T) − R·ln(p) (reference state at 1 K and 1 Pa). It lets the
//! real-fluid solvers be checked against the closed-form perfect-gas ones.

use anyhow::Context;

use super::{FlowState, RealFluid};
use crate::{
    PhysicsError,
    flow::PerfectGas,
    units::{
        MassDensity, Pressure, SpecificHeatCapacity, mass_density::kilogram_per_cubic_meter,
        pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
    },
};

impl RealFluid for PerfectGas {
    fn state_ps(&self, p: Pressure, s: SpecificHeatCapacity) -> anyhow::Result<FlowState> {
        let (p, s) = (p.get::<pascal>(), s.get::<joule_per_kilogram_kelvin>());
        check_pressure(p)?;
        let cp = self.cp();
        let t = ((s + self.gas_constant * p.ln()) / cp).exp();
        Ok(self.state_at(p, t))
    }

    fn state_prho(&self, p: Pressure, rho: MassDensity) -> anyhow::Result<FlowState> {
        let (p, rho) = (p.get::<pascal>(), rho.get::<kilogram_per_cubic_meter>());
        check_pressure(p)?;
        if !(rho > 0.0 && rho.is_finite()) {
            return Err(PhysicsError::InvalidState("density must be positive"))
                .with_context(|| format!("ρ = {rho} kg/m³"));
        }
        Ok(self.state_at(p, p / (rho * self.gas_constant)))
    }
}

impl PerfectGas {
    fn cp(&self) -> f64 {
        self.gamma * self.gas_constant / (self.gamma - 1.0)
    }

    fn state_at(&self, p: f64, t: f64) -> FlowState {
        let (r, cp) = (self.gas_constant, self.cp());
        FlowState::from_si(
            p,
            t,
            p / (r * t),
            cp * t,
            cp * t.ln() - r * p.ln(),
            (self.gamma * r * t).sqrt(),
            f64::NAN,
        )
    }
}

fn check_pressure(p: f64) -> anyhow::Result<()> {
    if !(p > 0.0 && p.is_finite()) {
        return Err(PhysicsError::InvalidState("pressure must be positive"))
            .with_context(|| format!("p = {p} Pa"));
    }
    Ok(())
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Air as a perfect gas.
    pub(crate) fn air() -> anyhow::Result<PerfectGas> {
        PerfectGas::new(
            1.4,
            SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(287.0),
        )
    }

    /// The perfect-gas state at (p, T).
    pub(crate) fn state(gas: &PerfectGas, p: f64, t: f64) -> anyhow::Result<FlowState> {
        gas.state_prho(
            Pressure::new::<pascal>(p),
            MassDensity::new::<kilogram_per_cubic_meter>(p / (gas.gas_constant * t)),
        )
    }

    #[test]
    fn test_perfect_gas_states_are_closed_form() -> anyhow::Result<()> {
        let gas = air()?;
        let x = state(&gas, 2e5, 400.0)?;
        anyhow::ensure!((x.t() - 400.0).abs() < 1e-12, "T = {}", x.t());
        let cp = 1.4 * 287.0 / 0.4;
        anyhow::ensure!((x.h() / (cp * 400.0) - 1.0).abs() < 1e-14, "h = {}", x.h());
        anyhow::ensure!(
            (x.w() / (1.4 * 287.0 * 400.0f64).sqrt() - 1.0).abs() < 1e-14,
            "w"
        );
        // Same entropy at a lower pressure: T/T₀ = (p/p₀)^((γ−1)/γ).
        let y = gas.state_ps(Pressure::new::<pascal>(1e5), x.entropy)?;
        let t = 400.0 * 0.5f64.powf(0.4 / 1.4);
        anyhow::ensure!((y.t() / t - 1.0).abs() < 1e-12, "T = {} vs {t}", y.t());
        anyhow::ensure!(y.vapour_quality.is_none(), "{y:?}");
        Ok(())
    }
}
