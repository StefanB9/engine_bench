//! The propellant as the flow solvers see it: the exact equations offline,
//! tables in real time.

use engine_bench_physics::{
    flow::real::{FlowState, RealFluid},
    fluid::{FrozenVapour, PureFluid, TabulatedFluid},
    units::{MassDensity, Pressure, SpecificHeatCapacity},
};

/// A real fluid, exact or tabulated.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Flow<E> {
    Exact(E),
    Table(TabulatedFluid),
}

impl<E: RealFluid> RealFluid for Flow<E> {
    fn state_ps(&self, p: Pressure, s: SpecificHeatCapacity) -> anyhow::Result<FlowState> {
        match self {
            Self::Exact(fluid) => fluid.state_ps(p, s),
            Self::Table(fluid) => fluid.state_ps(p, s),
        }
    }

    fn state_prho(&self, p: Pressure, rho: MassDensity) -> anyhow::Result<FlowState> {
        match self {
            Self::Exact(fluid) => fluid.state_prho(p, rho),
            Self::Table(fluid) => fluid.state_prho(p, rho),
        }
    }
}

/// Equilibrium for the tank and orifices, frozen vapour for the nozzle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Fluids {
    pub(super) equilibrium: Flow<PureFluid>,
    pub(super) frozen: Flow<FrozenVapour>,
}

impl Fluids {
    /// The exact equations of `fluid`.
    pub(super) fn exact(fluid: PureFluid) -> Self {
        Self {
            equilibrium: Flow::Exact(fluid),
            frozen: Flow::Exact(FrozenVapour::new(fluid)),
        }
    }

    /// `fluid`'s tables, built on first use in the process.
    pub(super) fn tabulated(fluid: PureFluid) -> anyhow::Result<Self> {
        Ok(Self {
            equilibrium: Flow::Table(TabulatedFluid::equilibrium(fluid)?),
            frozen: Flow::Table(TabulatedFluid::frozen(fluid)?),
        })
    }
}

#[cfg(test)]
mod tests {
    use engine_bench_physics::units::{
        ThermodynamicTemperature, pressure::pascal, thermodynamic_temperature::kelvin,
    };

    use super::*;

    #[test]
    fn test_exact_flow_passes_states_through() -> anyhow::Result<()> {
        let fluids = Fluids::exact(PureFluid::Argon);
        let x = PureFluid::Argon.properties(
            Pressure::new::<pascal>(1e6),
            ThermodynamicTemperature::new::<kelvin>(300.0),
        )?;
        let by_s = fluids.equilibrium.state_ps(x.pressure, x.entropy)?;
        let direct = PureFluid::Argon.state_ps(x.pressure, x.entropy)?;
        anyhow::ensure!(by_s == direct, "{by_s:?} vs {direct:?}");
        let by_rho = fluids.frozen.state_prho(x.pressure, x.density)?;
        let frozen = FrozenVapour::new(PureFluid::Argon).state_prho(x.pressure, x.density)?;
        anyhow::ensure!(by_rho == frozen, "{by_rho:?} vs {frozen:?}");
        Ok(())
    }
}
