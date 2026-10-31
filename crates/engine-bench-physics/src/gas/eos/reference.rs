//! Reference equations of state for a volume holding one pure fluid of a
//! species set: liquid, vapour, or both.

use std::{fmt, marker::PhantomData};

use anyhow::Context;

use super::{EquationOfState, FluidState};
use crate::{
    PhysicsError,
    flow::real::{FlowState, RealFluid},
    fluid::{Molar, Phase, PureFluid, RawState, density_si, molar, saturation_si, state_uv_si},
    thermo::{Mixture, SpeciesSet},
    units::{
        MassDensity, Pressure, SpecificHeatCapacity, mass_density::kilogram_per_cubic_meter,
        pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
    },
};

/// The reference equation of `fluid`, for volumes of the species set `S`
/// that hold that fluid alone. Any other species makes the state invalid.
#[derive(Clone, Copy)]
pub struct Reference<S: SpeciesSet<N>, const N: usize> {
    fluid: PureFluid,
    /// Index of the fluid's species in `S`.
    index: usize,
    set: PhantomData<S>,
}

// By hand: derives would demand Debug and PartialEq of the marker type S.
impl<S: SpeciesSet<N>, const N: usize> fmt::Debug for Reference<S, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Reference")
            .field("fluid", &self.fluid)
            .field("index", &self.index)
            .finish_non_exhaustive()
    }
}

impl<S: SpeciesSet<N>, const N: usize> PartialEq for Reference<S, N> {
    fn eq(&self, other: &Self) -> bool {
        self.fluid == other.fluid
    }
}

impl<S: SpeciesSet<N>, const N: usize> Reference<S, N> {
    /// The reference equation of `fluid`. Fails if `S` lacks its species.
    pub fn new(fluid: PureFluid) -> anyhow::Result<Self> {
        let index = S::index_of(fluid.species())
            .ok_or(PhysicsError::InvalidComposition(
                "the species set lacks the fluid",
            ))
            .with_context(|| format!("{fluid:?}"))?;
        Ok(Self {
            fluid,
            index,
            set: PhantomData,
        })
    }

    /// Fails unless `mixture` is this fluid alone.
    fn require_pure(&self, mixture: &Mixture<S, N>) -> anyhow::Result<()> {
        if self.is_pure(mixture) {
            return Ok(());
        }
        Err(PhysicsError::InvalidComposition(
            "a reference equation needs a pure fluid",
        ))
        .with_context(|| format!("{:?}", self.fluid))
    }

    /// Whether `mixture` is this fluid alone (to round-off).
    fn is_pure(&self, mixture: &Mixture<S, N>) -> bool {
        mixture.mass_fractions()[self.index] >= 1.0 - PURITY
    }

    /// Per-mass state from a molar equilibrium state.
    fn per_mass(&self, state: &RawState) -> FluidState {
        let m = self.fluid.data().data.molar_mass;
        let v = 1.0 / (state.rho * m);
        let port = state.phases.map_or((state.bulk, state.rho), |(_, vapour)| {
            (vapour.properties, vapour.rho)
        });
        let (outflow, rho_out) = port;
        let liquid_fraction = match state.phases {
            Some((liquid, _)) => (1.0 - state.quality) * state.rho / liquid.rho,
            None if state.rho >= self.fluid.data().data.rho_crit => 1.0,
            None => 0.0,
        };
        FluidState {
            t: state.t,
            p: state.bulk.p,
            v,
            u: state.bulk.u / m,
            h: outflow.h / m,
            cp: outflow.cp / m,
            cv: outflow.cv / m,
            gas_constant: outflow.p / (rho_out * m * state.t),
            quality: state.quality,
            liquid_fraction,
        }
    }
}

/// Mass fraction of other species tolerated as round-off.
const PURITY: f64 = 1e-9;

impl<S: SpeciesSet<N>, const N: usize> EquationOfState<S, N> for Reference<S, N> {
    fn state_pt(&self, mixture: &Mixture<S, N>, p: f64, t: f64) -> anyhow::Result<FluidState> {
        if !self.is_pure(mixture) {
            return Err(PhysicsError::InvalidComposition(
                "a reference equation needs a pure fluid",
            ))
            .with_context(|| format!("{:?}", self.fluid));
        }
        let fluid = self.fluid.data();
        let rho = density_si(fluid, p, t)?;
        Ok(self.per_mass(&RawState {
            t,
            rho,
            quality: f64::NAN,
            bulk: molar(fluid, rho, t),
            phases: None,
        }))
    }

    fn state_uv(
        &self,
        mixture: &Mixture<S, N>,
        u: f64,
        v: f64,
        previous: &FluidState,
    ) -> FluidState {
        if !self.is_pure(mixture) {
            return FluidState::INVALID;
        }
        let fluid = self.fluid.data();
        let m = fluid.data.molar_mass;
        // The solver needs only the previous temperature and whether it was
        // two-phase.
        let hint = previous.t.is_finite().then(|| {
            let nan = f64::NAN;
            let dummy = Molar {
                p: nan,
                dp_drho: nan,
                dp_dt: nan,
                u: nan,
                h: nan,
                s: nan,
                cv: nan,
                cp: nan,
                w: nan,
            };
            let phase = Phase {
                rho: nan,
                properties: dummy,
            };
            RawState {
                t: previous.t,
                rho: 1.0 / (previous.v * m),
                quality: previous.quality,
                bulk: dummy,
                phases: (!previous.quality.is_nan()).then_some((phase, phase)),
            }
        });
        state_uv_si(fluid, u * m, 1.0 / (v * m), hint.as_ref())
            .map_or(FluidState::INVALID, |state| self.per_mass(&state))
    }

    fn flow_state_ps(&self, mixture: &Mixture<S, N>, p: f64, s: f64) -> anyhow::Result<FlowState> {
        self.require_pure(mixture)?;
        self.fluid.state_ps(
            Pressure::new::<pascal>(p),
            SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(s),
        )
    }

    fn flow_state_prho(
        &self,
        mixture: &Mixture<S, N>,
        p: f64,
        rho: f64,
    ) -> anyhow::Result<FlowState> {
        self.require_pure(mixture)?;
        self.fluid.state_prho(
            Pressure::new::<pascal>(p),
            MassDensity::new::<kilogram_per_cubic_meter>(rho),
        )
    }

    fn entropy(&self, mixture: &Mixture<S, N>, state: &FluidState) -> f64 {
        let _ = mixture;
        let fluid = self.fluid.data();
        let m = fluid.data.molar_mass;
        if state.quality.is_nan() {
            return molar(fluid, 1.0 / (state.v * m), state.t).s / m;
        }
        saturation_si(fluid, state.t).map_or(f64::NAN, |(_, rho_l, rho_v)| {
            let (l, v) = (molar(fluid, rho_l, state.t), molar(fluid, rho_v, state.t));
            state.quality.mul_add(v.s - l.s, l.s) / m
        })
    }

    #[expect(
        clippy::many_single_char_names,
        reason = "thermodynamic symbols p, T, u, v, h, s, x"
    )]
    fn state_saturated(
        &self,
        mixture: &Mixture<S, N>,
        t: f64,
        liquid_fraction: f64,
    ) -> anyhow::Result<FluidState> {
        if !self.is_pure(mixture) {
            return Err(PhysicsError::InvalidComposition(
                "a reference equation needs a pure fluid",
            ))
            .with_context(|| format!("{:?}", self.fluid));
        }
        let fluid = self.fluid.data();
        let (p, rho_l, rho_v) = saturation_si(fluid, t)?;
        let (l, v) = (molar(fluid, rho_l, t), molar(fluid, rho_v, t));
        // Moles per unit volume in each phase.
        let (n_l, n_v) = (liquid_fraction * rho_l, (1.0 - liquid_fraction) * rho_v);
        let rho = n_l + n_v;
        let x = n_v / rho;
        let mix = |a: f64, b: f64| x.mul_add(b - a, a);
        Ok(self.per_mass(&RawState {
            t,
            rho,
            quality: x,
            bulk: Molar {
                p,
                u: mix(l.u, v.u),
                h: mix(l.h, v.h),
                s: mix(l.s, v.s),
                ..l
            },
            phases: Some((
                Phase {
                    rho: rho_l,
                    properties: l,
                },
                Phase {
                    rho: rho_v,
                    properties: v,
                },
            )),
        }))
    }
}
