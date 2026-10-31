//! Peng–Robinson cubic equation of state for gas mixtures: real-gas
//! departures added to the NASA ideal-gas properties.
//!
//! Single phase only: a state inside the spinodal, (∂p/∂v)_T ≥ 0, has no
//! stable homogeneous state and is reported as invalid. Pure fluids that can
//! condense use the reference equations instead.

mod critical;
mod kernel;

use std::{fmt, marker::PhantomData};

use anyhow::Context;
pub use critical::{Critical, critical_constants};
use kernel::Parameters;

use super::{EquationOfState, FluidState, IdealGas, check_pt};
use crate::{
    PhysicsError,
    flow::real::FlowState,
    thermo::{MOLAR_GAS_CONSTANT, Mixture, SpeciesSet},
    units::thermodynamic_temperature::kelvin,
};

/// Peng–Robinson (1976) with van der Waals one-fluid mixing (`k_ij` = 0) for
/// the species set `S`.
#[derive(Clone, Copy)]
pub struct PengRobinson<S: SpeciesSet<N>, const N: usize> {
    parameters: Parameters<N>,
    set: PhantomData<S>,
}

// By hand: derives would demand Debug and PartialEq of the marker type S.
impl<S: SpeciesSet<N>, const N: usize> fmt::Debug for PengRobinson<S, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PengRobinson")
            .field("parameters", &self.parameters)
            .finish()
    }
}

impl<S: SpeciesSet<N>, const N: usize> PartialEq for PengRobinson<S, N> {
    fn eq(&self, other: &Self) -> bool {
        self.parameters == other.parameters
    }
}

impl<S: SpeciesSet<N>, const N: usize> PengRobinson<S, N> {
    /// Peng–Robinson for `S`. Fails if a species has no critical point (a
    /// radical such as OH).
    pub fn new() -> anyhow::Result<Self> {
        let mut critical = [Critical {
            temperature: 0.0,
            pressure: 0.0,
            acentric_factor: 0.0,
        }; N];
        for (slot, species) in critical.iter_mut().zip(S::SPECIES) {
            *slot = critical_constants(species.name())
                .ok_or(PhysicsError::InvalidState("species has no critical point"))
                .with_context(|| format!("Peng-Robinson for {}", species.name()))?;
        }
        Ok(Self {
            parameters: Parameters::new(&critical),
            set: PhantomData,
        })
    }

    /// Per-mass state at `t` K and molar volume `vm` m³/mol, with the ideal-gas
    /// parts from the NASA polynomials.
    #[inline]
    #[expect(
        clippy::many_single_char_names,
        reason = "symbols of the Peng-Robinson equations"
    )]
    fn state_at(&self, mixture: &Mixture<S, N>, x: &[f64; N], t: f64, vm: f64) -> FluidState {
        let r = mixture.gas_constant_si();
        // Molar quantities per kilogram: divide by M = R_molar / R.
        let per_kg = r / MOLAR_GAS_CONSTANT;
        let m = self.parameters.evaluate(x, t, vm);
        let (u_ideal, cv_ideal) = mixture.internal_energy_and_cv_si(t);
        let v = vm * per_kg;
        let u = m.u_dep.mul_add(per_kg, u_ideal);
        let cv = m.cv_dep.mul_add(per_kg, cv_ideal);
        // c_p − c_v = −T (∂p/∂T)² / (∂p/∂v), molar.
        let cp = cv - t * m.dp_dt * m.dp_dt / m.dp_dv * per_kg;
        FluidState {
            t,
            p: m.p,
            v,
            u,
            h: m.p.mul_add(v, u),
            cp,
            cv,
            gas_constant: m.p * v / t,
            quality: f64::NAN,
            liquid_fraction: 0.0,
        }
    }
}

/// Newton iterations allowed in the flow-state solves.
const MAX_ITERATIONS: u32 = 60;

impl<S: SpeciesSet<N>, const N: usize> PengRobinson<S, N> {
    /// Flow state at `t` K and molar volume `vm`: w² = −v²·(cp/cv)·(∂p/∂v) at
    /// constant T.
    fn flow_at(&self, mixture: &Mixture<S, N>, x: &[f64; N], t: f64, vm: f64) -> FlowState {
        let state = self.state_at(mixture, x, t, vm);
        let s = self.entropy(mixture, &state);
        let dp_dv = self.parameters.evaluate(x, t, vm).dp_dv;
        let molar_mass = MOLAR_GAS_CONSTANT / mixture.gas_constant_si();
        let w = (-vm * vm * state.cp / state.cv * dp_dv / molar_mass).sqrt();
        FlowState::from_si(state.p, t, 1.0 / state.v, state.h, s, w, f64::NAN)
    }
}

impl<S: SpeciesSet<N>, const N: usize> EquationOfState<S, N> for PengRobinson<S, N> {
    fn state_pt(&self, mixture: &Mixture<S, N>, p: f64, t: f64) -> anyhow::Result<FluidState> {
        check_pt(p, t).context("Peng-Robinson state")?;
        let x = mixture.mole_fractions();
        let vm = self.parameters.molar_volume(&x, t, p);
        if !(vm.is_finite() && vm > self.parameters.covolume(&x)) {
            return Err(PhysicsError::InvalidState("no Peng-Robinson volume"))
                .with_context(|| format!("p = {p} Pa, T = {t} K"));
        }
        Ok(self.state_at(mixture, &x, t, vm))
    }

    fn state_uv(
        &self,
        mixture: &Mixture<S, N>,
        u: f64,
        v: f64,
        previous: &FluidState,
    ) -> FluidState {
        const MAX_ITERATIONS: u32 = 50;
        let x = mixture.mole_fractions();
        let vm = v * MOLAR_GAS_CONSTANT / mixture.gas_constant_si();
        if !(vm > self.parameters.covolume(&x) && u.is_finite()) {
            return FluidState::INVALID;
        }
        // Newton on u(T) at fixed v: c_v > 0, so u rises monotonically with T.
        let mut t = if previous.t.is_finite() && previous.t > 0.0 {
            previous.t
        } else {
            mixture.temperature_from_internal_energy_si(u, 300.0)
        };
        for _ in 0..MAX_ITERATIONS {
            let state = self.state_at(mixture, &x, t, vm);
            let step = (state.u - u) / state.cv;
            let next = (t - step).max(0.5 * t);
            if (next - t).abs() <= 1e-12 * t {
                let state = self.state_at(mixture, &x, next, vm);
                // Inside the spinodal no homogeneous state is stable.
                let stable = self.parameters.evaluate(&x, next, vm).dp_dv < 0.0;
                return if stable { state } else { FluidState::INVALID };
            }
            t = next;
        }
        FluidState::INVALID
    }

    fn flow_state_ps(&self, mixture: &Mixture<S, N>, p: f64, s: f64) -> anyhow::Result<FlowState> {
        check_pt(p, 1.0).context("isentropic Peng-Robinson state")?;
        let x = mixture.mole_fractions();
        // Newton on T at constant p from the ideal-gas answer: ∂s/∂T = c_p/T.
        let mut t = IdealGas
            .flow_state_ps(mixture, p, s)?
            .temperature
            .get::<kelvin>();
        for _ in 0..MAX_ITERATIONS {
            let vm = self.parameters.molar_volume(&x, t, p);
            let state = self.state_at(mixture, &x, t, vm);
            let step = (self.entropy(mixture, &state) - s) * t / state.cp;
            let next = (t - step).max(0.5 * t);
            if (next - t).abs() <= 1e-12 * t {
                let vm = self.parameters.molar_volume(&x, next, p);
                return Ok(self.flow_at(mixture, &x, next, vm));
            }
            t = next;
        }
        Err(PhysicsError::NoConvergence {
            iterations: MAX_ITERATIONS,
            residual: f64::NAN,
        })
        .with_context(|| format!("Peng-Robinson at {p} Pa, s = {s} J/(kg K)"))
    }

    fn flow_state_prho(
        &self,
        mixture: &Mixture<S, N>,
        p: f64,
        rho: f64,
    ) -> anyhow::Result<FlowState> {
        check_pt(p, rho).context("Peng-Robinson state at (p, ρ)")?;
        let x = mixture.mole_fractions();
        let vm = MOLAR_GAS_CONSTANT / (mixture.gas_constant_si() * rho);
        // Newton on T at constant v: ∂p/∂T > 0.
        let mut t = p * vm / MOLAR_GAS_CONSTANT;
        for _ in 0..MAX_ITERATIONS {
            let m = self.parameters.evaluate(&x, t, vm);
            let next = (t - (m.p - p) / m.dp_dt).max(0.5 * t);
            if (next - t).abs() <= 1e-12 * t {
                if self.parameters.evaluate(&x, next, vm).dp_dv >= 0.0 {
                    return Err(PhysicsError::InvalidState("state inside the spinodal"))
                        .with_context(|| format!("Peng-Robinson at {p} Pa, {rho} kg/m³"));
                }
                return Ok(self.flow_at(mixture, &x, next, vm));
            }
            t = next;
        }
        Err(PhysicsError::NoConvergence {
            iterations: MAX_ITERATIONS,
            residual: f64::NAN,
        })
        .with_context(|| format!("Peng-Robinson at {p} Pa, {rho} kg/m³"))
    }

    fn entropy(&self, mixture: &Mixture<S, N>, state: &FluidState) -> f64 {
        let r = mixture.gas_constant_si();
        let x = mixture.mole_fractions();
        let vm = state.v * MOLAR_GAS_CONSTANT / r;
        // Ideal gas at the same (T, v), plus the departure.
        let s_dep = self.parameters.evaluate(&x, state.t, vm).s_dep;
        mixture.entropy_si(state.t, r * state.t / state.v) + s_dep * r / MOLAR_GAS_CONSTANT
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::{
        gas::eos::IdealGas,
        thermo::{Air3, ColdGas, Species, gases},
    };

    /// A set holding a radical, which has no critical point.
    #[derive(Debug, Clone, Copy, PartialEq)]
    struct Hydroxyl;

    impl SpeciesSet<1> for Hydroxyl {
        const SPECIES: [&'static Species; 1] = [&gases::OH];
    }

    fn nitrogen() -> anyhow::Result<Mixture<ColdGas, 4>> {
        Mixture::from_mass_fractions([1.0, 0.0, 0.0, 0.0])
    }

    #[test]
    fn test_peng_robinson_needs_critical_constants_for_every_species() -> anyhow::Result<()> {
        PengRobinson::<Air3, 3>::new()?;
        PengRobinson::<ColdGas, 4>::new()?;
        anyhow::ensure!(PengRobinson::<Hydroxyl, 1>::new().is_err(), "OH accepted");
        Ok(())
    }

    #[test]
    fn test_peng_robinson_dense_nitrogen_is_compressed_below_ideal() -> anyhow::Result<()> {
        // 200 K, 100 bar: N2 is ~16 % denser than ideal (reference Z = 0.845).
        let pr = PengRobinson::<ColdGas, 4>::new()?;
        let n2 = nitrogen()?;
        let state = pr.state_pt(&n2, 1e7, 200.0)?;
        let z = state.gas_constant / n2.gas_constant_si();
        anyhow::ensure!((z - 0.824).abs() < 5e-3, "Z = {z}");
        // Attraction lowers u and h below the ideal gas at the same T.
        let ideal = IdealGas.state_pt(&n2, 1e7, 200.0)?;
        anyhow::ensure!(
            state.u < ideal.u && state.h < ideal.h,
            "{state:?} vs {ideal:?}"
        );
        anyhow::ensure!(state.cp > ideal.cp, "cp {} vs {}", state.cp, ideal.cp);
        Ok(())
    }

    #[test]
    fn test_peng_robinson_reports_states_inside_the_spinodal() -> anyhow::Result<()> {
        // CO2 at 250 K with a volume between its liquid and vapour.
        let pr = PengRobinson::<ColdGas, 4>::new()?;
        let co2 = Mixture::<ColdGas, 4>::from_mass_fractions([0.0, 0.0, 0.0, 1.0])?;
        let liquid = pr.state_pt(&co2, 30e5, 250.0)?;
        let state = pr.state_uv(&co2, liquid.u, 4.0 * liquid.v, &liquid);
        anyhow::ensure!(state.t.is_nan(), "{state:?}");
        Ok(())
    }

    #[test]
    fn test_peng_robinson_rejects_invalid_pressure_or_temperature() -> anyhow::Result<()> {
        let pr = PengRobinson::<ColdGas, 4>::new()?;
        let n2 = nitrogen()?;
        for (p, t) in [(0.0, 300.0), (1e5, 0.0), (f64::NAN, 300.0)] {
            anyhow::ensure!(pr.state_pt(&n2, p, t).is_err(), "p {p}, T {t} accepted");
        }
        Ok(())
    }

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_peng_robinson_state_uv_inverts_state_pt(
            t in 150.0f64..900.0,
            log_p in 4.0f64..7.6,
            y in 0.0f64..1.0,
        ) {
            // N2–He mixtures, supercritical everywhere in this range.
            let pr = PengRobinson::<ColdGas, 4>::new().map_err(|e| fail(&e))?;
            let mix = Mixture::<ColdGas, 4>::from_mass_fractions([y, 0.0, 1.0 - y, 0.0]).map_err(|e| fail(&e))?;
            let p = 10f64.powf(log_p);
            let state = pr.state_pt(&mix, p, t).map_err(|e| fail(&e))?;
            let start = FluidState { t: t * 1.2, ..state };
            let back = pr.state_uv(&mix, state.u, state.v, &start);
            prop_assert!((back.t / t - 1.0).abs() < 1e-10, "T {} vs {}", back.t, t);
            prop_assert!((back.p / p - 1.0).abs() < 1e-9, "p {} vs {}", back.p, p);
            // c_p is ∂h/∂T at constant p.
            let dt = 1e-4 * t;
            let hot = pr.state_pt(&mix, p, t + dt).map_err(|e| fail(&e))?;
            let cold = pr.state_pt(&mix, p, t - dt).map_err(|e| fail(&e))?;
            let cp = (hot.h - cold.h) / (2.0 * dt);
            prop_assert!((state.cp / cp - 1.0).abs() < 1e-5, "cp {} vs {}", state.cp, cp);
        }

        #[test]
        fn test_peng_robinson_reduces_to_the_ideal_gas_at_low_pressure(t in 150.0f64..900.0, y in 0.0f64..1.0) {
            let pr = PengRobinson::<ColdGas, 4>::new().map_err(|e| fail(&e))?;
            let mix = Mixture::<ColdGas, 4>::from_mass_fractions([y, 1.0 - y, 0.0, 0.0]).map_err(|e| fail(&e))?;
            let real = pr.state_pt(&mix, 10.0, t).map_err(|e| fail(&e))?;
            let ideal = IdealGas.state_pt(&mix, 10.0, t).map_err(|e| fail(&e))?;
            prop_assert!((real.v / ideal.v - 1.0).abs() < 1e-6);
            prop_assert!((real.h - ideal.h).abs() < 1e-3 * ideal.cp);
            let (s_real, s_ideal) = (pr.entropy(&mix, &real), IdealGas.entropy(&mix, &ideal));
            prop_assert!((s_real - s_ideal).abs() < 1e-3, "s {} vs {}", s_real, s_ideal);
        }
    }
}
