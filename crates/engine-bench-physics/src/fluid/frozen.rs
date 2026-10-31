//! Frozen (supersaturated) vapour expansion: states on the vapour branch of
//! the reference equation, metastable inside the two-phase dome and continued
//! below the triple point.
//!
//! Small, fast nozzles expand vapour far past the saturation line before
//! droplets can nucleate; this is that frozen limit. The ideal part of the
//! equation is analytic at any temperature and the residual part is small at
//! the low densities of such states. Below the triple point the equation is
//! extrapolated, only for nearly ideal vapour where its second-virial term
//! dominates: |Z − 1| ≤ 0.1 down to half the triple temperature, and dilute
//! vapour, |Z − 1| ≤ 0.02, down to a fifth of it (a vacuum nozzle's exit).
//! Beyond the vapour spinodal no vapour state exists (spontaneous
//! condensation) and the solvers fail.

use anyhow::Context;

use super::{
    Fluid, PureFluid,
    flow::single_phase,
    helmholtz::{Molar, molar},
    saturation::ancillary,
    state::max_density,
};
use crate::{
    PhysicsError,
    flow::real::{FlowState, RealFluid},
    numerics::roots::{RootOptions, newton_bracketed},
    units::{
        MassDensity, Pressure, SpecificHeatCapacity, mass_density::kilogram_per_cubic_meter,
        pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
    },
};

/// Frozen vapour reaches down to this fraction of the triple temperature:
/// each equation's second virial coefficient stays smooth that far.
const LOWEST_FRACTION_OF_TRIPLE: f64 = 0.2;
/// Below the triple point the vapour must stay this close to ideal, |Z − 1|,
/// so that the second-virial term dominates the extrapolated residual part
/// (the higher terms diverge there): down to half the triple temperature,
/// where each equation's second virial coefficient is still smooth, ...
const MAX_NON_IDEALITY: f64 = 0.1;
/// ... and below that, where only dilute vapour (a vacuum nozzle's exit)
/// is trusted, its residual part a small correction whatever its
/// extrapolation.
const MAX_NON_IDEALITY_COLD: f64 = 0.02;
/// Where the colder limit starts, as a fraction of the triple temperature.
const COLD_FRACTION_OF_TRIPLE: f64 = 0.5;
/// Density growth factor while bracketing the vapour root: small enough not
/// to step over the unstable part of a van der Waals loop.
const GROWTH: f64 = 1.1;

/// A pure fluid's vapour, frozen: it never condenses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FrozenVapour {
    fluid: PureFluid,
}

impl FrozenVapour {
    /// The frozen vapour of `fluid`.
    pub const fn new(fluid: PureFluid) -> Self {
        Self { fluid }
    }

    /// The fluid.
    pub const fn fluid(self) -> PureFluid {
        self.fluid
    }
}

impl RealFluid for FrozenVapour {
    #[expect(
        clippy::many_single_char_names,
        reason = "the thermodynamic symbols p, s, T, x"
    )]
    fn state_ps(&self, p: Pressure, s: SpecificHeatCapacity) -> anyhow::Result<FlowState> {
        let fluid = self.fluid.data();
        let m = fluid.data.molar_mass;
        let (p, s) = (p.get::<pascal>(), s.get::<joule_per_kilogram_kelvin>() * m);
        check_pressure(fluid, p)?;
        let (lo, hi) = (lowest(fluid), fluid.data.t_max);
        // s rises with T along the isobar, nearly linearly in ln T (∂s/∂ln T =
        // c_p). Where no vapour exists the temperature is too low for this
        // pressure: a negative residual steers the search up.
        let mut failure = None;
        let ln_t = newton_bracketed(
            |ln_t| {
                let t = ln_t.exp().clamp(lo, hi);
                match vapour_density(fluid, p, t) {
                    Ok(rho) => {
                        let x = molar(fluid, rho, t);
                        (x.s - s, x.cp)
                    }
                    Err(e) => {
                        failure = Some(e);
                        (-f64::MAX, f64::NAN)
                    }
                }
            },
            lo.ln(),
            hi.ln(),
            f64::midpoint(lo.ln(), hi.ln()),
            RootOptions::new(1e-13, 200)?,
        );
        let t = match ln_t {
            Ok(ln_t) => ln_t.exp().clamp(lo, hi),
            Err(e) => {
                let e = e.context(format!(
                    "frozen {:?} at {p} Pa, s = {s} J/(mol K)",
                    self.fluid
                ));
                return Err(match failure {
                    Some(cause) => e.context(format!("{cause:#}")),
                    None => e,
                });
            }
        };
        let rho = vapour_density(fluid, p, t)
            .with_context(|| format!("frozen {:?} at {p} Pa, {t} K", self.fluid))?;
        let x = molar(fluid, rho, t);
        // The search ends at the edge of the vapour's range when no vapour has
        // this entropy (a liquid-like state, or one too cold to trust).
        if (x.s - s).abs() > 1e-9 * s.abs().max(fluid.data.gas_constant) {
            return Err(PhysicsError::InvalidState(
                "no vapour state has this entropy",
            ))
            .with_context(|| format!("frozen {:?} at {p} Pa, s = {s} J/(mol K)", self.fluid));
        }
        Ok(single_phase(fluid, &x, rho, t))
    }

    fn state_prho(&self, p: Pressure, rho: MassDensity) -> anyhow::Result<FlowState> {
        let fluid = self.fluid.data();
        let m = fluid.data.molar_mass;
        let (p, rho) = (p.get::<pascal>(), rho.get::<kilogram_per_cubic_meter>() / m);
        check_pressure(fluid, p)?;
        if !(rho > 0.0 && rho.is_finite()) {
            return Err(PhysicsError::InvalidState("density must be positive"))
                .with_context(|| format!("ρ = {rho} mol/m³"));
        }
        let (lo, hi) = (lowest(fluid), fluid.data.t_max);
        // Along the isobar the vapour density falls as T rises, (∂ρ/∂T)_p =
        // −(∂p/∂T)_ρ / (∂p/∂ρ)_T; solving through the vapour root keeps the
        // answer on the vapour branch (the equation has spurious loops inside
        // the dome). Where no vapour exists T is too low: a positive residual
        // steers the search up.
        let ln_t = newton_bracketed(
            |ln_t| {
                let t = ln_t.exp().clamp(lo, hi);
                match vapour_density(fluid, p, t) {
                    Ok(root) => {
                        let x = molar(fluid, root, t);
                        (root - rho, -t * x.dp_dt / x.dp_drho)
                    }
                    Err(_) => (f64::MAX, f64::NAN),
                }
            },
            lo.ln(),
            hi.ln(),
            f64::midpoint(lo.ln(), hi.ln()),
            RootOptions::new(1e-13, 200)?,
        )
        .with_context(|| format!("frozen {:?} at {p} Pa, {rho} mol/m³", self.fluid))?;
        let t = ln_t.exp().clamp(lo, hi);
        let x = molar(fluid, rho, t);
        if (x.p - p).abs() > 1e-9 * p || !on_vapour_branch(fluid, rho, t, &x) {
            return Err(PhysicsError::InvalidState("not a vapour state"))
                .with_context(|| format!("frozen {:?} at {p} Pa, {rho} mol/m³", self.fluid));
        }
        Ok(single_phase(fluid, &x, rho, t))
    }
}

/// Lowest temperature of frozen vapour, K.
pub(super) fn lowest(fluid: &Fluid) -> f64 {
    LOWEST_FRACTION_OF_TRIPLE * fluid.data.t_triple
}

fn check_pressure(fluid: &Fluid, p: f64) -> anyhow::Result<()> {
    if !(p > 0.0 && p <= fluid.data.p_max) {
        return Err(PhysicsError::PressureOutOfRange {
            pressure: p,
            min: 0.0,
            max: fluid.data.p_max,
        })
        .with_context(|| fluid.data.species.to_string());
    }
    Ok(())
}

/// Below the critical temperature the vapour branch is the stable part of
/// p(ρ) under the critical density; above it the one fluid. Below the triple
/// point only nearly ideal vapour is trusted.
fn on_vapour_branch(fluid: &Fluid, rho: f64, t: f64, x: &Molar) -> bool {
    x.dp_drho > 0.0 && rho <= vapour_ceiling(fluid, t) && trusted(fluid, rho, t, x.p)
}

/// Whether the extrapolated equation holds at `rho` mol/m³, `t` K and `p` Pa.
fn trusted(fluid: &Fluid, rho: f64, t: f64, p: f64) -> bool {
    let t_triple = fluid.data.t_triple;
    if t >= t_triple {
        return true;
    }
    let limit = if t >= COLD_FRACTION_OF_TRIPLE * t_triple {
        MAX_NON_IDEALITY
    } else {
        MAX_NON_IDEALITY_COLD
    };
    (p / (rho * fluid.data.gas_constant * t) - 1.0).abs() <= limit
}

/// Densest vapour searched at `t` K, mol/m³: under the critical temperature
/// the vapour branch ends below the critical density.
fn vapour_ceiling(fluid: &Fluid, t: f64) -> f64 {
    if t < fluid.data.t_crit {
        fluid.data.rho_crit
    } else {
        max_density(fluid)
    }
}

/// Density, mol/m³, of the vapour at `p` Pa and `t` K, metastable or not.
///
/// p(ρ) rises from zero along the vapour branch up to the spinodal, so the
/// root is bracketed by growing the density from a stable vapour state while
/// the slope stays positive (and, below the critical temperature, keeps
/// falling).
pub(super) fn vapour_density(fluid: &Fluid, p: f64, t: f64) -> anyhow::Result<f64> {
    vapour_root(fluid, p, t).map_err(|end| {
        PhysicsError::InvalidState(match end {
            NoVapour::Spinodal => "beyond the vapour spinodal",
            NoVapour::Untrusted => "outside the trusted vapour range",
        })
        .into()
    })
}

/// Why no vapour state exists at some (p, T).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum NoVapour {
    /// Past the vapour spinodal: the vapour would condense spontaneously.
    Spinodal,
    /// Outside the range the equation is trusted in (below the triple point),
    /// or outside its pressure range.
    Untrusted,
}

/// [`vapour_density`] with the reason when there is none.
pub(super) fn vapour_root(fluid: &Fluid, p: f64, t: f64) -> Result<f64, NoVapour> {
    // Below the triple point a failing slope test is the extrapolation
    // breaking down, not a spinodal.
    let spinodal = if t < fluid.data.t_triple {
        NoVapour::Untrusted
    } else {
        NoVapour::Spinodal
    };
    let ceiling = vapour_ceiling(fluid, t);
    let f = |rho: f64| molar(fluid, rho, t);
    // Under T_c every density below saturated vapour's is stable vapour;
    // from there the march checks each step's slope (the ancillary is good to
    // ~1 %, hence the margin).
    let mut start = (p / (fluid.data.gas_constant * t)).min(ceiling);
    if t < fluid.data.t_crit {
        start = start.min(0.99 * ancillary(&fluid.data.vapour_density, t));
    }
    let (mut lo, mut hi) = (0.0, start);
    // Under T_c the vapour branch is concave, its slope falling to zero at the
    // spinodal; a rising slope means a step crossed onto another branch (the
    // equation has spurious loops inside the dome).
    let subcritical = t < fluid.data.t_crit;
    let mut slope = f64::INFINITY;
    loop {
        let x = f(hi);
        // At low density the slope is flat to round-off (∂p/∂ρ → RT).
        if x.dp_drho <= 0.0 || (subcritical && x.dp_drho > slope * (1.0 + 1e-9)) {
            return Err(spinodal);
        }
        slope = x.dp_drho;
        if x.p >= p {
            break;
        }
        if !trusted(fluid, hi, t, x.p) {
            return Err(NoVapour::Untrusted);
        }
        if hi >= ceiling {
            return Err(spinodal);
        }
        lo = hi;
        hi = (hi * GROWTH).min(ceiling);
    }
    let options = RootOptions::new(1e-14 * hi, 200).map_err(|_| NoVapour::Untrusted)?;
    let rho = newton_bracketed(
        |rho| {
            if rho <= 0.0 {
                return (-p, f64::NAN);
            }
            let x = f(rho);
            (x.p - p, x.dp_drho)
        },
        lo,
        hi,
        hi,
        options,
    )
    .map_err(|_| NoVapour::Untrusted)?;
    if !trusted(fluid, rho, t, p) {
        return Err(NoVapour::Untrusted);
    }
    Ok(rho)
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::units::{
        ThermodynamicTemperature, mass_density::kilogram_per_cubic_meter, pressure::pascal,
        specific_heat_capacity::joule_per_kilogram_kelvin, thermodynamic_temperature::kelvin,
    };

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
    }

    fn pa(p: f64) -> Pressure {
        Pressure::new::<pascal>(p)
    }

    fn entropy(s: f64) -> SpecificHeatCapacity {
        SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(s)
    }

    fn si(state: &FlowState) -> (f64, f64, f64, f64, f64) {
        (
            state.temperature.get::<kelvin>(),
            state.density.get::<kilogram_per_cubic_meter>(),
            state
                .enthalpy
                .get::<crate::units::available_energy::joule_per_kilogram>(),
            state.entropy.get::<joule_per_kilogram_kelvin>(),
            state
                .speed_of_sound
                .get::<crate::units::velocity::meter_per_second>(),
        )
    }

    /// Saturated vapour entropy of `fluid` at `t` K, J/(kg·K).
    fn vapour_entropy(fluid: PureFluid, t: f64) -> anyhow::Result<(f64, f64)> {
        let sat = fluid.saturation(ThermodynamicTemperature::new::<kelvin>(t))?;
        let v = fluid.state_prho(sat.pressure, sat.vapour_density)?;
        Ok((
            sat.pressure.get::<pascal>(),
            v.entropy.get::<joule_per_kilogram_kelvin>(),
        ))
    }

    #[test]
    fn test_frozen_vapour_equals_equilibrium_outside_the_dome() -> anyhow::Result<()> {
        for (fluid, p, t) in [
            (PureFluid::Nitrogen, 2e7, 300.0),
            (PureFluid::Nitrogen, 1e5, 90.0),
            (PureFluid::Argon, 5e6, 200.0),
            (PureFluid::Helium, 3e7, 20.0),
            (PureFluid::CarbonDioxide, 1e6, 260.0),
        ] {
            let props = fluid.properties(pa(p), ThermodynamicTemperature::new::<kelvin>(t))?;
            let frozen = FrozenVapour::new(fluid);
            let by_s = frozen.state_ps(props.pressure, props.entropy)?;
            let by_rho = frozen.state_prho(props.pressure, props.density)?;
            let exact = si(&fluid.state_ps(props.pressure, props.entropy)?);
            for (label, got) in [("(p, s)", si(&by_s)), ("(p, ρ)", si(&by_rho))] {
                let pairs = [
                    (got.0, exact.0),
                    (got.1, exact.1),
                    (got.2, exact.2),
                    (got.4, exact.4),
                ];
                for (a, b) in pairs {
                    anyhow::ensure!(
                        (a - b).abs() <= 1e-9 * b.abs().max(1.0),
                        "{fluid:?} {label} at {p} Pa, {t} K: {got:?} vs {exact:?}"
                    );
                }
            }
            anyhow::ensure!(by_s.vapour_quality.is_none(), "{by_s:?}");
        }
        Ok(())
    }

    #[test]
    fn test_frozen_vapour_supersaturates_inside_the_dome() -> anyhow::Result<()> {
        // Saturated N2 vapour at 100 K expanded a little: equilibrium would
        // condense, frozen vapour cools below the saturation temperature.
        let fluid = PureFluid::Nitrogen;
        let (p_sat, s_v) = vapour_entropy(fluid, 100.0)?;
        let p = 0.7 * p_sat;
        let frozen = FrozenVapour::new(fluid).state_ps(pa(p), entropy(s_v))?;
        let equilibrium = fluid.state_ps(pa(p), entropy(s_v))?;
        anyhow::ensure!(equilibrium.vapour_quality.is_some(), "{equilibrium:?}");
        anyhow::ensure!(frozen.vapour_quality.is_none(), "{frozen:?}");
        let t = frozen.temperature.get::<kelvin>();
        let t_sat = equilibrium.temperature.get::<kelvin>();
        anyhow::ensure!(t < t_sat - 1.0, "frozen {t} K vs saturation {t_sat} K");
        // A single phase: its speed of sound is above the equilibrium
        // mixture's, which drops on entering the dome.
        anyhow::ensure!(
            frozen.speed_of_sound > equilibrium.speed_of_sound,
            "{frozen:?} vs {equilibrium:?}"
        );
        Ok(())
    }

    #[test]
    fn test_frozen_vapour_expands_below_the_triple_point() -> anyhow::Result<()> {
        // N2 from 20 bar and 300 K down to 2 kPa: the perfect-gas isentrope
        // ends near 37 K, well below the 63.15 K triple point.
        let fluid = PureFluid::Nitrogen;
        let start = fluid.properties(pa(2e6), ThermodynamicTemperature::new::<kelvin>(300.0))?;
        let state = FrozenVapour::new(fluid).state_ps(pa(2e3), start.entropy)?;
        let t = state.temperature.get::<kelvin>();
        anyhow::ensure!(t < 63.0 && t > 30.0, "T = {t} K");
        // Nearly ideal at this density: Z within 1 %.
        let r = 8.314_462_618 / 0.028_013_4;
        let z = 2e3 / (state.density.get::<kilogram_per_cubic_meter>() * r * t);
        anyhow::ensure!((z - 1.0).abs() < 0.01, "Z = {z}");
        Ok(())
    }

    #[test]
    fn test_frozen_co2_expands_below_its_triple_point_at_low_density() -> anyhow::Result<()> {
        // From 20 bar and 320 K to 1.5 bar: ~170 K, below the 216.6 K triple
        // point; the vapour is dilute (Z ≈ 0.95), its ∂p/∂ρ flat to
        // round-off at the start of the bracket.
        let fluid = PureFluid::CarbonDioxide;
        let start = fluid.properties(pa(2e6), ThermodynamicTemperature::new::<kelvin>(320.0))?;
        let state = FrozenVapour::new(fluid).state_ps(pa(1.5e5), start.entropy)?;
        let t = state.temperature.get::<kelvin>();
        anyhow::ensure!(t < 216.0 && t > 150.0, "T = {t} K");
        Ok(())
    }

    #[test]
    fn test_frozen_vapour_fails_beyond_the_spinodal() -> anyhow::Result<()> {
        // Liquid entropy at 100 K but at a low pressure: no vapour has it.
        let fluid = PureFluid::Nitrogen;
        let sat = fluid.saturation(ThermodynamicTemperature::new::<kelvin>(100.0))?;
        let liquid = fluid.state_prho(sat.pressure, sat.liquid_density)?;
        let result = FrozenVapour::new(fluid).state_ps(pa(1e5), liquid.entropy);
        anyhow::ensure!(result.is_err(), "{result:?}");
        // A liquid density is not a vapour state either.
        let result = FrozenVapour::new(fluid).state_prho(sat.pressure, sat.liquid_density);
        anyhow::ensure!(result.is_err(), "{result:?}");
        Ok(())
    }

    #[test]
    fn test_frozen_vapour_keeps_its_fluid() {
        assert_eq!(
            FrozenVapour::new(PureFluid::Argon).fluid(),
            PureFluid::Argon
        );
    }

    fn fluids() -> impl Strategy<Value = PureFluid> {
        prop_oneof![
            Just(PureFluid::Nitrogen),
            Just(PureFluid::Argon),
            Just(PureFluid::Helium),
            Just(PureFluid::CarbonDioxide),
        ]
    }

    proptest! {
        // Most random expansions end deep inside the dome, beyond what frozen
        // vapour covers: many cases are rejected.
        #![proptest_config(ProptestConfig {
            max_global_rejects: 50_000,
            ..ProptestConfig::with_cases(1000)
        })]

        /// From a gas at 0.2–2 times the critical pressure and 1.2–3 times the critical temperature, expanded by up to 10³ in
        /// pressure but staying above the triple point and mildly
        /// supersaturated at most: (p, s) and back through (p, ρ) give the
        /// same state, s rises with T along the isobar, and the state is
        /// mechanically stable.
        #[test]
        fn test_frozen_vapour_isentropes_are_consistent(
            fluid in fluids(),
            pr in 0.2f64..2.0,
            tr in 1.2f64..3.0,
            log_ratio in 0.0f64..3.0,
        ) {
            let data = fluid.data().data;
            let critical = fluid
                .homogeneous(
                    ThermodynamicTemperature::new::<kelvin>(data.t_crit),
                    MassDensity::new::<kilogram_per_cubic_meter>(data.rho_crit * data.molar_mass),
                )
                .map_err(|e| fail(&e))?;
            let p0 = pr * critical.pressure.get::<pascal>();
            let t0 = (tr * data.t_crit).min(1500.0);
            // Ends above the triple point even for a monatomic gas, T ∝ p^0.4;
            // below it the explicit test covers the trusted range.
            prop_assume!(t0 * 10f64.powf(-0.4 * log_ratio) > 1.1 * data.t_triple);
            let start = fluid
                .properties(pa(p0), ThermodynamicTemperature::new::<kelvin>(t0))
                .map_err(|e| fail(&e))?;
            // A vapour-like start: above the critical entropy, else the
            // isentrope enters the dome from the liquid side.
            prop_assume!(start.entropy > critical.entropy);
            let frozen = FrozenVapour::new(fluid);
            let p = p0 / 10f64.powf(log_ratio);
            // Mildly supersaturated at most: in equilibrium it would condense
            // no more than 10 %, below half the critical pressure (deeper, or
            // nearer the critical point, it reaches the vapour spinodal).
            let equilibrium = fluid.state_ps(pa(p), start.entropy).map_err(|e| fail(&e))?;
            let p_c = critical.pressure.get::<pascal>();
            prop_assume!(equilibrium.vapour_quality.is_none_or(|q| q > 0.9 && p < 0.5 * p_c));
            let state = frozen.state_ps(pa(p), start.entropy).map_err(|e| fail(&e))?;
            let (t, rho, h, s, w) = si(&state);
            let s0 = start.entropy.get::<joule_per_kilogram_kelvin>();
            prop_assert!((s - s0).abs() < 1e-9 * s0.abs(), "s {} vs {}", s, s0);
            prop_assert!(w > 0.0 && w.is_finite(), "w = {}", w);
            let back = frozen.state_prho(pa(p), state.density).map_err(|e| fail(&e))?;
            let (t2, rho2, h2, _, _) = si(&back);
            prop_assert!((t2 / t - 1.0).abs() < 1e-9, "T {} vs {}", t2, t);
            prop_assert!((rho2 / rho - 1.0).abs() < 1e-9, "ρ {} vs {}", rho2, rho);
            prop_assert!((h2 - h).abs() < 1e-9 * h.abs().max(1e4), "h {} vs {}", h2, h);
            // Entropy rises with temperature along the isobar.
            let warmer = frozen.state_ps(pa(p), entropy(s0 + 10.0)).map_err(|e| fail(&e))?;
            prop_assert!(warmer.temperature.get::<kelvin>() > t, "not monotone at {} Pa", p);
        }
    }
}
