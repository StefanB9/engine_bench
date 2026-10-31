//! Tabulated reference fluids: real-fluid flow states fast enough for
//! real-time simulation, interpolated from tables built once per fluid.
//!
//! - Two branch tables over (ln p, s), vapour and liquid, each continued
//!   metastably across the saturation line so that no stencil straddles a kink;
//!   above the critical pressure they coincide. Nodes hold ln T, ln ρ, h and ln
//!   w; cubic Lagrange interpolation.
//! - The saturation line as a table in θ = (1 − T/Tc)^(1/3); inside the dome
//!   the equilibrium state is mixed exactly from the saturated phases.
//! - Equilibrium and frozen-vapour modes from the same tables.
//! - Wherever a stencil touches an invalid node, leaves the tables, or the
//!   phase is uncertain (the dome's ends), and near the critical point, the
//!   exact equation answers.

mod branch;
mod build;
mod grid;
mod saturation;

use std::{fmt, sync::OnceLock};

use anyhow::Context;

use self::{
    build::Envelope,
    grid::Grid2,
    saturation::{Saturated, SaturationTable},
};
use super::{FrozenVapour, PureFluid};
use crate::{
    PhysicsError,
    flow::real::{FlowState, RealFluid},
    units::{
        MassDensity, Pressure, SpecificHeatCapacity, mass_density::kilogram_per_cubic_meter,
        pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
    },
};

/// Node quantities.
const LN_T: usize = 0;
const LN_RHO: usize = 1;
const H: usize = 2;
const LN_W: usize = 3;
/// Newton iterations for (p, ρ) through the table.
const DENSITY_ITERATIONS: usize = 30;
/// Around the critical point properties change too fast for the tables (the
/// speed of sound by up to ~1e-2): states within these fractions of the
/// critical temperature and density come from the exact equation.
const CRITICAL_BAND_T: f64 = 0.05;
const CRITICAL_BAND_RHO: f64 = 0.5;

/// One fluid's tables.
pub(super) struct Tables {
    fluid: PureFluid,
    vapour: Grid2<4>,
    liquid: Grid2<4>,
    saturation: SaturationTable,
    /// Pressures of the dome, Pa: triple point to its top.
    dome: (f64, f64),
    /// Entropy, J/(kg·K), and density, kg/m³, at the critical point.
    s_crit: f64,
    rho_crit: f64,
    /// Lowest temperature of an equilibrium state, K.
    t_triple: f64,
    /// Critical temperature, K.
    t_crit: f64,
    /// Per vapour row, the entropy where the vapour branch ends (NaN if it
    /// does not end within the row), J/(kg·K).
    vapour_end: Vec<f64>,
    /// A gas state for ideal-gas estimates: s J/(kg·K) at 300 K and 1 bar,
    /// cp J/(kg·K) there, and R J/(kg·K).
    reference: (f64, f64, f64),
}

/// What the tables say about a state.
enum Lookup {
    /// The state.
    Hit(FlowState),
    /// Ask the exact equation.
    Exact,
    /// There is no such state (the exact equation would fail too).
    Reject(&'static str),
}

/// A reference fluid evaluated through its tables: the same states as the
/// exact equation to ~1e-6 (see the tests), at a fraction of the cost.
#[derive(Clone, Copy)]
pub struct TabulatedFluid {
    tables: &'static Tables,
    frozen: bool,
}

static CACHE: [OnceLock<Tables>; 4] = [const { OnceLock::new() }; 4];

impl TabulatedFluid {
    /// `fluid` in equilibrium (liquid and vapour inside the dome). The first
    /// call per fluid and process builds its tables (about a second).
    pub fn equilibrium(fluid: PureFluid) -> anyhow::Result<Self> {
        Ok(Self {
            tables: cached(fluid)?,
            frozen: false,
        })
    }

    /// `fluid` as frozen vapour, as [`FrozenVapour`].
    pub fn frozen(fluid: PureFluid) -> anyhow::Result<Self> {
        Ok(Self {
            tables: cached(fluid)?,
            frozen: true,
        })
    }

    /// The fluid.
    pub fn fluid(self) -> PureFluid {
        self.tables.fluid
    }

    /// Whether this is the frozen-vapour mode.
    pub fn is_frozen(self) -> bool {
        self.frozen
    }

    /// What the tables say about (p, s).
    fn lookup_ps(&self, p: f64, s: f64) -> Lookup {
        let t = self.tables;
        let below_triple_pressure = !self.frozen && p < t.dome.0;
        if let Some(state) = self.table_ps(p, s) {
            if below_triple_pressure && raw_t(&state) < t.t_triple * (1.0 - 1e-6) {
                return Lookup::Reject("colder than the triple point (solid)");
            }
            return if t.away_from_critical(&state) {
                Lookup::Hit(state)
            } else {
                Lookup::Exact
            };
        }
        if (self.frozen || below_triple_pressure) && t.past_vapour_end(p, s) {
            return Lookup::Reject("no vapour state has this entropy");
        }
        Lookup::Exact
    }

    /// The table answer for (p, s), or `None` where the exact equation must
    /// answer.
    fn table_ps(&self, p: f64, s: f64) -> Option<FlowState> {
        let t = self.tables;
        if self.frozen {
            return single(&t.vapour, p, s);
        }
        let (p_triple, p_top) = t.dome;
        if p < p_triple {
            // Below the triple pressure: vapour, unless colder than the triple
            // point (solid, which the exact equation reports).
            return single(&t.vapour, p, s);
        }
        if p <= p_top {
            let sat = t.saturation.at_pressure(p)?;
            return if s >= sat.s_v {
                single(&t.vapour, p, s)
            } else if s <= sat.s_l {
                single(&t.liquid, p, s)
            } else {
                Some(mixture(&sat, p, s))
            };
        }
        let (first, second) = if s >= t.s_crit {
            (&t.vapour, &t.liquid)
        } else {
            (&t.liquid, &t.vapour)
        };
        single(first, p, s).or_else(|| single(second, p, s))
    }

    /// The table answer for (p, ρ), or `None` where the exact equation must
    /// answer: single-phase states only, as [`RealFluid::state_prho`].
    fn table_prho(&self, p: f64, rho: f64) -> Option<FlowState> {
        let t = self.tables;
        if self.frozen {
            return by_density(&t.vapour, p, rho, t.gas_entropy(p, rho));
        }
        let (p_triple, p_top) = t.dome;
        if p < p_triple {
            return by_density(&t.vapour, p, rho, t.gas_entropy(p, rho))
                .filter(|x| raw_t(x) >= t.t_triple);
        }
        if p <= p_top {
            let sat = t.saturation.at_pressure(p)?;
            let v = 1.0 / rho;
            return if v >= sat.v_v {
                by_density(&t.vapour, p, rho, t.gas_entropy(p, rho).max(sat.s_v))
            } else if v <= sat.v_l {
                by_density(&t.liquid, p, rho, sat.s_l)
            } else {
                None
            };
        }
        let (first, second) = if rho < t.rho_crit {
            (&t.vapour, &t.liquid)
        } else {
            (&t.liquid, &t.vapour)
        };
        let start = if rho < t.rho_crit {
            t.gas_entropy(p, rho)
        } else {
            t.s_crit
        };
        by_density(first, p, rho, start).or_else(|| by_density(second, p, rho, start))
    }
}

impl RealFluid for TabulatedFluid {
    fn state_ps(&self, p: Pressure, s: SpecificHeatCapacity) -> anyhow::Result<FlowState> {
        let (p_si, s_si) = (p.get::<pascal>(), s.get::<joule_per_kilogram_kelvin>());
        match self.lookup_ps(p_si, s_si) {
            Lookup::Hit(state) => return Ok(state),
            Lookup::Reject(why) => {
                return Err(PhysicsError::InvalidState(why)).with_context(|| {
                    format!("{:?} at {p_si} Pa, s = {s_si} J/(kg K)", self.tables.fluid)
                });
            }
            Lookup::Exact => {}
        }
        if self.frozen {
            FrozenVapour::new(self.tables.fluid).state_ps(p, s)
        } else {
            self.tables.fluid.state_ps(p, s)
        }
        .context("tabulated fluid, exact fallback")
    }

    fn state_prho(&self, p: Pressure, rho: MassDensity) -> anyhow::Result<FlowState> {
        let (p_si, rho_si) = (p.get::<pascal>(), rho.get::<kilogram_per_cubic_meter>());
        if let Some(state) = self
            .table_prho(p_si, rho_si)
            .filter(|x| self.tables.away_from_critical(x))
        {
            return Ok(state);
        }
        if self.frozen {
            FrozenVapour::new(self.tables.fluid).state_prho(p, rho)
        } else {
            self.tables.fluid.state_prho(p, rho)
        }
        .context("tabulated fluid, exact fallback")
    }
}

impl Tables {
    /// Whether (p, s) lies clearly beyond where the vapour branch ends: below
    /// the end entropy of both neighbouring rows by two steps.
    fn past_vapour_end(&self, p: f64, s: f64) -> bool {
        let Some(i) = self.vapour.x.below(p.ln()) else {
            return false;
        };
        let ends = (self.vapour_end.get(i), self.vapour_end.get(i + 1));
        let (Some(&low), Some(&high)) = ends else {
            return false;
        };
        low.is_finite() && high.is_finite() && s < low.min(high) - 2.0 * self.vapour.y.step()
    }

    /// Entropy, J/(kg·K), of the ideal gas at `p` Pa and `rho` kg/m³: a start
    /// for solving (p, ρ) through the tables.
    fn gas_entropy(&self, p: f64, rho: f64) -> f64 {
        let (s_ref, cp, r) = self.reference;
        let t = p / (rho * r);
        s_ref + cp * (t / 300.0).ln() - r * (p / 1e5).ln()
    }

    /// Whether `state` lies outside the band around the critical point.
    fn away_from_critical(&self, state: &FlowState) -> bool {
        let rho = state.density.get::<kilogram_per_cubic_meter>();
        (raw_t(state) / self.t_crit - 1.0).abs() >= CRITICAL_BAND_T
            || (rho / self.rho_crit - 1.0).abs() >= CRITICAL_BAND_RHO
    }
}

impl fmt::Debug for TabulatedFluid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TabulatedFluid")
            .field("fluid", &self.tables.fluid)
            .field("frozen", &self.frozen)
            .finish_non_exhaustive()
    }
}

impl PartialEq for TabulatedFluid {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.tables, other.tables) && self.frozen == other.frozen
    }
}

/// The fluid's tables over the default envelope, built on first use.
fn cached(fluid: PureFluid) -> anyhow::Result<&'static Tables> {
    let index = PureFluid::ALL
        .iter()
        .position(|&f| f == fluid)
        .context("unknown fluid")?;
    let cell = CACHE.get(index).context("table cache")?;
    if let Some(tables) = cell.get() {
        return Ok(tables);
    }
    let tables = build::tables(fluid, &Envelope::DEFAULT)
        .with_context(|| format!("building the {fluid:?} tables"))?;
    // Another thread may have built them meanwhile; either copy is the same.
    let _ = cell.set(tables);
    cell.get().context("table cache")
}

fn raw_t(state: &FlowState) -> f64 {
    state
        .temperature
        .get::<crate::units::thermodynamic_temperature::kelvin>()
}

/// A single-phase state from `grid` at `p` Pa and `s` J/(kg·K).
fn single(grid: &Grid2<4>, p: f64, s: f64) -> Option<FlowState> {
    let x = grid.sample(p.ln(), s)?.value;
    Some(FlowState::from_si(
        p,
        x[LN_T].exp(),
        x[LN_RHO].exp(),
        x[H],
        s,
        x[LN_W].exp(),
        f64::NAN,
    ))
}

/// The equilibrium mixture at `p` Pa and `s` J/(kg·K) inside the dome.
fn mixture(sat: &Saturated, p: f64, s: f64) -> FlowState {
    let x = sat.mixture(s);
    FlowState::from_si(p, x.t, x.rho, x.h, s, x.w, x.quality)
}

/// The state of `grid` at `p` Pa with density `rho` kg/m³: Newton on s from
/// `s_start` (ln ρ falls as s rises along the isobar).
fn by_density(grid: &Grid2<4>, p: f64, rho: f64, s_start: f64) -> Option<FlowState> {
    let (ln_p, ln_rho) = (p.ln(), rho.ln());
    let mut s = s_start;
    for _ in 0..DENSITY_ITERATIONS {
        let x = grid.sample(ln_p, s)?;
        let slope = x.dy[LN_RHO];
        if slope.is_nan() || slope >= 0.0 {
            return None;
        }
        let step = (x.value[LN_RHO] - ln_rho) / slope;
        s -= step;
        if step.abs() <= 1e-12 * s.abs() {
            let state = single(grid, p, s)?;
            return Some(state);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::units::{
        ThermodynamicTemperature, available_energy::joule_per_kilogram,
        thermodynamic_temperature::kelvin, velocity::meter_per_second,
    };

    /// Nitrogen from 0.5 bar to 100 bar and 70 K to 400 K, at the default
    /// node spacing.
    const N2_ENVELOPE: Envelope = Envelope {
        p_min: 5e4,
        p_max: 1e7,
        t_min: 70.0,
        t_max: 400.0,
    };

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
    }

    /// The N2 tables over [`N2_ENVELOPE`], built once per test process.
    fn nitrogen(frozen: bool) -> anyhow::Result<TabulatedFluid> {
        static TABLES: OnceLock<Tables> = OnceLock::new();
        let tables = if let Some(t) = TABLES.get() {
            t
        } else {
            let built = build::tables(PureFluid::Nitrogen, &N2_ENVELOPE)?;
            let _ = TABLES.set(built);
            TABLES.get().context("tables")?
        };
        Ok(TabulatedFluid { tables, frozen })
    }

    /// Relative differences in T, ρ, h (against 10⁵ J/kg) and w.
    fn difference(a: &FlowState, b: &FlowState) -> [f64; 4] {
        let r = |x: f64, y: f64| (x / y - 1.0).abs();
        [
            r(raw_t(a), raw_t(b)),
            r(
                a.density.get::<kilogram_per_cubic_meter>(),
                b.density.get::<kilogram_per_cubic_meter>(),
            ),
            (a.enthalpy.get::<joule_per_kilogram>() - b.enthalpy.get::<joule_per_kilogram>()).abs()
                / 1e5,
            r(
                a.speed_of_sound.get::<meter_per_second>(),
                b.speed_of_sound.get::<meter_per_second>(),
            ),
        ]
    }

    fn pa(p: f64) -> Pressure {
        Pressure::new::<pascal>(p)
    }

    #[test]
    fn test_tabulated_fluid_reports_its_fluid_and_mode() -> anyhow::Result<()> {
        let (equilibrium, frozen) = (nitrogen(false)?, nitrogen(true)?);
        anyhow::ensure!(equilibrium.fluid() == PureFluid::Nitrogen, "fluid");
        anyhow::ensure!(!equilibrium.is_frozen() && frozen.is_frozen(), "modes");
        anyhow::ensure!(
            equilibrium != frozen && equilibrium == nitrogen(false)?,
            "equality"
        );
        let text = format!("{frozen:?}");
        anyhow::ensure!(
            text.contains("Nitrogen") && text.contains("frozen: true"),
            "{text}"
        );
        Ok(())
    }

    #[test]
    fn test_tabulated_fluid_mixes_inside_the_dome() -> anyhow::Result<()> {
        let tab = nitrogen(false)?;
        let sat = PureFluid::Nitrogen.saturation(ThermodynamicTemperature::new::<kelvin>(90.0))?;
        let l = PureFluid::Nitrogen.state_prho(sat.pressure, sat.liquid_density)?;
        let v = PureFluid::Nitrogen.state_prho(sat.pressure, sat.vapour_density)?;
        let s = l.entropy + (v.entropy - l.entropy) * 0.3;
        let (got, exact) = (
            tab.state_ps(sat.pressure, s)?,
            PureFluid::Nitrogen.state_ps(sat.pressure, s)?,
        );
        anyhow::ensure!(
            got.vapour_quality
                .zip(exact.vapour_quality)
                .is_some_and(|(a, b)| (a - b).abs() < 1e-7),
            "{got:?} vs {exact:?}"
        );
        let d = difference(&got, &exact);
        anyhow::ensure!(d[0] < 1e-8 && d[1] < 1e-7 && d[3] < 1e-4, "{d:?}");
        Ok(())
    }

    #[test]
    fn test_tabulated_fluid_falls_back_outside_its_tables() -> anyhow::Result<()> {
        let tab = nitrogen(false)?;
        // Beyond the envelope (300 bar, 1 kPa) and near the critical point.
        for (p, t) in [(3e7, 300.0), (1e3, 300.0), (3.5e6, 127.0)] {
            let x = PureFluid::Nitrogen
                .properties(pa(p), ThermodynamicTemperature::new::<kelvin>(t))?;
            let got = tab.state_ps(x.pressure, x.entropy)?;
            let exact = PureFluid::Nitrogen.state_ps(x.pressure, x.entropy)?;
            anyhow::ensure!(
                difference(&got, &exact).iter().all(|d| *d < 1e-12),
                "{p} Pa, {t} K"
            );
            let back = tab.state_prho(x.pressure, x.density)?;
            anyhow::ensure!(
                (raw_t(&back) / t - 1.0).abs() < 1e-9,
                "(p, ρ) at {p} Pa, {t} K"
            );
        }
        // Colder than the triple point: the exact equation's error.
        let x = PureFluid::Nitrogen
            .properties(pa(1e7), ThermodynamicTemperature::new::<kelvin>(300.0))?;
        anyhow::ensure!(
            tab.state_ps(pa(1e3), x.entropy).is_err(),
            "below the triple point"
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        /// Single-phase states and equilibrium mixtures over the envelope,
        /// away from the critical point: the table agrees with the exact
        /// equation, through (p, s) and back through (p, ρ).
        #[test]
        fn test_tabulated_fluid_matches_the_exact_equation(
            log_p in 4.8f64..7.0,
            t in 75.0f64..390.0,
            quality in 0.0f64..1.0,
        ) {
            prop_assume!((t / 126.19 - 1.0).abs() > 0.1);
            let fluid = PureFluid::Nitrogen;
            let tab = nitrogen(false).map_err(|e| fail(&e))?;
            let p = 10f64.powf(log_p);
            let exact = if let Some(sat) = fluid.saturation(ThermodynamicTemperature::new::<kelvin>(t)).ok().filter(|_| quality < 0.3) {
                // A mixture at the saturation pressure of T.
                let l = fluid.state_prho(sat.pressure, sat.liquid_density).map_err(|e| fail(&e))?;
                let v = fluid.state_prho(sat.pressure, sat.vapour_density).map_err(|e| fail(&e))?;
                let s = l.entropy + (v.entropy - l.entropy) * (quality / 0.3);
                fluid.state_ps(sat.pressure, s)
            } else {
                fluid.properties(pa(p), ThermodynamicTemperature::new::<kelvin>(t))
                    .and_then(|x| fluid.state_ps(x.pressure, x.entropy))
            }.map_err(|e| fail(&e))?;
            prop_assume!(exact.pressure.get::<pascal>() > 5.5e4);
            let got = tab.state_ps(exact.pressure, exact.entropy).map_err(|e| fail(&e))?;
            let d = difference(&got, &exact);
            prop_assert!(d[0] < 5e-6 && d[1] < 1e-5 && d[2] < 5e-6 && d[3] < 5e-5, "{:?}: {:?} vs {:?}", d, got, exact);
            if exact.vapour_quality.is_none() {
                let back = tab.state_prho(exact.pressure, exact.density).map_err(|e| fail(&e))?;
                let d = difference(&back, &exact);
                prop_assert!(d[0] < 5e-6 && d[2] < 5e-6 && d[3] < 5e-5, "(p, ρ) {:?}", d);
            }
        }

        /// Frozen expansions of N2 from 100 bar (as far as the envelope goes)
        /// agree with the exact frozen vapour.
        #[test]
        fn test_tabulated_frozen_vapour_matches_the_exact_one(
            t0 in 160.0f64..390.0,
            log_ratio in 0.0f64..2.2,
        ) {
            let fluid = PureFluid::Nitrogen;
            let tab = nitrogen(true).map_err(|e| fail(&e))?;
            let start = fluid.properties(pa(1e7), ThermodynamicTemperature::new::<kelvin>(t0)).map_err(|e| fail(&e))?;
            let p = pa(1e7 / 10f64.powf(log_ratio));
            let exact = FrozenVapour::new(fluid).state_ps(p, start.entropy);
            prop_assume!(exact.is_ok());
            let exact = exact.map_err(|e| fail(&e))?;
            let got = tab.state_ps(p, start.entropy).map_err(|e| fail(&e))?;
            let d = difference(&got, &exact);
            // Most states agree to ~1e-7. Deep in the metastable range, near
            // the vapour spinodal where w vanishes, the fields steepen: over
            // 5000 samples the worst were T 7e-6, ρ 2.4e-5, w 2.3e-4.
            prop_assert!(d[0] < 2e-5 && d[1] < 5e-5 && d[2] < 5e-6 && d[3] < 5e-4, "{:?}", d);
            prop_assert!(got.vapour_quality.is_none(), "{:?}", got);
        }
    }
}
