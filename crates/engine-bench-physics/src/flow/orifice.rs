//! Compressible flow through an orifice between two reservoirs.
//!
//! Isentropic nozzle flow from the upstream stagnation state to the downstream
//! pressure, scaled by a discharge coefficient; choked below the critical
//! pressure ratio. Close to pressure equality the exact √Δp behaviour has an
//! infinite slope, which makes ODE systems stiff, so within [`LINEAR_BAND`] of
//! equality a C¹ quadratic in (1 − `p_down/p_up`) takes over.

use anyhow::Context;

use crate::{
    PhysicsError,
    flow::{PerfectGas, isentropic},
    units::{
        Area, MassRate, Pressure, ThermodynamicTemperature, area::square_meter,
        mass_rate::kilogram_per_second, pressure::pascal, thermodynamic_temperature::kelvin,
    },
};

/// Width of the near-equilibrium band, as 1 − `p_down/p_up`, where flow is
/// quadratic in Δp.
pub(crate) const LINEAR_BAND: f64 = 1e-3;

/// Stagnation state on one side of an orifice.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Reservoir {
    /// Pa.
    pub(crate) pressure: f64,
    /// K.
    pub(crate) temperature: f64,
    pub(crate) gas: PerfectGas,
}

impl Reservoir {
    /// A reservoir at positive, finite pressure and temperature.
    pub fn new(p: Pressure, t: ThermodynamicTemperature, gas: PerfectGas) -> anyhow::Result<Self> {
        let (pressure, temperature) = (p.get::<pascal>(), t.get::<kelvin>());
        if !(pressure > 0.0 && pressure.is_finite()) {
            return Err(PhysicsError::InvalidState("pressure must be positive"))
                .with_context(|| format!("reservoir p = {pressure} Pa"));
        }
        if !(temperature > 0.0 && temperature.is_finite()) {
            return Err(PhysicsError::InvalidState("temperature must be positive"))
                .with_context(|| format!("reservoir T = {temperature} K"));
        }
        Ok(Self {
            pressure,
            temperature,
            gas,
        })
    }
}

/// A sharp-edged or nozzle-like opening with flow area and discharge
/// coefficient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Orifice {
    /// Cd·A in m².
    effective_area: f64,
}

impl Orifice {
    /// An orifice of geometric `area` and discharge coefficient
    /// `discharge_coefficient` in (0, 1].
    pub fn new(area: Area, discharge_coefficient: f64) -> anyhow::Result<Self> {
        let area = area.get::<square_meter>();
        if !(area > 0.0 && area.is_finite()) {
            return Err(PhysicsError::InvalidGeometry(
                "orifice area must be positive",
            ))
            .with_context(|| format!("A = {area} m²"));
        }
        if !(discharge_coefficient > 0.0 && discharge_coefficient <= 1.0) {
            return Err(PhysicsError::InvalidGeometry(
                "discharge coefficient must be in (0, 1]",
            ))
            .with_context(|| format!("Cd = {discharge_coefficient}"));
        }
        Ok(Self {
            effective_area: discharge_coefficient * area,
        })
    }

    /// Mass flow from `a` to `b` (negative when it runs from `b` to `a`), using
    /// the upstream side's gas and temperature.
    pub fn mass_flow(&self, a: &Reservoir, b: &Reservoir) -> MassRate {
        MassRate::new::<kilogram_per_second>(self.mass_flow_si(a, b))
    }

    /// Kernel form of [`mass_flow`](Self::mass_flow), in kg/s.
    #[inline]
    pub(crate) fn mass_flow_si(self, a: &Reservoir, b: &Reservoir) -> f64 {
        if a.pressure >= b.pressure {
            self.downstream_flow(a, b.pressure)
        } else {
            -self.downstream_flow(b, a.pressure)
        }
    }

    /// Flow out of `up` into pressure `p_down ≤ up.pressure`, in kg/s.
    #[inline]
    fn downstream_flow(self, up: &Reservoir, p_down: f64) -> f64 {
        let gamma = up.gas.gamma;
        // ṁ = K·√φ(r), r = p_down/p_up (Saint-Venant–Wantzel).
        let k = self.effective_area
            * up.pressure
            * (2.0 * gamma / ((gamma - 1.0) * up.gas.gas_constant * up.temperature)).sqrt();
        let x = 1.0 - p_down / up.pressure;
        if x >= LINEAR_BAND {
            let r = (1.0 - x).max(isentropic::critical_pressure_ratio(gamma));
            return k * flow_function(gamma, r).sqrt();
        }
        // f(x) = a·x + b·x² matching value and slope of K·√φ at x =
        // LINEAR_BAND.
        let r0 = 1.0 - LINEAR_BAND;
        let phi0 = flow_function(gamma, r0);
        let edge = k * phi0.sqrt();
        let slope = -k * flow_function_slope(gamma, r0) / (2.0 * phi0.sqrt());
        let b = (slope * LINEAR_BAND - edge) / (LINEAR_BAND * LINEAR_BAND);
        let a = slope - 2.0 * b * LINEAR_BAND;
        x * (a + b * x)
    }
}

/// φ(r) = r^(2/γ) − r^((γ+1)/γ).
#[inline]
fn flow_function(gamma: f64, r: f64) -> f64 {
    r.powf(2.0 / gamma) - r.powf((gamma + 1.0) / gamma)
}

/// dφ/dr.
#[inline]
fn flow_function_slope(gamma: f64, r: f64) -> f64 {
    2.0 / gamma * r.powf(2.0 / gamma - 1.0) - (gamma + 1.0) / gamma * r.powf(1.0 / gamma)
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::units::{SpecificHeatCapacity, specific_heat_capacity::joule_per_kilogram_kelvin};

    const P0: f64 = 1.0e6;
    const T0: f64 = 300.0;

    fn air() -> anyhow::Result<PerfectGas> {
        PerfectGas::new(
            1.4,
            SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(287.0),
        )
    }

    fn reservoir(p: f64, t: f64) -> anyhow::Result<Reservoir> {
        Reservoir::new(
            Pressure::new::<pascal>(p),
            ThermodynamicTemperature::new::<kelvin>(t),
            air()?,
        )
    }

    fn orifice() -> anyhow::Result<Orifice> {
        Orifice::new(Area::new::<square_meter>(1e-4), 1.0)
    }

    /// Isentropic flow to back pressure `p` through the full area, via the exit
    /// Mach number: an independent route through
    /// `isentropic::mass_flux_parameter`.
    fn reference_flow(p: f64) -> f64 {
        let gamma = 1.4;
        let r = (p / P0).max(isentropic::critical_pressure_ratio(gamma));
        let mach = ((r.powf(-(gamma - 1.0) / gamma) - 1.0) * 2.0 / (gamma - 1.0)).sqrt();
        1e-4 * P0 / (287.0 * T0).sqrt() * isentropic::mass_flux_parameter(gamma, mach)
    }

    #[test]
    fn test_reservoir_validates_state() -> anyhow::Result<()> {
        for (p, t) in [
            (0.0, 300.0),
            (1e5, 0.0),
            (f64::NAN, 300.0),
            (1e5, f64::INFINITY),
        ] {
            anyhow::ensure!(reservoir(p, t).is_err(), "accepted p {p}, T {t}");
        }
        Ok(())
    }

    #[test]
    fn test_orifice_validates_geometry() -> anyhow::Result<()> {
        for (area, cd) in [(0.0, 0.6), (1e-4, 0.0), (1e-4, 1.1), (f64::NAN, 0.6)] {
            anyhow::ensure!(
                Orifice::new(Area::new::<square_meter>(area), cd).is_err(),
                "accepted A {area}, Cd {cd}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_orifice_choked_flow_matches_analytic_and_ignores_back_pressure() -> anyhow::Result<()> {
        // A·p0/√(R·T0)·0.6847 for γ = 1.4.
        let expected = 1e-4 * P0 / (287.0 * T0).sqrt() * isentropic::mass_flux_parameter(1.4, 1.0);
        let up = reservoir(P0, T0)?;
        for p_down in [1e3, 1e5, 0.5e6] {
            let flow = orifice()?.mass_flow(&up, &reservoir(p_down, T0)?);
            let got = flow.get::<kilogram_per_second>();
            anyhow::ensure!(
                (got / expected - 1.0).abs() < 1e-12,
                "p_down {p_down}: {got}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_orifice_unchoked_flow_matches_isentropic_reference() -> anyhow::Result<()> {
        let up = reservoir(P0, T0)?;
        for p_down in [0.6e6, 0.8e6, 0.95e6] {
            let got = orifice()?.mass_flow_si(&up, &reservoir(p_down, T0)?);
            let want = reference_flow(p_down);
            anyhow::ensure!(
                (got / want - 1.0).abs() < 1e-12,
                "p_down {p_down}: {got} vs {want}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_orifice_flow_reverses_with_pressure_difference() -> anyhow::Result<()> {
        let (a, b) = (reservoir(P0, T0)?, reservoir(0.7e6, T0)?);
        let o = orifice()?;
        let forward = o.mass_flow_si(&a, &b);
        anyhow::ensure!(forward > 0.0, "forward {forward}");
        anyhow::ensure!(
            (o.mass_flow_si(&b, &a) + forward).abs() < 1e-15,
            "not antisymmetric"
        );
        anyhow::ensure!(o.mass_flow_si(&a, &a) == 0.0, "flow at equal pressure");
        Ok(())
    }

    #[test]
    fn test_orifice_slope_is_finite_at_equal_pressure() -> anyhow::Result<()> {
        let o = orifice()?;
        let up = reservoir(P0, T0)?;
        let dp = 1e-3;
        let slope = o.mass_flow_si(&up, &reservoir(P0 - dp, T0)?) / dp;
        let wider = o.mass_flow_si(&up, &reservoir(P0 - 10.0 * dp, T0)?) / (10.0 * dp);
        // Quadratic band: slope at Δp → 0 is the band's linear coefficient, not
        // ∞.
        anyhow::ensure!(
            slope.is_finite() && (slope / wider - 1.0).abs() < 1e-3,
            "{slope} {wider}"
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_orifice_flow_decreases_with_back_pressure(
            a in 0.0..1.0f64, b in 0.0..1.0f64,
        ) {
            // Ordering instead of rejecting: rejecting half the cases trips proptest's
            // global reject limit.
            let (lo, hi) = (a.min(b), a.max(b));
            let up = reservoir(P0, T0).map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let o = orifice().map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let flow = |r: f64| -> Result<f64, TestCaseError> {
                let down = reservoir((r * P0).max(1.0), T0)
                    .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
                Ok(o.mass_flow_si(&up, &down))
            };
            prop_assert!(flow(lo)? >= flow(hi)?, "not monotone between {lo} and {hi}");
        }

        #[test]
        fn test_orifice_is_c1_across_the_linear_band_edge(gamma in 1.05..1.67f64) {
            let gas = PerfectGas::new(gamma, SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(300.0))
                .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let side = |p: f64| Reservoir {pressure: p, temperature: T0, gas};
            let o = orifice().map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let up = side(P0);
            let flow = |x: f64| o.mass_flow_si(&up, &side(P0 * (1.0 - x)));
            let (x0, h) = (LINEAR_BAND, 1e-7);
            let (inside, edge, outside) = (flow(x0 - h), flow(x0), flow(x0 + h));
            prop_assert!((f64::midpoint(inside, outside) / edge - 1.0).abs() < 1e-6, "value jump");
            let slope_in = (edge - inside) / h;
            let slope_out = (outside - edge) / h;
            prop_assert!((slope_in / slope_out - 1.0).abs() < 1e-3,
                "slope {slope_in} vs {slope_out}");
        }
    }
}
