//! Quasi-1D converging–diverging nozzle operating point for every back-pressure
//! regime.
//!
//! `mod.rs` classifies the regime and assembles the exit state;
//! `shock_position` solves for a shock inside the divergent section.

mod shock_position;

use anyhow::Context;

use crate::{
    PhysicsError,
    flow::{
        Reservoir,
        isentropic::{self, Branch},
        shock,
    },
    units::{
        Area, Force, MassRate, Pressure, ThermodynamicTemperature, Time, Velocity,
        area::square_meter, force::newton, mass_rate::kilogram_per_second, pressure::pascal,
        thermodynamic_temperature::kelvin, time::second, velocity::meter_per_second,
    },
};

/// Standard gravity for specific impulse, m/s².
const STANDARD_GRAVITY: f64 = 9.806_65;

/// Back-pressure ratios within this distance (relative to chamber pressure) of
/// a regime boundary are reported as the boundary regime.
pub(crate) const BOUNDARY_TOLERANCE: f64 = 1e-9;

/// Summerfield criterion: an over-expanded jet separates from the wall once the
/// exit pressure drops below this fraction of the back pressure.
pub(crate) const SEPARATION_PRESSURE_RATIO: f64 = 0.4;

/// Flow regime set by the back pressure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NozzleRegime {
    /// Back pressure at or above chamber pressure.
    NoFlow,
    /// Subsonic throughout; the throat is not choked.
    Subsonic,
    /// Choked, with a normal shock in the divergent section.
    ShockInNozzle {
        /// Shock-plane area over throat area.
        shock_area_ratio: f64,
    },
    /// Choked, normal shock standing in the exit plane.
    ShockAtExit,
    /// Supersonic exit below back pressure (oblique shocks outside).
    OverExpanded,
    /// Supersonic exit at exactly back pressure.
    Ideal,
    /// Supersonic exit above back pressure (expansion fans outside).
    UnderExpanded,
}

/// Throat and exit areas of a converging–diverging (or converging, when equal)
/// nozzle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Nozzle {
    /// m².
    pub(crate) throat_area: f64,
    /// m².
    pub(crate) exit_area: f64,
}

/// Steady flow through a nozzle for given chamber and back pressure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NozzleOperatingPoint {
    /// Flow regime.
    pub regime: NozzleRegime,
    /// Mass flow.
    pub mass_flow: MassRate,
    /// Exit-plane Mach number (after the shock in
    /// [`NozzleRegime::ShockInNozzle`]).
    pub exit_mach: f64,
    /// Exit-plane static pressure.
    pub exit_pressure: Pressure,
    /// Exit-plane static temperature.
    pub exit_temperature: ThermodynamicTemperature,
    /// Exit-plane velocity.
    pub exit_velocity: Velocity,
    /// Thrust ṁ·vₑ + (pₑ − `p_b`)·Aₑ.
    pub thrust: Force,
    /// Over-expanded enough that the jet separates from the wall (Summerfield
    /// criterion); real thrust then differs from the attached-flow value.
    pub flow_separation: bool,
}

impl NozzleOperatingPoint {
    /// Thrust per weight flow F/(ṁ·g₀); `None` without flow.
    pub fn specific_impulse(&self) -> Option<Time> {
        let mass_flow = self.mass_flow.get::<kilogram_per_second>();
        (mass_flow > 0.0).then(|| {
            Time::new::<second>(self.thrust.get::<newton>() / (mass_flow * STANDARD_GRAVITY))
        })
    }
}

impl Nozzle {
    /// A nozzle with positive `throat` area and `exit` ≥ `throat`.
    pub fn new(throat: Area, exit: Area) -> anyhow::Result<Self> {
        let (throat_area, exit_area) = (throat.get::<square_meter>(), exit.get::<square_meter>());
        if !(throat_area > 0.0 && throat_area.is_finite()) {
            return Err(PhysicsError::InvalidGeometry(
                "throat area must be positive",
            ))
            .with_context(|| format!("A* = {throat_area} m²"));
        }
        if !(exit_area >= throat_area && exit_area.is_finite()) {
            return Err(PhysicsError::InvalidGeometry("exit area below throat area"))
                .with_context(|| format!("A* = {throat_area} m², Ae = {exit_area} m²"));
        }
        Ok(Self {
            throat_area,
            exit_area,
        })
    }

    /// Exit area over throat area, ε.
    pub fn expansion_ratio(&self) -> f64 {
        self.exit_area / self.throat_area
    }

    /// Steady operating point for the `chamber` stagnation state discharging
    /// into `back_pressure` (≥ 0; 0 is vacuum).
    pub fn operating_point(
        &self,
        chamber: &Reservoir,
        back_pressure: Pressure,
    ) -> anyhow::Result<NozzleOperatingPoint> {
        let pb = back_pressure.get::<pascal>();
        if !(pb >= 0.0 && pb.is_finite()) {
            return Err(PhysicsError::InvalidState(
                "back pressure must be non-negative",
            ))
            .with_context(|| format!("p_b = {pb} Pa"));
        }
        let (p0, t0) = (chamber.pressure, chamber.temperature);
        let gamma = chamber.gas.gamma;
        let ratio = pb / p0;
        if ratio >= 1.0 {
            return Ok(Self::no_flow(p0, t0));
        }
        // ṁ = A·flux_scale·mass_flux_parameter(M).
        let flux_scale = p0 / (chamber.gas.gas_constant * t0).sqrt();
        let epsilon = self.expansion_ratio();
        let subsonic_exit = isentropic::mach_from_area_ratio(gamma, epsilon, Branch::Subsonic)
            .context("subsonic design exit")?;
        let supersonic_exit = isentropic::mach_from_area_ratio(gamma, epsilon, Branch::Supersonic)
            .context("supersonic design exit")?;
        let exit_shock = shock::normal_shock(gamma, supersonic_exit);
        let subsonic_limit = isentropic::pressure_ratio(gamma, subsonic_exit);
        let design = isentropic::pressure_ratio(gamma, supersonic_exit);
        let shock_at_exit = design * exit_shock.pressure_ratio;
        let choked = self.throat_area * flux_scale * isentropic::mass_flux_parameter(gamma, 1.0);
        let near = |boundary: f64| (ratio - boundary).abs() <= BOUNDARY_TOLERANCE;

        let (regime, exit_mach, exit_ratio, mass_flow) = if ratio > subsonic_limit {
            let mach = (2.0 / (gamma - 1.0) * (ratio.powf(-(gamma - 1.0) / gamma) - 1.0)).sqrt();
            let flow = self.exit_area * flux_scale * isentropic::mass_flux_parameter(gamma, mach);
            (NozzleRegime::Subsonic, mach, ratio, flow)
        } else if near(design) {
            (NozzleRegime::Ideal, supersonic_exit, design, choked)
        } else if near(shock_at_exit) {
            let mach = exit_shock.downstream_mach;
            (NozzleRegime::ShockAtExit, mach, shock_at_exit, choked)
        } else if ratio > shock_at_exit {
            let (exit_mach, shock_area_ratio) =
                shock_position::shock_in_nozzle(gamma, epsilon, ratio, supersonic_exit)?;
            let regime = NozzleRegime::ShockInNozzle { shock_area_ratio };
            (regime, exit_mach, ratio, choked)
        } else if ratio > design {
            (NozzleRegime::OverExpanded, supersonic_exit, design, choked)
        } else {
            (NozzleRegime::UnderExpanded, supersonic_exit, design, choked)
        };

        let exit_temperature = t0 * isentropic::temperature_ratio(gamma, exit_mach);
        let exit_velocity =
            exit_mach * (gamma * chamber.gas.gas_constant * exit_temperature).sqrt();
        let exit_pressure = exit_ratio * p0;
        let thrust = mass_flow * exit_velocity + (exit_pressure - pb) * self.exit_area;
        Ok(NozzleOperatingPoint {
            regime,
            mass_flow: MassRate::new::<kilogram_per_second>(mass_flow),
            exit_mach,
            exit_pressure: Pressure::new::<pascal>(exit_pressure),
            exit_temperature: ThermodynamicTemperature::new::<kelvin>(exit_temperature),
            exit_velocity: Velocity::new::<meter_per_second>(exit_velocity),
            thrust: Force::new::<newton>(thrust),
            flow_separation: regime == NozzleRegime::OverExpanded
                && exit_pressure < SEPARATION_PRESSURE_RATIO * pb,
        })
    }

    pub(crate) fn no_flow(p0: f64, t0: f64) -> NozzleOperatingPoint {
        NozzleOperatingPoint {
            regime: NozzleRegime::NoFlow,
            mass_flow: MassRate::new::<kilogram_per_second>(0.0),
            exit_mach: 0.0,
            exit_pressure: Pressure::new::<pascal>(p0),
            exit_temperature: ThermodynamicTemperature::new::<kelvin>(t0),
            exit_velocity: Velocity::new::<meter_per_second>(0.0),
            thrust: Force::new::<newton>(0.0),
            flow_separation: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::{
        flow::PerfectGas,
        units::{SpecificHeatCapacity, specific_heat_capacity::joule_per_kilogram_kelvin},
    };

    const P0: f64 = 1.0e6;
    const T0: f64 = 300.0;
    const AT: f64 = 1e-4;

    fn chamber(gamma: f64) -> anyhow::Result<Reservoir> {
        let gas = PerfectGas::new(
            gamma,
            SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(287.0),
        )?;
        Reservoir::new(
            Pressure::new::<pascal>(P0),
            ThermodynamicTemperature::new::<kelvin>(T0),
            gas,
        )
    }

    fn nozzle(epsilon: f64) -> anyhow::Result<Nozzle> {
        Nozzle::new(
            Area::new::<square_meter>(AT),
            Area::new::<square_meter>(AT * epsilon),
        )
    }

    fn at(n: &Nozzle, gamma: f64, ratio: f64) -> anyhow::Result<NozzleOperatingPoint> {
        n.operating_point(&chamber(gamma)?, Pressure::new::<pascal>(ratio * P0))
    }

    /// Back-pressure ratios bounding the regimes: (subsonic limit, shock at
    /// exit, design).
    fn boundaries(gamma: f64, epsilon: f64) -> anyhow::Result<(f64, f64, f64)> {
        let sub = isentropic::mach_from_area_ratio(gamma, epsilon, Branch::Subsonic)?;
        let sup = isentropic::mach_from_area_ratio(gamma, epsilon, Branch::Supersonic)?;
        let design = isentropic::pressure_ratio(gamma, sup);
        Ok((
            isentropic::pressure_ratio(gamma, sub),
            design * shock::normal_shock(gamma, sup).pressure_ratio,
            design,
        ))
    }

    fn choked_flow(gamma: f64) -> f64 {
        AT * P0 / (287.0 * T0).sqrt() * isentropic::mass_flux_parameter(gamma, 1.0)
    }

    #[test]
    fn test_nozzle_validates_geometry() -> anyhow::Result<()> {
        for (throat, exit) in [(0.0, 1e-4), (1e-4, 0.5e-4), (f64::NAN, 1e-4)] {
            anyhow::ensure!(
                Nozzle::new(
                    Area::new::<square_meter>(throat),
                    Area::new::<square_meter>(exit)
                )
                .is_err(),
                "accepted throat {throat}, exit {exit}"
            );
        }
        anyhow::ensure!((nozzle(2.5)?.expansion_ratio() - 2.5).abs() < 1e-12, "ε");
        Ok(())
    }

    #[test]
    fn test_nozzle_rejects_invalid_back_pressure() -> anyhow::Result<()> {
        let n = nozzle(2.0)?;
        for pb in [-1.0, f64::NAN] {
            anyhow::ensure!(
                n.operating_point(&chamber(1.4)?, Pressure::new::<pascal>(pb))
                    .is_err(),
                "accepted p_b {pb}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_nozzle_classifies_every_regime() -> anyhow::Result<()> {
        let n = nozzle(1.6875)?;
        let (r1, r2, r3) = boundaries(1.4, 1.6875)?;
        anyhow::ensure!(
            (r3 - 0.1278).abs() < 1e-4 && (r2 - 0.5752).abs() < 1e-3,
            "{r2} {r3}"
        );
        let cases = [
            (1.0, "NoFlow"),
            (f64::midpoint(1.0, r1), "Subsonic"),
            (f64::midpoint(r1, r2), "ShockInNozzle"),
            (r2, "ShockAtExit"),
            (f64::midpoint(r2, r3), "OverExpanded"),
            (r3, "Ideal"),
            (0.5 * r3, "UnderExpanded"),
            (0.0, "UnderExpanded"),
        ];
        for (ratio, expected) in cases {
            let regime = at(&n, 1.4, ratio)?.regime;
            let name = format!("{regime:?}");
            anyhow::ensure!(name.starts_with(expected), "p_b/p0 = {ratio}: {name}");
        }
        Ok(())
    }

    #[test]
    fn test_nozzle_no_flow_at_or_above_chamber_pressure() -> anyhow::Result<()> {
        for ratio in [1.0, 1.5] {
            let point = at(&nozzle(2.0)?, 1.4, ratio)?;
            anyhow::ensure!(
                point.mass_flow.get::<kilogram_per_second>() == 0.0
                    && point.thrust.get::<newton>() == 0.0
                    && point.specific_impulse().is_none(),
                "p_b/p0 = {ratio}: {point:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_nozzle_ideal_expansion_matches_energy_equation() -> anyhow::Result<()> {
        let (gamma, epsilon) = (1.4, 4.0);
        let (_, _, r3) = boundaries(gamma, epsilon)?;
        let point = at(&nozzle(epsilon)?, gamma, r3)?;
        anyhow::ensure!(point.regime == NozzleRegime::Ideal, "{:?}", point.regime);
        let ve =
            (2.0 * gamma / (gamma - 1.0) * 287.0 * T0 * (1.0 - r3.powf((gamma - 1.0) / gamma)))
                .sqrt();
        let v = point.exit_velocity.get::<meter_per_second>();
        anyhow::ensure!((v / ve - 1.0).abs() < 1e-9, "vₑ {v} vs {ve}");
        let mdot = point.mass_flow.get::<kilogram_per_second>();
        anyhow::ensure!((mdot / choked_flow(gamma) - 1.0).abs() < 1e-12, "ṁ {mdot}");
        let thrust = point.thrust.get::<newton>();
        anyhow::ensure!((thrust / (mdot * ve) - 1.0).abs() < 1e-9, "F {thrust}");
        let isp = point.specific_impulse().map(|t| t.get::<second>());
        anyhow::ensure!(
            isp.is_some_and(|isp| (isp / (ve / STANDARD_GRAVITY) - 1.0).abs() < 1e-9),
            "Isp {isp:?}"
        );
        Ok(())
    }

    #[test]
    fn test_nozzle_shock_location_reproduces_back_pressure() -> anyhow::Result<()> {
        let (gamma, epsilon) = (1.4, 3.0);
        let (r1, r2, _) = boundaries(gamma, epsilon)?;
        let ratio = 0.3 * r1 + 0.7 * r2;
        let point = at(&nozzle(epsilon)?, gamma, ratio)?;
        let NozzleRegime::ShockInNozzle { shock_area_ratio } = point.regime else {
            anyhow::bail!("regime {:?}", point.regime);
        };
        anyhow::ensure!(
            shock_area_ratio > 1.0 && shock_area_ratio < epsilon,
            "{shock_area_ratio}"
        );
        // Independent route: shock Mach from its area, then the subsonic exit
        // behind it.
        let shock_mach =
            isentropic::mach_from_area_ratio(gamma, shock_area_ratio, Branch::Supersonic)?;
        let loss = shock::normal_shock(gamma, shock_mach).stagnation_pressure_ratio;
        let exit_mach = isentropic::mach_from_area_ratio(gamma, epsilon * loss, Branch::Subsonic)?;
        let exit_ratio = loss * isentropic::pressure_ratio(gamma, exit_mach);
        anyhow::ensure!(
            (exit_ratio / ratio - 1.0).abs() < 1e-9,
            "pₑ/p0 {exit_ratio} vs {ratio}"
        );
        anyhow::ensure!(
            (point.exit_mach / exit_mach - 1.0).abs() < 1e-9,
            "Mₑ {}",
            point.exit_mach
        );
        Ok(())
    }

    #[test]
    fn test_nozzle_subsonic_flow_is_below_choked() -> anyhow::Result<()> {
        let (r1, _, _) = boundaries(1.4, 2.0)?;
        let point = at(&nozzle(2.0)?, 1.4, f64::midpoint(1.0, r1))?;
        let mdot = point.mass_flow.get::<kilogram_per_second>();
        anyhow::ensure!(mdot > 0.0 && mdot < choked_flow(1.4), "ṁ {mdot}");
        anyhow::ensure!(point.exit_mach < 1.0, "Mₑ {}", point.exit_mach);
        let pe = point.exit_pressure.get::<pascal>();
        anyhow::ensure!(
            (pe / (f64::midpoint(1.0, r1) * P0) - 1.0).abs() < 1e-12,
            "pₑ {pe}"
        );
        Ok(())
    }

    #[test]
    fn test_nozzle_flags_separation_when_strongly_over_expanded() -> anyhow::Result<()> {
        // ε = 25 over-expands for 0.0019 < p_b/p0 < 0.055; at 0.02, pₑ ≈
        // 0.1·p_b.
        let point = at(&nozzle(25.0)?, 1.4, 0.02)?;
        anyhow::ensure!(
            point.regime == NozzleRegime::OverExpanded,
            "{:?}",
            point.regime
        );
        anyhow::ensure!(point.flow_separation, "not flagged");
        let mild = at(&nozzle(2.0)?, 1.4, 0.15)?;
        anyhow::ensure!(
            !mild.flow_separation,
            "mild over-expansion flagged: {mild:?}"
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_nozzle_thrust_is_continuous_across_regime_boundaries(
            gamma in 1.1..1.67f64, epsilon in 1.05..30.0f64,
        ) {
            let fail = |e: anyhow::Error| TestCaseError::fail(format!("{e:#}"));
            let n = nozzle(epsilon).map_err(fail)?;
            let (r1, r2, r3) = boundaries(gamma, epsilon).map_err(fail)?;
            // Slopes near the boundaries can be steep, so a fixed tolerance cannot tell a
            // slope from a step. A continuous function's jump shrinks with the probe width;
            // a discontinuity's does not.
            let jumps = |boundary: f64, width: f64| -> Result<(f64, f64), TestCaseError> {
                let below = at(&n, gamma, boundary * (1.0 - width)).map_err(fail)?;
                let above = at(&n, gamma, boundary * (1.0 + width)).map_err(fail)?;
                Ok((
                    (below.thrust - above.thrust).get::<newton>().abs(),
                    (below.mass_flow - above.mass_flow).get::<kilogram_per_second>().abs(),
                ))
            };
            let floor = (1e-12 * AT * P0, 1e-12 * choked_flow(gamma));
            for boundary in [r1, r2, r3] {
                let (wide_thrust, wide_flow) = jumps(boundary, 1e-6)?;
                let (thrust, flow) = jumps(boundary, 1e-8)?;
                prop_assert!(thrust <= 0.2 * wide_thrust + floor.0,
                    "thrust jump {thrust} N does not shrink ({wide_thrust}) at {boundary}");
                prop_assert!(flow <= 0.2 * wide_flow + floor.1,
                    "ṁ jump {flow} does not shrink ({wide_flow}) at {boundary}");
            }
        }

        #[test]
        fn test_nozzle_mass_flow_never_exceeds_choked(
            gamma in 1.1..1.67f64, epsilon in 1.0..30.0f64, ratio in 0.0..1.0f64,
        ) {
            let fail = |e: anyhow::Error| TestCaseError::fail(format!("{e:#}"));
            let point = at(&nozzle(epsilon).map_err(fail)?, gamma, ratio).map_err(fail)?;
            let mdot = point.mass_flow.get::<kilogram_per_second>();
            prop_assert!(mdot >= 0.0 && mdot <= choked_flow(gamma) * (1.0 + 1e-12), "ṁ {mdot}");
            prop_assert!(point.thrust.get::<newton>().is_finite(), "{point:?}");
        }
    }
}
