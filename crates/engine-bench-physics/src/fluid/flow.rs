//! Pure fluids as [`RealFluid`]s: states at (p, s) for isentropic flow (the
//! homogeneous equilibrium model inside the dome) and at (p, ρ) for shocks.

use anyhow::Context;

use super::{
    Fluid, PureFluid, density_si,
    helmholtz::{Molar, molar},
    saturation::{dome_top, saturation_si},
    state::dome_edge,
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

impl RealFluid for PureFluid {
    #[expect(
        clippy::many_single_char_names,
        reason = "the thermodynamic symbols p, s, T, h, w"
    )]
    fn state_ps(&self, p: Pressure, s: SpecificHeatCapacity) -> anyhow::Result<FlowState> {
        let fluid = self.data();
        let m = fluid.data.molar_mass;
        let (p, s) = (p.get::<pascal>(), s.get::<joule_per_kilogram_kelvin>() * m);
        check_pressure(fluid, p)?;
        let (lo, hi) = match saturation_temperature(fluid, p)? {
            Some(t_sat) => {
                let (_, rho_l, rho_v) = saturation_si(fluid, t_sat)?;
                let (l, v) = (molar(fluid, rho_l, t_sat), molar(fluid, rho_v, t_sat));
                if s >= l.s && s <= v.s {
                    return mixture_state(fluid, p, s, t_sat, (rho_l, l), (rho_v, v));
                }
                if s < l.s {
                    (fluid.data.t_triple, t_sat)
                } else {
                    (t_sat, fluid.data.t_max)
                }
            }
            None => (fluid.data.t_triple, fluid.data.t_max),
        };
        // Along an isobar s rises with T, nearly linearly in ln T: ∂s/∂ln T =
        // c_p. Next to the dome a bracket end can sit on the far side by
        // round-off; that saturated state is then the answer.
        let entropy_at = |t: f64| -> anyhow::Result<f64> {
            let rho = density_si(fluid, p, t)?;
            Ok(molar(fluid, rho, t).s - s)
        };
        let t = if entropy_at(lo)? >= 0.0 {
            lo
        } else if entropy_at(hi)? <= 0.0 {
            hi
        } else {
            let mut failure = None;
            let ln_t = newton_bracketed(
                |ln_t| {
                    // exp(ln T) can round past the bracket ends.
                    let t = ln_t.exp().clamp(lo, hi);
                    match density_si(fluid, p, t) {
                        Ok(rho) => {
                            let x = molar(fluid, rho, t);
                            (x.s - s, x.cp)
                        }
                        Err(e) => {
                            failure.get_or_insert(e);
                            (f64::NAN, f64::NAN)
                        }
                    }
                },
                lo.ln(),
                hi.ln(),
                f64::midpoint(lo.ln(), hi.ln()),
                RootOptions::new(1e-12, 200)?,
            );
            if let Some(e) = failure {
                return Err(e).context("isentropic temperature");
            }
            ln_t.with_context(|| format!("{self:?} at {p} Pa, s = {s} J/(mol K)"))?
                .exp()
                .clamp(lo, hi)
        };
        let rho = density_si(fluid, p, t)?;
        Ok(single_phase(fluid, &molar(fluid, rho, t), rho, t))
    }

    fn state_prho(&self, p: Pressure, rho: MassDensity) -> anyhow::Result<FlowState> {
        let fluid = self.data();
        let m = fluid.data.molar_mass;
        let (p, rho) = (p.get::<pascal>(), rho.get::<kilogram_per_cubic_meter>() / m);
        check_pressure(fluid, p)?;
        if !(rho > 0.0 && rho.is_finite()) {
            return Err(PhysicsError::InvalidState("density must be positive"))
                .with_context(|| format!("ρ = {rho} mol/m³"));
        }
        // At fixed density p rises with T: ∂p/∂T > 0. Inside the dome's density
        // range the homogeneous equation is unstable below the dome's edge and
        // can have other roots there, so the search starts at the edge.
        let lo = dome_edge(fluid, rho)?.unwrap_or(fluid.data.t_triple);
        let hi = fluid.data.t_max;
        let t = newton_bracketed(
            |t| {
                let x = molar(fluid, rho, t);
                (x.p - p, x.dp_dt)
            },
            lo,
            hi,
            f64::midpoint(lo, hi),
            RootOptions::new(1e-12 * hi, 200)?,
        )
        .with_context(|| format!("{self:?} at {p} Pa, {rho} mol/m³"))?;
        Ok(single_phase(fluid, &molar(fluid, rho, t), rho, t))
    }
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

fn single_phase(fluid: &Fluid, x: &Molar, rho: f64, t: f64) -> FlowState {
    let m = fluid.data.molar_mass;
    FlowState::from_si(x.p, t, rho * m, x.h / m, x.s / m, x.w, f64::NAN)
}

/// Saturation temperature, K, at `p` Pa, or `None` outside the dome's
/// pressure range. Newton with the Clausius–Clapeyron slope Δs/Δv.
fn saturation_temperature(fluid: &Fluid, p: f64) -> anyhow::Result<Option<f64>> {
    let (lo, hi) = (fluid.data.t_triple, dome_top(fluid));
    let (p_lo, _, _) = saturation_si(fluid, lo)?;
    let (p_hi, _, _) = saturation_si(fluid, hi)?;
    if !(p_lo..=p_hi).contains(&p) {
        return Ok(None);
    }
    // ln p is nearly linear in 1/T: a good first guess.
    let fraction = (p.ln() - p_lo.ln()) / (p_hi.ln() - p_lo.ln());
    let guess = 1.0 / (1.0 / lo + fraction * (1.0 / hi - 1.0 / lo));
    let mut failure = None;
    let t = newton_bracketed(
        |t| match saturation_si(fluid, t) {
            Ok((p_sat, rho_l, rho_v)) => {
                let (l, v) = (molar(fluid, rho_l, t), molar(fluid, rho_v, t));
                (p_sat - p, (v.s - l.s) / (1.0 / rho_v - 1.0 / rho_l))
            }
            Err(e) => {
                failure.get_or_insert(e);
                (f64::NAN, f64::NAN)
            }
        },
        lo,
        hi,
        guess,
        RootOptions::new(1e-12 * hi, 100)?,
    );
    if let Some(e) = failure {
        return Err(e).context("saturation temperature");
    }
    Ok(Some(t.context("saturation temperature")?))
}

/// Liquid and vapour in equilibrium at `p` Pa with molar entropy `s`: the
/// homogeneous equilibrium model. Its speed of sound, √(∂p/∂ρ)_s, comes from
/// the mixture density on the same isentrope a little above and below `p`.
#[expect(
    clippy::many_single_char_names,
    reason = "the thermodynamic symbols p, s, T, h, w"
)]
fn mixture_state(
    fluid: &Fluid,
    p: f64,
    s: f64,
    t: f64,
    (rho_l, l): (f64, Molar),
    (rho_v, v): (f64, Molar),
) -> anyhow::Result<FlowState> {
    let x = (s - l.s) / (v.s - l.s);
    let rho = 1.0 / ((1.0 - x) / rho_l + x / rho_v);
    let h = x.mul_add(v.h - l.h, l.h);
    let density_at = |p: f64| -> anyhow::Result<f64> {
        let t = saturation_temperature(fluid, p)?.context("mixture left the dome")?;
        let (_, rho_l, rho_v) = saturation_si(fluid, t)?;
        let (l, v) = (molar(fluid, rho_l, t), molar(fluid, rho_v, t));
        let x = (s - l.s) / (v.s - l.s);
        Ok(1.0 / ((1.0 - x) / rho_l + x / rho_v))
    };
    let dp = 1e-6 * p;
    let w = (2.0 * dp / (density_at(p + dp)? - density_at(p - dp)?)).sqrt();
    let m = fluid.data.molar_mass;
    Ok(FlowState::from_si(
        p,
        t,
        rho * m,
        h / m,
        s / m,
        w / m.sqrt(),
        x,
    ))
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::units::{
        ThermodynamicTemperature, thermodynamic_temperature::kelvin, velocity::meter_per_second,
    };

    fn fail(e: &anyhow::Error) -> TestCaseError {
        TestCaseError::fail(format!("{e:#}"))
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
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_pure_fluid_flow_states_round_trip_single_phase(fluid in fluids(), tr in 0.0f64..1.0, log_p in 4.0f64..7.3) {
            // Above the melting lines: N2 70 K, Ar 95 K, He 12 K, CO2 225 K.
            let t_low = match fluid {
                PureFluid::Nitrogen => 70.0,
                PureFluid::Argon => 95.0,
                PureFluid::Helium => 12.0,
                PureFluid::CarbonDioxide => 225.0,
            };
            let t = t_low + tr * (600.0 - t_low);
            let p = Pressure::new::<pascal>(10f64.powf(log_p));
            let x = fluid.properties(p, ThermodynamicTemperature::new::<kelvin>(t)).map_err(|e| fail(&e))?;
            let by_s = fluid.state_ps(p, x.entropy).map_err(|e| fail(&e))?;
            let by_rho = fluid.state_prho(p, x.density).map_err(|e| fail(&e))?;
            for y in [by_s, by_rho] {
                prop_assert!((y.temperature.get::<kelvin>() / t - 1.0).abs() < 1e-9, "{:?}: T {:?} vs {}", fluid, y.temperature, t);
                let (w, w0) = (y.speed_of_sound.get::<meter_per_second>(), x.speed_of_sound.get::<meter_per_second>());
                prop_assert!((w / w0 - 1.0).abs() < 1e-8, "w {} vs {}", w, w0);
            }
            prop_assert!(by_s.vapour_quality.is_none(), "{:?}", by_s);
        }
    }

    #[expect(
        clippy::many_single_char_names,
        reason = "the thermodynamic symbols p, s, T, h, w"
    )]
    #[test]
    fn test_pure_fluid_state_ps_inside_the_dome_is_an_equilibrium_mixture() -> anyhow::Result<()> {
        let co2 = PureFluid::CarbonDioxide;
        let t = ThermodynamicTemperature::new::<kelvin>(260.0);
        let sat = co2.saturation(t)?;
        let (l, v) = (
            co2.homogeneous(t, sat.liquid_density)?,
            co2.homogeneous(t, sat.vapour_density)?,
        );
        let s = l.entropy + (v.entropy - l.entropy) * 0.3;
        let x = co2.state_ps(sat.pressure, s)?;
        let quality = x
            .vapour_quality
            .ok_or_else(|| anyhow::anyhow!("single phase"))?;
        anyhow::ensure!((quality - 0.3).abs() < 1e-9, "x = {quality}");
        anyhow::ensure!(
            (x.temperature.get::<kelvin>() - 260.0).abs() < 1e-8,
            "{x:?}"
        );
        let v_mix = 0.7 / sat.liquid_density.get::<kilogram_per_cubic_meter>()
            + 0.3 / sat.vapour_density.get::<kilogram_per_cubic_meter>();
        anyhow::ensure!(
            (x.density.get::<kilogram_per_cubic_meter>() * v_mix - 1.0).abs() < 1e-9,
            "ρ"
        );
        // The equilibrium mixture is softer than either phase alone.
        let w = x.speed_of_sound.get::<meter_per_second>();
        anyhow::ensure!(
            w > 0.0 && w < v.speed_of_sound.get::<meter_per_second>(),
            "w = {w}"
        );
        Ok(())
    }
}
