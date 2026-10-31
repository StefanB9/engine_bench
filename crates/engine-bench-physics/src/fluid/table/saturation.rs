//! The saturation line as a table in θ = (1 − T/Tc)^(1/3): the phases'
//! properties are smooth in θ up to the critical point, where they are not in
//! T. Inside the dome a state is mixed exactly (homogeneous equilibrium) from
//! the interpolated saturated liquid and vapour.

use super::{
    super::{Fluid, helmholtz::molar, saturation::dome_top, saturation_si},
    grid::{Axis, Grid1},
};

/// Nodes along the saturation line.
const NODES: usize = 400;
/// Newton iterations inverting ln p(θ) within one cell.
const INVERSION_ITERATIONS: usize = 8;

/// Per node: ln p, `s_l`, `s_v` (J/(kg·K)), `h_l`, `h_v` (J/kg), ln `v_l`, ln
/// `v_v` (m³/kg).
const LN_P: usize = 0;
const S_L: usize = 1;
const S_V: usize = 2;
const H_L: usize = 3;
const H_V: usize = 4;
const LN_V_L: usize = 5;
const LN_V_V: usize = 6;

/// The saturation line of one fluid.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct SaturationTable {
    grid: Grid1<7>,
    t_crit: f64,
}

/// Saturated liquid and vapour at one pressure, with their derivatives along
/// the saturation line, d/dp.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Saturated {
    pub(super) t: f64,
    pub(super) s_l: f64,
    pub(super) s_v: f64,
    h_l: f64,
    h_v: f64,
    pub(super) v_l: f64,
    pub(super) v_v: f64,
    ds_l: f64,
    ds_v: f64,
    dv_l: f64,
    dv_v: f64,
}

/// A two-phase equilibrium state: T K, density kg/m³, h J/kg, speed of sound
/// m/s and vapour quality.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Mixture {
    pub(super) t: f64,
    pub(super) rho: f64,
    pub(super) h: f64,
    pub(super) w: f64,
    pub(super) quality: f64,
}

impl SaturationTable {
    /// The table from the triple point to the top of the dome.
    #[expect(
        clippy::many_single_char_names,
        reason = "the thermodynamic symbols p, T, l, v, m"
    )]
    pub(super) fn build(fluid: &Fluid) -> anyhow::Result<Self> {
        let data = fluid.data;
        let t_crit = data.t_crit;
        let theta = |t: f64| (1.0 - t / t_crit).cbrt();
        let axis = Axis::new(theta(dome_top(fluid)), theta(data.t_triple), NODES);
        let m = data.molar_mass;
        let mut values = Vec::with_capacity(NODES);
        for i in 0..NODES {
            let t = (t_crit * (1.0 - axis.node(i).powi(3))).clamp(data.t_triple, dome_top(fluid));
            let (p, rho_l, rho_v) = saturation_si(fluid, t)?;
            let (l, v) = (molar(fluid, rho_l, t), molar(fluid, rho_v, t));
            values.push([
                p.ln(),
                l.s / m,
                v.s / m,
                l.h / m,
                v.h / m,
                -(rho_l * m).ln(),
                -(rho_v * m).ln(),
            ]);
        }
        Ok(Self {
            grid: Grid1::new(axis, values),
            t_crit,
        })
    }

    /// The saturated phases at `p` Pa, or `None` outside the table's
    /// interior (near the triple point, near the critical point, or off the
    /// line).
    pub(super) fn at_pressure(&self, p: f64) -> Option<Saturated> {
        let theta = self.theta(p.ln())?;
        let x = self.grid.sample(theta)?;
        let v = x.value;
        // d/dp = (d/dθ) / (dp/dθ), with dp/dθ = p · d ln p/dθ.
        let dtheta_dp = 1.0 / (p * x.dy[LN_P]);
        let (v_l, v_v) = (v[LN_V_L].exp(), v[LN_V_V].exp());
        Some(Saturated {
            t: self.t_crit * (1.0 - theta.powi(3)),
            s_l: v[S_L],
            s_v: v[S_V],
            h_l: v[H_L],
            h_v: v[H_V],
            v_l,
            v_v,
            ds_l: x.dy[S_L] * dtheta_dp,
            ds_v: x.dy[S_V] * dtheta_dp,
            dv_l: v_l * x.dy[LN_V_L] * dtheta_dp,
            dv_v: v_v * x.dy[LN_V_V] * dtheta_dp,
        })
    }

    /// θ where ln p(θ) = `ln_p`: the cell by bisection over the nodes (ln p
    /// falls as θ rises), then Newton on the cubic.
    fn theta(&self, ln_p: f64) -> Option<f64> {
        let nodes = self.grid.nodes();
        let (first, last) = (nodes.first()?[LN_P], nodes.last()?[LN_P]);
        if !(ln_p <= first && ln_p >= last) {
            return None;
        }
        let i = nodes
            .partition_point(|n| n[LN_P] > ln_p)
            .clamp(1, nodes.len() - 1);
        let axis = self.grid.x;
        let (a, b) = (nodes[i - 1][LN_P], nodes[i][LN_P]);
        let mut theta = axis.node(i - 1) + (axis.node(i) - axis.node(i - 1)) * (a - ln_p) / (a - b);
        for _ in 0..INVERSION_ITERATIONS {
            let x = self.grid.sample(theta)?;
            let step = (x.value[LN_P] - ln_p) / x.dy[LN_P];
            theta -= step;
            if step.abs() < 1e-15 {
                break;
            }
        }
        Some(theta)
    }
}

impl Saturated {
    /// The equilibrium mixture with entropy `s` J/(kg·K), between `s_l` and
    /// `s_v`. Its speed of sound is v·√(−1/(dv/dp)) along the isentrope, where
    /// the quality shifts as the saturated phases do.
    pub(super) fn mixture(&self, s: f64) -> Mixture {
        let ds = self.s_v - self.s_l;
        let x = (s - self.s_l) / ds;
        let v = x.mul_add(self.v_v - self.v_l, self.v_l);
        let dx = -x.mul_add(self.ds_v - self.ds_l, self.ds_l) / ds;
        let dv = x.mul_add(self.dv_v - self.dv_l, self.dv_l) + dx * (self.v_v - self.v_l);
        Mixture {
            t: self.t,
            rho: 1.0 / v,
            h: x.mul_add(self.h_v - self.h_l, self.h_l),
            w: v * (-1.0 / dv).sqrt(),
            quality: x,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        flow::real::RealFluid,
        fluid::PureFluid,
        units::{
            SpecificHeatCapacity, ThermodynamicTemperature, mass_density::kilogram_per_cubic_meter,
            pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
            thermodynamic_temperature::kelvin, velocity::meter_per_second,
        },
    };

    #[test]
    #[expect(
        clippy::many_single_char_names,
        reason = "the thermodynamic symbols p, T, l, v, m"
    )]
    fn test_saturation_table_matches_the_exact_saturation_line() -> anyhow::Result<()> {
        for fluid in PureFluid::ALL {
            let data = fluid.data();
            let table = SaturationTable::build(data)?;
            let m = data.data.molar_mass;
            let (t_lo, t_hi) = (data.data.t_triple, data.data.t_crit);
            for fraction in [0.05, 0.3, 0.6, 0.9, 0.97] {
                let t = t_lo + fraction * (t_hi - t_lo);
                let (p, rho_l, rho_v) = saturation_si(data, t)?;
                let sat = table
                    .at_pressure(p)
                    .ok_or_else(|| anyhow::anyhow!("{fluid:?} at {t} K off the table"))?;
                let (l, v) = (molar(data, rho_l, t), molar(data, rho_v, t));
                let pairs = [
                    ("T", sat.t, t),
                    ("s_l", sat.s_l, l.s / m),
                    ("s_v", sat.s_v, v.s / m),
                    ("v_l", sat.v_l, 1.0 / (rho_l * m)),
                    ("v_v", sat.v_v, 1.0 / (rho_v * m)),
                ];
                for (name, got, want) in pairs {
                    anyhow::ensure!(
                        (got - want).abs() <= 1e-7 * want.abs().max(1.0),
                        "{fluid:?} {t} K {name}: {got} vs {want}"
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_saturation_table_mixture_matches_the_exact_hem_state() -> anyhow::Result<()> {
        for (fluid, t) in [
            (PureFluid::CarbonDioxide, 260.0),
            (PureFluid::Nitrogen, 90.0),
            (PureFluid::Argon, 120.0),
        ] {
            let data = fluid.data();
            let table = SaturationTable::build(data)?;
            let sat = fluid.saturation(ThermodynamicTemperature::new::<kelvin>(t))?;
            let p = sat.pressure.get::<pascal>();
            let s_bounds = table.at_pressure(p).ok_or_else(|| anyhow::anyhow!("off"))?;
            let s = s_bounds.s_l + 0.6 * (s_bounds.s_v - s_bounds.s_l);
            let mix = s_bounds.mixture(s);
            let exact = fluid.state_ps(
                sat.pressure,
                SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(s),
            )?;
            let rho = exact.density.get::<kilogram_per_cubic_meter>();
            let w = exact.speed_of_sound.get::<meter_per_second>();
            anyhow::ensure!(
                (mix.rho / rho - 1.0).abs() < 1e-7,
                "{fluid:?} ρ {} vs {rho}",
                mix.rho
            );
            // The exact speed of sound is a finite difference.
            anyhow::ensure!(
                (mix.w / w - 1.0).abs() < 1e-4,
                "{fluid:?} w {} vs {w}",
                mix.w
            );
            anyhow::ensure!(
                exact
                    .vapour_quality
                    .is_some_and(|q| (q - mix.quality).abs() < 1e-7),
                "{fluid:?} quality {} vs {:?}",
                mix.quality,
                exact.vapour_quality
            );
        }
        Ok(())
    }

    #[test]
    fn test_saturation_table_has_no_states_off_the_line() -> anyhow::Result<()> {
        let table = SaturationTable::build(PureFluid::Nitrogen.data())?;
        anyhow::ensure!(table.at_pressure(1.0).is_none(), "1 Pa");
        anyhow::ensure!(table.at_pressure(1e8).is_none(), "100 MPa");
        anyhow::ensure!(table.at_pressure(f64::NAN).is_none(), "NaN");
        Ok(())
    }
}
