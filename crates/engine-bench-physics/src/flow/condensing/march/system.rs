//! The march as an ODE in the axial position x: the impulse function and the
//! droplets' liquid change; everything else follows at each station.

use std::f64::consts::PI;

use super::{
    classes::{CLASSES, Classes, radius},
    state::{Conserved, Station, saturation, solve},
};
use crate::{
    flow::condensing::{
        CondensationModel, CondensingFluid, Contour,
        growth::{Growth, Surroundings},
        nucleation::{Supersaturated, nucleation},
    },
    numerics::ode::OdeSystem,
    units::{
        ThermodynamicTemperature, dynamic_viscosity::pascal_second,
        thermal_conductivity::watt_per_meter_kelvin, thermodynamic_temperature::kelvin,
    },
};

/// ODE state: impulse function / its scale, the open class's droplets per
/// kg / the number scale and liquid / `Y_SCALE`, then each closed class's
/// liquid / `Y_SCALE`.
pub(super) const N: usize = 3 + CLASSES;
pub(super) const IMPULSE: usize = 0;
pub(super) const OPEN_NUMBER: usize = 1;
pub(super) const OPEN_LIQUID: usize = 2;
pub(super) const FIRST_CLASS: usize = 3;
/// Liquid fractions are carried in units of 1 %.
pub(super) const Y_SCALE: f64 = 0.01;

/// Everything the march knows of the flow at one x.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Local {
    pub(super) station: Station,
    /// Nucleation rate, 1/(m³·s), and critical radius, m.
    pub(super) rate: f64,
    pub(super) critical_radius: f64,
    /// Droplets per kg, Σ n r², Σ n r³ over all classes.
    pub(super) number: f64,
    pub(super) area_moment: f64,
    pub(super) volume_moment: f64,
}

/// The march: fluid, nozzle, model and the state carried between steps.
pub(super) struct March<'a, F> {
    pub(super) fluid: &'a F,
    pub(super) contour: &'a Contour,
    pub(super) model: CondensationModel,
    pub(super) conserved: Conserved,
    /// Mass of a molecule, kg.
    pub(super) molecule: f64,
    pub(super) impulse_scale: f64,
    pub(super) number_scale: f64,
    pub(super) classes: Classes,
    /// The last station solved: the next solve's starting point.
    pub(super) last: Station,
}

impl<F: CondensingFluid> March<'_, F> {
    /// The flow at `x` with ODE state `y`, and the state's derivative into
    /// `dydx` if given; `None` if the station cannot be solved.
    #[expect(
        clippy::many_single_char_names,
        reason = "the march's symbols: x, y, s, v, d"
    )]
    pub(super) fn local(
        &mut self,
        x: f64,
        y: &[f64; N],
        dydx: Option<&mut [f64; N]>,
    ) -> Option<Local> {
        let (area, darea) = self.contour.area_si(x);
        let impulse = y[IMPULSE] * self.impulse_scale;
        let open_number = y[OPEN_NUMBER] * self.number_scale;
        let open_liquid = y[OPEN_LIQUID] * Y_SCALE;
        let count = self.classes.count();
        let class_liquid = |k: usize| y.get(FIRST_CLASS + k).copied().unwrap_or(0.0) * Y_SCALE;
        let liquid = open_liquid + (0..count).map(class_liquid).sum::<f64>();
        let station = solve(
            self.fluid,
            &self.conserved,
            area,
            impulse,
            liquid.max(0.0),
            &self.last,
        )?;
        self.last = station;
        let (s, v) = (&station, &station.vapour);

        let (rate, critical_radius) = self.nucleation(s)?;
        let has_droplets = count > 0 || open_number > 0.0;
        let growth = if has_droplets {
            Some(self.growth(s)?)
        } else {
            None
        };

        let rho_l = s.droplets.rho_liquid;
        let flux = s.rho * s.u;
        // Liquid gained per unit length by `number` droplets per kg of
        // radius r: n·4πr²ρ_l·(dr/dt)/u.
        let gain = |number: f64, r: f64| match growth {
            Some(g) if r > 0.0 => number * 4.0 * PI * r * r * rho_l * g.rate(r) / s.u,
            _ => 0.0,
        };
        let mut local = Local {
            station,
            rate,
            critical_radius,
            number: 0.0,
            area_moment: 0.0,
            volume_moment: 0.0,
        };
        let mut add = |number: f64, r: f64| {
            local.number += number;
            local.area_moment += number * r * r;
            local.volume_moment += number * r * r * r;
        };
        let open_radius = radius(open_number, open_liquid, rho_l);
        add(open_number, open_radius);
        for k in 0..count {
            add(
                self.classes.number(k),
                radius(self.classes.number(k), class_liquid(k), rho_l),
            );
        }

        if let Some(d) = dydx {
            d.fill(0.0);
            d[IMPULSE] = v.p * darea / self.impulse_scale;
            let born = if rate > 0.0 { rate / flux } else { 0.0 };
            let born_mass = if born > 0.0 {
                born * 4.0 / 3.0 * PI * critical_radius.powi(3) * rho_l
            } else {
                0.0
            };
            d[OPEN_NUMBER] = born / self.number_scale;
            d[OPEN_LIQUID] = (born_mass + gain(open_number, open_radius)) / Y_SCALE;
            for k in 0..count {
                let n = self.classes.number(k);
                let yk = class_liquid(k);
                if let Some(slot) = d.get_mut(FIRST_CLASS + k) {
                    *slot = if yk > 0.0 {
                        gain(n, radius(n, yk, rho_l)) / Y_SCALE
                    } else {
                        0.0
                    };
                }
            }
        }
        Some(local)
    }
}

impl<F: CondensingFluid> March<'_, F> {
    /// Nucleation rate, 1/(m³·s), and critical radius, m, at station `s`;
    /// `None` if the fluid's properties are out of range there.
    fn nucleation(&self, s: &Station) -> Option<(f64, f64)> {
        if !self.model.nucleation || s.t >= self.fluid.dome_top() {
            return Some((0.0, f64::INFINITY));
        }
        let v = &s.vapour;
        let at_t = saturation(self.fluid, s.t)?;
        Some(nucleation(
            &Supersaturated {
                p: v.p,
                t: s.t,
                rho: v.rho,
                gamma: v.cp / v.cv,
                p_sat: at_t.p,
                sigma: at_t.sigma,
                rho_liquid: at_t.rho_liquid,
                latent: at_t.latent,
                molecule: self.molecule,
            },
            self.model.condensation_coefficient,
            self.model.kantrowitz,
        ))
    }

    /// The growth law at station `s`; `None` above the dome or out of the
    /// transport correlations' range.
    fn growth(&self, s: &Station) -> Option<Growth> {
        if !s.t_sat.is_finite() {
            return None;
        }
        let v = &s.vapour;
        let transport = self
            .fluid
            .vapour_transport(ThermodynamicTemperature::new::<kelvin>(s.t))
            .ok()?;
        Some(Growth::new(
            &Surroundings {
                p: v.p,
                t: s.t,
                t_sat: s.t_sat,
                sigma: s.droplets.sigma,
                rho_liquid: s.droplets.rho_liquid,
                latent: s.droplets.latent,
                viscosity: transport.viscosity.get::<pascal_second>(),
                conductivity: transport
                    .thermal_conductivity
                    .get::<watt_per_meter_kelvin>(),
                cp: v.cp,
                gamma: v.cp / v.cv,
                gas_constant: self.conserved.gas_constant,
            },
            &self.model,
        ))
    }
}

impl<F: CondensingFluid> OdeSystem<N> for March<'_, F> {
    fn derivatives(&mut self, x: f64, y: &[f64; N], dydx: &mut [f64; N]) {
        if self.local(x, y, Some(dydx)).is_none() {
            dydx.fill(f64::NAN);
        }
    }
}
