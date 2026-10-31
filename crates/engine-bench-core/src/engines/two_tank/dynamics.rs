//! The two-tank ODE system.

use engine_bench_physics::{
    flow::Orifice,
    gas::orifice_transfer,
    numerics::ode::OdeSystem,
    units::{energy::joule, thermodynamic_temperature::kelvin},
};

use super::{INITIAL_TEMPERATURE, Tank, VOLUME_A, VOLUME_B};

/// The ODE: state `[m_N2, m_O2, m_Ar, U]` of tank A, then of tank B.
pub(super) struct Dynamics<'a> {
    pub(super) a: &'a mut Tank,
    pub(super) b: &'a mut Tank,
    pub(super) valve: Option<Orifice>,
    /// W into tank B.
    pub(super) heater: f64,
    /// Wall heat transfer coefficient, W/(m²·K).
    pub(super) wall_h: f64,
}

/// Wall area of a cube of `volume` m³, m².
pub(super) fn wall_area(volume: f64) -> f64 {
    6.0 * volume.powf(2.0 / 3.0)
}

/// Heat flow into `tank` (W) through a wall of `area` with coefficient `h`,
/// from air at the initial temperature.
pub(super) fn wall_heat(tank: &Tank, h: f64, area: f64) -> f64 {
    h * area * (INITIAL_TEMPERATURE - tank.temperature().get::<kelvin>())
}

impl Dynamics<'_> {
    pub(super) fn state(&self) -> [f64; 8] {
        let (ma, mb) = (self.a.species_masses(), self.b.species_masses());
        [
            ma[0],
            ma[1],
            ma[2],
            self.a.internal_energy().get::<joule>(),
            mb[0],
            mb[1],
            mb[2],
            self.b.internal_energy().get::<joule>(),
        ]
    }

    pub(super) fn unpack(&mut self, y: &[f64; 8]) {
        self.a.set_conserved(&[y[0], y[1], y[2]], y[3], VOLUME_A);
        self.b.set_conserved(&[y[4], y[5], y[6]], y[7], VOLUME_B);
    }
}

impl OdeSystem<8> for Dynamics<'_> {
    fn derivatives(&mut self, _t: f64, y: &[f64; 8], dydt: &mut [f64; 8]) {
        self.unpack(y);
        *dydt = [0.0; 8];
        if let Some(valve) = self.valve {
            let flow = orifice_transfer(valve, &*self.a, &*self.b);
            for i in 0..3 {
                dydt[i] = -flow.species[i];
                dydt[4 + i] = flow.species[i];
            }
            dydt[3] = -flow.energy;
            dydt[7] = flow.energy;
        }
        dydt[3] += wall_heat(self.a, self.wall_h, wall_area(VOLUME_A));
        dydt[7] += self.heater + wall_heat(self.b, self.wall_h, wall_area(VOLUME_B));
    }
}
