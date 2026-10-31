//! Orifice flow between volumes with a real equation of state: the upstream
//! port's stagnation state expands along its own isentrope, saturated vapour
//! when a two-phase volume is drawn from the top.

use anyhow::Context;

use super::{
    ControlVolume, Environment, GasPort, Transfer,
    eos::{EosFluid, EquationOfState, IdealGas},
};
use crate::{
    flow::{
        Orifice,
        real::{RealFluid, Stagnation},
    },
    thermo::SpeciesSet,
    units::{
        MassDensity, MassRate, Pressure, available_energy::joule_per_kilogram,
        mass_density::kilogram_per_cubic_meter, mass_rate::kilogram_per_second, pressure::pascal,
    },
};

/// A port real-fluid orifice flow can leave from.
pub trait RealPort<S: SpeciesSet<N>, const N: usize>: GasPort<S, N> {
    /// Pressure at the port.
    fn port_pressure(&self) -> Pressure;
    /// Stagnation state of what leaves through the port.
    fn port_stagnation(&self) -> anyhow::Result<Stagnation>;
    /// Mass flow out through `orifice` into `p_down` (≤ the port pressure).
    fn outflow(&self, orifice: &Orifice, p_down: Pressure) -> anyhow::Result<MassRate>;
    /// Stagnation state of what leaves through the port, any state the port
    /// cannot give directly coming from `fluid` (which must describe the
    /// port's contents).
    fn port_stagnation_in<F: RealFluid>(&self, fluid: &F) -> anyhow::Result<Stagnation> {
        let _ = fluid;
        self.port_stagnation()
    }
}

impl<S: SpeciesSet<N>, const N: usize, E: EquationOfState<S, N>> RealPort<S, N>
    for ControlVolume<S, N, E>
{
    fn port_pressure(&self) -> Pressure {
        self.pressure()
    }

    fn port_stagnation(&self) -> anyhow::Result<Stagnation> {
        self.port_stagnation_in(&EosFluid {
            eos: self.eos(),
            mixture: self.mixture(),
        })
    }

    fn port_stagnation_in<F: RealFluid>(&self, fluid: &F) -> anyhow::Result<Stagnation> {
        let state = self.fluid_state();
        if state.quality.is_nan() {
            return Ok(Stagnation {
                p: state.p,
                h: state.h,
                s: self.eos().entropy(self.mixture(), state),
            });
        }
        // Two-phase: the port draws saturated vapour, of density p/(Z·R·T).
        let rho = state.p / (state.gas_constant * state.t);
        let vapour = fluid
            .state_prho(
                Pressure::new::<pascal>(state.p),
                MassDensity::new::<kilogram_per_cubic_meter>(rho),
            )
            .context("saturated vapour at the port")?;
        Ok(Stagnation::at_rest(&vapour))
    }

    fn outflow(&self, orifice: &Orifice, p_down: Pressure) -> anyhow::Result<MassRate> {
        let fluid = EosFluid {
            eos: self.eos(),
            mixture: self.mixture(),
        };
        orifice.real_mass_flow(&fluid, &self.port_stagnation()?, p_down)
    }
}

impl<S: SpeciesSet<N>, const N: usize> RealPort<S, N> for Environment<S, N> {
    fn port_pressure(&self) -> Pressure {
        Pressure::new::<pascal>(self.reservoir().pressure)
    }

    fn port_stagnation(&self) -> anyhow::Result<Stagnation> {
        let r = self.reservoir();
        let rho = r.pressure / (self.mixture().gas_constant_si() * r.temperature);
        let state = IdealGas
            .flow_state_prho(self.mixture(), r.pressure, rho)
            .context("surroundings")?;
        Ok(Stagnation::at_rest(&state))
    }

    fn outflow(&self, orifice: &Orifice, p_down: Pressure) -> anyhow::Result<MassRate> {
        let fluid = EosFluid {
            eos: &IdealGas,
            mixture: self.mixture(),
        };
        orifice.real_mass_flow(&fluid, &self.port_stagnation()?, p_down)
    }
}

/// Flow through `orifice` from `a` to `b` (negative values flow from `b` to
/// `a`) along the upstream port's real-fluid isentrope, carrying its
/// composition and enthalpy. Subtract it from `a`'s rates and add it to `b`'s.
pub fn real_orifice_transfer<S: SpeciesSet<N>, const N: usize>(
    orifice: &Orifice,
    a: &impl RealPort<S, N>,
    b: &impl RealPort<S, N>,
) -> anyhow::Result<Transfer<N>> {
    transfer(
        a,
        b,
        |a, p_down| a.outflow(orifice, p_down),
        |b, p_down| b.outflow(orifice, p_down),
    )
}

/// [`real_orifice_transfer`] with the expansion computed in `fluid`, which
/// must describe what both ports hold: e.g. a
/// [`TabulatedFluid`](crate::fluid::TabulatedFluid) for real time.
pub fn real_orifice_transfer_in<S: SpeciesSet<N>, const N: usize, F: RealFluid>(
    fluid: &F,
    orifice: &Orifice,
    a: &impl RealPort<S, N>,
    b: &impl RealPort<S, N>,
) -> anyhow::Result<Transfer<N>> {
    transfer(
        a,
        b,
        |a, p_down| orifice.real_mass_flow(fluid, &a.port_stagnation_in(fluid)?, p_down),
        |b, p_down| orifice.real_mass_flow(fluid, &b.port_stagnation_in(fluid)?, p_down),
    )
}

/// The transfer between `a` and `b`, the mass flow out of the upstream port
/// from `out_a` or `out_b` given the downstream pressure.
fn transfer<S: SpeciesSet<N>, const N: usize, A: RealPort<S, N>, B: RealPort<S, N>>(
    a: &A,
    b: &B,
    out_a: impl FnOnce(&A, Pressure) -> anyhow::Result<MassRate>,
    out_b: impl FnOnce(&B, Pressure) -> anyhow::Result<MassRate>,
) -> anyhow::Result<Transfer<N>> {
    let (pa, pb) = (a.port_pressure(), b.port_pressure());
    let (mass_flow, fractions, enthalpy) = if pa >= pb {
        let flow = out_a(a, pb).context("flow out of the first port")?;
        (
            flow.get::<kilogram_per_second>(),
            a.mass_fractions(),
            a.specific_enthalpy(),
        )
    } else {
        let flow = out_b(b, pa).context("flow out of the second port")?;
        (
            -flow.get::<kilogram_per_second>(),
            b.mass_fractions(),
            b.specific_enthalpy(),
        )
    };
    Ok(Transfer {
        species: fractions.map(|y| mass_flow * y),
        energy: mass_flow * enthalpy.get::<joule_per_kilogram>(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        fluid::PureFluid,
        gas::{
            eos::{PengRobinson, Reference},
            orifice_transfer,
        },
        numerics::ode::{OdeSystem, rk4_step},
        thermo::{Air3, ColdGas, Mixture},
        units::{
            Area, Ratio, ThermodynamicTemperature, Volume, area::square_meter, energy::joule,
            pressure::pascal, ratio::ratio, thermodynamic_temperature::kelvin, volume::cubic_meter,
        },
    };

    fn air_tank(p: f64, t: f64) -> anyhow::Result<ControlVolume<Air3, 3>> {
        ControlVolume::new(
            Mixture::dry_air()?,
            Pressure::new::<pascal>(p),
            ThermodynamicTemperature::new::<kelvin>(t),
            Volume::new::<cubic_meter>(0.01),
        )
    }

    fn orifice() -> anyhow::Result<Orifice> {
        Orifice::new(Area::new::<square_meter>(1e-5), 0.8)
    }

    #[test]
    fn test_real_transfer_matches_the_perfect_gas_transfer_for_cold_ideal_air() -> anyhow::Result<()>
    {
        // Thermally perfect (c_p(T)) against frozen γ: close for cold air.
        let (a, b) = (air_tank(5e5, 300.0)?, air_tank(1e5, 300.0)?);
        let real = real_orifice_transfer(&orifice()?, &a, &b)?;
        let perfect = orifice_transfer(orifice()?, &a, &b);
        anyhow::ensure!(
            (real.mass_flow() / perfect.mass_flow() - 1.0).abs() < 2e-3,
            "{} vs {}",
            real.mass_flow(),
            perfect.mass_flow()
        );
        anyhow::ensure!((real.energy / perfect.energy - 1.0).abs() < 2e-3, "energy");
        Ok(())
    }

    #[test]
    fn test_real_transfer_is_antisymmetric() -> anyhow::Result<()> {
        let (a, b) = (air_tank(5e5, 400.0)?, air_tank(1e5, 300.0)?);
        let forward = real_orifice_transfer(&orifice()?, &a, &b)?;
        let backward = real_orifice_transfer(&orifice()?, &b, &a)?;
        anyhow::ensure!(
            forward.energy.to_bits() == (-backward.energy).to_bits(),
            "energy"
        );
        anyhow::ensure!(
            forward
                .species
                .iter()
                .zip(backward.species)
                .all(|(f, b)| f.to_bits() == (-b).to_bits()),
            "species"
        );
        Ok(())
    }

    #[test]
    fn test_real_transfer_in_a_given_fluid_matches_the_volumes_own() -> anyhow::Result<()> {
        // A CO2 bottle (liquid and vapour) feeding a CO2 plenum, both ways.
        let co2 = Mixture::<ColdGas, 4>::from_mass_fractions([0.0, 0.0, 0.0, 1.0])?;
        let bottle = ControlVolume::with_saturated(
            Reference::new(PureFluid::CarbonDioxide)?,
            co2,
            ThermodynamicTemperature::new::<kelvin>(290.0),
            Ratio::new::<ratio>(0.6),
            Volume::new::<cubic_meter>(1e-3),
        )?;
        let plenum = ControlVolume::with_eos(
            Reference::new(PureFluid::CarbonDioxide)?,
            co2,
            Pressure::new::<pascal>(1e6),
            ThermodynamicTemperature::new::<kelvin>(280.0),
            Volume::new::<cubic_meter>(1e-4),
        )?;
        let eos = Reference::<ColdGas, 4>::new(PureFluid::CarbonDioxide)?;
        let fluid = EosFluid {
            eos: &eos,
            mixture: &co2,
        };
        for (a, b) in [(&bottle, &plenum), (&plenum, &bottle)] {
            let given = real_orifice_transfer_in(&fluid, &orifice()?, a, b)?;
            let own = real_orifice_transfer(&orifice()?, a, b)?;
            anyhow::ensure!(
                given.energy.to_bits() == own.energy.to_bits()
                    && given.mass_flow().to_bits() == own.mass_flow().to_bits(),
                "{given:?} vs {own:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_real_transfer_vents_saturated_vapour_from_a_co2_bottle() -> anyhow::Result<()> {
        let bottle = ControlVolume::with_saturated(
            Reference::new(PureFluid::CarbonDioxide)?,
            Mixture::<ColdGas, 4>::from_mass_fractions([0.0, 0.0, 0.0, 1.0])?,
            ThermodynamicTemperature::new::<kelvin>(290.0),
            Ratio::new::<ratio>(0.6),
            Volume::new::<cubic_meter>(1e-3),
        )?;
        let room = Environment::new(
            Mixture::<ColdGas, 4>::from_mass_fractions([1.0, 0.0, 0.0, 0.0])?,
            Pressure::new::<pascal>(1e5),
            ThermodynamicTemperature::new::<kelvin>(290.0),
        )?;
        let t = real_orifice_transfer(&orifice()?, &bottle, &room)?;
        anyhow::ensure!(
            t.mass_flow() > 0.0 && t.species[3].to_bits() == t.mass_flow().to_bits(),
            "{t:?}"
        );
        let h = bottle
            .specific_enthalpy()
            .get::<crate::units::available_energy::joule_per_kilogram>();
        anyhow::ensure!(
            (t.energy / (t.mass_flow() * h) - 1.0).abs() < 1e-12,
            "energy"
        );
        Ok(())
    }

    /// Two Peng–Robinson nitrogen tanks joined by an orifice.
    struct Pair {
        tanks: [ControlVolume<ColdGas, 4, PengRobinson<ColdGas, 4>>; 2],
        orifice: Orifice,
    }

    impl OdeSystem<4> for Pair {
        fn derivatives(&mut self, _t: f64, y: &[f64; 4], dydt: &mut [f64; 4]) {
            for (i, tank) in self.tanks.iter_mut().enumerate() {
                tank.set_conserved(&[y[2 * i], 0.0, 0.0, 0.0], y[2 * i + 1], 0.01);
            }
            let flow = real_orifice_transfer(&self.orifice, &self.tanks[0], &self.tanks[1])
                .map_or([f64::NAN; 2], |t| [t.mass_flow(), t.energy]);
            *dydt = [-flow[0], -flow[1], flow[0], flow[1]];
        }
    }

    #[test]
    fn test_real_transfer_conserves_mass_and_energy_between_dense_tanks() -> anyhow::Result<()> {
        let n2 = Mixture::<ColdGas, 4>::from_mass_fractions([1.0, 0.0, 0.0, 0.0])?;
        let tank = |p: f64| -> anyhow::Result<_> {
            ControlVolume::with_eos(
                PengRobinson::new()?,
                n2,
                Pressure::new::<pascal>(p),
                ThermodynamicTemperature::new::<kelvin>(300.0),
                Volume::new::<cubic_meter>(0.01),
            )
        };
        let mut pair = Pair {
            tanks: [tank(2e7)?, tank(1e6)?],
            orifice: orifice()?,
        };
        let mut y = [0.0; 4];
        for (i, cv) in pair.tanks.iter().enumerate() {
            y[2 * i] = cv.species_masses()[0];
            y[2 * i + 1] = cv.internal_energy().get::<joule>();
        }
        let (mass, energy) = (y[0] + y[2], y[1] + y[3]);
        let mut t = 0.0;
        let start = y[0];
        for _ in 0..200 {
            rk4_step(&mut pair, &mut t, &mut y, 1e-3).context("RK4 step")?;
        }
        anyhow::ensure!(y[0] < 0.99 * start, "no flow: {} kg of {start}", y[0]);
        anyhow::ensure!(((y[0] + y[2]) / mass - 1.0).abs() < 1e-13, "mass drift");
        anyhow::ensure!(((y[1] + y[3]) / energy - 1.0).abs() < 1e-12, "energy drift");
        Ok(())
    }
}
