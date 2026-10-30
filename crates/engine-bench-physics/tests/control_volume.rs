//! Control volumes coupled through orifices, integrated in time: conservation
//! and analytic blowdown, through the public API only.

use anyhow::Context;
use engine_bench_physics::{
    flow::Orifice,
    gas::{ControlVolume, Environment, GasPort, orifice_transfer},
    numerics::ode::{AdaptiveOptions, DormandPrince45, NoEvents, OdeSystem, rk4_step},
    thermo::{Air3, Mixture},
    units::{
        Area, Pressure, ThermodynamicTemperature, Volume, area::square_meter, energy::joule,
        mass::kilogram, pressure::pascal, specific_heat_capacity::joule_per_kilogram_kelvin,
        thermodynamic_temperature::kelvin, volume::cubic_meter,
    },
};
use proptest::prelude::*;

type Tank = ControlVolume<Air3, 3>;

fn tank(mixture: Mixture<Air3, 3>, p: f64, t: f64, v: f64) -> anyhow::Result<Tank> {
    ControlVolume::new(
        mixture,
        Pressure::new::<pascal>(p),
        ThermodynamicTemperature::new::<kelvin>(t),
        Volume::new::<cubic_meter>(v),
    )
}

/// Packs a tank as `[m_N2, m_O2, m_Ar, U]`.
fn pack(cv: &Tank, out: &mut [f64; 4]) {
    out[..3].copy_from_slice(cv.species_masses());
    out[3] = cv.internal_energy().get::<joule>();
}

/// Two rigid tanks joined by an orifice. State: [tank a (4), tank b (4)].
struct TwoTanks {
    a: Tank,
    b: Tank,
    orifice: Orifice,
}

impl TwoTanks {
    fn unpack(&mut self, y: &[f64; 8]) {
        let va = self.a.volume().get::<cubic_meter>();
        let vb = self.b.volume().get::<cubic_meter>();
        self.a.set_conserved(&[y[0], y[1], y[2]], y[3], va);
        self.b.set_conserved(&[y[4], y[5], y[6]], y[7], vb);
    }

    fn state(&self) -> [f64; 8] {
        let (mut a, mut b) = ([0.0; 4], [0.0; 4]);
        pack(&self.a, &mut a);
        pack(&self.b, &mut b);
        [a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3]]
    }
}

impl OdeSystem<8> for TwoTanks {
    fn derivatives(&mut self, _t: f64, y: &[f64; 8], dydt: &mut [f64; 8]) {
        self.unpack(y);
        let flow = orifice_transfer(self.orifice, &self.a, &self.b);
        for i in 0..3 {
            dydt[i] = -flow.species[i];
            dydt[4 + i] = flow.species[i];
        }
        dydt[3] = -flow.energy;
        dydt[7] = flow.energy;
    }
}

fn totals(y: &[f64; 8]) -> ([f64; 3], f64) {
    ([y[0] + y[4], y[1] + y[5], y[2] + y[6]], y[3] + y[7])
}

fn two_tanks(pa: f64, ta: f64, pb: f64, tb: f64) -> anyhow::Result<TwoTanks> {
    Ok(TwoTanks {
        a: tank(Mixture::dry_air()?, pa, ta, 0.02)?,
        b: tank(Mixture::dry_air()?, pb, tb, 0.05)?,
        orifice: Orifice::new(Area::new::<square_meter>(2e-5), 0.8)?,
    })
}

#[test]
fn test_two_tanks_equalise_and_conserve_mass_and_energy() -> anyhow::Result<()> {
    let mut system = two_tanks(8e5, 500.0, 1e5, 300.0)?;
    let mut y = system.state();
    let (mass0, energy0) = totals(&y);
    let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-9, 1e-12)?);
    let mut t = 0.0;
    dp.integrate(&mut system, &mut NoEvents, &mut t, &mut y, 30.0)?;
    system.unpack(&y);
    system.a.validate().context("tank a")?;
    system.b.validate().context("tank b")?;

    let (mass1, energy1) = totals(&y);
    for (m0, m1) in mass0.iter().zip(mass1) {
        anyhow::ensure!((m1 / m0 - 1.0).abs() < 1e-12, "species mass {m0} → {m1}");
    }
    anyhow::ensure!(
        (energy1 / energy0 - 1.0).abs() < 1e-12,
        "energy {energy0} → {energy1}"
    );
    let (pa, pb) = (
        system.a.pressure().get::<pascal>(),
        system.b.pressure().get::<pascal>(),
    );
    anyhow::ensure!(
        (pa / pb - 1.0).abs() < 1e-6,
        "pressures did not equalise: {pa} vs {pb}"
    );
    Ok(())
}

#[test]
fn test_two_tanks_conserve_each_species_while_mixing() -> anyhow::Result<()> {
    let nitrogen = Mixture::<Air3, 3>::from_mass_fractions([1.0, 0.0, 0.0])?;
    let oxygen = Mixture::<Air3, 3>::from_mass_fractions([0.0, 1.0, 0.0])?;
    let mut system = TwoTanks {
        a: tank(nitrogen, 6e5, 300.0, 0.02)?,
        b: tank(oxygen, 2e5, 300.0, 0.05)?,
        orifice: Orifice::new(Area::new::<square_meter>(2e-5), 0.8)?,
    };
    let mut y = system.state();
    let (mass0, _) = totals(&y);
    let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-9, 1e-12)?);
    let mut t = 0.0;
    dp.integrate(&mut system, &mut NoEvents, &mut t, &mut y, 5.0)?;
    let (mass1, _) = totals(&y);
    anyhow::ensure!(mass1[2] == 0.0 && mass0[2] == 0.0, "argon appeared");
    for (m0, m1) in mass0.iter().zip(mass1).take(2) {
        anyhow::ensure!((m1 / m0 - 1.0).abs() < 1e-12, "species mass {m0} → {m1}");
    }
    // Mostly a → b, so b gains nitrogen. Temperatures keep changing after the
    // pressures meet, so the flow briefly reverses and a gains a trace of
    // oxygen; Runge–Kutta's negative weights can leave such an absent species a
    // roundoff-sized amount below zero.
    anyhow::ensure!(y[4] > 0.0, "no nitrogen reached b: {y:?}");
    let total = mass1[0] + mass1[1];
    anyhow::ensure!(
        y[1] >= -1e-8 * total && y[1] < 1e-3 * total,
        "oxygen in a: {}",
        y[1]
    );
    system.unpack(&y);
    system.a.validate().context("tank a after mixing")?;
    system.b.validate().context("tank b after mixing")
}

/// One tank emptying into an environment. State: `[m_N2, m_O2, m_Ar, U]`.
struct Blowdown {
    tank: Tank,
    outside: Environment<Air3, 3>,
    orifice: Orifice,
}

impl OdeSystem<4> for Blowdown {
    fn derivatives(&mut self, _t: f64, y: &[f64; 4], dydt: &mut [f64; 4]) {
        let volume = self.tank.volume().get::<cubic_meter>();
        self.tank.set_conserved(&[y[0], y[1], y[2]], y[3], volume);
        let flow = orifice_transfer(self.orifice, &self.tank, &self.outside);
        for i in 0..3 {
            dydt[i] = -flow.species[i];
        }
        dydt[3] = -flow.energy;
    }
}

/// Argon (constant γ = 5/3) emptying into vacuum through a choked orifice
/// follows the isentropic solution p/p₀ = (1 + (γ−1)/2·rate·t)^(−2γ/(γ−1)),
/// with rate = Cd·A·c₀·(2/(γ+1))^((γ+1)/(2(γ−1)))/V.
#[test]
fn test_argon_blowdown_into_vacuum_matches_isentropic_solution() -> anyhow::Result<()> {
    let argon = Mixture::<Air3, 3>::from_mass_fractions([0.0, 0.0, 1.0])?;
    let (p0, t0, volume, area, cd) = (2e6, 300.0, 0.01, 1e-5, 0.9);
    let mut system = Blowdown {
        tank: tank(argon, p0, t0, volume)?,
        outside: Environment::new(
            argon,
            Pressure::new::<pascal>(0.0),
            ThermodynamicTemperature::new::<kelvin>(300.0),
        )?,
        orifice: Orifice::new(Area::new::<square_meter>(area), cd)?,
    };
    let gamma: f64 = 5.0 / 3.0;
    let gas_constant = argon.gas_constant().get::<joule_per_kilogram_kelvin>();
    let choking = (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)));
    let rate = cd * area * (gamma * gas_constant * t0).sqrt() * choking / volume;
    let exact = |time: f64| {
        p0 * (1.0 + 0.5 * (gamma - 1.0) * rate * time).powf(-2.0 * gamma / (gamma - 1.0))
    };

    let mut y = [0.0; 4];
    pack(&system.tank, &mut y);
    let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-11, 1e-14)?);
    let mut t = 0.0;
    // Stop while T stays above the 200 K data floor (p/p₀ > 0.36).
    for t_end in [0.5, 1.0, 2.0, 3.0] {
        dp.integrate(&mut system, &mut NoEvents, &mut t, &mut y, t_end)?;
        let volume = system.tank.volume().get::<cubic_meter>();
        system.tank.set_conserved(&[y[0], y[1], y[2]], y[3], volume);
        system.tank.validate().with_context(|| format!("t = {t}"))?;
        let p = system.tank.pressure().get::<pascal>();
        anyhow::ensure!(
            (p / exact(t) - 1.0).abs() < 1e-6,
            "t = {t}: p {p} vs {}",
            exact(t)
        );
    }
    anyhow::ensure!(
        system.tank.pressure().get::<pascal>() < 0.6 * p0,
        "barely blew down"
    );
    Ok(())
}

/// A closed volume squeezed by a moving boundary. State: `[m_N2, m_O2, m_Ar,
/// U]`.
struct Compression {
    tank: Tank,
    initial_volume: f64,
}

impl Compression {
    /// Linear 5:1 compression over one second.
    fn volume(&self, t: f64) -> f64 {
        self.initial_volume * (1.0 - 0.8 * t)
    }
}

impl OdeSystem<4> for Compression {
    fn derivatives(&mut self, t: f64, y: &[f64; 4], dydt: &mut [f64; 4]) {
        self.tank
            .set_conserved(&[y[0], y[1], y[2]], y[3], self.volume(t));
        dydt[..3].fill(0.0);
        dydt[3] = self.tank.boundary_work_rate(-0.8 * self.initial_volume);
    }
}

#[test]
fn test_adiabatic_compression_conserves_entropy() -> anyhow::Result<()> {
    let initial_volume = 1e-3;
    let mut system = Compression {
        tank: tank(Mixture::dry_air()?, 1e5, 300.0, initial_volume)?,
        initial_volume,
    };
    let s0 = system
        .tank
        .specific_entropy()
        .get::<joule_per_kilogram_kelvin>();
    let mut y = [0.0; 4];
    pack(&system.tank, &mut y);
    let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-11, 1e-14)?);
    let mut t = 0.0;
    dp.integrate(&mut system, &mut NoEvents, &mut t, &mut y, 1.0)?;
    let final_volume = system.volume(1.0);
    system
        .tank
        .set_conserved(&[y[0], y[1], y[2]], y[3], final_volume);
    system.tank.validate()?;

    let s1 = system
        .tank
        .specific_entropy()
        .get::<joule_per_kilogram_kelvin>();
    anyhow::ensure!((s1 - s0).abs() < 1e-6 * s0.abs(), "entropy {s0} → {s1}");
    // Thermally perfect air heats less than constant-γ air would (γ falls as T
    // rises), but well past the 300 K start: 5^0.4 · 300 K = 571 K for γ = 1.4.
    let t1 = system.tank.temperature().get::<kelvin>();
    anyhow::ensure!(
        (540.0..571.0).contains(&t1),
        "T after 5:1 compression: {t1}"
    );
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1000))]

    #[test]
    fn test_two_tanks_conserve_for_random_states(
        pa in 2e4..5e6f64, ta in 220.0..2000.0f64, pb in 2e4..5e6f64, tb in 220.0..2000.0f64,
    ) {
        let fail = |e: anyhow::Error| TestCaseError::fail(format!("{e:#}"));
        let mut system = two_tanks(pa, ta, pb, tb).map_err(fail)?;
        let mut y = system.state();
        let (mass0, energy0) = totals(&y);
        let mut t = 0.0;
        for _ in 0..20 {
            rk4_step(&mut system, &mut t, &mut y, 1e-4).map_err(fail)?;
        }
        let (mass1, energy1) = totals(&y);
        for (m0, m1) in mass0.iter().zip(mass1) {
            prop_assert!((m1 - m0).abs() <= 1e-12 * m0, "species mass {m0} → {m1}");
        }
        prop_assert!((energy1 - energy0).abs() <= 1e-12 * energy0.abs().max(1.0),
            "energy {energy0} → {energy1}");
        system.unpack(&y);
        let total_mass = system.a.mass().get::<kilogram>() + system.b.mass().get::<kilogram>();
        prop_assert!(total_mass > 0.0);
        let _ = system.a.reservoir();
    }
}
