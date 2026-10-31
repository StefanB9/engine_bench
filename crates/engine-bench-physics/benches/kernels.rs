//! Hot-kernel benchmarks: a baseline for the real-time loop, not a pass/fail
//! gate.
//!
//! Run with `cargo bench -p engine-bench-physics --profile bench-unwind`.

use std::{hint::black_box, io::Write as _};

use criterion::{Criterion, criterion_group, criterion_main};
use engine_bench_physics::{
    flow::{Nozzle, Orifice, PerfectGas, Reservoir},
    gas::{ControlVolume, eos::PengRobinson, orifice_transfer},
    numerics::{
        ode::{AdaptiveOptions, DormandPrince45, NoEvents, OdeSystem, rk4_step},
        roots::{RootOptions, brent, newton_bracketed},
    },
    steam,
    thermo::{Air3, ColdGas, Major7, Mixture, TemperaturePowers, gases},
    units::{
        Area, Pressure, Ratio, SpecificHeatCapacity, ThermodynamicTemperature, Volume,
        area::square_meter, energy::joule, mass::kilogram, pressure::pascal, ratio::ratio,
        specific_heat_capacity::joule_per_kilogram_kelvin, thermodynamic_temperature::kelvin,
        volume::cubic_meter,
    },
};

/// Lorenz system: small, nonlinear and non-trivial for step-size control.
struct Lorenz;

impl OdeSystem<3> for Lorenz {
    fn derivatives(&mut self, _t: f64, y: &[f64; 3], dydt: &mut [f64; 3]) {
        dydt[0] = 10.0 * (y[1] - y[0]);
        dydt[1] = y[0] * (28.0 - y[2]) - y[1];
        dydt[2] = y[0] * y[1] - 8.0 / 3.0 * y[2];
    }
}

fn ode(c: &mut Criterion) {
    c.bench_function("rk4_step/lorenz", |b| {
        let (mut t, mut y) = (0.0, [1.0, 1.0, 1.0]);
        b.iter(|| rk4_step(&mut Lorenz, &mut t, black_box(&mut y), black_box(1e-4)));
    });

    c.bench_function("dormand_prince/lorenz_1s_rtol1e-8", |b| {
        b.iter(|| -> anyhow::Result<[f64; 3]> {
            let mut dp = DormandPrince45::new(AdaptiveOptions::new(1e-8, 1e-10)?);
            let (mut t, mut y) = (0.0, black_box([1.0, 1.0, 1.0]));
            dp.integrate(&mut Lorenz, &mut NoEvents, &mut t, &mut y, 1.0)?;
            Ok(y)
        });
    });
}

fn roots(c: &mut Criterion) {
    // Shaped like recovering temperature from enthalpy: monotone, mildly
    // curved.
    let h = |t: f64| 1000.0 * t + 0.1 * t * t - 1.0e6;
    let dh = |t: f64| 1000.0 + 0.2 * t;

    c.bench_function("brent/enthalpy_like", |b| {
        b.iter(|| {
            brent(
                h,
                black_box(200.0),
                black_box(6000.0),
                RootOptions::default(),
            )
        });
    });
    c.bench_function("newton_bracketed/enthalpy_like_warm_start", |b| {
        b.iter(|| {
            newton_bracketed(
                |t| (h(t), dh(t)),
                200.0,
                6000.0,
                black_box(950.0),
                RootOptions::default(),
            )
        });
    });
}

fn thermo(c: &mut Criterion) {
    let t = ThermodynamicTemperature::new::<kelvin>(1850.0);
    c.bench_function("species/co2_cp_h_s", |b| {
        b.iter(|| {
            let p = TemperaturePowers::new(black_box(t));
            (
                gases::CO2.cp_over_r(&p),
                gases::CO2.h_over_rt(&p),
                gases::CO2.s_over_r(&p),
            )
        });
    });

    // Lean hydrogen-air combustion products.
    let Ok(products) =
        Mixture::<Major7, 7>::from_mass_fractions([0.72, 0.12, 0.012, 0.148, 0.0, 0.0, 0.0])
    else {
        return;
    };
    c.bench_function("mixture_major7/cp", |b| {
        b.iter(|| products.cp(black_box(t)));
    });
    c.bench_function("mixture_major7/enthalpy", |b| {
        b.iter(|| products.enthalpy(black_box(t)));
    });
    let u = products.internal_energy(t);
    let warm = ThermodynamicTemperature::new::<kelvin>(1840.0);
    let cold = ThermodynamicTemperature::new::<kelvin>(300.0);
    c.bench_function("mixture_major7/t_from_u_warm_start", |b| {
        b.iter(|| products.temperature_from_internal_energy(black_box(u), black_box(warm)));
    });
    c.bench_function("mixture_major7/t_from_u_cold_start", |b| {
        b.iter(|| products.temperature_from_internal_energy(black_box(u), black_box(cold)));
    });
}

fn flow(c: &mut Criterion) -> anyhow::Result<()> {
    let air = PerfectGas::new(
        1.4,
        SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(287.0),
    )?;
    let reservoir = |p: f64| {
        Reservoir::new(
            Pressure::new::<pascal>(p),
            ThermodynamicTemperature::new::<kelvin>(300.0),
            air,
        )
    };
    let tank = reservoir(1.0e6)?;
    let orifice = Orifice::new(Area::new::<square_meter>(1e-5), 0.8)?;
    for (name, p_down) in [
        ("linear_band", 0.9995e6),
        ("unchoked", 0.8e6),
        ("choked", 1e5),
    ] {
        let down = reservoir(p_down)?;
        c.bench_function(&format!("orifice/{name}"), |b| {
            b.iter(|| orifice.mass_flow(black_box(&tank), black_box(&down)));
        });
    }

    let nozzle = Nozzle::new(
        Area::new::<square_meter>(1e-4),
        Area::new::<square_meter>(4e-4),
    )?;
    for (name, p_back) in [("under_expanded", 1e4), ("shock_in_nozzle", 6e5)] {
        let p_back = Pressure::new::<pascal>(p_back);
        c.bench_function(&format!("nozzle/{name}"), |b| {
            b.iter(|| nozzle.operating_point(black_box(&tank), black_box(p_back)));
        });
    }
    Ok(())
}

/// What one ODE derivative evaluation of a two-tank model costs.
fn gas(c: &mut Criterion) -> anyhow::Result<()> {
    let tank = |p: f64, t: f64| -> anyhow::Result<ControlVolume<Air3, 3>> {
        ControlVolume::new(
            Mixture::dry_air()?,
            Pressure::new::<pascal>(p),
            ThermodynamicTemperature::new::<kelvin>(t),
            Volume::new::<cubic_meter>(0.01),
        )
    };
    let (mut a, b) = (tank(8e5, 500.0)?, tank(1e5, 300.0)?);
    let masses = *a.species_masses();
    let energy = a.internal_energy().get::<joule>();
    c.bench_function("control_volume/set_conserved_air3", |bench| {
        bench.iter(|| {
            // A slightly perturbed state, as between Runge–Kutta stages.
            a.set_conserved(black_box(&masses), black_box(energy * 1.000_001), 0.01);
            a.pressure()
        });
    });
    // The same update for 200 bar of nitrogen through Peng-Robinson.
    let mut n2 = ControlVolume::with_eos(
        PengRobinson::new()?,
        Mixture::<ColdGas, 4>::from_mass_fractions([1.0, 0.0, 0.0, 0.0])?,
        Pressure::new::<pascal>(2e7),
        ThermodynamicTemperature::new::<kelvin>(300.0),
        Volume::new::<cubic_meter>(0.01),
    )?;
    let (n2_masses, n2_energy) = (*n2.species_masses(), n2.internal_energy().get::<joule>());
    c.bench_function("control_volume/set_conserved_peng_robinson_n2", |bench| {
        bench.iter(|| {
            n2.set_conserved(
                black_box(&n2_masses),
                black_box(n2_energy * 1.000_001),
                0.01,
            );
            n2.pressure()
        });
    });
    let orifice = Orifice::new(Area::new::<square_meter>(2e-5), 0.8)?;
    c.bench_function("gas/orifice_transfer_air3", |bench| {
        bench.iter(|| orifice_transfer(black_box(orifice), black_box(&a), black_box(&b)));
    });
    Ok(())
}

/// IF97 forward equations, the (u, v) state solve and the exact T(p, h)
/// inverse.
#[expect(
    clippy::many_single_char_names,
    reason = "thermodynamic symbols p, T, u, h, m"
)]
fn steam_kernels(c: &mut Criterion) -> anyhow::Result<()> {
    for (name, p, t) in [
        ("steam/properties_region1", 3e6, 300.0),
        ("steam/properties_region2", 1e6, 500.0),
        ("steam/properties_region5", 3e7, 1500.0),
    ] {
        let (p, t) = (
            Pressure::new::<pascal>(p),
            ThermodynamicTemperature::new::<kelvin>(t),
        );
        steam::properties(p, t)?;
        c.bench_function(name, |b| {
            b.iter(|| steam::properties(black_box(p), black_box(t)));
        });
    }
    let (p3, t3) = (
        Pressure::new::<pascal>(25e6),
        ThermodynamicTemperature::new::<kelvin>(650.0),
    );
    steam::properties(p3, t3)?;
    c.bench_function("steam/properties_region3", |b| {
        b.iter(|| steam::properties(black_box(p3), black_box(t3)));
    });
    // A boiler at 450 K, 30 % liquid by volume: the (u, v) solve from scratch,
    // and the warm-started update a simulation step makes.
    let mut boiler = steam::SteamVolume::saturated(
        ThermodynamicTemperature::new::<kelvin>(450.0),
        Ratio::new::<ratio>(0.3),
        Volume::new::<cubic_meter>(1e-3),
    )?;
    let state = boiler.state();
    c.bench_function("steam/state_from_uv_two_phase_cold", |b| {
        b.iter(|| {
            steam::state_from_uv(
                black_box(state.internal_energy),
                black_box(state.specific_volume),
            )
        });
    });
    let (m, u) = (
        boiler.mass().get::<kilogram>(),
        boiler.internal_energy().get::<joule>(),
    );
    c.bench_function("steam/volume_set_conserved_two_phase_warm", |b| {
        b.iter(|| {
            boiler.set_conserved(black_box(m), black_box(u * 1.000_001), 1e-3);
            boiler.pressure()
        });
    });
    let p = Pressure::new::<pascal>(1e6);
    let h = steam::properties(p, ThermodynamicTemperature::new::<kelvin>(500.0))?.enthalpy;
    c.bench_function("steam/temperature_from_enthalpy_region2", |b| {
        b.iter(|| steam::temperature_from_enthalpy(black_box(p), black_box(h)));
    });
    Ok(())
}

fn flow_benches(c: &mut Criterion) {
    // Criterion offers no error channel; report a setup failure instead of
    // skipping silently.
    if let Err(e) = flow(c) {
        let _ = writeln!(std::io::stderr(), "flow benchmarks skipped: {e:#}");
    }
    if let Err(e) = gas(c) {
        let _ = writeln!(std::io::stderr(), "gas benchmarks skipped: {e:#}");
    }
    if let Err(e) = steam_kernels(c) {
        let _ = writeln!(std::io::stderr(), "steam benchmarks skipped: {e:#}");
    }
}

criterion_group!(kernels, ode, roots, thermo, flow_benches);
criterion_main!(kernels);
