//! Hot-kernel benchmarks: a baseline for the real-time loop, not a pass/fail
//! gate.
//!
//! Run with `cargo bench -p engine-bench-physics --profile bench-unwind`.

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use engine_bench_physics::{
    numerics::{
        ode::{AdaptiveOptions, DormandPrince45, NoEvents, OdeSystem, rk4_step},
        roots::{RootOptions, brent, newton_bracketed},
    },
    thermo::{Major7, Mixture, TemperaturePowers, gases},
    units::{ThermodynamicTemperature, thermodynamic_temperature::kelvin},
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

criterion_group!(kernels, ode, roots, thermo);
criterion_main!(kernels);
