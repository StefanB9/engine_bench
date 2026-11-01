//! Condensing nozzle flow against measurements.
//!
//! Steam: Moses and Stein's nozzle, cases 252 and 257, as read from the
//! figures of Starzmann et al., "Results of the International Wet Steam
//! Modeling Project", Proc. IMechE Part A 232(5), 2018 (plan 0005). The
//! effective area is fitted to the dry pressures, as the paper recommends for
//! 1D validation (the small throat's boundary layers block part of it).
//!
//! The tolerances are the agreement measured with the default model
//! (Young's α = 9), kept as a regression bound; the paper's participants
//! scatter over a similar range, most predicting too small droplets at the
//! lower pressure (case 252).

use engine_bench_physics::{
    flow::condensing::{CondensationModel, CondensingFlow, CondensingNozzle, Contour},
    steam::Water,
    units::{
        Angle, Length, Pressure, ThermodynamicTemperature,
        angle::degree,
        length::{meter, millimeter},
        pressure::{kilopascal, pascal},
        thermodynamic_temperature::kelvin,
    },
};

/// One measured case: supply state, (x cm, p/p₀) and (x cm, Sauter radius
/// nm), and how many leading pressure points are dry.
struct Case {
    p0_kpa: f64,
    t0: f64,
    pressures: &'static [(f64, f64)],
    radii: &'static [(f64, f64)],
    dry: usize,
}

#[expect(
    clippy::approx_constant,
    reason = "3.141 cm is a measured position, not π"
)]
const CASE_252: Case = Case {
    p0_kpa: 40.05,
    t0: 374.3,
    pressures: &[
        (2.054, 0.3969),
        (2.259, 0.3829),
        (2.469, 0.3687),
        (2.682, 0.3544),
        (2.902, 0.3402),
        (3.141, 0.3259),
        (3.371, 0.3119),
        (3.559, 0.3004),
        (3.732, 0.2970),
        (3.911, 0.2971),
        (4.085, 0.2989),
        (4.259, 0.2985),
        (4.424, 0.2953),
        (4.560, 0.2896),
        (4.770, 0.2798),
        (5.068, 0.2661),
        (5.412, 0.2508),
        (5.691, 0.2381),
        (5.993, 0.2241),
    ],
    radii: &[
        (3.465, 10.57),
        (3.665, 11.90),
        (3.853, 12.92),
        (4.023, 13.75),
        (4.209, 14.53),
        (4.411, 15.24),
        (4.604, 15.81),
        (4.791, 16.35),
        (4.979, 16.78),
        (5.165, 17.13),
        (5.377, 17.44),
        (5.575, 17.75),
        (5.773, 18.03),
        (5.970, 18.19),
    ],
    dry: 6,
};

const CASE_257: Case = Case {
    p0_kpa: 67.66,
    t0: 376.7,
    pressures: &[
        (1.957, 0.4000),
        (2.212, 0.3897),
        (2.840, 0.3814),
        (3.446, 0.3841),
        (3.920, 0.3557),
        (4.480, 0.3236),
        (5.111, 0.2870),
    ],
    radii: &[(2.495, 16.83), (3.895, 29.24), (6.003, 34.34)],
    dry: 2,
};

/// Moses and Stein's nozzle (10 mm deep, 10 mm throat) as a throat arc of
/// radius `curvature` m to 7 cm.
fn moses_stein(curvature: f64) -> anyhow::Result<Contour> {
    let rise = curvature - (curvature * curvature - 0.07f64 * 0.07).sqrt();
    Contour::planar(
        Length::new::<millimeter>(10.0),
        Length::new::<millimeter>(10.0),
        Length::new::<meter>(0.01 + 2.0 * rise),
        Length::new::<meter>(curvature),
        Angle::new::<degree>(15.0),
    )
}

fn march(case: &Case, curvature: f64, model: CondensationModel) -> anyhow::Result<CondensingFlow> {
    CondensingNozzle::new(Water, moses_stein(curvature)?, model).expand(
        Pressure::new::<kilopascal>(case.p0_kpa),
        ThermodynamicTemperature::new::<kelvin>(case.t0),
    )
}

/// The flow's p/p₀ at `x_cm`, by linear interpolation.
fn pressure_at(flow: &CondensingFlow, p0: f64, x_cm: f64) -> f64 {
    sample_at(flow, x_cm, |s| s.pressure.get::<pascal>() / p0)
}

fn sample_at(
    flow: &CondensingFlow,
    x_cm: f64,
    f: impl Fn(&engine_bench_physics::flow::condensing::FlowSample) -> f64,
) -> f64 {
    let x = x_cm / 100.0;
    let s = flow.samples();
    s.windows(2)
        .find(|w| w[0].position.get::<meter>() <= x && x <= w[1].position.get::<meter>())
        .map_or(f64::NAN, |w| {
            let (x0, x1) = (w[0].position.get::<meter>(), w[1].position.get::<meter>());
            let k = if x1 > x0 { (x - x0) / (x1 - x0) } else { 0.0 };
            f(&w[0]) + k * (f(&w[1]) - f(&w[0]))
        })
}

/// Throat curvature, m, whose frozen expansion best fits the dry points
/// (golden-section search on the squared error).
fn fit_curvature(case: &Case) -> anyhow::Result<f64> {
    let p0 = case.p0_kpa * 1e3;
    let error = |r: f64| -> anyhow::Result<f64> {
        let flow = march(case, r, CondensationModel::frozen())?;
        Ok(case
            .pressures
            .iter()
            .take(case.dry)
            .map(|&(x, p)| (pressure_at(&flow, p0, x) - p).powi(2))
            .sum())
    };
    let golden = 0.5 * (5f64.sqrt() - 1.0);
    let (mut a, mut b) = (0.4, 2.0);
    let (mut c, mut d) = (b - golden * (b - a), a + golden * (b - a));
    let (mut fc, mut fd) = (error(c)?, error(d)?);
    while b - a > 1e-4 {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - golden * (b - a);
            fc = error(c)?;
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + golden * (b - a);
            fd = error(d)?;
        }
    }
    Ok(f64::midpoint(a, b))
}

fn fitted(case: &Case) -> anyhow::Result<(f64, CondensingFlow)> {
    let r = fit_curvature(case)?;
    Ok((r, march(case, r, CondensationModel::default())?))
}

/// The largest |Δ(p/p₀)| over `points`.
fn worst_pressure(flow: &CondensingFlow, case: &Case, points: &[(f64, f64)]) -> f64 {
    points
        .iter()
        .map(|&(x, p)| (pressure_at(flow, case.p0_kpa * 1e3, x) - p).abs())
        .fold(0.0, f64::max)
}

/// Model over measured Sauter radius at each measured point.
fn radius_ratios(flow: &CondensingFlow, case: &Case) -> Vec<f64> {
    case.radii
        .iter()
        .map(|&(x, r)| sample_at(flow, x, |s| s.sauter_radius.get::<meter>() * 1e9) / r)
        .collect()
}

#[test]
fn test_moses_stein_252_condensation_against_measurements() -> anyhow::Result<()> {
    let (r, flow) = fitted(&CASE_252)?;
    // The dry points fix a 0.82 m effective throat arc (the wall's is
    // 0.686 m: boundary layers slow the effective area's growth).
    anyhow::ensure!((0.80..0.84).contains(&r), "effective curvature {r} m");
    let dry = CASE_252.pressures.get(..CASE_252.dry).unwrap_or_default();
    let dry_error = worst_pressure(&flow, &CASE_252, dry);
    anyhow::ensure!(dry_error < 0.002, "dry pressures off by {dry_error}");
    // The condensation rise comes ~0.5 cm late and overshoots downstream.
    let error = worst_pressure(&flow, &CASE_252, CASE_252.pressures);
    anyhow::ensure!(error < 0.025, "pressures off by {error}");
    // Final droplets 30 % small (the measurement is ±20 %).
    let last = radius_ratios(&flow, &CASE_252)
        .last()
        .copied()
        .unwrap_or(f64::NAN);
    anyhow::ensure!((0.65..1.2).contains(&last), "final radius ratio {last}");
    let wilson = flow.wilson_point().map_or(f64::NAN, |x| x.get::<meter>());
    anyhow::ensure!((0.034..0.042).contains(&wilson), "Wilson point {wilson} m");
    Ok(())
}

#[test]
fn test_moses_stein_257_condensation_against_measurements() -> anyhow::Result<()> {
    let (r, flow) = fitted(&CASE_257)?;
    // The same effective nozzle as case 252's, fitted independently.
    anyhow::ensure!((0.80..0.84).contains(&r), "effective curvature {r} m");
    let error = worst_pressure(&flow, &CASE_257, CASE_257.pressures);
    anyhow::ensure!(error < 0.005, "pressures off by {error}");
    for ratio in radius_ratios(&flow, &CASE_257) {
        anyhow::ensure!((0.75..1.2).contains(&ratio), "radius ratio {ratio}");
    }
    // At the higher pressure condensation sets in sooner.
    let wilson = flow.wilson_point().map_or(f64::NAN, |x| x.get::<meter>());
    anyhow::ensure!((0.018..0.028).contains(&wilson), "Wilson point {wilson} m");
    Ok(())
}
