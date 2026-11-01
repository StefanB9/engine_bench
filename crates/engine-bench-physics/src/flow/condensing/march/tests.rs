use super::*;
use crate::{
    flow::{Nozzle, condensing::kernel::VapourKernel, real::Stagnation},
    fluid::{FrozenVapour, PureFluid},
    steam::Water,
    units::{
        Angle,
        angle::degree,
        force::newton,
        length::millimeter,
        mass_rate::kilogram_per_second,
        pressure::{kilopascal, torr},
    },
};

/// Yale nozzle B (Wegener et al. 1980): planar, 7 mm deep, 4 mm throat,
/// 20 mm throat curvature, 3.4° walls, 13.3 mm exit.
fn yale() -> anyhow::Result<Contour> {
    let mm = Length::new::<millimeter>;
    Contour::planar(
        mm(7.0),
        mm(4.0),
        mm(13.3),
        mm(20.0),
        Angle::new::<degree>(3.4),
    )
}

/// Moses and Stein's nozzle with its geometric 68.6 cm throat arc to 7 cm.
fn moses_stein() -> anyhow::Result<Contour> {
    let mm = Length::new::<millimeter>;
    let r = 0.686;
    let rise = r - (r * r - 0.07f64 * 0.07).sqrt();
    Contour::planar(
        mm(10.0),
        mm(10.0),
        Length::new::<meter>(0.01 + 2.0 * rise),
        Length::new::<meter>(r),
        Angle::new::<degree>(15.0),
    )
}

fn kelvin_t(t: f64) -> ThermodynamicTemperature {
    ThermodynamicTemperature::new::<kelvin>(t)
}

/// Moses and Stein case 252: 40.05 kPa, 374.3 K.
fn steam_252(model: CondensationModel) -> anyhow::Result<CondensingFlow> {
    CondensingNozzle::new(Water, moses_stein()?, model)
        .expand(Pressure::new::<kilopascal>(40.05), kelvin_t(374.3))
}

#[test]
fn test_condensing_frozen_limit_matches_the_real_fluid_nozzle() -> anyhow::Result<()> {
    // Without nucleation the march is the isentropic expansion of frozen
    // vapour, which the real-fluid nozzle solves independently.
    let contour = yale()?;
    let (p0, t0) = (Pressure::new::<torr>(463.0), kelvin_t(113.0));
    let flow = CondensingNozzle::new(PureFluid::Nitrogen, contour, CondensationModel::frozen())
        .expand(p0, t0)?;
    let fluid = PureFluid::Nitrogen;
    let a0 = fluid.vapour_variable(p0.get::<pascal>(), 113.0)?;
    let v0 = fluid.vapour(a0, 113.0);
    let chamber = Stagnation {
        p: p0.get::<pascal>(),
        h: v0.h,
        s: v0.s,
    };
    let point = Nozzle::new(contour.throat_area(), contour.exit_area())?.real_operating_point(
        &FrozenVapour::new(fluid),
        &chamber,
        Pressure::new::<pascal>(0.0),
    )?;
    let exit = flow.exit().ok_or(PhysicsError::InvalidState("no exit"))?;
    let close = |a: f64, b: f64, what: &str| {
        anyhow::ensure!((a / b - 1.0).abs() < 2e-5, "{what}: {a} vs {b}");
        Ok(())
    };
    close(
        flow.mass_flow().get::<kilogram_per_second>(),
        point.mass_flow.get::<kilogram_per_second>(),
        "mass flow",
    )?;
    close(
        exit.pressure.get::<pascal>(),
        point.exit_pressure.get::<pascal>(),
        "exit pressure",
    )?;
    close(
        exit.temperature.get::<kelvin>(),
        point.exit_temperature.get::<kelvin>(),
        "exit temperature",
    )?;
    let thrust = flow
        .thrust(Pressure::new::<pascal>(0.0))
        .ok_or(PhysicsError::InvalidState("thrust"))?;
    close(
        thrust.get::<newton>(),
        point.thrust.get::<newton>(),
        "thrust",
    )?;
    anyhow::ensure!(
        exit.liquid_fraction == 0.0 && exit.droplets_per_kilogram == 0.0,
        "dry"
    );
    Ok(())
}

#[test]
fn test_condensing_march_conserves_mass_energy_and_obeys_the_second_law() -> anyhow::Result<()> {
    let flow = steam_252(CondensationModel::default())?;
    let samples = flow.samples();
    let first = samples
        .first()
        .ok_or(PhysicsError::InvalidState("no samples"))?;
    let exit = flow.exit().ok_or(PhysicsError::InvalidState("no exit"))?;
    anyhow::ensure!(
        exit.liquid_fraction > 0.02,
        "condensed {}",
        exit.liquid_fraction
    );
    let contour = moses_stein()?;
    let m = flow.mass_flow().get::<kilogram_per_second>();
    let s0 = first.entropy.value;
    let mut previous = s0;
    for s in samples {
        let area = contour
            .area(s.position)
            .get::<crate::units::area::square_meter>();
        let flux = s.density.value * s.velocity.value * area;
        anyhow::ensure!(
            (flux / m - 1.0).abs() < 1e-9,
            "mass at {:?}: {flux}",
            s.position
        );
        // Entropy never falls: nucleation and growth across a finite
        // subcooling are irreversible.
        anyhow::ensure!(
            s.entropy.value >= previous - 1e-6,
            "entropy falls at {:?}",
            s.position
        );
        previous = s.entropy.value;
    }
    anyhow::ensure!(previous > s0 + 1.0, "no entropy produced");
    Ok(())
}

#[test]
fn test_condensing_faster_growth_approaches_the_reversible_limit() -> anyhow::Result<()> {
    // Faster growth relaxes the subcooling sooner: a lower peak, less entropy
    // produced, and a condensed fraction that converges. (Nucleation
    // quenches itself sooner too, so fewer, larger droplets carry it.)
    let run = |factor: f64| {
        steam_252(CondensationModel {
            growth_factor: factor,
            ..CondensationModel::default()
        })
    };
    let runs = [run(1.0)?, run(30.0)?, run(1000.0)?];
    let mut last: Option<(f64, f64, f64, f64)> = None;
    for flow in &runs {
        let exit = flow.exit().ok_or(PhysicsError::InvalidState("no exit"))?;
        let peak = flow
            .samples()
            .iter()
            .map(|s| s.subcooling.value)
            .fold(f64::MIN, f64::max);
        let now = (
            peak,
            exit.entropy.value,
            exit.liquid_fraction,
            exit.droplets_per_kilogram,
        );
        if let Some((peak0, s0, _, n0)) = last {
            anyhow::ensure!(
                now.0 < peak0 && now.1 < s0 && now.3 < n0,
                "{now:?} after {last:?}"
            );
        }
        last = Some(now);
    }
    let y = |i: usize| {
        runs.get(i)
            .and_then(CondensingFlow::exit)
            .map_or(f64::NAN, |e| e.liquid_fraction)
    };
    anyhow::ensure!(
        (y(2) / y(1) - 1.0).abs() < 0.01,
        "condensed {} vs {}",
        y(2),
        y(1)
    );
    Ok(())
}

#[test]
fn test_condensing_droplets_follow_the_classes() -> anyhow::Result<()> {
    // Nucleation ends before the exit: the number is then fixed, and the
    // Sauter radius holds the liquid.
    let flow = steam_252(CondensationModel::default())?;
    let exit = flow.exit().ok_or(PhysicsError::InvalidState("no exit"))?;
    anyhow::ensure!(exit.nucleation_rate < 1.0, "still nucleating");
    let r = exit.sauter_radius.get::<meter>();
    anyhow::ensure!(r > 5e-9 && r < 50e-9, "Sauter radius {r}");
    let wilson = flow
        .wilson_point()
        .ok_or(PhysicsError::InvalidState("no Wilson point"))?;
    let onset = flow
        .onset(1e-4)
        .ok_or(PhysicsError::InvalidState("no onset"))?;
    anyhow::ensure!(wilson > onset, "Wilson {wilson:?} before onset {onset:?}");
    Ok(())
}

#[test]
fn test_condensing_shock_conserves_and_evaporates() -> anyhow::Result<()> {
    let nozzle = CondensingNozzle::new(Water, moses_stein()?, CondensationModel::default());
    let flow = nozzle.expand_with_shock(
        Pressure::new::<kilopascal>(40.05),
        kelvin_t(374.3),
        Length::new::<meter>(0.05),
    )?;
    anyhow::ensure!(
        flow.shock()
            .is_some_and(|x| (x.get::<meter>() - 0.05).abs() < 1e-12),
        "shock"
    );
    let samples = flow.samples();
    let at = samples
        .windows(2)
        .position(|w| w[0].position == w[1].position)
        .ok_or(PhysicsError::InvalidState("no shock in the samples"))?;
    let (before, after) = (&samples[at], &samples[at + 1]);
    anyhow::ensure!(before.frozen_mach > 1.0 && after.frozen_mach < 1.0, "Mach");
    anyhow::ensure!(
        (after.impulse / before.impulse - 1.0).abs() < 1e-9,
        "impulse"
    );
    anyhow::ensure!(
        (after.liquid_fraction - before.liquid_fraction).abs() < f64::EPSILON,
        "droplets pass unchanged"
    );
    anyhow::ensure!(
        after.pressure > before.pressure && after.entropy > before.entropy,
        "jump"
    );
    // Behind the shock the vapour is superheated and the droplets evaporate.
    let exit = flow.exit().ok_or(PhysicsError::InvalidState("no exit"))?;
    anyhow::ensure!(
        after.subcooling.value < 0.0,
        "superheated {:?}",
        after.subcooling
    );
    anyhow::ensure!(exit.liquid_fraction < after.liquid_fraction, "evaporation");
    anyhow::ensure!(exit.frozen_mach < 1.0, "subsonic exit");
    Ok(())
}

#[test]
fn test_condensing_back_pressure_places_the_shock() -> anyhow::Result<()> {
    let nozzle = CondensingNozzle::new(Water, moses_stein()?, CondensationModel::default());
    let (p0, t0) = (Pressure::new::<kilopascal>(40.05), kelvin_t(374.3));
    let pb = Pressure::new::<kilopascal>(25.0);
    let flow = nozzle.expand_against(p0, t0, pb)?;
    let shock = flow.shock().ok_or(PhysicsError::InvalidState("no shock"))?;
    anyhow::ensure!(
        shock.get::<meter>() > 0.0 && shock.get::<meter>() < 0.07,
        "{shock:?}"
    );
    let exit = flow.exit().ok_or(PhysicsError::InvalidState("no exit"))?;
    let error = (exit.pressure / pb).value - 1.0;
    anyhow::ensure!(error.abs() < 1e-5, "exit pressure off by {error}");
    // Into a low back pressure the flow stays supersonic; too high a one
    // unchokes the nozzle, which the march does not model.
    let free = nozzle.expand_against(p0, t0, Pressure::new::<kilopascal>(1.0))?;
    anyhow::ensure!(free.shock().is_none(), "supersonic exit");
    anyhow::ensure!(
        nozzle
            .expand_against(p0, t0, Pressure::new::<kilopascal>(39.0))
            .is_err(),
        "unchoked"
    );
    Ok(())
}

#[test]
fn test_condensing_nozzle_rejects_what_it_cannot_march() -> anyhow::Result<()> {
    let nozzle = CondensingNozzle::new(Water, moses_stein()?, CondensationModel::default());
    let p0 = Pressure::new::<kilopascal>(40.05);
    // A reservoir of liquid water, and shocks outside the nozzle.
    anyhow::ensure!(
        nozzle.expand(p0, kelvin_t(300.0)).is_err()
            || nozzle
                .expand(Pressure::new::<kilopascal>(-1.0), kelvin_t(374.3))
                .is_err(),
        "liquid"
    );
    for x in [-0.01, 0.0, 0.07, 0.2] {
        anyhow::ensure!(
            nozzle
                .expand_with_shock(p0, kelvin_t(374.3), Length::new::<meter>(x))
                .is_err(),
            "shock at {x}"
        );
    }
    Ok(())
}
