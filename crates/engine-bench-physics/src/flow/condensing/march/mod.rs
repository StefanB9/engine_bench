//! The condensing nozzle: steady quasi-1D flow of a vapour that nucleates
//! and grows droplets, marched downstream from just past the sonic throat.
//!
//! Mass flow, stagnation enthalpy and the impulse function F = pA + ṁu
//! (dF/dx = p·dA/dx, inviscid) are conserved exactly at every station; the
//! ODE in x carries F and the droplet classes. A normal shock conserves the
//! same three, so the state behind it is the subsonic root at the same F,
//! with the droplets passing through unchanged; downstream they evaporate
//! by the same growth law.

mod classes;
mod output;
mod start;
mod state;
mod system;

use anyhow::Context;

pub use self::output::{CondensingFlow, FlowSample};
use self::{
    classes::{CLASSES, Classes},
    start::{Reservoir, sonic_throat, supersonic_station},
    state::{Conserved, Station, solve},
    system::{FIRST_CLASS, IMPULSE, Local, March, N, OPEN_LIQUID, OPEN_NUMBER},
};
use super::{CondensationModel, CondensingFluid, Contour, nucleation::BOLTZMANN};
use crate::{
    PhysicsError,
    numerics::{
        ode::{AdaptiveOptions, DormandPrince45, NoEvents, StepOutcome},
        roots::{RootOptions, brent},
    },
    units::{
        Length, Pressure, ThermodynamicTemperature, length::meter, mass::kilogram,
        pressure::pascal, thermodynamic_temperature::kelvin,
    },
};

/// The march starts where the area is this much above the throat's: close
/// enough that nothing nucleates before, far enough from the sonic point's
/// singularity.
const START_AREA_RATIO: f64 = 1.001;
/// Relative and absolute (scaled state) tolerances of the march.
const RTOL: f64 = 1e-7;
const ATOL: f64 = 1e-9;
/// Steps are at most this fraction of the nozzle.
const MAX_STEP_FRACTION: f64 = 0.01;
/// Droplets per kg below which a population is negligible.
const NUMBER_FLOOR: f64 = 1e6;
/// A march that stops within this of frozen Mach 1 has choked: heat
/// release (or, behind a shock, area growth) drove it to the sonic point.
const CHOKE_BAND: f64 = 0.05;
/// Bisection of a condensation shock's position stops at this fraction of
/// the nozzle length.
const SHOCK_TOLERANCE: f64 = 1e-8;
/// Offsets past a subsonic choke, as fractions of the nozzle length, at
/// which the flow is tried on its supersonic branch again: the first whose
/// supersonic march reaches the exit is taken.
const REEXPANSION_OFFSETS: [f64; 8] = [0.0, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2];

/// How a march ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MarchEnd {
    /// It reached where it was sent.
    Reached,
    /// It stopped at the sonic point.
    Choked,
}

/// A nozzle contour with a condensing fluid and its condensation model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CondensingNozzle<F> {
    fluid: F,
    contour: Contour,
    model: CondensationModel,
}

/// One accepted step's full state, to restart the march from.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Checkpoint {
    x: f64,
    y: [f64; N],
    classes: Classes,
    number_scale: f64,
    station: Station,
}

impl<F: CondensingFluid> CondensingNozzle<F> {
    /// `fluid` through `contour`, condensing by `model`.
    pub fn new(fluid: F, contour: Contour, model: CondensationModel) -> Self {
        Self {
            fluid,
            contour,
            model,
        }
    }

    /// Fully supersonic flow from vapour at rest at `p0`, `t0` to the exit.
    pub fn expand(
        &self,
        p0: Pressure,
        t0: ThermodynamicTemperature,
    ) -> anyhow::Result<CondensingFlow> {
        let mut run = self.begin(p0, t0)?;
        let end = self.contour.length_si();
        if run.march(end, true)? == MarchEnd::Choked {
            let (run, x) = run.condensation_shock()?;
            return Ok(run.finish(None, Some(x)));
        }
        Ok(run.finish(None, None))
    }

    /// Flow with a normal shock at `shock` downstream of the throat, the
    /// droplets passing through it and evaporating behind.
    pub fn expand_with_shock(
        &self,
        p0: Pressure,
        t0: ThermodynamicTemperature,
        shock: Length,
    ) -> anyhow::Result<CondensingFlow> {
        let x_shock = shock.get::<meter>();
        let mut run = self.begin(p0, t0)?;
        let end = self.contour.length_si();
        if !(x_shock > run.start && x_shock < end) {
            return Err(PhysicsError::InvalidGeometry(
                "shock outside the diverging section",
            ))
            .with_context(|| format!("shock at {x_shock} m, nozzle {end} m"));
        }
        run.march_through(x_shock)?;
        run.shock()?;
        run.march_through(end)?;
        Ok(run.finish(Some(x_shock), None))
    }

    /// Flow exhausting against `back_pressure`: fully supersonic if the exit
    /// pressure behind a shock at the exit would not reach it, else with the
    /// shock where the subsonic exit pressure meets it.
    pub fn expand_against(
        &self,
        p0: Pressure,
        t0: ThermodynamicTemperature,
        back_pressure: Pressure,
    ) -> anyhow::Result<CondensingFlow> {
        let pb = back_pressure.get::<pascal>();
        let mut run = self.begin(p0, t0)?;
        let end = self.contour.length_si();
        let mut condensation = None;
        let mut first = run.start * 1.5;
        if run.march(end, true)? == MarchEnd::Choked {
            let (shocked, x) = run.condensation_shock()?;
            run = shocked;
            condensation = Some(x);
            // A back-pressure shock can then stand only where the flow is
            // supersonic again, downstream of the sonic point.
            first = run
                .reexpansion()
                .context("back pressure behind a condensation shock")?;
        }
        let behind_exit = run.exit_pressure_with_shock_at(end)?;
        if pb <= behind_exit {
            return Ok(run.finish(None, condensation));
        }
        if condensation.is_some() {
            return Err(PhysicsError::InvalidState(
                "back pressure shock with a condensation shock upstream is not modelled",
            ))
            .with_context(|| format!("back pressure {pb} Pa, exit limit {behind_exit} Pa"));
        }
        let behind_first = run.exit_pressure_with_shock_at(first)?;
        if pb >= behind_first {
            return Err(PhysicsError::InvalidState(
                "back pressure too high for a choked nozzle (subsonic throughout)",
            ))
            .with_context(|| format!("back pressure {pb} Pa, limit {behind_first} Pa"));
        }
        let mut failure = None;
        let x_shock = brent(
            |x| match run.exit_pressure_with_shock_at(x) {
                Ok(p) => p - pb,
                Err(e) => {
                    failure.get_or_insert(e);
                    f64::NAN
                }
            },
            first,
            end,
            RootOptions::new(1e-7 * end, 100)?,
        );
        if let Some(e) = failure {
            return Err(e).context("shock position");
        }
        let x_shock = x_shock.context("shock position")?;
        self.expand_with_shock(p0, t0, Length::new::<meter>(x_shock))
    }

    /// The reservoir, the frozen throat and the march's first station.
    fn begin(&self, p0: Pressure, t0: ThermodynamicTemperature) -> anyhow::Result<Run<'_, F>> {
        let (p0, t0) = (p0.get::<pascal>(), t0.get::<kelvin>());
        let reservoir = Reservoir::new(&self.fluid, p0, t0)?;
        let molecule = self.fluid.molecular_mass().get::<kilogram>();
        let gas_constant = BOLTZMANN / molecule;
        let (_, flux) = sonic_throat(&self.fluid, &reservoir)?;
        let throat_area = self.contour.area_si(0.0).0;
        let mass_flow = flux * throat_area;
        let start = start_position(&self.contour, throat_area * START_AREA_RATIO);
        let station = supersonic_station(
            &self.fluid,
            &reservoir,
            flux,
            self.contour.area_si(start).0 / throat_area,
            gas_constant,
        )?;
        let area = self.contour.area_si(start).0;
        let impulse = station.vapour.p * area + mass_flow * station.u;
        let mut y = [0.0; N];
        y[IMPULSE] = 1.0;
        let march = March {
            fluid: &self.fluid,
            contour: &self.contour,
            model: self.model,
            conserved: Conserved {
                mass_flow,
                h0: reservoir.h0,
                gas_constant,
            },
            molecule,
            impulse_scale: impulse,
            number_scale: NUMBER_FLOOR,
            classes: Classes::default(),
            last: station,
        };
        let mut run = Run {
            march,
            start,
            checkpoints: Vec::with_capacity(512),
            samples: Vec::with_capacity(512),
            reexpanded: None,
        };
        run.checkpoints.push(Checkpoint {
            x: start,
            y,
            classes: Classes::default(),
            number_scale: NUMBER_FLOOR,
            station,
        });
        run.record(start, &y)?;
        Ok(run)
    }
}

/// The vapour's natural variable where its pressure is `p` Pa at `t` K, by
/// Newton's method on p(a, T) from `start`; `None` if it does not converge.
fn variable_at<F: CondensingFluid>(fluid: &F, p: f64, t: f64, start: f64) -> Option<f64> {
    let mut a = start;
    for _ in 0..60 {
        let v = fluid.vapour(a, t);
        let step = (v.p - p) / v.dp_da;
        if !step.is_finite() || v.dp_da <= 0.0 {
            return None;
        }
        let bounded = step.clamp(-0.5 * a.abs(), 0.5 * a.abs());
        a -= bounded;
        if step.abs() <= 1e-13 * a.abs() {
            return Some(a);
        }
    }
    None
}

/// Where the area reaches `area`, by bisection on the growing contour.
fn start_position(contour: &Contour, area: f64) -> f64 {
    let (mut lo, mut hi) = (0.0, contour.length_si());
    for _ in 0..200 {
        let mid = f64::midpoint(lo, hi);
        if contour.area_si(mid).0 < area {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    hi
}

/// A march in progress.
struct Run<'a, F> {
    march: March<'a, F>,
    start: f64,
    checkpoints: Vec<Checkpoint>,
    samples: Vec<FlowSample>,
    /// Where the flow re-expanded through the sonic point behind a
    /// condensation shock.
    reexpanded: Option<f64>,
}

impl<F: CondensingFluid> Run<'_, F> {
    fn current(&self) -> anyhow::Result<&Checkpoint> {
        self.checkpoints
            .last()
            .ok_or_else(|| PhysicsError::InvalidState("march has no state").into())
    }

    /// Marches to `x_end`; a choke is an error.
    fn march_through(&mut self, x_end: f64) -> anyhow::Result<()> {
        match self.march(x_end, true)? {
            MarchEnd::Reached => Ok(()),
            MarchEnd::Choked => Err(PhysicsError::InvalidState(
                "the flow chokes by heat addition before the requested shock",
            )
            .into()),
        }
    }

    /// Whether the last station is within [`CHOKE_BAND`] of frozen Mach 1.
    fn near_sonic(&self) -> bool {
        let s = &self.march.last;
        (s.u / s.frozen_sound_speed() - 1.0).abs() < CHOKE_BAND
    }

    /// Marches from the last checkpoint to `x_end`, recording every step if
    /// `record`; stops short where the flow chokes.
    fn march(&mut self, x_end: f64, record: bool) -> anyhow::Result<MarchEnd> {
        let mut at = *self.current()?;
        self.march.classes = at.classes;
        self.march.number_scale = at.number_scale;
        self.march.last = at.station;
        let length = self.march.contour.length_si();
        let options =
            AdaptiveOptions::new(RTOL, ATOL)?.with_max_step(MAX_STEP_FRACTION * length)?;
        let mut solver = DormandPrince45::<N>::new(options);
        let (mut x, mut y) = (at.x, at.y);
        while x < x_end {
            let stepped = solver.step(&mut self.march, &mut NoEvents, &mut x, &mut y, x_end);
            if stepped.is_err() && self.near_sonic() {
                self.march.last = at.station;
                return Ok(MarchEnd::Choked);
            }
            let outcome = stepped.with_context(|| {
                let s = &self.march.last;
                format!(
                    "condensing march at x = {x:.6} m (vapour at {:.2} K, {:.1} Pa, frozen Mach \
                     {:.4}, liquid {:.4})",
                    s.t,
                    s.vapour.p,
                    s.u / s.frozen_sound_speed(),
                    s.y
                )
            })?;
            if self.close_open_class(&mut y) {
                solver.forget_derivative();
            }
            // The station at the accepted point, for the next restart.
            let local = self.march.local(x, &y, None);
            if local.is_none() && self.near_sonic() {
                self.march.last = at.station;
                return Ok(MarchEnd::Choked);
            }
            let local = local.ok_or_else(|| {
                let s = &self.march.last;
                anyhow::Error::from(PhysicsError::InvalidState(
                    "no flow state or condensation properties: beyond the fluid's property range",
                ))
                .context(format!(
                    "condensing march at x = {x:.6} m (vapour at {:.2} K, {:.1} Pa)",
                    s.t, s.vapour.p
                ))
            })?;
            at = Checkpoint {
                x,
                y,
                classes: self.march.classes,
                number_scale: self.march.number_scale,
                station: local.station,
            };
            self.checkpoints.push(at);
            if record {
                self.samples
                    .push(FlowSample::from_local(x, &local, self.march.contour));
            }
            if outcome == StepOutcome::ReachedEnd {
                break;
            }
        }
        Ok(MarchEnd::Reached)
    }

    /// Closes the open class into the list (merging if full) and drops
    /// evaporated classes; whether the state changed.
    fn close_open_class(&mut self, y: &mut [f64; N]) -> bool {
        let number = y[OPEN_NUMBER] * self.march.number_scale;
        let open_liquid = y[OPEN_LIQUID];
        let mut fractions = [0.0; CLASSES];
        for (k, f) in fractions.iter_mut().enumerate() {
            *f = y.get(FIRST_CLASS + k).copied().unwrap_or(0.0);
        }
        let had = self.march.classes.count();
        let mut changed = false;
        if number > 0.0 && open_liquid > 0.0 {
            let rho_l = self.march.last.droplets.rho_liquid;
            self.march
                .classes
                .push(number, open_liquid, &mut fractions, rho_l);
            changed = true;
        }
        self.march.classes.drop_empty(&mut fractions);
        changed |= self.march.classes.count() != had + usize::from(changed);
        if changed {
            y[OPEN_NUMBER] = 0.0;
            y[OPEN_LIQUID] = 0.0;
            for (k, f) in fractions.iter().enumerate() {
                if let Some(slot) = y.get_mut(FIRST_CLASS + k) {
                    *slot = *f;
                }
            }
            self.march.number_scale = self.march.classes.total().max(NUMBER_FLOOR);
        }
        changed
    }

    fn record(&mut self, x: f64, y: &[f64; N]) -> anyhow::Result<()> {
        let local = self
            .march
            .local(x, y, None)
            .ok_or(PhysicsError::InvalidState(
                "station at the start of the march",
            ))?;
        self.samples
            .push(FlowSample::from_local(x, &local, self.march.contour));
        Ok(())
    }

    /// Jumps to the subsonic root at the current station: a normal shock.
    #[expect(
        clippy::many_single_char_names,
        reason = "the symbols of the jump: s, p, t, a, g"
    )]
    fn shock(&mut self) -> anyhow::Result<()> {
        let mut at = *self.current()?;
        let s = at.station;
        let mach = s.u / s.frozen_sound_speed();
        if mach.is_nan() || mach <= 1.0 {
            return Err(PhysicsError::InvalidState("shock in subsonic flow"))
                .with_context(|| format!("frozen Mach {mach}"));
        }
        // Perfect-gas jump of the vapour as the starting point.
        let g = s.vapour.cp / s.vapour.cv;
        let m2 = mach * mach;
        let p_ratio = 1.0 + 2.0 * g / (g + 1.0) * (m2 - 1.0);
        let t_ratio = p_ratio * (2.0 + (g - 1.0) * m2) / ((g + 1.0) * m2);
        let (p, t) = (s.vapour.p * p_ratio, s.t * t_ratio);
        // From the state ahead along p(a, T) (the density solve trusts less
        // than the march's states reach below the triple point).
        let a = variable_at(self.march.fluid, p, t, s.a).ok_or(PhysicsError::InvalidState(
            "no vapour state for the jump's guess",
        ))?;
        let area = self.march.contour.area_si(at.x).0;
        let vapour = self.march.fluid.vapour(a, t);
        let guess = Station {
            a,
            t,
            vapour,
            // Subsonic, as the root sought.
            u: self.march.conserved.mass_flow / (vapour.rho * area),
            ..s
        };
        let impulse = at.y[IMPULSE] * self.march.impulse_scale;
        let behind = solve(
            self.march.fluid,
            &self.march.conserved,
            area,
            impulse,
            s.y,
            &guess,
        )
        .ok_or(PhysicsError::NoConvergence {
            iterations: 60,
            residual: f64::NAN,
        })
        .context("state behind the shock")?;
        if behind.u.is_nan() || behind.u >= behind.frozen_sound_speed() {
            return Err(PhysicsError::InvalidState("no subsonic state behind the shock").into());
        }
        at.station = behind;
        self.march.last = behind;
        self.checkpoints.push(at);
        let local = self
            .march
            .local(at.x, &at.y, None)
            .ok_or(PhysicsError::InvalidState("station behind the shock"))?;
        self.samples
            .push(FlowSample::from_local(at.x, &local, self.march.contour));
        Ok(())
    }

    /// A copy of this run cut back to its last checkpoint at or before `x`,
    /// without samples.
    fn trial(&self, x: f64) -> anyhow::Result<Self> {
        let index = self
            .checkpoints
            .iter()
            .rposition(|c| c.x <= x)
            .ok_or(PhysicsError::InvalidState("shock upstream of the march"))?;
        Ok(Run {
            march: March {
                classes: Classes::default(),
                ..self.march.clone_shallow()
            },
            start: self.start,
            checkpoints: self
                .checkpoints
                .get(..=index)
                .map_or_else(Vec::new, <[_]>::to_vec),
            samples: Vec::new(),
            reexpanded: None,
        })
    }

    /// Exit pressure, Pa, with a shock at `x_shock`, from the supersonic
    /// march's checkpoints (the run itself is left unchanged).
    fn exit_pressure_with_shock_at(&self, x_shock: f64) -> anyhow::Result<f64> {
        let mut trial = self.trial(x_shock)?;
        trial.march_unrecorded(x_shock)?;
        trial.shock()?;
        let end = self.march.contour.length_si();
        trial.march_unrecorded(end)?;
        Ok(trial.current()?.station.vapour.p)
    }

    /// Marches without recording; a choke is an error.
    fn march_unrecorded(&mut self, x_end: f64) -> anyhow::Result<()> {
        match self.march(x_end, false)? {
            MarchEnd::Reached => Ok(()),
            MarchEnd::Choked => Err(PhysicsError::InvalidState("the flow chokes").into()),
        }
    }

    /// With the supersonic march choked by heat addition, the quasi-steady
    /// condensation shock: its position, bisected between shocks behind
    /// which the flow decelerates to a subsonic exit (upstream) and shocks
    /// behind which it chokes again (downstream), and the flow re-expanding
    /// from the sonic point it then reaches. Returns the recorded run and
    /// the shock's position, m.
    fn condensation_shock(self) -> anyhow::Result<(Self, f64)> {
        let end = self.march.contour.length_si();
        let choked_at = self.current()?.x;
        // Upstream bound: the last station still dry, where a shock stops
        // condensation altogether.
        let mut lo = self
            .checkpoints
            .iter()
            .rev()
            .find(|c| c.station.y <= 0.0 && c.x > self.start)
            .map_or(self.start * 1.5, |c| c.x);
        // Downstream bound: the last station clearly supersonic, so that the
        // shock's subsonic root stands apart from the supersonic one.
        let mut hi = self
            .checkpoints
            .iter()
            .rev()
            .find(|c| c.station.u > (1.0 + 0.4 * CHOKE_BAND) * c.station.frozen_sound_speed())
            .map_or(choked_at, |c| c.x);
        let outcome = |x: f64| -> anyhow::Result<MarchEnd> {
            let mut trial = self.trial(x)?;
            trial.march_unrecorded(x)?;
            trial.shock()?;
            trial.march(end, false)
        };
        if outcome(lo)? != MarchEnd::Reached || outcome(hi)? != MarchEnd::Choked {
            return Err(PhysicsError::InvalidState(
                "no quasi-steady condensation shock brackets the choke",
            ))
            .with_context(|| format!("choked at x = {choked_at} m"));
        }
        while hi - lo > SHOCK_TOLERANCE * end {
            let mid = f64::midpoint(lo, hi);
            if outcome(mid)? == MarchEnd::Choked {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        // The recorded flow: up to the shock, behind it to the sonic point,
        // across to the supersonic branch (as close past the saddle as lets
        // the supersonic flow reach the exit) and on.
        let recorded = || -> anyhow::Result<Self> {
            let mut run = self.trial(hi)?;
            let kept = run.current()?.x;
            run.samples = self
                .samples
                .iter()
                .copied()
                .filter(|s| s.position.get::<meter>() <= kept)
                .collect();
            if run.march(hi, true)? == MarchEnd::Choked {
                return Err(
                    PhysicsError::InvalidState("choked before the condensation shock").into(),
                );
            }
            run.shock()?;
            Ok(run)
        };
        let mut behind = recorded()?;
        if behind.march(end, true)? == MarchEnd::Reached {
            return Ok((behind, hi));
        }
        for offset in REEXPANSION_OFFSETS {
            let mut run = behind.clone_run();
            if run.supersonic_again(offset).is_ok() && run.march(end, true)? == MarchEnd::Reached {
                return Ok((run, hi));
            }
        }
        Err(PhysicsError::InvalidState(
            "the flow chokes again past the condensation shock",
        ))
        .context("re-expansion")
    }

    /// A full copy of this run.
    fn clone_run(&self) -> Self {
        Self {
            march: self.march.clone_shallow(),
            start: self.start,
            checkpoints: self.checkpoints.clone(),
            samples: self.samples.clone(),
            reexpanded: self.reexpanded,
        }
    }

    /// Where the flow passed back to supersonic after a condensation shock,
    /// m.
    fn reexpansion(&self) -> anyhow::Result<f64> {
        self.reexpanded
            .ok_or_else(|| PhysicsError::InvalidState("no re-expansion").into())
    }

    /// From a subsonic station at the sonic point, the supersonic root at
    /// the same x: past the saddle where heat release and area growth
    /// balance, the flow continues supersonic (the two roots meet there, so
    /// the jump is as small as the bisection is fine).
    #[expect(
        clippy::many_single_char_names,
        reason = "the symbols of the jump: s, g, p, t, a, x"
    )]
    fn supersonic_again(&mut self, offset: f64) -> anyhow::Result<()> {
        let at = *self.current()?;
        let s = at.station;
        let below = s.u / s.frozen_sound_speed();
        let above = 2.0 - below.min(0.995);
        let g = s.vapour.cp / s.vapour.cv;
        let stagnation = |m: f64| 1.0 + 0.5 * (g - 1.0) * m * m;
        let p = s.vapour.p * (stagnation(below) / stagnation(above)).powf(g / (g - 1.0));
        let t = s.t * stagnation(below) / stagnation(above);
        // From the state ahead along p(a, T) (the density solve trusts less
        // than the march's states reach below the triple point).
        let a = variable_at(self.march.fluid, p, t, s.a).ok_or(PhysicsError::InvalidState(
            "no vapour state for the jump's guess",
        ))?;
        let vapour = self.march.fluid.vapour(a, t);
        let length = self.march.contour.length_si();
        let area = self.march.contour.area_si(at.x).0;
        let impulse = at.y[IMPULSE] * self.march.impulse_scale;
        // Where the subsonic flow choked the two roots merge; the supersonic
        // one exists just past the saddle, where area growth takes over.
        // Across that sliver the impulse gains p·dA; the droplets are held.
        {
            let x = at.x + offset * length;
            let new_area = self.march.contour.area_si(x).0;
            let shifted = impulse + s.vapour.p * (new_area - area);
            let guess = Station {
                a,
                t,
                vapour,
                u: self.march.conserved.mass_flow / (vapour.rho * new_area),
                ..s
            };
            let ahead = solve(
                self.march.fluid,
                &self.march.conserved,
                new_area,
                shifted,
                s.y,
                &guess,
            )
            .ok_or(PhysicsError::NoConvergence {
                iterations: 60,
                residual: f64::NAN,
            })
            .context("supersonic state past the sonic point")?;
            let mut next = at;
            next.x = x;
            next.y[IMPULSE] = shifted / self.march.impulse_scale;
            next.station = ahead;
            self.march.last = ahead;
            self.checkpoints.push(next);
            self.reexpanded = Some(x);
            let local = self
                .march
                .local(x, &next.y, None)
                .ok_or(PhysicsError::InvalidState("station past the sonic point"))?;
            self.samples
                .push(FlowSample::from_local(x, &local, self.march.contour));
        }
        Ok(())
    }

    fn finish(self, shock: Option<f64>, condensation_shock: Option<f64>) -> CondensingFlow {
        CondensingFlow::new(
            self.samples,
            self.march.conserved.mass_flow,
            self.march.contour.area_si(self.march.contour.length_si()).0,
            shock,
            condensation_shock,
        )
    }
}

impl<F> March<'_, F> {
    fn clone_shallow(&self) -> Self {
        Self {
            fluid: self.fluid,
            contour: self.contour,
            model: self.model,
            conserved: self.conserved,
            molecule: self.molecule,
            impulse_scale: self.impulse_scale,
            number_scale: self.number_scale,
            classes: self.classes,
            last: self.last,
        }
    }
}

impl Local {
    /// Sauter mean radius, m, of all droplets (zero without).
    fn sauter_radius(&self) -> f64 {
        if self.area_moment > 0.0 {
            self.volume_moment / self.area_moment
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests;
