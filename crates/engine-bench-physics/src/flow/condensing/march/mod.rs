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
        run.march(end, true)?;
        Ok(run.finish(None))
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
        run.march(x_shock, true)?;
        run.shock()?;
        run.march(end, true)?;
        Ok(run.finish(Some(x_shock)))
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
        run.march(end, true)?;
        let behind_exit = run.exit_pressure_with_shock_at(end)?;
        if pb <= behind_exit {
            return Ok(run.finish(None));
        }
        let first = run.start * 1.5;
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
}

impl<F: CondensingFluid> Run<'_, F> {
    fn current(&self) -> anyhow::Result<&Checkpoint> {
        self.checkpoints
            .last()
            .ok_or_else(|| PhysicsError::InvalidState("march has no state").into())
    }

    /// Marches from the last checkpoint to `x_end`, recording every step if
    /// `record`.
    fn march(&mut self, x_end: f64, record: bool) -> anyhow::Result<()> {
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
            let outcome = solver
                .step(&mut self.march, &mut NoEvents, &mut x, &mut y, x_end)
                .with_context(|| {
                    let s = &self.march.last;
                    format!(
                        "condensing march at x = {x:.6} m (frozen Mach {:.4}, liquid {:.4})",
                        s.u / s.frozen_sound_speed(),
                        s.y
                    )
                })?;
            if self.close_open_class(&mut y) {
                solver.forget_derivative();
            }
            // The station at the accepted point, for the next restart.
            let local = self
                .march
                .local(x, &y, None)
                .ok_or(PhysicsError::InvalidState("station after an accepted step"))?;
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
        Ok(())
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
        let a = self
            .march
            .fluid
            .vapour_variable(p, t)
            .context("vapour behind the shock")?;
        let guess = Station {
            a,
            t,
            vapour: self.march.fluid.vapour(a, t),
            ..s
        };
        let area = self.march.contour.area_si(at.x).0;
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

    /// Exit pressure, Pa, with a shock at `x_shock`, from the supersonic
    /// march's checkpoints (the run itself is left unchanged).
    fn exit_pressure_with_shock_at(&self, x_shock: f64) -> anyhow::Result<f64> {
        let index = self
            .checkpoints
            .iter()
            .rposition(|c| c.x <= x_shock)
            .ok_or(PhysicsError::InvalidState("shock upstream of the march"))?;
        let mut trial = Run {
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
        };
        trial.march(x_shock, false)?;
        trial.shock()?;
        let end = self.march.contour.length_si();
        trial.march(end, false)?;
        Ok(trial.current()?.station.vapour.p)
    }

    fn finish(self, shock: Option<f64>) -> CondensingFlow {
        CondensingFlow::new(
            self.samples,
            self.march.conserved.mass_flow,
            self.march.contour.area_si(self.march.contour.length_si()).0,
            shock,
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
