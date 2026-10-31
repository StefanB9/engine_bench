//! The live simulation as a pure state machine: commands in, engine stepped by
//! a pacer, snapshots out. Time comes from a [`Clock`], so tests drive it with
//! a mock.

use anyhow::Context;
use engine_bench_core::{
    clock::Clock,
    engine::{AnyEngine, Engine, EngineId},
    pacer::{PacerSettings, RealTimePacer},
    spec::{ControlId, StepContext},
};
use engine_bench_physics::units::{Time, time::second};

/// Physics step of the live view, simulated seconds.
pub(crate) const PHYSICS_STEP: f64 = 1e-3;
/// At most this many physics steps per tick; beyond that the simulation falls
/// behind real time instead of freezing the machine.
const MAX_STEPS_PER_TICK: u32 = 250;

/// A request from the UI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Command {
    /// Switch to a fresh instance of another engine.
    SelectEngine(EngineId),
    /// Set a control (SI value).
    SetControl(ControlId, f64),
    /// Simulated seconds per wall second.
    SetTimeScale(f64),
    /// Stop advancing.
    Pause,
    /// Continue advancing.
    Resume,
    /// Advance one physics step (while paused).
    Step,
    /// Back to the initial state.
    Reset,
    /// Stop the simulation thread.
    Shutdown,
}

/// What the UI shows of the live simulation. Reused between ticks: refilling it
/// does not allocate once its vectors have reached the engine's sizes.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Snapshot {
    /// Engine being simulated.
    pub(crate) engine: EngineId,
    /// Simulation time, s.
    pub(crate) time: f64,
    /// Latest channel values (SI).
    pub(crate) sample: Vec<f64>,
    /// Current control values (SI).
    pub(crate) controls: Vec<f64>,
    /// Whether paused.
    pub(crate) paused: bool,
    /// Requested time scale.
    pub(crate) time_scale: f64,
    /// Achieved time scale, if measurable yet.
    pub(crate) achieved_time_scale: Option<f64>,
    /// Most recent error (rejected control, failed step), until the next reset.
    pub(crate) error: Option<String>,
    /// Plot rows the UI was too slow to take.
    pub(crate) dropped_rows: u64,
}

impl Snapshot {
    /// An empty snapshot of `engine`, before the first tick.
    pub(crate) fn empty(engine: EngineId) -> Self {
        Self {
            engine,
            time: 0.0,
            sample: Vec::new(),
            controls: Vec::new(),
            paused: false,
            time_scale: 1.0,
            achieved_time_scale: None,
            error: None,
            dropped_rows: 0,
        }
    }
}

/// Engine, pacer and clock.
pub(crate) struct Simulator<C: Clock> {
    clock: C,
    pacer: RealTimePacer,
    engine: AnyEngine,
    sample: Vec<f64>,
    error: Option<String>,
}

impl<C: Clock> Simulator<C> {
    /// A running simulation of `engine` at default controls.
    pub(crate) fn new(clock: C, engine: EngineId) -> anyhow::Result<Self> {
        let settings = PacerSettings::new(PHYSICS_STEP, MAX_STEPS_PER_TICK)?;
        let engine = AnyEngine::new(engine)?;
        let mut sim = Self {
            clock,
            pacer: RealTimePacer::new(settings),
            sample: vec![0.0; engine.channels().len()],
            engine,
            error: None,
        };
        sim.refresh_sample();
        Ok(sim)
    }

    /// Applies `command`; returns `false` once the simulation should stop.
    pub(crate) fn handle(&mut self, command: Command) -> bool {
        let result = match command {
            Command::SelectEngine(id) => self.select(id),
            Command::SetControl(id, value) => self.engine.set_control(id, value),
            Command::SetTimeScale(scale) => self.pacer.set_time_scale(scale),
            Command::Pause => {
                self.pacer.pause();
                Ok(())
            }
            Command::Resume => {
                self.pacer.resume();
                Ok(())
            }
            Command::Step => {
                self.pacer.request_step();
                Ok(())
            }
            Command::Reset => {
                self.error = None;
                self.engine.reset()
            }
            Command::Shutdown => return false,
        };
        if let Err(error) = result {
            self.error = Some(format!("{error:#}"));
        }
        self.refresh_sample();
        true
    }

    fn select(&mut self, id: EngineId) -> anyhow::Result<()> {
        self.engine = AnyEngine::new(id)?;
        self.sample.resize(self.engine.channels().len(), 0.0);
        self.error = None;
        Ok(())
    }

    /// Advances by the wall time since the last tick. Returns the new
    /// simulation time and sample when the engine moved.
    pub(crate) fn tick(&mut self) -> Option<(f64, &[f64])> {
        let plan = self.pacer.plan(self.clock.now());
        if plan.steps == 0 {
            return None;
        }
        let dt = Time::new::<second>(PHYSICS_STEP);
        let context = StepContext::real_time();
        for _ in 0..plan.steps {
            if let Err(error) = self.engine.advance(dt, &context) {
                // A failed step leaves the engine mid-state: stop until the
                // user resets.
                self.error = Some(format!("{error:#}"));
                self.pacer.pause();
                break;
            }
        }
        self.refresh_sample();
        Some((self.engine.time().get::<second>(), &self.sample))
    }

    fn refresh_sample(&mut self) {
        if let Err(error) = self.engine.sample(&mut self.sample).context("sampling") {
            self.error = Some(format!("{error:#}"));
        }
    }

    /// Copies the current state into `snapshot`.
    pub(crate) fn write_snapshot(&self, snapshot: &mut Snapshot) {
        snapshot.engine = self.engine.descriptor().id;
        snapshot.time = self.engine.time().get::<second>();
        snapshot.sample.clear();
        snapshot.sample.extend_from_slice(&self.sample);
        snapshot.controls.clear();
        snapshot
            .controls
            .extend_from_slice(self.engine.control_values());
        snapshot.paused = self.pacer.is_paused();
        snapshot.time_scale = self.pacer.time_scale();
        snapshot.achieved_time_scale = self.pacer.achieved_time_scale();
        snapshot.error.clone_from(&self.error);
    }
}

#[cfg(test)]
mod tests {
    use std::time::Duration;

    use engine_bench_core::{clock::MockClock, spec::find_control};

    use super::*;

    const DEMO: EngineId = EngineId::TwoTankDemo;

    fn simulator() -> anyhow::Result<Simulator<MockClock>> {
        let mut sim = Simulator::new(MockClock::new(), DEMO)?;
        // The first tick only starts the pacer's clock.
        sim.tick();
        Ok(sim)
    }

    fn snapshot(sim: &Simulator<MockClock>) -> Snapshot {
        let mut snapshot = Snapshot::empty(DEMO);
        sim.write_snapshot(&mut snapshot);
        snapshot
    }

    fn wait(sim: &mut Simulator<MockClock>, ms: u64) -> Option<f64> {
        sim.clock.advance(Duration::from_millis(ms));
        sim.tick().map(|(t, _)| t)
    }

    fn valve() -> anyhow::Result<ControlId> {
        find_control(DEMO.controls(), "valve").context("valve")
    }

    #[test]
    fn test_simulator_starts_at_rest_with_full_snapshot() -> anyhow::Result<()> {
        let sim = simulator()?;
        let s = snapshot(&sim);
        anyhow::ensure!(s.engine == DEMO && s.time == 0.0 && !s.paused, "{s:?}");
        anyhow::ensure!(s.sample.len() == DEMO.channels().len(), "sample length");
        anyhow::ensure!(s.controls.len() == DEMO.controls().len(), "controls length");
        anyhow::ensure!(
            s.error.is_none() && (s.time_scale - 1.0).abs() < 1e-15,
            "{s:?}"
        );
        Ok(())
    }

    #[test]
    fn test_simulator_follows_wall_time() -> anyhow::Result<()> {
        let mut sim = simulator()?;
        let t = wait(&mut sim, 100).context("no progress")?;
        anyhow::ensure!((t - 0.1).abs() < 1e-9, "t = {t}");
        anyhow::ensure!(wait(&mut sim, 0).is_none(), "moved without wall time");
        Ok(())
    }

    #[test]
    fn test_simulator_pause_step_resume_and_scale() -> anyhow::Result<()> {
        let mut sim = simulator()?;
        sim.handle(Command::Pause);
        anyhow::ensure!(
            wait(&mut sim, 100).is_none() && snapshot(&sim).paused,
            "moved while paused"
        );
        sim.handle(Command::Step);
        let stepped = wait(&mut sim, 10).context("single step")?;
        anyhow::ensure!((stepped - PHYSICS_STEP).abs() < 1e-12, "step {stepped}");
        sim.handle(Command::Resume);
        sim.handle(Command::SetTimeScale(0.5));
        let t = wait(&mut sim, 100).context("resumed")?;
        anyhow::ensure!((t - (PHYSICS_STEP + 0.05)).abs() < 1e-9, "t = {t}");
        anyhow::ensure!((snapshot(&sim).time_scale - 0.5).abs() < 1e-15, "scale");
        Ok(())
    }

    #[test]
    fn test_simulator_applies_controls_and_reports_rejections() -> anyhow::Result<()> {
        let mut sim = simulator()?;
        sim.handle(Command::SetControl(valve()?, 0.8));
        let s = snapshot(&sim);
        anyhow::ensure!(
            (s.controls[valve()?.index()] - 0.8).abs() < 1e-15,
            "{:?}",
            s.controls
        );
        anyhow::ensure!(s.error.is_none(), "{:?}", s.error);
        sim.handle(Command::SetControl(valve()?, 7.0));
        let s = snapshot(&sim);
        anyhow::ensure!(
            s.error.as_deref().is_some_and(|e| e.contains("valve")),
            "{:?}",
            s.error
        );
        anyhow::ensure!(
            (s.controls[valve()?.index()] - 0.8).abs() < 1e-15,
            "value changed"
        );
        sim.handle(Command::SetTimeScale(-1.0));
        anyhow::ensure!(
            snapshot(&sim).error.is_some(),
            "bad scale accepted silently"
        );
        Ok(())
    }

    #[test]
    fn test_simulator_reset_and_engine_selection_restart_time() -> anyhow::Result<()> {
        let mut sim = simulator()?;
        sim.handle(Command::SetControl(valve()?, 7.0));
        wait(&mut sim, 200);
        sim.handle(Command::Reset);
        let s = snapshot(&sim);
        anyhow::ensure!(s.time == 0.0 && s.error.is_none(), "{s:?}");
        wait(&mut sim, 50);
        sim.handle(Command::SelectEngine(DEMO));
        anyhow::ensure!(snapshot(&sim).time == 0.0, "selection did not restart");
        Ok(())
    }

    #[test]
    fn test_simulator_shutdown_stops() -> anyhow::Result<()> {
        let mut sim = simulator()?;
        anyhow::ensure!(sim.handle(Command::Pause), "pause stopped the loop");
        anyhow::ensure!(!sim.handle(Command::Shutdown), "shutdown ignored");
        Ok(())
    }

    #[test]
    fn test_simulator_caps_work_per_tick() -> anyhow::Result<()> {
        let mut sim = simulator()?;
        // Ten seconds of wall time at once must not run ten thousand steps.
        let t = wait(&mut sim, 10_000).context("no progress")?;
        let cap = f64::from(MAX_STEPS_PER_TICK) * PHYSICS_STEP;
        anyhow::ensure!((t - cap).abs() < 1e-9, "t = {t}, cap {cap}");
        Ok(())
    }
}
