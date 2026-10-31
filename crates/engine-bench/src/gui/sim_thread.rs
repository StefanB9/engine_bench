//! Runs a [`Simulator`] on its own thread.
//!
//! Commands arrive over a channel; the latest [`Snapshot`] is handed over
//! through a triple buffer (the UI always reads the newest, nothing blocks);
//! plot rows go through a bounded channel of fixed-size rows, dropped and
//! counted when the UI falls behind.

use std::{
    sync::mpsc::{self, Receiver, Sender, SyncSender},
    thread::JoinHandle,
    time::Duration,
};

use anyhow::Context;
use engine_bench_core::{clock::SystemClock, engine::EngineId};

use super::simulator::{Command, Simulator, Snapshot};

/// Most channels a plot row carries; engines beyond this plot their first
/// channels only.
pub(crate) const MAX_PLOT_CHANNELS: usize = 64;
/// Longest wait between ticks when no command arrives.
const TICK: Duration = Duration::from_millis(2);
/// Plot rows buffered for the UI (about four seconds of ticks).
const ROW_QUEUE: usize = 2048;

/// One row of channel values for plotting; `Copy`, so sending it never
/// allocates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PlotRow {
    /// Simulation time, s.
    pub(crate) time: f64,
    len: usize,
    values: [f64; MAX_PLOT_CHANNELS],
}

impl PlotRow {
    /// A row holding the first [`MAX_PLOT_CHANNELS`] values of `sample`.
    pub(crate) fn new(time: f64, sample: &[f64]) -> Self {
        let len = sample.len().min(MAX_PLOT_CHANNELS);
        let mut values = [0.0; MAX_PLOT_CHANNELS];
        values[..len].copy_from_slice(&sample[..len]);
        Self { time, len, values }
    }

    /// The channel values.
    pub(crate) fn values(&self) -> &[f64] {
        &self.values[..self.len]
    }
}

/// The UI's side of a running simulation. Dropping it stops and joins the
/// thread.
pub(crate) struct SimulationHandle {
    commands: Sender<Command>,
    snapshots: triple_buffer::Output<Snapshot>,
    rows: Receiver<PlotRow>,
    thread: Option<JoinHandle<()>>,
}

impl SimulationHandle {
    /// Starts simulating `engine` in real time.
    pub(crate) fn spawn(engine: EngineId) -> anyhow::Result<Self> {
        let sim = Simulator::new(SystemClock::new(), engine)?;
        let mut initial = Snapshot::empty(engine);
        sim.write_snapshot(&mut initial);
        let (snapshot_in, snapshots) = triple_buffer::triple_buffer(&initial);
        let (commands, command_rx) = mpsc::channel();
        let (row_tx, rows) = mpsc::sync_channel(ROW_QUEUE);
        let thread = std::thread::Builder::new()
            .name("simulation".to_owned())
            .spawn(move || run(sim, &command_rx, &row_tx, snapshot_in))
            .context("starting the simulation thread")?;
        Ok(Self {
            commands,
            snapshots,
            rows,
            thread: Some(thread),
        })
    }

    /// Sends a command (ignored once the thread has stopped).
    pub(crate) fn send(&self, command: Command) {
        let _ = self.commands.send(command);
    }

    /// The newest snapshot.
    pub(crate) fn snapshot(&mut self) -> &Snapshot {
        self.snapshots.read()
    }

    /// Plot rows produced since the last call.
    pub(crate) fn drain_rows(&self) -> impl Iterator<Item = PlotRow> + '_ {
        self.rows.try_iter()
    }
}

/// The simulation thread: commands, tick, publish; until shutdown or the UI
/// goes away.
fn run(
    mut sim: Simulator<SystemClock>,
    commands: &Receiver<Command>,
    rows: &SyncSender<PlotRow>,
    mut snapshots: triple_buffer::Input<Snapshot>,
) {
    let mut dropped = 0;
    loop {
        // Waiting on the channel instead of sleeping applies commands
        // immediately.
        match commands.recv_timeout(TICK) {
            Ok(command) => {
                if !sim.handle(command) {
                    return;
                }
                for command in commands.try_iter() {
                    if !sim.handle(command) {
                        return;
                    }
                }
            }
            Err(mpsc::RecvTimeoutError::Timeout) => {}
            Err(mpsc::RecvTimeoutError::Disconnected) => return,
        }
        if let Some((time, sample)) = sim.tick()
            && rows.try_send(PlotRow::new(time, sample)).is_err()
        {
            dropped += 1;
        }
        let snapshot = snapshots.input_buffer_mut();
        sim.write_snapshot(snapshot);
        snapshot.dropped_rows = dropped;
        snapshots.publish();
    }
}

impl Drop for SimulationHandle {
    fn drop(&mut self) {
        let _ = self.commands.send(Command::Shutdown);
        if let Some(thread) = self.thread.take() {
            let _ = thread.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;

    /// Polls until `done` holds for the newest snapshot, for at most two
    /// seconds.
    fn wait_for(
        handle: &mut SimulationHandle,
        done: impl Fn(&Snapshot) -> bool,
    ) -> anyhow::Result<Snapshot> {
        let start = Instant::now();
        while start.elapsed() < Duration::from_secs(2) {
            let snapshot = handle.snapshot().clone();
            if done(&snapshot) {
                return Ok(snapshot);
            }
            std::thread::sleep(Duration::from_millis(5));
        }
        anyhow::bail!("timed out; last snapshot {:?}", handle.snapshot())
    }

    #[test]
    fn test_plot_row_keeps_values_up_to_the_limit() -> anyhow::Result<()> {
        let row = PlotRow::new(1.5, &[1.0, 2.0, 3.0]);
        anyhow::ensure!(
            row.time.to_bits() == 1.5f64.to_bits() && row.values() == [1.0, 2.0, 3.0],
            "{row:?}"
        );
        let wide = vec![7.0; MAX_PLOT_CHANNELS + 5];
        anyhow::ensure!(
            PlotRow::new(0.0, &wide).values().len() == MAX_PLOT_CHANNELS,
            "truncation"
        );
        Ok(())
    }

    #[test]
    fn test_simulation_thread_runs_streams_and_obeys_commands() -> anyhow::Result<()> {
        let mut handle = SimulationHandle::spawn(EngineId::TwoTankDemo)?;
        wait_for(&mut handle, |s| s.time > 0.02 && !s.sample.is_empty())?;
        anyhow::ensure!(handle.drain_rows().count() > 0, "no plot rows");

        handle.send(Command::Pause);
        let paused = wait_for(&mut handle, |s| s.paused)?;
        std::thread::sleep(Duration::from_millis(30));
        let later = wait_for(&mut handle, |_| true)?;
        anyhow::ensure!(
            later.time.to_bits() == paused.time.to_bits(),
            "advanced while paused"
        );

        handle.send(Command::Reset);
        wait_for(&mut handle, |s| s.time == 0.0)?;
        // Dropping the handle must stop and join the thread (the test would
        // hang otherwise).
        drop(handle);
        Ok(())
    }
}
