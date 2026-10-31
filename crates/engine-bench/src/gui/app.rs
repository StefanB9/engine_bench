//! The application window: live simulation or replay, scene, plots and status.

use std::path::Path;

use eframe::egui::{self, Sense};
use engine_bench_core::{engine::EngineId, scene::Scene};

use super::{
    controls,
    history::PlotHistory,
    painter::{self, ViewTransform},
    plots::{self, PlotSelection},
    replay::Replay,
    sim_thread::SimulationHandle,
    simulator::{Command, Snapshot},
};
use crate::{
    recording_io::read_parquet,
    units::{Preset, UnitSelection},
};

/// Plot rows kept for live plots (about 20 s at the tick rate).
const HISTORY_ROWS: usize = 10_000;

/// Opens the window and runs until it is closed.
pub(crate) fn run() -> anyhow::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Engine Bench")
            .with_inner_size([1280.0, 800.0]),
        ..Default::default()
    };
    eframe::run_native(
        "Engine Bench",
        options,
        Box::new(|_| Ok(Box::new(App::new(EngineId::TwoTankDemo)?))),
    )
    .map_err(|e| anyhow::anyhow!("GUI failed: {e}"))
}

enum Mode {
    Live {
        sim: SimulationHandle,
        snapshot: Snapshot,
    },
    Replay {
        replay: Option<Replay>,
        path: String,
    },
}

/// Live simulation or replay, with everything the window shows.
pub(crate) struct App {
    mode: Mode,
    engine: EngineId,
    units: UnitSelection,
    scene: Scene,
    history: PlotHistory,
    plotted: PlotSelection,
    commands: Vec<Command>,
    error: Option<String>,
}

impl App {
    /// Starts a live simulation of `engine`.
    pub(crate) fn new(engine: EngineId) -> anyhow::Result<Self> {
        Ok(Self {
            mode: live(engine)?,
            engine,
            units: UnitSelection::new(Preset::Engineering),
            scene: Scene::new(),
            history: PlotHistory::new(engine.channels().len(), HISTORY_ROWS),
            plotted: PlotSelection::new(engine.channels()),
            commands: Vec::new(),
            error: None,
        })
    }

    /// Draws one frame into `ui`.
    pub(crate) fn show(&mut self, ui: &mut egui::Ui) {
        self.top_bar(ui);
        match &mut self.mode {
            Mode::Live { sim, snapshot } => {
                snapshot.clone_from(sim.snapshot());
                for row in sim.drain_rows() {
                    self.history.push(row.time, row.values());
                }
            }
            Mode::Replay { replay, .. } => {
                if let Some(replay) = replay {
                    replay.advance(f64::from(ui.input(|i| i.stable_dt)));
                }
            }
        }
        egui::Panel::left("controls").show(ui, |ui| self.side_panel(ui));
        egui::Panel::bottom("plots")
            .resizable(true)
            .default_size(220.0)
            .show(ui, |ui| self.plot_panel(ui));
        egui::CentralPanel::default().show(ui, |ui| self.scene_view(ui));

        if let Mode::Live { sim, .. } = &self.mode {
            for command in self.commands.drain(..) {
                sim.send(command);
            }
        }
        ui.ctx().request_repaint();
    }

    fn top_bar(&mut self, ui: &mut egui::Ui) {
        egui::Panel::top("mode").show(ui, |ui| {
            ui.horizontal(|ui| {
                let is_live = matches!(self.mode, Mode::Live { .. });
                if ui.selectable_label(is_live, "Live").clicked() && !is_live {
                    self.switch_to_live();
                }
                if ui.selectable_label(!is_live, "Replay").clicked() && is_live {
                    self.mode = Mode::Replay {
                        replay: None,
                        path: String::new(),
                    };
                }
                ui.separator();
                if is_live {
                    egui::ComboBox::from_id_salt("engine")
                        .selected_text(self.engine.descriptor().name)
                        .show_ui(ui, |ui| {
                            for id in EngineId::ALL {
                                if ui
                                    .selectable_label(id == self.engine, id.descriptor().name)
                                    .clicked()
                                {
                                    self.select_engine(id);
                                }
                            }
                        });
                }
                ui.menu_button("Units", |ui| controls::unit_picker(ui, &mut self.units));
            });
        });
    }

    fn switch_to_live(&mut self) {
        match live(self.engine) {
            Ok(mode) => {
                self.mode = mode;
                self.reset_plots(self.engine);
                self.error = None;
            }
            Err(error) => self.error = Some(format!("error: {error:#}")),
        }
    }

    fn select_engine(&mut self, id: EngineId) {
        self.engine = id;
        self.commands.push(Command::SelectEngine(id));
        self.reset_plots(id);
    }

    fn reset_plots(&mut self, id: EngineId) {
        self.history = PlotHistory::new(id.channels().len(), HISTORY_ROWS);
        self.plotted = PlotSelection::new(id.channels());
    }

    fn side_panel(&mut self, ui: &mut egui::Ui) {
        match &mut self.mode {
            Mode::Live { snapshot, .. } => {
                let engine = snapshot.engine;
                ui.heading(engine.descriptor().name);
                ui.label(engine.descriptor().summary);
                ui.separator();
                controls::transport(ui, snapshot, &mut self.commands);
                ui.separator();
                controls::engine_controls(
                    ui,
                    engine.controls(),
                    &snapshot.controls,
                    &self.units,
                    &mut self.commands,
                );
                if let Some(error) = &snapshot.error {
                    ui.colored_label(ui.visuals().error_fg_color, format!("error: {error}"));
                }
            }
            Mode::Replay { replay, path } => {
                if let Some(message) = replay_panel(ui, replay, path, &mut self.history) {
                    self.error = Some(message);
                }
            }
        }
        if let Some(error) = &self.error {
            ui.colored_label(ui.visuals().error_fg_color, error);
        }
    }

    fn plot_panel(&mut self, ui: &mut egui::Ui) {
        let channels = self.current_engine().channels();
        self.status_line(ui);
        self.plotted.picker(ui, channels);
        match &self.mode {
            Mode::Live { .. } => {
                plots::plots(
                    ui,
                    channels,
                    &self.plotted,
                    &self.units,
                    |i| self.history.points(i),
                    None,
                );
            }
            Mode::Replay {
                replay: Some(replay),
                ..
            } => {
                let recording = &replay.file().runs[replay.run()];
                let points = |i: usize| {
                    recording
                        .time()
                        .iter()
                        .zip(recording.column(i).unwrap_or_default())
                        .map(|(t, v)| [*t, *v])
                        .collect()
                };
                plots::plots(
                    ui,
                    channels,
                    &self.plotted,
                    &self.units,
                    points,
                    Some(replay.cursor()),
                );
            }
            Mode::Replay { replay: None, .. } => {}
        }
    }

    fn status_line(&self, ui: &mut egui::Ui) {
        let text = match &self.mode {
            Mode::Live { snapshot, .. } => {
                let achieved = snapshot
                    .achieved_time_scale
                    .map_or_else(|| "–".to_owned(), |f| format!("{f:.2}"));
                format!(
                    "t = {:.3} s · achieved ×{achieved} of ×{:.3} · dropped rows {}",
                    snapshot.time, snapshot.time_scale, snapshot.dropped_rows
                )
            }
            Mode::Replay {
                replay: Some(replay),
                ..
            } => {
                format!("t = {:.3} s of {:.3} s", replay.cursor(), replay.duration())
            }
            Mode::Replay { replay: None, .. } => "no recording loaded".to_owned(),
        };
        ui.label(text);
    }

    fn current_engine(&self) -> EngineId {
        match &self.mode {
            Mode::Live { snapshot, .. } => snapshot.engine,
            Mode::Replay {
                replay: Some(replay),
                ..
            } => replay.engine(),
            Mode::Replay { replay: None, .. } => self.engine,
        }
    }

    fn scene_view(&mut self, ui: &mut egui::Ui) {
        let drawn = match &mut self.mode {
            Mode::Live { snapshot, .. } => snapshot.engine.draw(&snapshot.sample, &mut self.scene),
            Mode::Replay {
                replay: Some(replay),
                ..
            } => {
                let engine = replay.engine();
                engine.draw(replay.sample(), &mut self.scene)
            }
            Mode::Replay { replay: None, .. } => {
                ui.centered_and_justified(|ui| ui.label("Open or drop a Parquet recording."));
                return;
            }
        };
        if let Err(error) = drawn {
            ui.colored_label(ui.visuals().error_fg_color, format!("error: {error:#}"));
            return;
        }
        let (response, painter) = ui.allocate_painter(ui.available_size(), Sense::hover());
        let Some(bounds) = self.scene.view_or_bounds() else {
            return;
        };
        let view = ViewTransform::fit(bounds, response.rect, 40.0);
        painter.extend(painter::to_shapes(
            &self.scene,
            &view,
            &self.units,
            &painter,
        ));
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.show(ui);
    }
}

fn live(engine: EngineId) -> anyhow::Result<Mode> {
    Ok(Mode::Live {
        sim: SimulationHandle::spawn(engine)?,
        snapshot: Snapshot::empty(engine),
    })
}

/// Replay controls; returns an error message for the user, if any.
fn replay_panel(
    ui: &mut egui::Ui,
    replay: &mut Option<Replay>,
    path: &mut String,
    history: &mut PlotHistory,
) -> Option<String> {
    let mut message = None;
    ui.heading("Replay");
    ui.horizontal(|ui| {
        ui.text_edit_singleline(path);
        let dropped = ui.input(|i| i.raw.dropped_files.first().map(|f| f.path().to_path_buf()));
        if let Some(dropped) = &dropped {
            *path = dropped.display().to_string();
        }
        if ui.button("Open").clicked() || dropped.is_some() {
            match load(Path::new(path.trim())) {
                Ok(loaded) => {
                    *replay = Some(loaded);
                    history.clear();
                }
                Err(error) => message = Some(format!("error: {error:#}")),
            }
        }
    });
    let Some(replay) = replay else {
        return message;
    };
    ui.label(replay.engine().descriptor().name);
    // Cloned: selecting a run below needs `replay` mutably.
    if let Some(sweep) = replay.file().sweep.clone() {
        let current = replay.run();
        egui::ComboBox::from_label(format!("run ({})", sweep.control))
            .selected_text(format!("{current}: {}", sweep.values[current]))
            .show_ui(ui, |ui| {
                for (index, value) in sweep.values.iter().enumerate() {
                    if ui
                        .selectable_label(index == current, format!("{index}: {value}"))
                        .clicked()
                    {
                        let _ = replay.select_run(index);
                    }
                }
            });
    }
    ui.horizontal(|ui| {
        let label = if replay.is_playing() { "Pause" } else { "Play" };
        if ui.button(label).clicked() {
            replay.set_playing(!replay.is_playing());
        }
        let mut speed = replay.speed();
        let slider = egui::Slider::new(&mut speed, 0.01..=100.0)
            .logarithmic(true)
            .text("speed")
            .suffix("×");
        if ui.add(slider).changed() {
            replay.set_speed(speed);
        }
    });
    let mut cursor = replay.cursor();
    let scrub = egui::Slider::new(&mut cursor, 0.0..=replay.duration())
        .text("time")
        .suffix(" s");
    if ui.add(scrub).changed() {
        replay.seek(cursor);
    }
    message
}

fn load(path: &Path) -> anyhow::Result<Replay> {
    anyhow::ensure!(
        !path.as_os_str().is_empty(),
        "enter the path of a Parquet recording"
    );
    Replay::new(read_parquet(path)?)
}

#[cfg(test)]
mod tests {
    use std::time::{Duration, Instant};

    use egui_kittest::{Harness, kittest::Queryable};

    use super::*;

    /// Runs frames until `done` holds, for at most two seconds.
    fn frames_until(
        h: &mut Harness<'_, App>,
        done: impl Fn(&Harness<'_, App>) -> bool,
    ) -> anyhow::Result<()> {
        let start = Instant::now();
        while start.elapsed() < Duration::from_secs(2) {
            h.step();
            if done(h) {
                return Ok(());
            }
            std::thread::sleep(Duration::from_millis(10));
        }
        anyhow::bail!("condition not reached in time")
    }

    #[test]
    fn test_app_shows_live_demo_and_pauses_it() -> anyhow::Result<()> {
        let app = App::new(EngineId::TwoTankDemo)?;
        let mut h = Harness::new_ui_state(|ui, app: &mut App| app.show(ui), app);
        h.step();
        anyhow::ensure!(
            h.query_all_by_label_contains("Two-tank demo")
                .next()
                .is_some(),
            "engine name missing"
        );
        anyhow::ensure!(
            h.query_all_by_label_contains("Valve opening")
                .next()
                .is_some(),
            "controls missing"
        );
        frames_until(&mut h, |h| {
            h.query_all_by_label_contains("t = 0.0").next().is_none()
        })?;
        h.get_by_label("Pause").click();
        frames_until(&mut h, |h| h.query_all_by_label("Resume").next().is_some())?;
        Ok(())
    }

    #[test]
    fn test_app_switches_to_replay_and_reports_bad_paths() -> anyhow::Result<()> {
        let app = App::new(EngineId::TwoTankDemo)?;
        let mut h = Harness::new_ui_state(|ui, app: &mut App| app.show(ui), app);
        h.step();
        h.get_by_label("Replay").click();
        h.step();
        h.get_by_label("Open").click();
        h.step();
        anyhow::ensure!(
            h.query_all_by_label_contains("error").next().is_some(),
            "no error for an empty path"
        );
        Ok(())
    }
}
