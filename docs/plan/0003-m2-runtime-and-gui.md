# Plan 0003 — M2 Runtime, GUI Shell and Batch Mode

Parent plan: `0001-architecture-and-roadmap.md`, milestone M2. Builds on M1 (`engine-bench-physics`).

## Decisions (from planning Q&A)

| Topic | Decision |
|---|---|
| Content | A non-catalogue **two-tank demo** engine proves the pipeline end to end |
| Time stepping | `Engine::advance(dt, &StepContext)`: each engine integrates itself with physics helpers; the loop only paces |
| Renderer | eframe with **wgpu** (shares the device with a future GPU CFD backend) |
| Recordings | **Polars**: Parquet files, CSV export, DataFrame API for sweeps |
| Platforms | **Windows and Linux** (X11 and Wayland). Web later |
| Display units | **Per-quantity unit picker**: SI default, presets (SI / engineering / imperial), per-quantity override; recordings stored in SI with the quantity kind in metadata |
| Replay | **Animated via scene-state channels**: everything a scene needs is a recorded channel, and one `draw(sample)` function serves live view and replay |

## Crate placement

- **core** (`engine-bench-core`, now depends on physics): `Engine` trait, control and channel
  descriptors, `Scene`, in-memory `Recording`, `Simulation` runner and real-time pacer with a
  `Clock` trait, `AnyEngine` enum, the two-tank demo. No GUI, no file formats.
- **binary** (`engine-bench`): CLI (clap), recording I/O (polars), display-unit system, eframe
  app (simulation thread, engine picker, controls, scene painter, plots, replay).

New dependencies, all with `default-features = false`:

| Crate | Crate placement | Features | Why |
|---|---|---|---|
| eframe 0.37 | binary | `wgpu`, `default_fonts`, `x11`, `wayland` | Native window + immediate-mode UI; decided above |
| egui_plot 0.38 | binary | none | Live plots; same family as egui, no alternative needed |
| clap 4.6 | binary | `std`, `derive`, `help`, `usage`, `error-context` | Standard CLI parser; derive keeps arguments typed |
| polars 0.55 | binary | `parquet`, `csv` | Decided above |
| triple_buffer 9.0 | binary | none | Lock-free latest-value handoff of the scene to the UI; a `Mutex<Snapshot>` would make the sim thread wait on rendering |
| egui_kittest 0.37 (dev) | binary | none | Headless UI tests, so the GUI can be test-driven |

## Core API (sketch; final signatures settle during TDD)

```rust
pub trait Engine {
    fn descriptor(&self) -> &EngineDescriptor;          // id, name, catalogue no., summary
    fn controls(&self) -> &[ControlSpec];               // slider/toggle/button, range, quantity kind
    fn channels(&self) -> &[ChannelSpec];               // observable name + quantity kind
    fn set_control(&mut self, id: ControlId, value: f64) -> anyhow::Result<()>;
    fn advance(&mut self, dt: Time, context: &StepContext) -> anyhow::Result<()>;
    fn time(&self) -> Time;
    fn sample(&self, out: &mut [f64]);                  // channel values, SI, no allocation
    fn reset(&mut self);
}

impl EngineId {
    /// Builds the scene from one sample of the engine's channels (refills a reused Scene).
    fn draw(self, sample: &[f64], scene: &mut Scene);
}
```

- **Scene from channels.** The scene is a pure function of the channel sample: the live view
  calls `sample()` then `draw()`, and replay calls `draw()` on recorded rows. Live and replay
  run the same code on the same numbers, so they cannot drift apart. Every value a scene
  needs, constant geometry included, is a channel. Constant columns cost almost nothing in
  Parquet (run-length encoded).
- `ChannelSpec` has a role: `Observable` (offered in plots) or `SceneState` (recorded for
  drawing, hidden from plot pickers by default).

- `StepContext { accuracy: RealTime | Offline { rtol } }`: engines pick RK4 sub-steps or DP45.
- `AnyEngine` is an enum over all engines, implementing `Engine` by `match` (static
  dispatch, per plan 0001). An `EngineId` enum plus `AnyEngine::new(EngineId)` act as the registry.
- `QuantityKind` (Pressure, Temperature, MassFlow, Force, AngularVelocity, Power, …) tags
  controls and channels, so units and plots need no per-engine code.
- **Scene**: a renderer-agnostic list reused between frames. World coordinates in metres,
  y up; bounds for fitting the view. Primitives: polygon (fill + stroke), polyline, circle,
  arrow (vector quantity), text label. Colours are semantic (metal, liquid, gas coloured by
  T or p through a colormap). Scalar-field images arrive with M4.
- **Recording**: SI columns (time + one per channel) with names, quantity kinds and roles,
  plus the engine id and a channel-schema version so replay can pick the right `draw`;
  pre-allocated by duration × sample rate.

## Runtime

- **Pacing** (core, `Clock` trait + mock): a fixed physics step `dt`, a time-scale factor
  (slow motion, pause, single step) and a per-frame catch-up cap. It reports the achieved
  real-time factor when the machine can't keep up. Pure logic, tested with a mock clock.
- **Real-time** (binary): the simulation runs on its own `std::thread`.
  - Scene and latest values go to the UI through `triple_buffer` (latest wins, no allocation
    after warm-up).
  - Channel samples go through a bounded `sync_channel` at the plot sample rate; when it's
    full, samples are dropped and counted rather than blocking the simulation.
  - UI commands (select engine, set control, scale, pause, step, reset) go through an
    `mpsc` channel.
- **Offline** (core runner + binary CLI):
  - `engine-bench run <engine> --duration 10s --rate 1kHz --set valve=0.5 -o run.parquet`
  - `engine-bench sweep <engine> --vary valve=0.1:1.0:10 ...`: runs in parallel on scoped
    threads, one Parquet with a `run` column.
  - `engine-bench export run.parquet --csv --units engineering`.
  - Running `engine-bench` with no subcommand opens the GUI.

## GUI (binary)

Engine picker → control panel generated from `ControlSpec`s → central scene view (painter
fits scene bounds, pan/zoom) → plot panel (channels selectable, history window) → status
bar (sim time, real-time factor, dropped samples). Unit picker: preset plus per-quantity
override. Replay: open a Parquet recording, then scrub or play it back; the scene animates
through `EngineId::draw` on the recorded rows, and the plots show a time cursor. A recording
whose engine id or schema version is unknown opens plots-only, with a warning.

## Two-tank demo

Tank A (high pressure) and tank B joined by a pipe with a valve. Controls: valve opening
(scales orifice area), initial pressures, heater power on tank B, reset. Observable channels:
p, T and m of each tank, mass flow. Scene-state channels: tank volumes, valve opening. Scene: tanks coloured by temperature, valve symbol, flow arrow
scaled by ṁ. Uses M1 `ControlVolume`, `orifice_transfer`, RK4 / DP45.

## Testing

- Core: unit tests for every public item; a proptest that pacing never exceeds the catch-up
  cap and that its time accounting is exact; tests that the demo conserves mass and energy
  through `advance` in both accuracy modes; scene bounds/contents tests; `draw` on a recorded
  row equals the live scene at the same time.
- Binary:
  - unit conversions (golden values, and proptest round trips SI → display → SI);
  - CLI parsing;
  - integration test: `run` → Parquet → read back equals the in-memory recording;
  - CSV export in chosen units;
  - `egui_kittest` (no GPU snapshots in this container): app starts, picker lists engines,
    a slider sends a control command, pause/step works.
- Visual appearance: checked by you on Windows; I can only test headless.

## Delivery: 4 sequential PRs

1. `feature/m2-engine-core`: Engine trait, specs, `QuantityKind`, `Scene`, `Recording`,
   `AnyEngine`, two-tank demo.
2. `feature/m2-simulation`: `Clock`, pacer, offline runner.
3. `feature/m2-cli`: clap commands, polars Parquet/CSV I/O, display-unit system.
4. `feature/m2-gui`: eframe app, sim thread, picker, controls, scene painter, plots,
   unit picker, replay.

Each PR updates CLAUDE.md for the modules it adds.
