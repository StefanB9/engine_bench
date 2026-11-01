# Plan 0001 — Architecture and Roadmap

## Goal

Physically simulate and visually present 38 engines (Aeolipile → HCCI/RCCI), with
credible numbers (thrust, torque, power, Isp, efficiency, p-V / T-s diagrams) and animated
2D cutaways. Two modes: real-time interactive exploration and offline batch runs.

## Decisions (from planning Q&A)

| Topic | Decision |
|---|---|
| Fidelity | Lumped 0D control volumes + 1D gas dynamics everywhere; **2D CFD where needed** (nozzles, plumes, inlets, detonation, scramjet) |
| Visuals | 2D animated cutaways + live plots, egui/eframe (WASM-capable later) |
| Modes | Real-time (live controls, slow-mo, pause) **and** offline batch (sweeps, high-res, replay) |
| CFD compute | CPU first (rayon), solver behind a trait so a wgpu backend can follow |
| Chemistry | NASA 7-coefficient polynomials + chemical equilibrium (CEA-style) |
| Units | `uom` at API boundaries (zero runtime cost); raw `f64` struct-of-arrays inside kernels |
| Crates | Split now: physics / core / binary |
| Priority | **Runtime performance over compile time, in every trade-off** (added to CLAUDE.md) |
| Milestone 1 | Framework + engines 1–3 |

## Crate layout

```
crates/
├── engine-bench-physics/   Units, thermo, chemistry, numerics, 0D/1D/2D solvers. No engine knowledge.
├── engine-bench-core/      Engine trait, engine models, scenarios, recording, renderer-agnostic Scene.
└── engine-bench/           Binary: eframe GUI, CLI (batch mode), logging, Scene → egui painter.
```

Dependency boundary: physics ← core ← engine-bench. Never reversed. No GUI types in
core: engines describe what to draw as a `Scene` (polygons, labeled parts, colormapped
scalar fields, particles, arrows), and the binary renders it.

## Core abstractions (sketch, final signatures decided during TDD)

- `Engine` trait: `controls() -> &[ControlSpec]`, `step(&mut self, dt, &Inputs) -> Result<()>`,
  `observables() -> Observables` (thrust, torque, rpm, efficiency, named channels),
  `scene(&self, &mut Scene)`, `diagrams()` (p-V, T-s, cycle traces).
  Statically dispatched inside the sim loop (enum of engines, not `Box<dyn>`).
- `GasModel` / `Mixture`: species set, NASA polynomials (cp, h, s, g), mixing, equilibrium
  solver (Gibbs minimisation, element-potential method) for products.
- Numerics: fixed-step RK4, adaptive RK45 (Dormand–Prince), later implicit for stiff
  problems. Root finding (Brent, Newton with bracketing).
- 0D building blocks: control volume (mass/energy), orifice and choked/unchoked nozzle
  flow, valves/ports, heat transfer (Woschni-style later), rigid-body rotation and
  translation, crank-slider kinematics (piston milestone).
- 1D solver: finite-volume Euler on ducts with area variation, HLLC flux, MUSCL
  reconstruction, boundary conditions coupling to 0D volumes.
- 2D solver: structured-grid finite-volume Euler (HLLC + MUSCL-Hancock), immersed/masked
  geometry, axisymmetric option for nozzles and plumes, rayon-parallel. Coarse grid in real
  time, fine grid offline. Behind a `FlowSolver2d` trait for a future GPU backend.
- Clock: simulation time is owned by the sim loop; a `Clock` trait abstracts wall time for
  the real-time pacer (mockable in tests).

## Runtime model

- **Real-time:** sim on its own thread with a fixed physics timestep (substeps per frame
  as needed for stability), time scale slider (slow-mo / pause / step). Snapshots go to
  the UI over a lock-free single-slot handoff (latest value wins), controls come back
  over a channel. Zero allocation per step after warm-up.
- **Offline:** headless run from the CLI (`engine-bench run <engine> --config ...`,
  `engine-bench sweep ...`), records time series to disk; the GUI can load and scrub a
  recording.

## Validation strategy

Each model gets golden-value tests against analytic solutions or published reference data:
Sod shock tube and Shu–Osher (1D), oblique shock / Prandtl–Meyer / quasi-1D nozzle (2D),
NASA CEA flame temperatures and Isp (chemistry), ideal-cycle limits (Otto, Diesel,
Brayton, Carnot/Stirling) as upper bounds. Proptest for conservation (mass, energy,
momentum) and thermodynamic consistency (h = ∫cp dT, Gibbs relation).

## Milestones

| # | Content | Physics added |
|---|---|---|
| M0 | Crate split, `uom`, CLAUDE.md update (structure, perf-over-compile-time rule) | — |
| M1 | Physics foundations: NASA polys for air/N2/O2/H2O/CO2/CO/H2/Ar, 0D control volume, isentropic and choked nozzle flow, RK4/RK45 | thermo, 0D |
| M2 | Sim loop (real-time + offline), `Engine` trait, `Scene`, eframe shell with engine picker, live plots, CLI batch mode | runtime |
| M3 | **Engines 1–3**: Aeolipile (boiler + steam tables*, tangential jets, rotor with friction/drag), cold gas thruster (tank blowdown, adiabatic/isothermal, thrust/Isp), water rocket (water expulsion → air blowdown → ballistic flight with drag) | water/steam, flight |
| M4 | 2D Euler solver + validation; first consumer: cold-gas-thruster plume view | 2D CFD |
| M5 | Equilibrium chemistry; engines 4–7 (Newcomen, Watt, ideal Stirling, solid rocket with grain burnback) | equilibrium, condensation |
| M6+ | Piston family 8–18 (crank kinematics, Wiebe combustion, 1D ducts, turbo maps), then air-breathing and turbomachinery 19–25, 32–34, then advanced rockets 21, 26, 29–31, 36, then detonation/electric/HCCI 27–28, 35, 37–38 | per milestone plan |

Each milestone after M0 gets its own detailed plan in `docs/plan/` before implementation.

**Visuals are finished per engine**: an engine ships with its final animated scene, since
everything a scene draws must be a recorded channel and a later visual pass would change
channel lists (schema bumps, old recordings unreadable). First a shared visual language in
the painter (flow colouring, droplet and exhaust particles, common parts), piloted on the
cold-gas thruster after its condensation work (plan 0005).

\* Steam: IAPWS-IF97 implemented in physics (no crate dependency) — decided in the M3 plan.

## Planned dependencies (each justified in its milestone plan)

- `uom` (physics, core): `default-features = false`, features `f64`, `si`, `std`.
- `rayon` (physics, M4): data-parallel 2D CFD.
- `eframe`/`egui`, `egui_plot` (binary, M2): GUI and plots.
- `clap` (binary, M2): CLI for batch mode.
- dev: `proptest`, `approx`.

## Known risks and open questions (deferred to later milestone plans)

1. **2D CFD inside piston cylinders** needs moving boundaries. Proposal: piston engines stay
   0D/1D, with 2D CFD reserved for open flows. Revisit at M6.
2. **Auto-ignition and detonation** (#28, #35, #38) cannot be captured by equilibrium alone.
   Options: calibrated ignition-delay correlations (Livengood–Wu) or reduced finite-rate
   kinetics. Decide at that milestone.
3. **Real-time 2D CFD budget**: grid size versus frame rate on CPU; GPU backend if insufficient.
4. **Ion/Hall thrusters** (#37) need plasma models (0D power balance + particle optics),
   a different physics family. Separate plan.
