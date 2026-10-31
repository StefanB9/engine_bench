//! The interactive application: live simulation and animated replay.
//!
//! `simulator` is the pure engine-plus-pacer state machine, `sim_thread` runs
//! it on its own thread, `history` keeps plot data, `painter` turns a
//! [`Scene`](engine_bench_core::scene) into egui shapes, `replay` steps through
//! recordings and `ui` draws the panels.

mod app;
mod controls;
mod history;
mod painter;
mod plots;
mod replay;
mod sim_thread;
mod simulator;

pub(crate) use self::app::run;
