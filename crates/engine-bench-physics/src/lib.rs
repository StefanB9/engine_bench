//! Engine Bench physics: units, thermodynamics, chemistry, numerics and flow
//! solvers.
//!
//! Knows nothing about specific engines; `engine-bench-core` builds engine
//! models on top.

pub mod error;
pub mod numerics;
pub mod units;

pub use error::PhysicsError;
