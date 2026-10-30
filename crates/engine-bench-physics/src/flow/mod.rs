//! Compressible flow of a calorically perfect gas: isentropic relations, normal
//! shocks, orifices and converging–diverging nozzles.
//!
//! γ is frozen at the upstream (stagnation) state. That is exact for cold gas
//! and a small error for warm gas; a thermally perfect variant follows with
//! equilibrium chemistry (M5).

pub mod isentropic;
mod orifice;
mod perfect_gas;
pub mod shock;

pub use orifice::{Orifice, Reservoir};
pub use perfect_gas::PerfectGas;
