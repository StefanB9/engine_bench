//! Compressible flow of a calorically perfect gas: isentropic relations, normal
//! shocks, orifices and converging–diverging nozzles.
//!
//! γ is frozen at the upstream (stagnation) state. That is exact for cold gas
//! and a small error for warm gas; a thermally perfect variant follows with
//! equilibrium chemistry (M5).

pub mod isentropic;
pub(crate) mod nozzle;
pub(crate) mod orifice;
mod perfect_gas;
pub mod real;
pub mod shock;

pub use nozzle::{Nozzle, NozzleOperatingPoint, NozzleRegime};
pub use orifice::{Orifice, Reservoir};
pub use perfect_gas::PerfectGas;
