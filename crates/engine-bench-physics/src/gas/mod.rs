//! Lumped (0D) gas volumes and the flows between them.
//!
//! A [`ControlVolume`] stores conserved quantities (species masses, internal
//! energy) and derives p and T from them. Anything that can feed an orifice
//! implements [`GasPort`]; [`orifice_transfer`] returns what leaves one port
//! and enters the other, so mass and energy are conserved by construction.

mod control_volume;
pub mod eos;
mod port;

pub use control_volume::{ControlVolume, NEGATIVE_MASS_TOLERANCE};
pub use port::{Environment, GasPort, Transfer, orifice_transfer};
