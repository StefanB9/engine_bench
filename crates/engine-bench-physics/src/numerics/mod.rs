//! Generic numerical methods on raw `f64`.
//!
//! These work on dimensionless or mixed-unit state (an ODE state vector holds
//! masses, energies and angles side by side), so they take plain `f64`; callers
//! wrap results in `units` types.

pub mod ode;
pub mod roots;
