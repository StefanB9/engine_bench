//! Position of a normal shock inside the divergent section.

use anyhow::Context;

use crate::{
    flow::{isentropic, shock},
    numerics::roots::{RootOptions, brent},
};

/// Exit Mach number and shock area ratio for a normal shock inside the
/// divergent section that brings the exit pressure to `back_ratio`·p₀.
///
/// Behind the shock, pₑ·Aₑ = p₀₂·A₂*·(p/p₀)(A/A*)(Mₑ) and p₀·A* = p₀₂·A₂* (same
/// choked mass flow), so (p/p₀)(A/A*)(Mₑ) = `back_ratio`·ε. That product equals
/// C/(M·√(1 + (γ−1)/2·M²)) with C = (2/(γ+1))^((γ+1)/(2(γ−1))): a quadratic in
/// M². The stagnation-pressure loss then fixes the shock Mach number.
pub(super) fn shock_in_nozzle(
    gamma: f64,
    epsilon: f64,
    back_ratio: f64,
    supersonic_exit: f64,
) -> anyhow::Result<(f64, f64)> {
    let a = 0.5 * (gamma - 1.0);
    let c = (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)));
    let k = back_ratio * epsilon;
    let mach_squared = ((1.0 + 4.0 * a * c * c / (k * k)).sqrt() - 1.0) / (2.0 * a);
    let exit_mach = mach_squared.sqrt();

    let loss = back_ratio / isentropic::pressure_ratio(gamma, exit_mach);
    let options = RootOptions::new(1e-13 * supersonic_exit, 200).context("root options")?;
    let shock_mach = brent(
        |m| shock::normal_shock(gamma, m).stagnation_pressure_ratio - loss,
        1.0,
        supersonic_exit,
        options,
    )
    .with_context(|| format!("shock Mach for stagnation pressure ratio {loss}"))?;
    Ok((exit_mach, isentropic::area_ratio(gamma, shock_mach)))
}
