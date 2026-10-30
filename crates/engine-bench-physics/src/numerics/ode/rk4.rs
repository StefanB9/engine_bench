//! Classic fixed-step fourth-order Runge–Kutta.

use anyhow::Context;

use super::{OdeSystem, all_finite, combine};
use crate::PhysicsError;

/// Advances `y` from `t` to `t + h` with one classic fourth-order Runge–Kutta
/// step.
///
/// On a non-finite result `t` and `y` are left unchanged and
/// [`PhysicsError::NonFiniteState`] is returned.
pub fn rk4_step<S: OdeSystem<N>, const N: usize>(
    system: &mut S,
    t: &mut f64,
    y: &mut [f64; N],
    h: f64,
) -> anyhow::Result<()> {
    let t0 = *t;
    let mut k1 = [0.0; N];
    let mut k2 = [0.0; N];
    let mut k3 = [0.0; N];
    let mut k4 = [0.0; N];
    let mut stage = [0.0; N];

    system.derivatives(t0, y, &mut k1);
    combine(&mut stage, y, &[(0.5 * h, &k1)]);
    system.derivatives(t0 + 0.5 * h, &stage, &mut k2);
    combine(&mut stage, y, &[(0.5 * h, &k2)]);
    system.derivatives(t0 + 0.5 * h, &stage, &mut k3);
    combine(&mut stage, y, &[(h, &k3)]);
    system.derivatives(t0 + h, &stage, &mut k4);
    let sixth = h / 6.0;
    combine(
        &mut stage,
        y,
        &[
            (sixth, &k1),
            (2.0 * sixth, &k2),
            (2.0 * sixth, &k3),
            (sixth, &k4),
        ],
    );

    if !all_finite(&stage) {
        return Err(PhysicsError::NonFiniteState { t: t0 }).context("rk4 step");
    }
    *y = stage;
    *t = t0 + h;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::f64::consts::TAU;

    use proptest::prelude::*;

    use super::*;
    use crate::numerics::ode::test_systems::{AlwaysNan, Decay, Oscillator, physics_error};

    fn rk4_decay_error(steps: u32) -> anyhow::Result<f64> {
        let (mut t, mut y) = (0.0, [1.0]);
        let h = 1.0 / f64::from(steps);
        for _ in 0..steps {
            rk4_step(&mut Decay { rate: 1.0 }, &mut t, &mut y, h)?;
        }
        Ok((y[0] - (-t).exp()).abs())
    }

    #[test]
    fn test_rk4_step_is_fourth_order() -> anyhow::Result<()> {
        let ratio = rk4_decay_error(10)? / rk4_decay_error(20)?;
        anyhow::ensure!(
            (14.0..18.0).contains(&ratio),
            "error ratio {ratio}, expected ≈ 16"
        );
        Ok(())
    }

    #[test]
    fn test_rk4_step_advances_oscillator_one_period() -> anyhow::Result<()> {
        let (mut t, mut y) = (0.0, [1.0, 0.0]);
        let h = TAU / 1000.0;
        for _ in 0..1000 {
            rk4_step(&mut Oscillator, &mut t, &mut y, h)?;
        }
        anyhow::ensure!((t - TAU).abs() < 1e-12, "t = {t}");
        anyhow::ensure!((y[0] - 1.0).abs() < 1e-9 && y[1].abs() < 1e-9, "y = {y:?}");
        Ok(())
    }

    #[test]
    fn test_rk4_step_rejects_non_finite_and_keeps_state() -> anyhow::Result<()> {
        let (mut t, mut y) = (0.5, [2.0]);
        let Err(e) = rk4_step(&mut AlwaysNan, &mut t, &mut y, 0.1) else {
            anyhow::bail!("accepted NaN derivative");
        };
        anyhow::ensure!(
            matches!(physics_error(&e), Some(PhysicsError::NonFiniteState { .. })),
            "wrong error: {e:#}"
        );
        anyhow::ensure!(
            t.to_bits() == 0.5f64.to_bits() && y[0].to_bits() == 2f64.to_bits(),
            "state changed"
        );
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_rk4_step_conserves_oscillator_energy_to_order(
            x0 in -10.0..10.0f64, v0 in -10.0..10.0f64,
        ) {
            let (mut t, mut y) = (0.0, [x0, v0]);
            for _ in 0..100 {
                rk4_step(&mut Oscillator, &mut t, &mut y, 0.01)
                    .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            }
            let e0 = x0 * x0 + v0 * v0;
            let e1 = y[0] * y[0] + y[1] * y[1];
            prop_assert!((e1 - e0).abs() <= 1e-10 * e0.max(1.0), "energy {e0} → {e1}");
        }
    }
}
