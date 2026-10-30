//! Scalar root finding: Brent's method and bracketed (safeguarded) Newton.

use anyhow::Context;

use crate::PhysicsError;

/// Convergence settings for the root finders.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RootOptions {
    x_tol: f64,
    max_iterations: u32,
}

impl RootOptions {
    /// Converge when the bracket (or Newton step) is smaller than `x_tol`,
    /// giving up after `max_iterations`.
    pub fn new(x_tol: f64, max_iterations: u32) -> anyhow::Result<Self> {
        if !(x_tol > 0.0 && x_tol.is_finite()) {
            return Err(PhysicsError::InvalidSetting(
                "x_tol must be positive and finite",
            ))
            .with_context(|| format!("x_tol = {x_tol}"));
        }
        if max_iterations == 0 {
            return Err(PhysicsError::InvalidSetting("max_iterations must be at least 1").into());
        }
        Ok(Self {
            x_tol,
            max_iterations,
        })
    }

    /// Absolute tolerance on the root location.
    pub fn x_tol(&self) -> f64 {
        self.x_tol
    }

    /// Iteration limit.
    pub fn max_iterations(&self) -> u32 {
        self.max_iterations
    }
}

impl Default for RootOptions {
    fn default() -> Self {
        Self {
            x_tol: 1e-12,
            max_iterations: 100,
        }
    }
}

/// Finds a root of `f` in `[a, b]` with Brent's method. `f(a)` and `f(b)` must
/// differ in sign.
#[expect(
    clippy::many_single_char_names,
    reason = "names follow Brent (1973) and Numerical Recipes so the code can be checked against \
              them"
)]
pub fn brent(
    mut f: impl FnMut(f64) -> f64,
    a: f64,
    b: f64,
    options: RootOptions,
) -> anyhow::Result<f64> {
    // Brent (1973) as in Numerical Recipes `zbrent`: inverse quadratic
    // interpolation or secant, falling back to bisection whenever the
    // interpolated step would not shrink the bracket.
    let (mut a, mut b) = (a, b);
    let (mut fa, mut fb) = (f(a), f(b));
    if fa == 0.0 {
        return Ok(a);
    }
    if fb == 0.0 {
        return Ok(b);
    }
    ensure_bracket(a, fa, b, fb)?;

    let (mut c, mut fc) = (a, fa);
    let mut d = b - a;
    let mut e = d;
    for _ in 0..options.max_iterations {
        if fb.signum() == fc.signum() {
            (c, fc) = (a, fa);
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            (a, fa) = (b, fb);
            (b, fb) = (c, fc);
            (c, fc) = (a, fa);
        }
        let tol = 2.0 * f64::EPSILON * b.abs() + 0.5 * options.x_tol;
        let half = 0.5 * (c - b);
        if half.abs() <= tol || fb == 0.0 {
            return Ok(b);
        }
        if e.abs() >= tol && fa.abs() > fb.abs() {
            let s = fb / fa;
            #[expect(
                clippy::float_cmp,
                reason = "c is assigned from a: exact identity selects secant over inverse \
                          quadratic"
            )]
            let only_two_points = a == c;
            let (mut p, mut q) = if only_two_points {
                (2.0 * half * s, 1.0 - s)
            } else {
                let qa = fa / fc;
                let r = fb / fc;
                (
                    s * (2.0 * half * qa * (qa - r) - (b - a) * (r - 1.0)),
                    (qa - 1.0) * (r - 1.0) * (s - 1.0),
                )
            };
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();
            let interpolation_limit = (3.0 * half * q - (tol * q).abs()).min((e * q).abs());
            if 2.0 * p < interpolation_limit {
                e = d;
                d = p / q;
            } else {
                d = half;
                e = d;
            }
        } else {
            d = half;
            e = d;
        }
        (a, fa) = (b, fb);
        b += if d.abs() > tol { d } else { tol.copysign(half) };
        fb = f(b);
    }
    Err(PhysicsError::NoConvergence {
        iterations: options.max_iterations,
        residual: fb.abs(),
    })
    .context("brent")
}

/// Finds a root of `f` in `[lo, hi]` starting from `x0` with Newton's method,
/// falling back to bisection whenever a Newton step leaves the bracket or fails
/// to shrink it fast enough. `f_df` returns the function value and its
/// derivative.
pub fn newton_bracketed(
    mut f_df: impl FnMut(f64) -> (f64, f64),
    lo: f64,
    hi: f64,
    x0: f64,
    options: RootOptions,
) -> anyhow::Result<f64> {
    // Numerical Recipes `rtsafe`: keep a bracket [x_neg, x_pos] with f(x_neg) <
    // 0 < f(x_pos) and bisect when Newton would leave it or converge slower
    // than bisection.
    let (f_lo, _) = f_df(lo);
    let (f_hi, _) = f_df(hi);
    if f_lo == 0.0 {
        return Ok(lo);
    }
    if f_hi == 0.0 {
        return Ok(hi);
    }
    ensure_bracket(lo, f_lo, hi, f_hi)?;

    let (mut x_neg, mut x_pos) = if f_lo < 0.0 { (lo, hi) } else { (hi, lo) };
    let mut x = x0.clamp(lo.min(hi), lo.max(hi));
    let mut dx_old = (hi - lo).abs();
    let mut dx = dx_old;
    let (mut fx, mut dfx) = f_df(x);
    for _ in 0..options.max_iterations {
        let newton_dx = fx / dfx;
        let leaves_bracket = ((x - x_pos) * dfx - fx) * ((x - x_neg) * dfx - fx) > 0.0;
        let too_slow = (2.0 * fx).abs() > (dx_old * dfx).abs();
        dx_old = dx;
        if leaves_bracket || too_slow || !newton_dx.is_finite() {
            dx = 0.5 * (x_pos - x_neg);
            x = x_neg + dx;
        } else {
            dx = newton_dx;
            x -= dx;
        }
        if dx.abs() < options.x_tol {
            return Ok(x);
        }
        (fx, dfx) = f_df(x);
        if fx == 0.0 {
            return Ok(x);
        }
        if fx < 0.0 {
            x_neg = x;
        } else {
            x_pos = x;
        }
    }
    Err(PhysicsError::NoConvergence {
        iterations: options.max_iterations,
        residual: fx.abs(),
    })
    .context("bracketed newton")
}

fn ensure_bracket(a: f64, fa: f64, b: f64, fb: f64) -> anyhow::Result<()> {
    if fa.is_finite() && fb.is_finite() && fa.signum() != fb.signum() {
        Ok(())
    } else {
        Err(PhysicsError::NotBracketed { a, fa, b, fb }.into())
    }
}

#[cfg(test)]
mod tests {
    use approx::relative_eq;
    use proptest::prelude::*;

    use super::*;

    fn opts() -> anyhow::Result<RootOptions> {
        RootOptions::new(1e-13, 200).context("valid options")
    }

    #[test]
    fn test_root_options_rejects_nonpositive_tolerance() -> anyhow::Result<()> {
        for tol in [0.0, -1.0, f64::NAN] {
            anyhow::ensure!(RootOptions::new(tol, 10).is_err(), "accepted tol {tol}");
        }
        anyhow::ensure!(RootOptions::new(1e-9, 0).is_err(), "accepted 0 iterations");
        let o = opts()?;
        anyhow::ensure!(
            relative_eq!(o.x_tol(), 1e-13) && o.max_iterations() == 200,
            "{o:?}"
        );
        Ok(())
    }

    #[test]
    fn test_brent_finds_cube_root_of_two() -> anyhow::Result<()> {
        let x = brent(|x| x * x * x - 2.0, 0.0, 2.0, opts()?)?;
        anyhow::ensure!((x - 2f64.cbrt()).abs() < 1e-12, "x = {x}");
        Ok(())
    }

    #[test]
    fn test_brent_accepts_reversed_interval_and_exact_endpoint() -> anyhow::Result<()> {
        let x = brent(|x| x.cos() - x, 1.0, 0.0, opts()?)?;
        anyhow::ensure!((x - 0.739_085_133_215_160_6).abs() < 1e-12, "x = {x}");
        let at_end = brent(|x| x - 3.0, 0.0, 3.0, opts()?)?;
        anyhow::ensure!(relative_eq!(at_end, 3.0), "endpoint root {at_end}");
        Ok(())
    }

    #[test]
    fn test_brent_rejects_unbracketed_interval() -> anyhow::Result<()> {
        let Err(e) = brent(|x| x * x + 1.0, -1.0, 1.0, opts()?) else {
            anyhow::bail!("found a root of x² + 1");
        };
        anyhow::ensure!(
            matches!(
                e.downcast_ref::<PhysicsError>(),
                Some(PhysicsError::NotBracketed { .. })
            ),
            "wrong error: {e:#}"
        );
        Ok(())
    }

    #[test]
    fn test_brent_reports_no_convergence_when_iterations_run_out() -> anyhow::Result<()> {
        let tight = RootOptions::new(1e-15, 2)?;
        let Err(e) = brent(|x| x.exp() - 10.0, 0.0, 10.0, tight) else {
            anyhow::bail!("converged in 2 iterations");
        };
        anyhow::ensure!(
            matches!(
                e.downcast_ref::<PhysicsError>(),
                Some(PhysicsError::NoConvergence { iterations: 2, .. })
            ),
            "wrong error: {e:#}"
        );
        Ok(())
    }

    #[test]
    fn test_newton_bracketed_converges_from_poor_start() -> anyhow::Result<()> {
        // Plain Newton from x0 = 0.1 on atan overshoots and diverges; the
        // bracket saves it.
        let x = newton_bracketed(
            |x| (x.atan(), 1.0 / (1.0 + x * x)),
            -20.0,
            30.0,
            25.0,
            opts()?,
        )?;
        anyhow::ensure!(x.abs() < 1e-12, "x = {x}");
        Ok(())
    }

    #[test]
    fn test_newton_bracketed_survives_zero_derivative() -> anyhow::Result<()> {
        // f'(0) = 0 at the start point; must bisect instead of dividing by
        // zero.
        let x = newton_bracketed(|x| (x * x * x - 8.0, 3.0 * x * x), -1.0, 5.0, 0.0, opts()?)?;
        anyhow::ensure!((x - 2.0).abs() < 1e-12, "x = {x}");
        Ok(())
    }

    #[test]
    fn test_newton_bracketed_rejects_unbracketed_interval() -> anyhow::Result<()> {
        let result = newton_bracketed(|x| (x * x + 1.0, 2.0 * x), -1.0, 1.0, 0.5, opts()?);
        anyhow::ensure!(result.is_err(), "found a root of x² + 1: {result:?}");
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_brent_finds_root_of_random_monotone_cubic(
            root in -100.0..100.0f64, a in 0.01..10.0f64, c in 0.0..5.0f64,
            below in 0.1..50.0f64, above in 0.1..50.0f64,
        ) {
            let f = |x: f64| { let d = x - root; a * d + c * d * d * d };
            let x = brent(f, root - below, root + above, RootOptions::default())
                .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            prop_assert!((x - root).abs() < 1e-9, "root {root}, found {x}");
        }

        #[test]
        fn test_newton_bracketed_matches_brent(
            target in 0.5..1000.0f64, x0 in 0.0..40.0f64,
        ) {
            // ln-type function similar to T-from-h inversion: monotone, curved.
            let f = |x: f64| x * x + x.ln_1p() - target;
            let df = |x: f64| 2.0 * x + 1.0 / (1.0 + x);
            let hi = 40.0;
            let n = newton_bracketed(|x| (f(x), df(x)), 0.0, hi, x0, RootOptions::default())
                .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            let b = brent(f, 0.0, hi, RootOptions::default())
                .map_err(|e| TestCaseError::fail(format!("{e:#}")))?;
            prop_assert!((n - b).abs() < 1e-9, "newton {n}, brent {b}");
        }
    }
}
