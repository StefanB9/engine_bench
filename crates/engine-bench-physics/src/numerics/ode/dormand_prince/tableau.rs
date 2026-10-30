//! Dormand–Prince 5(4) stages, embedded error estimate and dense-output
//! coefficients.

use super::DormandPrince45;
use crate::numerics::ode::{OdeSystem, all_finite, combine};

// Dormand & Prince (1980) tableau; error weights are b − b̂; dense output from
// Hairer's DOPRI5.
const C2: f64 = 1.0 / 5.0;
const C3: f64 = 3.0 / 10.0;
const C4: f64 = 4.0 / 5.0;
const C5: f64 = 8.0 / 9.0;
const A21: f64 = 1.0 / 5.0;
const A31: f64 = 3.0 / 40.0;
const A32: f64 = 9.0 / 40.0;
const A41: f64 = 44.0 / 45.0;
const A42: f64 = -56.0 / 15.0;
const A43: f64 = 32.0 / 9.0;
const A51: f64 = 19372.0 / 6561.0;
const A52: f64 = -25360.0 / 2187.0;
const A53: f64 = 64448.0 / 6561.0;
const A54: f64 = -212.0 / 729.0;
const A61: f64 = 9017.0 / 3168.0;
const A62: f64 = -355.0 / 33.0;
const A63: f64 = 46732.0 / 5247.0;
const A64: f64 = 49.0 / 176.0;
const A65: f64 = -5103.0 / 18656.0;
const A71: f64 = 35.0 / 384.0;
const A73: f64 = 500.0 / 1113.0;
const A74: f64 = 125.0 / 192.0;
const A75: f64 = -2187.0 / 6784.0;
const A76: f64 = 11.0 / 84.0;
const E1: f64 = 71.0 / 57600.0;
const E3: f64 = -71.0 / 16695.0;
const E4: f64 = 71.0 / 1920.0;
const E5: f64 = -17253.0 / 339_200.0;
const E6: f64 = 22.0 / 525.0;
const E7: f64 = -1.0 / 40.0;
const D1: f64 = -12_715_105_075.0 / 11_282_082_432.0;
const D3: f64 = 87_487_479_700.0 / 32_700_410_799.0;
const D4: f64 = -10_690_763_975.0 / 1_880_347_072.0;
const D5: f64 = 701_980_252_875.0 / 199_316_789_632.0;
const D6: f64 = -1_453_857_185.0 / 822_651_844.0;
const D7: f64 = 69_997_945.0 / 29_380_423.0;

/// Dense-output polynomial of the last accepted step, valid on `[t0, t0 + h]`.
#[derive(Debug, Clone, Copy)]
pub(super) struct Dense<const N: usize> {
    pub(super) t0: f64,
    pub(super) h: f64,
    pub(super) r: [[f64; N]; 5],
}

impl<const N: usize> Dense<N> {
    #[inline]
    pub(super) fn eval(&self, t: f64, out: &mut [f64; N]) {
        let theta = (t - self.t0) / self.h;
        let theta1 = 1.0 - theta;
        let [r1, r2, r3, r4, r5] = &self.r;
        for i in 0..N {
            out[i] = r1[i] + theta * (r2[i] + theta1 * (r3[i] + theta * (r4[i] + theta1 * r5[i])));
        }
    }
}

pub(super) struct Trial<const N: usize> {
    pub(super) y: [f64; N],
    pub(super) k7: [f64; N],
    /// RMS of the scaled error estimate; NaN if any stage was non-finite.
    pub(super) error: f64,
    pub(super) dense: Dense<N>,
}

/// `(err / scale)²`, with an exact zero error on a zero-scale component
/// counting as zero.
#[inline]
fn scaled_square(err: f64, scale: f64) -> f64 {
    if err == 0.0 {
        0.0
    } else {
        (err / scale).powi(2)
    }
}

impl<const N: usize> DormandPrince45<N> {
    pub(super) fn trial_step<S: OdeSystem<N>>(
        &mut self,
        system: &mut S,
        t: f64,
        y: &[f64; N],
        k1: &[f64; N],
        h: f64,
    ) -> Trial<N> {
        let mut k2 = [0.0; N];
        let mut k3 = [0.0; N];
        let mut k4 = [0.0; N];
        let mut k5 = [0.0; N];
        let mut k6 = [0.0; N];
        let mut k7 = [0.0; N];
        let mut stage = [0.0; N];

        combine(&mut stage, y, &[(h * A21, k1)]);
        self.evaluate(system, t + C2 * h, &stage, &mut k2);
        combine(&mut stage, y, &[(h * A31, k1), (h * A32, &k2)]);
        self.evaluate(system, t + C3 * h, &stage, &mut k3);
        combine(
            &mut stage,
            y,
            &[(h * A41, k1), (h * A42, &k2), (h * A43, &k3)],
        );
        self.evaluate(system, t + C4 * h, &stage, &mut k4);
        combine(
            &mut stage,
            y,
            &[
                (h * A51, k1),
                (h * A52, &k2),
                (h * A53, &k3),
                (h * A54, &k4),
            ],
        );
        self.evaluate(system, t + C5 * h, &stage, &mut k5);
        combine(
            &mut stage,
            y,
            &[
                (h * A61, k1),
                (h * A62, &k2),
                (h * A63, &k3),
                (h * A64, &k4),
                (h * A65, &k5),
            ],
        );
        self.evaluate(system, t + h, &stage, &mut k6);
        let mut y1 = [0.0; N];
        combine(
            &mut y1,
            y,
            &[
                (h * A71, k1),
                (h * A73, &k3),
                (h * A74, &k4),
                (h * A75, &k5),
                (h * A76, &k6),
            ],
        );
        self.evaluate(system, t + h, &y1, &mut k7);

        let mut sum = 0.0;
        let mut dense = Dense {
            t0: t,
            h,
            r: [[0.0; N]; 5],
        };
        for i in 0..N {
            let err =
                h * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i] + E5 * k5[i] + E6 * k6[i] + E7 * k7[i]);
            sum += scaled_square(err, self.scale(y[i], y1[i]));

            let diff = y1[i] - y[i];
            let b_spline = h * k1[i] - diff;
            dense.r[0][i] = y[i];
            dense.r[1][i] = diff;
            dense.r[2][i] = b_spline;
            dense.r[3][i] = diff - h * k7[i] - b_spline;
            dense.r[4][i] =
                h * (D1 * k1[i] + D3 * k3[i] + D4 * k4[i] + D5 * k5[i] + D6 * k6[i] + D7 * k7[i]);
        }
        #[expect(clippy::cast_precision_loss, reason = "N is a small state dimension")]
        let error = (sum / N.max(1) as f64).sqrt();
        // Non-finite stages must never be accepted.
        let error = if all_finite(&y1) && all_finite(&k7) {
            error
        } else {
            f64::NAN
        };
        Trial {
            y: y1,
            k7,
            error,
            dense,
        }
    }

    pub(super) fn dense_at(&self, t: f64, out: &mut [f64; N]) {
        if let Some(dense) = &self.dense {
            dense.eval(t, out);
        }
    }
}
