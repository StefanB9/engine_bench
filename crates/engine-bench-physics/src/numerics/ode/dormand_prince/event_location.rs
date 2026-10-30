//! Locating event crossings on the dense output of an accepted step.

use super::DormandPrince45;
use crate::numerics::ode::{EventSet, events::crosses};

impl<const N: usize> DormandPrince45<N> {
    /// Earliest event crossing in `(t0, t1]`, located on the dense output to
    /// the far side of the crossing so that resuming from it cannot fire the
    /// same event again.
    pub(super) fn locate_event<E: EventSet<N, M>, const M: usize>(
        &self,
        events: &mut E,
        t0: f64,
        y0: &[f64; N],
        t1: f64,
        y1: &[f64; N],
    ) -> Option<(usize, f64)> {
        if M == 0 {
            return None;
        }
        let directions = events.directions();
        let mut g0 = [0.0; M];
        let mut g1 = [0.0; M];
        events.values(t0, y0, &mut g0);
        events.values(t1, y1, &mut g1);

        let mut earliest: Option<(usize, f64)> = None;
        for i in 0..M {
            if !crosses(directions[i], g0[i], g1[i]) {
                continue;
            }
            let t_event = self.bracket_crossing(events, i, (t0, g0[i]), (t1, g1[i]));
            if earliest.is_none_or(|(_, t_best)| t_event < t_best) {
                earliest = Some((i, t_event));
            }
        }
        earliest
    }

    /// Illinois-modified regula falsi on `g_index` along the dense output.
    /// Returns the end of the final bracket on the post-crossing side.
    fn bracket_crossing<E: EventSet<N, M>, const M: usize>(
        &self,
        events: &mut E,
        index: usize,
        (mut t_pre, mut g_pre): (f64, f64),
        (mut t_post, mut g_post): (f64, f64),
    ) -> f64 {
        let pre_positive = g_pre > 0.0;
        let tolerance = 4.0 * f64::EPSILON * t_post.abs().max(1.0);
        let mut y = [0.0; N];
        let mut g = [0.0; M];
        let mut retained_side = 0_i8;
        for _ in 0..200 {
            if t_post - t_pre <= tolerance {
                break;
            }
            let secant = t_post - g_post * (t_post - t_pre) / (g_post - g_pre);
            let t_mid = if secant > t_pre && secant < t_post {
                secant
            } else {
                f64::midpoint(t_pre, t_post)
            };
            self.dense_at(t_mid, &mut y);
            events.values(t_mid, &y, &mut g);
            let g_mid = g[index];
            if g_mid != 0.0 && (g_mid > 0.0) == pre_positive {
                (t_pre, g_pre) = (t_mid, g_mid);
                if retained_side == 1 {
                    g_post *= 0.5;
                }
                retained_side = 1;
            } else {
                (t_post, g_post) = (t_mid, g_mid);
                if g_mid == 0.0 {
                    break;
                }
                if retained_side == -1 {
                    g_pre *= 0.5;
                }
                retained_side = -1;
            }
        }
        t_post
    }
}
