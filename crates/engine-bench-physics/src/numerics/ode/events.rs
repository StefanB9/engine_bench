//! Event functions whose zero crossings stop the adaptive integrator.

/// Which zero crossings of an event function trigger it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// `g` goes from negative to non-negative.
    Rising,
    /// `g` goes from positive to non-positive.
    Falling,
    /// Either crossing.
    Either,
}

/// `M` event functions `g_i(t, y)` whose zero crossings stop the integration.
///
/// Crossings are detected by sign change across an accepted step and located on
/// the dense output, so a function that crosses zero twice within one step is
/// missed; cap the step with
/// [`AdaptiveOptions::with_max_step`](super::AdaptiveOptions::with_max_step) if
/// that can happen.
pub trait EventSet<const N: usize, const M: usize> {
    /// Crossing direction for each event.
    fn directions(&self) -> [Direction; M];
    /// Writes `g_i(t, y)` for every event into `g`.
    fn values(&mut self, t: f64, y: &[f64; N], g: &mut [f64; M]);
}

/// The empty event set.
#[derive(Debug, Clone, Copy, Default)]
pub struct NoEvents;

impl<const N: usize> EventSet<N, 0> for NoEvents {
    fn directions(&self) -> [Direction; 0] {
        []
    }

    fn values(&mut self, _t: f64, _y: &[f64; N], _g: &mut [f64; 0]) {}
}

/// Whether `g` crossed zero from `g0` to `g1` in `direction`.
pub(super) fn crosses(direction: Direction, g0: f64, g1: f64) -> bool {
    let rising = g0 < 0.0 && g1 >= 0.0;
    let falling = g0 > 0.0 && g1 <= 0.0;
    match direction {
        Direction::Rising => rising,
        Direction::Falling => falling,
        Direction::Either => rising || falling,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crosses_respects_direction() -> anyhow::Result<()> {
        use Direction::{Either, Falling, Rising};
        let cases = [
            (Rising, -1.0, 1.0, true),
            (Rising, -1.0, 0.0, true),
            (Rising, 1.0, -1.0, false),
            (Rising, 0.0, 1.0, false),
            (Falling, 1.0, -1.0, true),
            (Falling, 1.0, 0.0, true),
            (Falling, -1.0, 1.0, false),
            (Falling, 0.0, -1.0, false),
            (Either, -1.0, 1.0, true),
            (Either, 1.0, -1.0, true),
            (Either, 1.0, 2.0, false),
        ];
        for (direction, g0, g1, expected) in cases {
            let got = crosses(direction, g0, g1);
            anyhow::ensure!(got == expected, "{direction:?} {g0} → {g1}: {got}");
        }
        Ok(())
    }

    #[test]
    fn test_no_events_is_empty() -> anyhow::Result<()> {
        let mut g = [];
        EventSet::<1, 0>::values(&mut NoEvents, 0.0, &[1.0], &mut g);
        anyhow::ensure!(
            EventSet::<1, 0>::directions(&NoEvents).is_empty(),
            "has directions"
        );
        Ok(())
    }
}
