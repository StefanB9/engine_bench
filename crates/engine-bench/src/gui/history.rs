//! Recent channel values for live plots, bounded in length.

use std::collections::VecDeque;

/// The last `capacity` rows of every channel.
#[derive(Debug, Clone)]
pub(crate) struct PlotHistory {
    capacity: usize,
    time: VecDeque<f64>,
    columns: Vec<VecDeque<f64>>,
}

impl PlotHistory {
    /// Empty history of `channels` channels keeping at most `capacity` rows.
    pub(crate) fn new(channels: usize, capacity: usize) -> Self {
        Self {
            capacity,
            time: VecDeque::with_capacity(capacity),
            columns: (0..channels)
                .map(|_| VecDeque::with_capacity(capacity))
                .collect(),
        }
    }

    /// Appends a row. Time going backwards (reset, engine switch) starts a new
    /// history; a different channel count re-shapes it.
    pub(crate) fn push(&mut self, time: f64, sample: &[f64]) {
        let went_back = self.time.back().is_some_and(|last| time < *last);
        if went_back || sample.len() != self.columns.len() {
            *self = Self::new(sample.len(), self.capacity);
        }
        if self.time.len() == self.capacity {
            self.time.pop_front();
            for column in &mut self.columns {
                column.pop_front();
            }
        }
        self.time.push_back(time);
        for (column, value) in self.columns.iter_mut().zip(sample) {
            column.push_back(*value);
        }
    }

    /// Rows held.
    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.time.len()
    }

    /// `[time, value]` pairs of `channel`, oldest first, ready for plotting.
    pub(crate) fn points(&self, channel: usize) -> Vec<[f64; 2]> {
        let Some(column) = self.columns.get(channel) else {
            return Vec::new();
        };
        self.time
            .iter()
            .zip(column)
            .map(|(t, v)| [*t, *v])
            .collect()
    }

    /// Forgets every row.
    pub(crate) fn clear(&mut self) {
        self.time.clear();
        for column in &mut self.columns {
            column.clear();
        }
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    #[test]
    fn test_history_keeps_rows_in_order_per_channel() -> anyhow::Result<()> {
        let mut history = PlotHistory::new(2, 10);
        history.push(0.0, &[1.0, 10.0]);
        history.push(0.1, &[2.0, 20.0]);
        anyhow::ensure!(history.len() == 2, "len");
        anyhow::ensure!(
            history.points(1) == [[0.0, 10.0], [0.1, 20.0]],
            "{:?}",
            history.points(1)
        );
        anyhow::ensure!(history.points(5).is_empty(), "missing channel has points");
        Ok(())
    }

    #[test]
    fn test_history_restarts_when_time_goes_back() -> anyhow::Result<()> {
        let mut history = PlotHistory::new(1, 10);
        history.push(1.0, &[1.0]);
        history.push(2.0, &[2.0]);
        history.push(0.0, &[5.0]);
        anyhow::ensure!(history.points(0) == [[0.0, 5.0]], "{:?}", history.points(0));
        Ok(())
    }

    #[test]
    fn test_history_reshapes_for_new_channel_count() -> anyhow::Result<()> {
        let mut history = PlotHistory::new(1, 10);
        history.push(0.0, &[1.0]);
        history.push(0.1, &[1.0, 2.0, 3.0]);
        anyhow::ensure!(
            history.len() == 1 && history.points(2) == [[0.1, 3.0]],
            "{:?}",
            history.points(2)
        );
        history.clear();
        anyhow::ensure!(history.len() == 0, "clear");
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        #[test]
        fn test_history_never_exceeds_capacity_and_keeps_newest(
            capacity in 1..50usize, rows in 0..200u32,
        ) {
            let mut history = PlotHistory::new(1, capacity);
            for row in 0..rows {
                history.push(f64::from(row), &[f64::from(row)]);
            }
            prop_assert!(history.len() <= capacity);
            let points = history.points(0);
            if rows > 0 {
                prop_assert_eq!(points.last().map(|p| p[0]), Some(f64::from(rows - 1)));
            }
            prop_assert!(points.windows(2).all(|w| w[0][0] < w[1][0]));
        }
    }
}
