//! Wall-clock time behind a trait, so pacing logic can be tested with a mock.

use std::{cell::Cell, time::Duration};

/// A monotonic clock.
pub trait Clock {
    /// Time since the clock's origin; never decreases.
    fn now(&self) -> Duration;
}

/// The real monotonic clock, measured from its creation.
#[derive(Debug, Clone, Copy)]
pub struct SystemClock {
    origin: std::time::Instant,
}

impl SystemClock {
    /// A clock starting at zero now.
    pub fn new() -> Self {
        Self {
            origin: std::time::Instant::now(),
        }
    }
}

impl Default for SystemClock {
    fn default() -> Self {
        Self::new()
    }
}

impl Clock for SystemClock {
    fn now(&self) -> Duration {
        self.origin.elapsed()
    }
}

/// A clock that only moves when told to (tests).
#[derive(Debug, Clone, Default)]
pub struct MockClock {
    now: Cell<Duration>,
}

impl MockClock {
    /// A mock clock at zero.
    pub fn new() -> Self {
        Self::default()
    }

    /// Jumps to `time`.
    pub fn set_time(&self, time: Duration) {
        self.now.set(time);
    }

    /// Moves forward by `step`.
    pub fn advance(&self, step: Duration) {
        self.now.set(self.now.get() + step);
    }
}

impl Clock for MockClock {
    fn now(&self) -> Duration {
        self.now.get()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mock_clock_set_and_advance() -> anyhow::Result<()> {
        let clock = MockClock::new();
        anyhow::ensure!(clock.now() == Duration::ZERO, "start");
        clock.advance(Duration::from_millis(16));
        clock.advance(Duration::from_millis(4));
        anyhow::ensure!(
            clock.now() == Duration::from_millis(20),
            "{:?}",
            clock.now()
        );
        clock.set_time(Duration::from_secs(3));
        anyhow::ensure!(clock.now() == Duration::from_secs(3), "{:?}", clock.now());
        Ok(())
    }

    #[test]
    fn test_system_clock_moves_forward() -> anyhow::Result<()> {
        let clock = SystemClock::new();
        let first = clock.now();
        std::thread::sleep(Duration::from_millis(2));
        let second = clock.now();
        anyhow::ensure!(
            second >= first + Duration::from_millis(2),
            "{first:?} → {second:?}"
        );
        Ok(())
    }
}
