//! The solenoid fire valve: commanded open or shut by the fire mode, it
//! travels between closed and open in finite times.

/// Time to travel from closed to open, s.
pub(super) const OPENING_TIME: f64 = 5e-3;
/// Time to travel from open to closed, s.
pub(super) const CLOSING_TIME: f64 = 3e-3;
/// Open flow area as a multiple of the throat: the valve must not choke.
pub(super) const AREA_PER_THROAT: f64 = 4.0;
/// Discharge coefficient when open.
pub(super) const DISCHARGE_COEFFICIENT: f64 = 0.8;

/// How the valve is commanded.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum FireMode {
    /// Shut.
    Off,
    /// Open.
    Continuous,
    /// Open for a pulse width at the start of every period.
    Pulsed,
}

impl FireMode {
    /// The mode at option `index` of the fire-mode control.
    pub(super) fn from_index(index: usize) -> Self {
        match index {
            1 => Self::Continuous,
            2 => Self::Pulsed,
            _ => Self::Off,
        }
    }

    /// Whether the valve is commanded open at time `t` s, pulses lasting
    /// `width` s every `period` s from t = 0.
    pub(super) fn open(self, t: f64, width: f64, period: f64) -> bool {
        match self {
            Self::Off => false,
            Self::Continuous => true,
            Self::Pulsed => t.rem_euclid(period) < width,
        }
    }
}

/// Rate of change of the valve position (0 shut, 1 open), 1/s, travelling
/// towards the command at constant speed and resting at the ends.
pub(super) fn travel(position: f64, open: bool) -> f64 {
    if open && position < 1.0 {
        1.0 / OPENING_TIME
    } else if !open && position > 0.0 {
        -1.0 / CLOSING_TIME
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fire_mode_follows_the_control_options() {
        assert_eq!(FireMode::from_index(0), FireMode::Off);
        assert_eq!(FireMode::from_index(1), FireMode::Continuous);
        assert_eq!(FireMode::from_index(2), FireMode::Pulsed);
    }

    #[test]
    fn test_fire_mode_commands_pulses_at_the_start_of_each_period() {
        let pulsed = FireMode::Pulsed;
        assert!(pulsed.open(0.0, 0.05, 0.5));
        assert!(pulsed.open(0.049, 0.05, 0.5));
        assert!(!pulsed.open(0.051, 0.05, 0.5));
        assert!(pulsed.open(1.02, 0.05, 0.5));
        assert!(!FireMode::Off.open(0.0, 0.05, 0.5));
        assert!(FireMode::Continuous.open(0.3, 0.05, 0.5));
    }

    #[test]
    fn test_valve_travels_in_its_opening_and_closing_times() {
        assert!((travel(0.0, true) * OPENING_TIME - 1.0).abs() < 1e-12);
        assert!((travel(1.0, false) * CLOSING_TIME + 1.0).abs() < 1e-12);
        assert!(travel(1.0, true).abs() < f64::EPSILON);
        assert!(travel(0.0, false).abs() < f64::EPSILON);
    }
}
