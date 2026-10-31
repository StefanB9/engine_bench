//! The pressure regulator: a spring-loaded, unbalanced poppet.
//!
//! The main spring, preloaded to the setpoint, pushes the poppet open; outlet
//! pressure under the diaphragm and inlet pressure on the seated poppet push
//! it closed. Lift opens a curtain area around the seat. The result shows
//! droop (outlet pressure falls with flow, as the spring relaxes), lock-up
//! (closed at a little above the setpoint without flow) and the supply
//! pressure effect (outlet pressure rises as the tank empties).

use std::f64::consts::PI;

/// Diaphragm area, m² (25 mm diameter).
const DIAPHRAGM_AREA: f64 = PI * 0.025 * 0.025 / 4.0;
/// Seat diameter, m.
const SEAT_DIAMETER: f64 = 1e-3;
/// Seat area, m².
const SEAT_AREA: f64 = PI * SEAT_DIAMETER * SEAT_DIAMETER / 4.0;
/// Spring rate, N/m: about 4 % droop at full flow.
const SPRING_RATE: f64 = 8e4;
/// Moving mass (poppet, stem, diaphragm), kg.
const MASS: f64 = 0.03;
/// Dashpot, N·s/m. Lift, flow and plenum pressure form a loop, m·s³ + c·s² +
/// k·s + K = 0 with K the diaphragm area × γRT/V × dṁ/dx, stable while
/// c·k > m·K (Routh); this keeps it so for every propellant up to full tank
/// pressure.
const DAMPING: f64 = 200.0;
/// Full lift, m.
pub(super) const FULL_LIFT: f64 = 5e-4;
/// Discharge coefficient of the curtain.
pub(super) const DISCHARGE_COEFFICIENT: f64 = 0.7;

/// Poppet acceleration, m/s², at lift `x` m and velocity `v` m/s with inlet
/// pressure `p_in` and outlet pressure `p_out` Pa, for setpoint `p_set` Pa.
pub(super) fn acceleration(x: f64, v: f64, p_in: f64, p_out: f64, p_set: f64) -> f64 {
    let force = (p_set - p_out) * DIAPHRAGM_AREA
        - SPRING_RATE * x
        - (p_in - p_out) * SEAT_AREA
        - DAMPING * v;
    force / MASS
}

/// Lift and lift rate rates of change, with the poppet held at its stops:
/// it neither moves nor accelerates past closed or full lift.
pub(super) fn motion(x: f64, v: f64, p_in: f64, p_out: f64, p_set: f64) -> (f64, f64) {
    let a = acceleration(x, v, p_in, p_out, p_set);
    let closed = x <= 0.0 && v <= 0.0 && a <= 0.0;
    let full = x >= FULL_LIFT && v >= 0.0 && a >= 0.0;
    if closed || full { (0.0, 0.0) } else { (v, a) }
}

/// The lift and velocity after a step, put back within the stops (the poppet
/// hits them inelastically).
pub(super) fn at_stops(x: f64, v: f64) -> (f64, f64) {
    if x <= 0.0 {
        (0.0, v.max(0.0))
    } else if x >= FULL_LIFT {
        (FULL_LIFT, v.min(0.0))
    } else {
        (x, v)
    }
}

/// Geometric flow area at lift `x`, m²: the curtain around the seat, at most
/// the seat itself.
pub(super) fn flow_area(x: f64) -> f64 {
    (PI * SEAT_DIAMETER * x.clamp(0.0, FULL_LIFT)).min(SEAT_AREA)
}

/// Outlet pressure at which the poppet closes with inlet pressure `p_in`
/// (lock-up), Pa.
pub(super) fn lock_up(p_in: f64, p_set: f64) -> f64 {
    (p_set * DIAPHRAGM_AREA - p_in * SEAT_AREA) / (DIAPHRAGM_AREA - SEAT_AREA)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_regulator_holds_still_at_lock_up() -> anyhow::Result<()> {
        let p_out = lock_up(2e7, 1e6);
        let a = acceleration(0.0, 0.0, 2e7, p_out, 1e6);
        anyhow::ensure!(a.abs() < 1e-6, "a = {a}");
        // The supply pressure effect: a fuller tank closes it at a lower
        // outlet pressure.
        anyhow::ensure!(lock_up(2e7, 1e6) < lock_up(2e6, 1e6), "supply effect");
        anyhow::ensure!((lock_up(1e6, 1e6) - 1e6).abs() < 1e-6, "no differential");
        Ok(())
    }

    #[test]
    fn test_regulator_opens_below_and_closes_above_the_setpoint() {
        assert!(acceleration(1e-4, 0.0, 2e7, 5e5, 1e6) > 0.0);
        assert!(acceleration(1e-4, 0.0, 2e7, 1.5e6, 1e6) < 0.0);
    }

    #[test]
    fn test_regulator_stops_hold_the_poppet() {
        // Closed and pushed shut, or fully open and pushed open: no motion.
        assert_eq!(motion(0.0, 0.0, 2e7, 2e6, 1e6), (0.0, 0.0));
        assert_eq!(motion(FULL_LIFT, 0.0, 2e7, 0.0, 1e6), (0.0, 0.0));
        // Closed but pushed open: it starts to lift.
        let (dx, dv) = motion(0.0, 0.0, 2e7, 5e5, 1e6);
        assert!(dx.abs() < f64::EPSILON && dv > 0.0);
        assert_eq!(at_stops(-1e-6, -0.1), (0.0, 0.0));
        assert_eq!(at_stops(1e-3, 0.2), (FULL_LIFT, 0.0));
        assert_eq!(at_stops(1e-4, 0.2), (1e-4, 0.2));
    }

    #[test]
    fn test_regulator_curtain_area_saturates_at_the_seat() {
        assert!(flow_area(0.0).abs() < f64::EPSILON);
        assert!((flow_area(1e-4) - PI * SEAT_DIAMETER * 1e-4).abs() < 1e-18);
        assert!((flow_area(FULL_LIFT) - SEAT_AREA).abs() < 1e-18);
        assert!((flow_area(-1.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn test_regulator_pressure_loop_is_stable_at_full_tank_pressure() {
        // Worst case: helium (high γRT, so the plenum pressure moves fastest)
        // from a full 300 bar tank at 293 K, perfect-gas estimates. Choked
        // mass flux G* = p·√γ·(2/(γ+1))^((γ+1)/(2(γ−1))) / √(RT).
        for (gamma, r) in [
            (5.0_f64 / 3.0, 2077.3_f64),
            (1.4, 296.8),
            (5.0 / 3.0, 208.1),
        ] {
            let t = 293.0;
            let critical: f64 = (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)));
            let flux = 3e7 * gamma.sqrt() * critical / (r * t).sqrt();
            let flow_gain = DISCHARGE_COEFFICIENT * PI * SEAT_DIAMETER * flux;
            let pressure_gain = gamma * r * t / super::super::hardware::PLENUM_VOLUME;
            let loop_gain = DIAPHRAGM_AREA * pressure_gain * flow_gain;
            let margin = DAMPING * SPRING_RATE / (MASS * loop_gain);
            assert!(margin > 1.5, "R = {r}: Routh margin {margin}");
        }
    }
}
