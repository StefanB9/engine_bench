//! Hardware sized from the controls: the tank and its wall, the plumbing
//! volumes and the nozzle.

use std::f64::consts::PI;

/// The tank is a cylinder this many diameters long.
const TANK_ASPECT: f64 = 3.0;
/// Aluminium wall: yield stress, Pa, safety factor, density, kg/m³, and
/// specific heat, J/(kg·K).
const YIELD_STRESS: f64 = 2.75e8;
const SAFETY_FACTOR: f64 = 2.0;
const WALL_DENSITY: f64 = 2700.0;
const WALL_SPECIFIC_HEAT: f64 = 900.0;
/// Thinnest wall, m: below this a tank is not made.
const MIN_WALL: f64 = 1e-3;
/// Line and manifold from the regulator to the valve, m³.
pub(super) const PLENUM_VOLUME: f64 = 1e-4;
/// Between the valve seat and the throat (outlet passage and nozzle inlet),
/// m³.
pub(super) const OUTLET_VOLUME: f64 = 2e-6;

/// The tank's shell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Tank {
    /// Inner volume, m³.
    pub(super) volume: f64,
    /// Inner diameter, m.
    pub(super) diameter: f64,
    /// Wall area, m² (inner, also used outside: the wall is thin).
    pub(super) area: f64,
    /// Wall heat capacity, J/K.
    pub(super) heat_capacity: f64,
}

impl Tank {
    /// A cylindrical tank of `volume` m³ whose wall carries `pressure` Pa: hoop
    /// stress p·r/t at the yield stress over the safety factor.
    pub(super) fn sized(volume: f64, pressure: f64) -> Self {
        let diameter = (4.0 * volume / (PI * TANK_ASPECT)).cbrt();
        let length = TANK_ASPECT * diameter;
        let area = PI * diameter * length + 2.0 * PI * diameter * diameter / 4.0;
        let thickness = (pressure * 0.5 * diameter * SAFETY_FACTOR / YIELD_STRESS).max(MIN_WALL);
        Self {
            volume,
            diameter,
            area,
            heat_capacity: WALL_DENSITY * area * thickness * WALL_SPECIFIC_HEAT,
        }
    }

    /// Length, m.
    #[cfg(test)]
    pub(super) fn length(&self) -> f64 {
        TANK_ASPECT * self.diameter
    }
}

/// Area of a circle of `diameter` m, m².
pub(super) fn circle_area(diameter: f64) -> f64 {
    PI * diameter * diameter / 4.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tank_holds_its_volume_in_a_cylinder() {
        let tank = Tank::sized(2e-3, 2e7);
        let volume = circle_area(tank.diameter) * tank.length();
        assert!((volume / 2e-3 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_tank_wall_carries_the_pressure() {
        // 2 L at 200 bar: D ≈ 9.5 cm, wall ≈ 6.9 mm, about 1.3 kg of aluminium.
        let tank = Tank::sized(2e-3, 2e7);
        let thickness = 2e7 * 0.5 * tank.diameter * 2.0 / 2.75e8;
        let mass = tank.heat_capacity / 900.0;
        assert!((mass - 2700.0 * tank.area * thickness).abs() < 1e-9);
        assert!(mass > 1.0 && mass < 2.0, "{mass} kg");
        // A low-pressure tank still has a minimum wall.
        let thin = Tank::sized(2e-3, 1e5);
        assert!((thin.heat_capacity - 2700.0 * thin.area * 1e-3 * 900.0).abs() < 1e-9);
    }
}
