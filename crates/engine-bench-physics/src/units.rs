//! SI quantities for public APIs.
//!
//! `uom` quantities are a single `f64` at runtime (no conversion cost, verified
//! by test), so they are used at every API boundary. Hot kernels unwrap to raw
//! `f64` struct-of-arrays and re-wrap results on the way out.

pub use uom::si::{
    Quantity, acceleration, angle, angular_velocity, area, available_energy, energy, f64::*, force,
    heat_capacity, length, mass, mass_density, mass_rate, molar_energy, molar_heat_capacity,
    molar_mass, power, pressure, ratio, specific_heat_capacity, thermodynamic_temperature, time,
    torque, velocity, volume,
};

#[cfg(test)]
mod tests {
    use approx::relative_eq;

    use super::*;

    #[test]
    fn test_units_pressure_converts_bar_to_pascal() -> anyhow::Result<()> {
        let pa = Pressure::new::<pressure::bar>(1.0).get::<pressure::pascal>();
        anyhow::ensure!(relative_eq!(pa, 100_000.0), "1 bar = {pa} Pa");
        Ok(())
    }

    #[test]
    fn test_units_temperature_converts_celsius_to_kelvin() -> anyhow::Result<()> {
        let k = ThermodynamicTemperature::new::<thermodynamic_temperature::degree_celsius>(25.0)
            .get::<thermodynamic_temperature::kelvin>();
        anyhow::ensure!(relative_eq!(k, 298.15), "25 °C = {k} K");
        Ok(())
    }

    // Guards the zero-cost claim the module relies on: if uom ever adds a
    // field, this fails.
    #[test]
    fn test_units_quantities_are_plain_f64_in_memory() -> anyhow::Result<()> {
        for (name, size, align) in [
            ("Pressure", size_of::<Pressure>(), align_of::<Pressure>()),
            ("MassRate", size_of::<MassRate>(), align_of::<MassRate>()),
            ("Velocity", size_of::<Velocity>(), align_of::<Velocity>()),
        ] {
            anyhow::ensure!(
                size == size_of::<f64>() && align == align_of::<f64>(),
                "{name} is {size} bytes / align {align}, expected f64 layout"
            );
        }
        Ok(())
    }
}
