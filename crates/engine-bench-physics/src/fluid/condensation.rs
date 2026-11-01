//! Condensation properties of the reference fluids: surface tension from
//! Mulero et al. (2012), latent heat, liquid density and vapour pressure from
//! the reference equations; below the triple point the liquid is supercooled.

use anyhow::Context;

use super::{
    PureFluid, helmholtz::molar, saturation::dome_top, saturation_si, transport::dilute_transport,
};
use crate::{
    PhysicsError,
    flow::condensing::{Condensable, CondensationProperties, VapourTransport},
    units::{
        AvailableEnergy, DynamicViscosity, Mass, MassDensity, SurfaceTension, ThermalConductivity,
        ThermodynamicTemperature, available_energy::joule_per_kilogram,
        dynamic_viscosity::pascal_second, mass::kilogram, mass_density::kilogram_per_cubic_meter,
        surface_tension::newton_per_meter, thermal_conductivity::watt_per_meter_kelvin,
        thermodynamic_temperature::kelvin,
    },
};

/// Avogadro's number, 1/mol.
const AVOGADRO: f64 = 6.022_140_76e23;

/// Mulero, Cachadiña and Cuadrado (2012), J. Phys. Chem. Ref. Data 41, 043105:
/// σ = Σ aᵢ·(1 − T/Tc)^nᵢ N/m, with their critical temperatures.
fn mulero(fluid: PureFluid) -> (f64, &'static [(f64, f64)]) {
    match fluid {
        PureFluid::Nitrogen => (126.192, &[(0.02898, 1.246)]),
        PureFluid::Argon => (150.687, &[(0.037, 1.25)]),
        PureFluid::Helium => (
            5.1953,
            &[(0.000_465_6, 1.04), (0.001_889, 2.468), (-0.002_006, 2.661)],
        ),
        PureFluid::CarbonDioxide => (304.128, &[(0.078_63, 1.254)]),
    }
}

/// Surface tension, N/m, at `t` K (below the critical temperature).
fn surface_tension(fluid: PureFluid, t: f64) -> f64 {
    let (t_crit, terms) = mulero(fluid);
    let tau = (1.0 - t / t_crit).max(0.0);
    terms.iter().map(|&(a, n)| a * tau.powf(n)).sum()
}

impl PureFluid {
    /// Saturated liquid density, kg/m³, latent heat, J/kg, and vapour
    /// pressure, Pa, at `t` K from the triple point to the dome's top.
    #[expect(
        clippy::many_single_char_names,
        reason = "the thermodynamic symbols p, T, m and liquid l, vapour v"
    )]
    fn saturated_phases(self, t: f64) -> anyhow::Result<(f64, f64, f64)> {
        let fluid = self.data();
        let m = fluid.data.molar_mass;
        let (p, rho_l, rho_v) = saturation_si(fluid, t)?;
        let (l, v) = (molar(fluid, rho_l, t), molar(fluid, rho_v, t));
        Ok((rho_l * m, (v.h - l.h) / m, p))
    }
}

impl Condensable for PureFluid {
    /// Below the triple point the liquid is supercooled: its density follows
    /// the saturated liquid's slope at the triple point, the latent heat stays
    /// at its triple-point value and the vapour pressure is the
    /// Clausius–Clapeyron extrapolation of [`PureFluid::vapour_pressure`] —
    /// the three consistent with each other.
    fn condensation(&self, t: ThermodynamicTemperature) -> anyhow::Result<CondensationProperties> {
        let fluid = self.data();
        let t = t.get::<kelvin>();
        let (t_triple, top) = (fluid.data.t_triple, dome_top(fluid));
        if !(t > 0.0 && t <= top) {
            return Err(PhysicsError::TemperatureOutOfRange {
                temperature: t,
                min: 0.0,
                max: top,
            })
            .with_context(|| format!("{self:?} condensation properties"));
        }
        let (liquid_density, latent_heat, vapour_pressure) = if t >= t_triple {
            self.saturated_phases(t)?
        } else {
            let (rho, latent, _) = self.saturated_phases(t_triple)?;
            let t_next = t_triple * (1.0 + 1e-4);
            let (rho_next, _, _) = self.saturated_phases(t_next)?;
            let slope = (rho_next - rho) / (t_next - t_triple);
            let p = self.vapour_pressure(ThermodynamicTemperature::new::<kelvin>(t))?;
            (
                slope.mul_add(t - t_triple, rho),
                latent,
                p.get::<crate::units::pressure::pascal>(),
            )
        };
        Ok(CondensationProperties {
            vapour_pressure: crate::units::Pressure::new::<crate::units::pressure::pascal>(
                vapour_pressure,
            ),
            liquid_density: MassDensity::new::<kilogram_per_cubic_meter>(liquid_density),
            latent_heat: AvailableEnergy::new::<joule_per_kilogram>(latent_heat),
            surface_tension: SurfaceTension::new::<newton_per_meter>(surface_tension(*self, t)),
        })
    }

    fn vapour_transport(&self, t: ThermodynamicTemperature) -> anyhow::Result<VapourTransport> {
        let kelvin_t = t.get::<kelvin>();
        if !(kelvin_t > 0.0 && kelvin_t.is_finite()) {
            return Err(PhysicsError::InvalidState("temperature must be positive"))
                .with_context(|| format!("{self:?} vapour transport at {kelvin_t} K"));
        }
        let (mu, k) = dilute_transport(*self, kelvin_t);
        Ok(VapourTransport {
            viscosity: DynamicViscosity::new::<pascal_second>(mu),
            thermal_conductivity: ThermalConductivity::new::<watt_per_meter_kelvin>(k),
        })
    }

    fn molecular_mass(&self) -> Mass {
        Mass::new::<kilogram>(self.data().data.molar_mass / AVOGADRO)
    }
}

#[cfg(test)]
#[expect(
    clippy::unreadable_literal,
    reason = "golden values exactly as tools/reference/condensation_golden.py prints them"
)]
mod tests {
    use proptest::prelude::*;

    use super::*;
    use crate::units::{
        available_energy::joule_per_kilogram, mass::kilogram,
        mass_density::kilogram_per_cubic_meter, pressure::pascal,
        surface_tension::newton_per_meter, thermodynamic_temperature::kelvin,
    };

    fn at(fluid: PureFluid, t: f64) -> anyhow::Result<[f64; 4]> {
        let x = fluid.condensation(ThermodynamicTemperature::new::<kelvin>(t))?;
        Ok([
            x.surface_tension.get::<newton_per_meter>(),
            x.latent_heat.get::<joule_per_kilogram>(),
            x.liquid_density.get::<kilogram_per_cubic_meter>(),
            x.vapour_pressure.get::<pascal>(),
        ])
    }

    /// (T K, σ N/m, latent heat J/kg, liquid kg/m³, vapour pressure Pa).
    type Row = (f64, f64, f64, f64, f64);
    type Golden = (PureFluid, [Row; 3]);

    // (T K, σ N/m, latent heat J/kg, liquid density kg/m³, p_sat Pa) from
    // tools/reference/condensation_golden.py (CoolProp 8.0).
    const GOLDEN: [Golden; 4] = [
        (
            PureFluid::Nitrogen,
            [
                (
                    65.98756999998032,
                    0.011524703911867041,
                    212507.83463074057,
                    855.4880578457547,
                    20579.954785127382,
                ),
                (
                    85.84355999984251,
                    0.006999594114270946,
                    187222.28591682157,
                    765.9975404860468,
                    248050.5166290526,
                ),
                (
                    108.53611999968503,
                    0.0024994064265441624,
                    138864.501653862,
                    632.5611183744085,
                    1345168.1757670508,
                ),
            ],
        ),
        (
            PureFluid::Argon,
            [
                (
                    86.7733325000005,
                    0.012664824162628414,
                    161536.19915306487,
                    1398.6520849974424,
                    95786.54158838271,
                ),
                (
                    107.54466000000403,
                    0.0077488385471208635,
                    142759.84745987057,
                    1260.8794966274966,
                    564472.025259583,
                ),
                (
                    131.28332000000805,
                    0.0028540535695923777,
                    106810.10037968487,
                    1054.364651608827,
                    2152107.400372201,
                ),
            ],
        ),
        (
            PureFluid::Helium,
            [
                (
                    2.314736750647708,
                    0.000275187223997737,
                    22884.63562573213,
                    145.7123307600881,
                    6949.802114358763,
                ),
                (
                    3.2802940051816614,
                    0.00018487364958343446,
                    23402.88504486816,
                    138.42926751714722,
                    35666.304586577426,
                ),
                (
                    4.383788010363323,
                    7.250667421703676e-05,
                    19518.169970732055,
                    121.28350417348425,
                    117358.9782596124,
                ),
            ],
        ),
        (
            PureFluid::CarbonDioxide,
            [
                (
                    220.20848950014158,
                    0.015644427386065217,
                    344570.7445432449,
                    1165.3791670922942,
                    604388.9684784702,
                ),
                (
                    245.52391600113268,
                    0.009972799937691063,
                    298699.8886238156,
                    1065.6134517573255,
                    1544987.3397560755,
                ),
                (
                    274.4558320022653,
                    0.0042477753340345115,
                    226894.3904065031,
                    919.5125251897678,
                    3607154.6647481667,
                ),
            ],
        ),
    ];

    #[test]
    fn test_condensation_properties_match_coolprop() -> anyhow::Result<()> {
        for (fluid, rows) in GOLDEN {
            for (t, sigma, latent, rho_l, p) in rows {
                let got = at(fluid, t)?;
                for (name, value, want) in [
                    ("σ", got[0], sigma),
                    ("latent heat", got[1], latent),
                    ("liquid density", got[2], rho_l),
                    ("p_sat", got[3], p),
                ] {
                    // CoolProp's N2 reducing density is 11183.9014646 mol/m³
                    // against the published 11183.9: 1.3e-7 in density.
                    let tolerance = if fluid == PureFluid::Nitrogen {
                        1e-6
                    } else {
                        1e-7
                    };
                    anyhow::ensure!(
                        (value / want - 1.0).abs() < tolerance,
                        "{fluid:?} at {t} K: {name} {value} vs {want}"
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_condensation_properties_continue_below_the_triple_point() -> anyhow::Result<()> {
        for fluid in [
            PureFluid::Nitrogen,
            PureFluid::Argon,
            PureFluid::CarbonDioxide,
        ] {
            let t_triple = fluid.data().data.t_triple;
            let (above, below) = (at(fluid, t_triple)?, at(fluid, t_triple * (1.0 - 1e-9))?);
            for k in 0..4 {
                anyhow::ensure!(
                    (below[k] / above[k] - 1.0).abs() < 1e-6,
                    "{fluid:?} property {k}: {} vs {}",
                    below[k],
                    above[k]
                );
            }
            // Colder: denser liquid, higher surface tension, lower vapour
            // pressure; the latent heat stays at its triple-point value.
            let cold = at(fluid, 0.5 * t_triple)?;
            anyhow::ensure!(
                cold[0] > above[0] && cold[2] > above[2] && cold[3] < above[3],
                "{fluid:?}: {cold:?} vs {above:?}"
            );
            anyhow::ensure!(
                (cold[1] / above[1] - 1.0).abs() < 1e-12,
                "{fluid:?} latent heat"
            );
        }
        Ok(())
    }

    #[test]
    fn test_condensation_properties_end_at_the_critical_point() -> anyhow::Result<()> {
        for fluid in PureFluid::ALL {
            let t_crit = fluid.data().data.t_crit;
            anyhow::ensure!(at(fluid, t_crit).is_err(), "{fluid:?} at T_c");
            anyhow::ensure!(
                at(fluid, 0.0).is_err() && at(fluid, f64::NAN).is_err(),
                "{fluid:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_molecular_mass_is_molar_mass_over_avogadro() {
        let m = PureFluid::Nitrogen.molecular_mass().get::<kilogram>();
        assert!(
            (m / (0.028_013_4 / 6.022_140_76e23) - 1.0).abs() < 1e-5,
            "{m}"
        );
    }

    #[test]
    fn test_vapour_transport_wraps_the_dilute_correlations() -> anyhow::Result<()> {
        use crate::units::{
            dynamic_viscosity::pascal_second, thermal_conductivity::watt_per_meter_kelvin,
        };
        let x =
            PureFluid::Nitrogen.vapour_transport(ThermodynamicTemperature::new::<kelvin>(300.0))?;
        let (mu, k) = dilute_transport(PureFluid::Nitrogen, 300.0);
        anyhow::ensure!(
            (x.viscosity.get::<pascal_second>() - mu).abs() < 1e-18
                && (x.thermal_conductivity.get::<watt_per_meter_kelvin>() - k).abs() < 1e-15,
            "{x:?}"
        );
        for t in [0.0, -1.0, f64::NAN] {
            anyhow::ensure!(
                PureFluid::Argon
                    .vapour_transport(ThermodynamicTemperature::new::<kelvin>(t))
                    .is_err(),
                "accepted {t} K"
            );
        }
        Ok(())
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(1000))]

        /// Along the saturation line the properties are smooth and physical:
        /// the surface tension and latent heat fall towards the critical point,
        /// the liquid thins, the vapour pressure rises.
        #[test]
        fn test_condensation_properties_vary_monotonically(fluid_index in 0usize..4, f in 0.0f64..0.94) {
            let fluid = PureFluid::ALL[fluid_index];
            let data = fluid.data().data;
            let t = data.t_triple + f * (data.t_crit - data.t_triple);
            let (a, b) = (at(fluid, t), at(fluid, t + 0.005 * (data.t_crit - data.t_triple)));
            let (a, b) = (a.map_err(|e| TestCaseError::fail(format!("{e:#}")))?, b.map_err(|e| TestCaseError::fail(format!("{e:#}")))?);
            prop_assert!(b[0] < a[0] && b[2] < a[2] && b[3] > a[3], "{:?} at {} K: {:?} then {:?}", fluid, t, a, b);
            // Helium's latent heat peaks inside the range.
            if fluid != PureFluid::Helium {
                prop_assert!(b[1] < a[1], "{:?} latent heat at {} K", fluid, t);
            }
        }
    }
}
