//! Dilute-gas viscosity and thermal conductivity of the reference fluids, the
//! zero-density terms of the correlations CoolProp uses:
//!
//! - N₂, Ar: Lemmon and Jacobsen, Int. J. Thermophys. 25 (2004) 21: viscosity
//!   from a collision integral, conductivity from it plus a power series.
//! - He: Arp, McCarty and Friend, NIST TN 1334 (1998) viscosity; Hands and Arp,
//!   Cryogenics 21 (1981) 697 conductivity.
//! - CO₂: Laesecke and Muzny, JPCRD 46 (2017) 013107 viscosity; Huber et al.,
//!   JPCRD 45 (2016) 013102 conductivity.

use super::PureFluid;

/// Lemmon–Jacobsen collision-integral coefficients bᵢ of ln Ω = Σ bᵢ (ln T*)ⁱ.
const COLLISION: [f64; 5] = [0.431, -0.4623, 0.084_06, 0.005_341, -0.003_31];
/// Prefactor of η⁰ = C √(M T) / (σ² Ω), M in kg/kmol, σ in nm, η in Pa·s.
const LEMMON_C: f64 = 2.669_58e-8;

/// Lemmon–Jacobsen parameters: molar mass kg/kmol, σ nm, ε/k K, critical
/// temperature K and the conductivity terms (N₁; (Nᵢ, tᵢ) of Nᵢ τ^tᵢ, W/(m·K)).
struct Lemmon {
    molar_mass: f64,
    sigma: f64,
    epsilon: f64,
    t_crit: f64,
    n1: f64,
    terms: &'static [(f64, f64)],
}

const NITROGEN: Lemmon = Lemmon {
    molar_mass: 28.013_48,
    sigma: 0.3656,
    epsilon: 98.94,
    t_crit: 126.192,
    n1: 1.511e-3,
    terms: &[(2.117e-3, -1.0), (-3.332e-3, -0.7)],
};

const ARGON: Lemmon = Lemmon {
    molar_mass: 39.948,
    sigma: 0.335,
    epsilon: 143.2,
    t_crit: 150.687,
    n1: 8.158e-4,
    terms: &[(-4.32e-4, -0.77)],
};

impl Lemmon {
    fn viscosity(&self, t: f64) -> f64 {
        let ln_t = (t / self.epsilon).ln();
        let ln_omega = COLLISION
            .iter()
            .rev()
            .fold(0.0_f64, |acc, b| acc.mul_add(ln_t, *b));
        LEMMON_C * (self.molar_mass * t).sqrt() / (self.sigma * self.sigma * ln_omega.exp())
    }

    fn conductivity(&self, t: f64, viscosity: f64) -> f64 {
        let tau = self.t_crit / t;
        self.terms
            .iter()
            .fold(self.n1 * viscosity * 1e6, |acc, &(n, e)| {
                acc + n * tau.powf(e)
            })
    }
}

/// Helium viscosity, Pa·s (Arp et al.: zero-density terms in x = ln T, T
/// capped at 300 K in x, the high-temperature form above 100 K).
fn helium_viscosity(t: f64) -> f64 {
    let x = t.min(300.0).ln();
    let ln_eta = -0.135_311_743 / x + 1.003_478_41 + 1.206_546_49 * x - 0.149_564_551 * x * x
        + 0.012_520_841 * x * x * x;
    let micro_poise = if t <= 100.0 {
        ln_eta.exp()
    } else {
        196.0 * t.powf(0.719_38) * (12.451 / t - 295.67 / (t * t) - 4.1249).exp()
    };
    micro_poise * 1e-7
}

/// Helium conductivity, W/(m·K) (Hands and Arp's dilute term).
fn helium_conductivity(t: f64) -> f64 {
    let inv = 1.0 / t;
    let sum = inv
        * (3.739_232_544 + inv * (-26.203_169_69 + inv * (59.822_522_46 - inv * 49.263_976_34)));
    2.787_003_4e-3 * t.powf(0.703_400_705_7) * sum.exp()
}

/// CO₂ viscosity, Pa·s (Laesecke and Muzny, eq. 4).
fn carbon_dioxide_viscosity(t: f64) -> f64 {
    const A: [f64; 7] = [
        1_749.354_893_188_35,
        -369.069_300_007_128,
        5_423_856.348_876_91,
        -2.212_838_521_683_56,
        -269_503.247_933_569,
        73_145.021_531_826,
        5.343_686_495_092_78,
    ];
    let (cube, sixth) = (t.cbrt(), t.powf(1.0 / 6.0));
    let den = A[0]
        + A[1] * sixth
        + A[2] * (A[3] * cube).exp()
        + (A[4] + A[5] * cube) / cube.exp()
        + A[6] * t.sqrt();
    1.0055e-3 * t.sqrt() / den
}

/// CO₂ conductivity, W/(m·K) (Huber et al., eq. 3).
fn carbon_dioxide_conductivity(t: f64) -> f64 {
    const L: [f64; 4] = [
        0.015_187_430_7,
        0.028_067_404,
        0.022_856_419,
        -0.007_416_242_1,
    ];
    let tau = 304.1282 / t;
    let den = L.iter().rev().fold(0.0_f64, |acc, l| acc.mul_add(tau, *l));
    1e-3 / (tau.sqrt() * den)
}

/// Lennard-Jones σ, m, of the viscosity correlations (Lemmon and Jacobsen
/// for N₂ and Ar, Laesecke and Muzny for CO₂; for He Hirschfelder, Curtiss
/// and Bird's 2.576 Å).
pub(crate) fn lennard_jones_diameter(fluid: PureFluid) -> f64 {
    match fluid {
        PureFluid::Nitrogen => NITROGEN.sigma * 1e-9,
        PureFluid::Argon => ARGON.sigma * 1e-9,
        PureFluid::Helium => 2.576e-10,
        PureFluid::CarbonDioxide => 3.784_21e-10,
    }
}

/// Dilute-gas viscosity, Pa·s, and thermal conductivity, W/(m·K), of
/// `fluid` at `t` K.
pub(crate) fn dilute_transport(fluid: PureFluid, t: f64) -> (f64, f64) {
    match fluid {
        PureFluid::Nitrogen | PureFluid::Argon => {
            let law = if fluid == PureFluid::Nitrogen {
                &NITROGEN
            } else {
                &ARGON
            };
            let mu = law.viscosity(t);
            (mu, law.conductivity(t, mu))
        }
        PureFluid::Helium => (helium_viscosity(t), helium_conductivity(t)),
        PureFluid::CarbonDioxide => (carbon_dioxide_viscosity(t), carbon_dioxide_conductivity(t)),
    }
}

#[cfg(test)]
#[expect(
    clippy::unreadable_literal,
    reason = "golden values exactly as tools/reference/transport_golden.py prints them"
)]
mod tests {
    use super::*;

    /// One fluid's rows of (T K, viscosity Pa·s, conductivity W/(m·K)).
    type Golden = (PureFluid, &'static [(f64, f64, f64)]);

    /// The rows from
    /// `tools/reference/transport_golden.py` (CoolProp 8.0, 1 µmol/m³).
    const GOLDEN: [Golden; 4] = [
        (
            PureFluid::Nitrogen,
            &[
                (70.0, 4.850163495274111e-06, 0.006297193817972508),
                (100.0, 6.903488327132431e-06, 0.009277493929510156),
                (200.0, 1.2888865785908165e-05, 0.018230854150803777),
                (300.0, 1.7877064146328245e-05, 0.025936086672010222),
            ],
        ),
        (
            PureFluid::Argon,
            &[
                (90.0, 7.3488579707189445e-06, 0.005704707706294764),
                (150.0, 1.2223660721864426e-05, 0.009541579759324978),
                (300.0, 2.2724103711893546e-05, 0.01780423473967644),
            ],
        ),
        (
            PureFluid::Helium,
            &[
                (3.0, 7.70542496432322e-07, 0.00569808935611798),
                (5.0, 1.2504002840428083e-06, 0.009549838856351476),
                (20.0, 3.5416174406973773e-06, 0.026070144799544337),
                (100.0, 9.764967249693512e-06, 0.07363189335523569),
            ],
        ),
        (
            PureFluid::CarbonDioxide,
            &[
                (220.0, 1.1117727298291105e-05, 0.010893607966549142),
                (300.0, 1.4993786441644809e-05, 0.016719181705923833),
                (500.0, 2.3910303022709982e-05, 0.03283743403011052),
            ],
        ),
    ];

    #[test]
    fn test_dilute_transport_matches_coolprop() -> anyhow::Result<()> {
        // Helium's critical-enhancement term (3.5 K to 12 K) still adds ~1e-8
        // at the goldens' 1 µmol/m³: it vanishes only as √ρ.
        for (fluid, rows) in GOLDEN {
            for &(t, mu, k) in rows {
                let (got_mu, got_k) = dilute_transport(fluid, t);
                anyhow::ensure!(
                    (got_mu / mu - 1.0).abs() < 1e-9,
                    "{fluid:?} {t} K viscosity {got_mu} vs {mu}"
                );
                anyhow::ensure!(
                    (got_k / k - 1.0).abs() < 1e-7,
                    "{fluid:?} {t} K conductivity {got_k} vs {k}"
                );
            }
        }
        Ok(())
    }
}
