//! Dilute-gas viscosity and thermal conductivity of water vapour: the
//! zero-density terms of IAPWS R12-08 (viscosity) and R15-11 (conductivity).

/// Critical temperature of the IAPWS transport releases, K.
const T_STAR: f64 = 647.096;

/// Dilute-gas viscosity, Pa·s, and thermal conductivity, W/(m·K), at `t` K.
pub(crate) fn dilute_transport(t: f64) -> (f64, f64) {
    const H: [f64; 4] = [1.677_52, 2.204_62, 0.636_656_4, -0.241_605];
    const L: [f64; 5] = [
        2.443_221e-3,
        1.323_095e-2,
        6.770_357e-3,
        -3.454_586e-3,
        4.096_266e-4,
    ];
    let tbar = t / T_STAR;
    let inv = 1.0 / tbar;
    let series = |c: &[f64]| c.iter().rev().fold(0.0_f64, |acc, x| acc.mul_add(inv, *x));
    (
        100e-6 * tbar.sqrt() / series(&H),
        1e-3 * tbar.sqrt() / series(&L),
    )
}

#[cfg(test)]
#[expect(
    clippy::unreadable_literal,
    reason = "golden values exactly as tools/reference/transport_golden.py prints them"
)]
mod tests {
    use super::*;

    #[test]
    fn test_water_dilute_transport_matches_coolprop() -> anyhow::Result<()> {
        // tools/reference/transport_golden.py (CoolProp 8.0, 1 µmol/m³).
        for (t, mu, k) in [
            (280.0, 9.147960470716126e-06, 0.017208632759857068),
            (373.15, 1.2337030686429913e-05, 0.024155842561257836),
            (500.0, 1.732628304228029e-05, 0.035780357253938355),
        ] {
            let (got_mu, got_k) = dilute_transport(t);
            anyhow::ensure!((got_mu / mu - 1.0).abs() < 1e-9, "{t} K viscosity {got_mu}");
            anyhow::ensure!((got_k / k - 1.0).abs() < 1e-9, "{t} K conductivity {got_k}");
        }
        Ok(())
    }
}
