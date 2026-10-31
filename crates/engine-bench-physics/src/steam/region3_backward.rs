//! Region 3 backward equations v(p, T) of IAPWS SR5-05(2016): 26 subregions
//! (3a … 3z) separated by boundary curves T(p).
//!
//! They agree with the Helmholtz equation to about 1e-6 away from the critical
//! point, so [`super::region3`] uses them only as starting points for a Newton
//! solve on p(ρ, T). Infallible kernels on raw SI values.

use super::{
    coefficients::{
        Subregion,
        region3::{
            LINE_AB, LINE_CD, LINE_GH, LINE_IJ, LINE_JK, LINE_MN, LINE_OP, LINE_QU, LINE_RX,
            LINE_UV, LINE_WX, SUBREGIONS,
        },
    },
    saturation::saturation_temperature_si,
};

/// Which side of the saturation line a state on it belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Side {
    /// Not on the saturation line: pick the subregion from (p, T) alone.
    Either,
    /// Saturated liquid.
    Liquid,
    /// Saturated vapour.
    Vapour,
}

/// Estimated specific volume, m³/kg, at p Pa and T K in region 3, or `None`
/// where the subregion boundaries leave a gap (only in slivers next to the
/// critical point).
pub(crate) fn specific_volume_estimate(p: f64, t: f64, side: Side) -> Option<f64> {
    let letter = adjust_for_side(side, p, subregion(p, t)?);
    let v = evaluate(&SUBREGIONS[usize::from(letter - b'a')], p, t);
    v.is_finite().then_some(v)
}

/// v(p, T) of one subregion.
fn evaluate(sub: &Subregion, p: f64, t: f64) -> f64 {
    let x = p / sub.p_star - sub.a;
    let y = t / sub.t_star - sub.b;
    if sub.exponential {
        let sum: f64 = sub
            .terms
            .iter()
            .map(|term| term.n * x.powi(term.i) * y.powi(term.j))
            .sum();
        return sum.exp() * sub.v_star;
    }
    let (xc, yd) = (x.powf(sub.c), y.powf(sub.d));
    let sum: f64 = sub
        .terms
        .iter()
        .map(|term| term.n * xc.powi(term.i) * yd.powi(term.j))
        .sum();
    sum.powf(sub.e) * sub.v_star
}

/// Σ n x^I.
fn line(terms: &[(i32, f64)], x: f64) -> f64 {
    terms.iter().map(|&(i, n)| n * x.powi(i)).sum()
}

/// Boundary temperatures, K, at p Pa (SR5-05 Eqs. 2–4).
fn t_ab(p: f64) -> f64 {
    line(&LINE_AB, (p * 1e-6).ln())
}
fn t_cd(p: f64) -> f64 {
    line(&LINE_CD, p * 1e-6)
}
fn t_ef(p: f64) -> f64 {
    3.727_888_004f64.mul_add(p * 1e-6 - 22.064, 647.096)
}
fn t_gh(p: f64) -> f64 {
    line(&LINE_GH, p * 1e-6)
}
fn t_ij(p: f64) -> f64 {
    line(&LINE_IJ, p * 1e-6)
}
fn t_jk(p: f64) -> f64 {
    line(&LINE_JK, p * 1e-6)
}
fn t_mn(p: f64) -> f64 {
    line(&LINE_MN, p * 1e-6)
}
fn t_op(p: f64) -> f64 {
    line(&LINE_OP, (p * 1e-6).ln())
}
fn t_qu(p: f64) -> f64 {
    line(&LINE_QU, p * 1e-6)
}
fn t_rx(p: f64) -> f64 {
    line(&LINE_RX, p * 1e-6)
}
fn t_uv(p: f64) -> f64 {
    line(&LINE_UV, p * 1e-6)
}
fn t_wx(p: f64) -> f64 {
    line(&LINE_WX, (p * 1e-6).ln())
}

/// Subregion letter of (p, T) (SR5-05 Table 2, Section 3.4).
fn subregion(p: f64, t: f64) -> Option<u8> {
    let mpa = p * 1e-6;
    let letter = if mpa > 40.0 {
        if t <= t_ab(p) { b'a' } else { b'b' }
    } else if mpa > 25.0 {
        if t <= t_cd(p) {
            b'c'
        } else if t <= t_ab(p) {
            b'd'
        } else if t <= t_ef(p) {
            b'e'
        } else {
            b'f'
        }
    } else if mpa > 23.5 {
        if t <= t_cd(p) {
            b'c'
        } else if t <= t_gh(p) {
            b'g'
        } else if t <= t_ef(p) {
            b'h'
        } else if t <= t_ij(p) {
            b'i'
        } else if t <= t_jk(p) {
            b'j'
        } else {
            b'k'
        }
    } else if mpa > 23.0 {
        if t <= t_cd(p) {
            b'c'
        } else if t <= t_gh(p) {
            b'l'
        } else if t <= t_ef(p) {
            b'h'
        } else if t <= t_ij(p) {
            b'i'
        } else if t <= t_jk(p) {
            b'j'
        } else {
            b'k'
        }
    } else if mpa > 22.5 {
        if t <= t_cd(p) {
            b'c'
        } else if t <= t_gh(p) {
            b'l'
        } else if t <= t_mn(p) {
            b'm'
        } else if t <= t_ef(p) {
            b'n'
        } else if t <= t_op(p) {
            b'o'
        } else if t <= t_ij(p) {
            b'p'
        } else if t <= t_jk(p) {
            b'j'
        } else {
            b'k'
        }
    } else if mpa > 21.043_367_32 {
        if t <= t_cd(p) {
            b'c'
        } else if t <= t_qu(p) {
            b'q'
        } else if t > t_jk(p) {
            b'k'
        } else if t > t_rx(p) {
            b'r'
        } else {
            return near_critical(p, t);
        }
    } else if mpa > 20.5 {
        if t <= t_cd(p) {
            b'c'
        } else if t <= saturation_temperature_si(p) {
            b's'
        } else if t <= t_jk(p) {
            b'r'
        } else {
            b'k'
        }
    } else if mpa > 19.008_811_891_739_29 {
        if t <= t_cd(p) {
            b'c'
        } else if t <= saturation_temperature_si(p) {
            b's'
        } else {
            b't'
        }
    } else if t <= saturation_temperature_si(p) {
        b'c'
    } else {
        b't'
    };
    Some(letter)
}

/// Subregions 3u … 3z between the qu and rx boundaries, 21.04–22.5 MPa
/// (SR5-05 Section 4.3).
fn near_critical(p: f64, t: f64) -> Option<u8> {
    let mpa = p * 1e-6;
    let between = |lo: f64, hi: f64| lo < t && t <= hi;
    if mpa > 22.11 {
        if between(t_qu(p), t_uv(p)) {
            Some(b'u')
        } else if between(t_uv(p), t_ef(p)) {
            Some(b'v')
        } else if between(t_ef(p), t_wx(p)) {
            Some(b'w')
        } else if between(t_wx(p), t_rx(p)) {
            Some(b'x')
        } else {
            None
        }
    } else if mpa > 22.064 {
        if between(t_qu(p), t_uv(p)) {
            Some(b'u')
        } else if between(t_uv(p), t_ef(p)) {
            Some(b'y')
        } else if between(t_ef(p), t_wx(p)) {
            Some(b'z')
        } else if between(t_wx(p), t_rx(p)) {
            Some(b'x')
        } else {
            None
        }
    } else if t <= saturation_temperature_si(p) {
        if mpa > 21.931_615_51 {
            if between(t_qu(p), t_uv(p)) {
                Some(b'u')
            } else if t > t_uv(p) {
                Some(b'y')
            } else {
                None
            }
        } else {
            Some(b'u')
        }
    } else if mpa > 21.900_962_65 {
        if t <= t_wx(p) {
            Some(b'z')
        } else if t <= t_rx(p) {
            Some(b'x')
        } else {
            None
        }
    } else {
        Some(b'x')
    }
}

/// On the saturation line (p, T) alone cannot tell liquid from vapour; move
/// to the subregion on the requested side (as CoolProp/IF97 does).
fn adjust_for_side(side: Side, p: f64, letter: u8) -> u8 {
    let mpa = p * 1e-6;
    match (side, letter) {
        (Side::Vapour, b'c') => b't',
        (Side::Vapour, b's') => {
            if mpa < 20.5 {
                b't'
            } else {
                b'r'
            }
        }
        (Side::Vapour, b'u') => {
            if mpa < 21.900_962_65 {
                b'x'
            } else {
                b'z'
            }
        }
        (Side::Vapour, b'y') => b'z',
        (Side::Liquid, b'z') => {
            if mpa > 21.931_615_51 {
                b'y'
            } else {
                b'u'
            }
        }
        (Side::Liquid, b'x') => b'u',
        (Side::Liquid, b'r' | b'k') => b's',
        (Side::Liquid, b't') => {
            if mpa > 19.008_811_891_739_29 {
                b's'
            } else {
                b'c'
            }
        }
        _ => letter,
    }
}

#[cfg(test)]
#[expect(
    clippy::unreadable_literal,
    reason = "verification values exactly as IAPWS SR5-05(2016) prints them"
)]
mod tests {
    use super::*;

    /// SR5-05(2016) Tables 5 and 13 as reproduced in CoolProp/IF97: p MPa,
    /// T K, v m³/kg; two points per subregion 3a … 3z (3a: one).
    const TABLE: [(f64, f64, f64); 51] = [
        (80.0, 670.0, 1.503831359e-3),
        (50.0, 710.0, 2.204728587e-3),
        (80.0, 750.0, 1.973692940e-3),
        (20.0, 630.0, 1.761696406e-3),
        (30.0, 650.0, 1.819560617e-3),
        (26.0, 656.0, 2.245587720e-3),
        (30.0, 670.0, 2.506897702e-3),
        (26.0, 661.0, 2.970225962e-3),
        (30.0, 675.0, 3.004627086e-3),
        (26.0, 671.0, 5.019029401e-3),
        (30.0, 690.0, 4.656470142e-3),
        (23.6, 649.0, 2.163198378e-3),
        (24.0, 650.0, 2.166044161e-3),
        (23.6, 652.0, 2.651081407e-3),
        (24.0, 654.0, 2.967802335e-3),
        (23.6, 653.0, 3.273916816e-3),
        (24.0, 655.0, 3.550329864e-3),
        (23.5, 655.0, 4.545001142e-3),
        (24.0, 660.0, 5.100267704e-3),
        (23.0, 660.0, 6.109525997e-3),
        (24.0, 670.0, 6.427325645e-3),
        (22.6, 646.0, 2.117860851e-3),
        (23.0, 646.0, 2.062374674e-3),
        (22.6, 648.6, 2.533063780e-3),
        (22.8, 649.3, 2.572971781e-3),
        (22.6, 649.0, 2.923432711e-3),
        (22.8, 649.7, 2.913311494e-3),
        (22.6, 649.1, 3.131208996e-3),
        (22.8, 649.9, 3.221160278e-3),
        (22.6, 649.4, 3.715596186e-3),
        (22.8, 650.2, 3.664754790e-3),
        (21.1, 640.0, 1.970999272e-3),
        (21.8, 643.0, 2.043919161e-3),
        (21.1, 644.0, 5.251009921e-3),
        (21.8, 648.0, 5.256844741e-3),
        (19.1, 635.0, 1.932829079e-3),
        (20.0, 638.0, 1.985387227e-3),
        (17.0, 626.0, 8.483262001e-3),
        (20.0, 640.0, 6.227528101e-3),
        (21.5, 644.6, 2.268366647e-3),
        (22.0, 646.1, 2.296350553e-3),
        (22.5, 648.6, 2.832373260e-3),
        (22.3, 647.9, 2.811424405e-3),
        (22.15, 647.5, 3.694032281e-3),
        (22.3, 648.1, 3.622226305e-3),
        (22.11, 648.0, 4.528072649e-3),
        (22.3, 649.0, 4.556905799e-3),
        (22.0, 646.84, 2.698354719e-3),
        (22.064, 647.05, 2.717655648e-3),
        (22.0, 646.89, 3.798732962e-3),
        (22.064, 647.15, 3.701940010e-3),
    ];

    #[test]
    fn test_region3_backward_matches_iapws_sr5_tables_5_and_13() -> anyhow::Result<()> {
        for (p_mpa, t, v) in TABLE {
            let got = specific_volume_estimate(p_mpa * 1e6, t, Side::Either)
                .ok_or_else(|| anyhow::anyhow!("no subregion at {p_mpa} MPa, {t} K"))?;
            let error = (got / v - 1.0).abs();
            anyhow::ensure!(error < 1e-9, "{p_mpa} MPa, {t} K: {got} vs {v} ({error:e})");
        }
        Ok(())
    }

    #[test]
    fn test_region3_backward_puts_saturated_states_on_their_side() -> anyhow::Result<()> {
        // Just below the critical point the liquid is denser than the vapour at
        // the same (p_sat, T_sat), whichever subregion (p, T) alone would pick.
        for p in [17e6, 20e6, 21e6, 21.95e6, 22.05e6] {
            let t = saturation_temperature_si(p);
            let liquid = specific_volume_estimate(p, t, Side::Liquid)
                .ok_or_else(|| anyhow::anyhow!("no liquid estimate at {p} Pa"))?;
            let vapour = specific_volume_estimate(p, t, Side::Vapour)
                .ok_or_else(|| anyhow::anyhow!("no vapour estimate at {p} Pa"))?;
            anyhow::ensure!(
                liquid < 1.0 / 322.0 && vapour > 1.0 / 322.0,
                "{p} Pa: v' = {liquid}, v'' = {vapour}"
            );
        }
        Ok(())
    }
}
