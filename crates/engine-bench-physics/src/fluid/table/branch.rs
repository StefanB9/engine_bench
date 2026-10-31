//! Single-phase states on one branch of the reference equation, for table
//! nodes: the vapour branch (frozen vapour, metastable inside the dome) and
//! the liquid branch (superheated liquid inside the dome). Above the critical
//! temperature both are the one fluid.

use super::super::{
    Fluid,
    frozen::{NoVapour, lowest, vapour_root},
    helmholtz::molar,
    saturation::ancillary,
    state::{density_si, max_density},
};
use crate::{
    PhysicsError,
    numerics::roots::{RootOptions, newton_bracketed},
};

/// Density step while bracketing the superheated liquid root.
const SHRINK: f64 = 0.995;
/// Warm-started Newton iterations before falling back to a bracketed solve.
const WARM_ITERATIONS: usize = 12;

/// Which branch of p(ρ) a state is on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Branch {
    Vapour,
    Liquid,
}

/// A table node: ln T, ln ρ (kg/m³), h (J/kg), ln w (m/s).
pub(super) type Node = [f64; 4];

/// Every quantity invalid.
pub(super) const INVALID: Node = [f64::NAN; 4];

/// Why a branch has no state at some node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum End {
    /// Its spinodal, where the fields turn singular (w → 0).
    Spinodal,
    /// Its temperature or trusted range.
    Range,
}

impl Branch {
    /// Temperatures the branch covers, K.
    fn range(self, fluid: &Fluid) -> (f64, f64) {
        match self {
            Self::Vapour => (lowest(fluid), fluid.data.t_max),
            Self::Liquid => (fluid.data.t_triple, fluid.data.t_max),
        }
    }

    /// Density, mol/m³, on this branch at `p` Pa and `t` K.
    pub(super) fn density(self, fluid: &Fluid, p: f64, t: f64) -> Result<f64, End> {
        match self {
            Self::Vapour => vapour_root(fluid, p, t).map_err(|e| match e {
                NoVapour::Spinodal => End::Spinodal,
                NoVapour::Untrusted => End::Range,
            }),
            Self::Liquid => liquid_root(fluid, p, t),
        }
    }

    /// Residual sign where no state of this branch exists at T: vapour needs a
    /// higher temperature, liquid a lower one.
    fn missing(self) -> f64 {
        match self {
            Self::Vapour => -f64::MAX,
            Self::Liquid => f64::MAX,
        }
    }

    /// The node at `p` Pa and molar entropy `s`, its temperature found from
    /// `guess` K if given, else by a bracketed search; where no state of this
    /// branch has that entropy, how the branch ends.
    #[expect(
        clippy::many_single_char_names,
        reason = "the thermodynamic symbols p, s, T, ρ, x"
    )]
    pub(super) fn node(
        self,
        fluid: &Fluid,
        p: f64,
        s: f64,
        guess: Option<f64>,
    ) -> Result<Node, End> {
        let t = guess
            .and_then(|t| self.warm(fluid, p, s, t))
            .or_else(|| self.bracketed(fluid, p, s))
            .ok_or(End::Range)?;
        let rho = self.density(fluid, p, t)?;
        let x = molar(fluid, rho, t);
        let m = fluid.data.molar_mass;
        let close = (x.s - s).abs() <= 1e-10 * s.abs().max(fluid.data.gas_constant);
        if close && x.w > 0.0 {
            return Ok([t.ln(), (rho * m).ln(), x.h / m, x.w.ln()]);
        }
        // The search stopped at the branch's edge: why it ends just beyond.
        let beyond = match self {
            Self::Vapour => t * (1.0 - 1e-6),
            Self::Liquid => t * (1.0 + 1e-6),
        };
        Err(self.density(fluid, p, beyond).err().unwrap_or(End::Range))
    }

    /// Newton on ln T from `t` (∂s/∂ln T = cp).
    fn warm(self, fluid: &Fluid, p: f64, s: f64, mut t: f64) -> Option<f64> {
        let (lo, hi) = self.range(fluid);
        for _ in 0..WARM_ITERATIONS {
            let rho = self.density(fluid, p, t).ok()?;
            let x = molar(fluid, rho, t);
            let step = ((x.s - s) / x.cp).clamp(-0.2, 0.2);
            let next = (t * (-step).exp()).clamp(lo, hi);
            if (next / t - 1.0).abs() < 1e-13 {
                return Some(next);
            }
            t = next;
        }
        None
    }

    /// Bracketed Newton on ln T over the branch's range.
    fn bracketed(self, fluid: &Fluid, p: f64, s: f64) -> Option<f64> {
        let (lo, hi) = self.range(fluid);
        let missing = self.missing();
        newton_bracketed(
            |ln_t| {
                let t = ln_t.exp().clamp(lo, hi);
                match self.density(fluid, p, t) {
                    Ok(rho) => {
                        let x = molar(fluid, rho, t);
                        (x.s - s, x.cp)
                    }
                    Err(_) => (missing, f64::NAN),
                }
            },
            lo.ln(),
            hi.ln(),
            f64::midpoint(lo.ln(), hi.ln()),
            RootOptions::new(1e-13, 200).ok()?,
        )
        .ok()
        .map(|ln_t| ln_t.exp().clamp(lo, hi))
    }
}

/// Density, mol/m³, of the liquid at `p` Pa and `t` K, superheated
/// (metastable) inside the dome.
///
/// Under the critical temperature p(ρ) falls along the liquid branch from the
/// saturated liquid down to the spinodal, its slope falling to zero there; the
/// root is bracketed by stepping down from compressed liquid while the slope
/// stays positive and falling. Above it the one fluid.
pub(super) fn liquid_density(fluid: &Fluid, p: f64, t: f64) -> anyhow::Result<f64> {
    liquid_root(fluid, p, t).map_err(|end| {
        PhysicsError::InvalidState(match end {
            End::Spinodal => "beyond the liquid spinodal",
            End::Range => "outside the liquid's range",
        })
        .into()
    })
}

/// [`liquid_density`] with how the branch ends when there is no root.
fn liquid_root(fluid: &Fluid, p: f64, t: f64) -> Result<f64, End> {
    let data = fluid.data;
    if t >= data.t_crit {
        return density_si(fluid, p, t).map_err(|_| End::Range);
    }
    if !(t >= data.t_triple && p > 0.0 && p <= data.p_max) {
        return Err(End::Range);
    }
    let f = |rho: f64| molar(fluid, rho, t);
    // Denser than saturated liquid (the ancillary is good to ~1 %).
    let start = (1.01 * ancillary(&data.liquid_density, t)).max(1.001 * data.rho_crit);
    let (lo, hi) = if f(start).p < p {
        // Compressed liquid: p rises monotonically with ρ.
        let mut hi = start;
        while f(hi).p < p {
            if hi >= max_density(fluid) {
                return Err(End::Range);
            }
            hi = (hi * 1.05).min(max_density(fluid));
        }
        (start, hi)
    } else {
        let (mut lo, mut hi) = (start, start);
        let mut slope = f64::INFINITY;
        loop {
            let x = f(lo);
            if x.dp_drho <= 0.0 || x.dp_drho > slope || lo <= data.rho_crit {
                return Err(End::Spinodal);
            }
            if x.p <= p {
                break;
            }
            slope = x.dp_drho;
            hi = lo;
            lo *= SHRINK;
        }
        (lo, hi)
    };
    let options = RootOptions::new(1e-14 * hi, 200).map_err(|_| End::Range)?;
    newton_bracketed(
        |rho| {
            let x = f(rho);
            (x.p - p, x.dp_drho)
        },
        lo,
        hi,
        hi,
        options,
    )
    .map_err(|_| End::Range)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid::{PureFluid, saturation_si};

    #[test]
    fn test_liquid_density_matches_the_stable_liquid() -> anyhow::Result<()> {
        let fluid = PureFluid::Nitrogen.data();
        // Compressed liquid at 80 K, 50 bar.
        let (stable, branch) = (
            density_si(fluid, 5e6, 80.0)?,
            liquid_density(fluid, 5e6, 80.0)?,
        );
        anyhow::ensure!(
            (branch / stable - 1.0).abs() < 1e-12,
            "{branch} vs {stable}"
        );
        // Supercritical, both branches the one fluid.
        let (stable, branch) = (
            density_si(fluid, 5e6, 200.0)?,
            liquid_density(fluid, 5e6, 200.0)?,
        );
        anyhow::ensure!(
            (branch / stable - 1.0).abs() < 1e-12,
            "{branch} vs {stable}"
        );
        Ok(())
    }

    #[test]
    fn test_liquid_density_superheats_inside_the_dome() -> anyhow::Result<()> {
        let fluid = PureFluid::Nitrogen.data();
        let (p_sat, rho_l, _) = saturation_si(fluid, 100.0)?;
        let rho = liquid_density(fluid, 0.8 * p_sat, 100.0)?;
        anyhow::ensure!(rho < rho_l && rho > fluid.data.rho_crit, "{rho} vs {rho_l}");
        let x = molar(fluid, rho, 100.0);
        anyhow::ensure!(
            (x.p / (0.8 * p_sat) - 1.0).abs() < 1e-10 && x.dp_drho > 0.0,
            "{x:?}"
        );
        // Far below saturation the liquid has passed its spinodal.
        anyhow::ensure!(liquid_density(fluid, 1e2, 120.0).is_err(), "120 K, 100 Pa");
        Ok(())
    }

    #[test]
    fn test_branch_nodes_match_the_exact_states() -> anyhow::Result<()> {
        let fluid = PureFluid::Nitrogen.data();
        let m = fluid.data.molar_mass;
        for (branch, p, t) in [
            (Branch::Vapour, 1e5, 250.0),
            (Branch::Liquid, 5e6, 80.0),
            (Branch::Vapour, 2e7, 300.0),
        ] {
            let rho = branch
                .density(fluid, p, t)
                .map_err(|e| anyhow::anyhow!("{branch:?} at {p} Pa, {t} K: {e:?}"))?;
            let x = molar(fluid, rho, t);
            for guess in [None, Some(0.9 * t)] {
                let node = branch
                    .node(fluid, p, x.s, guess)
                    .map_err(|e| anyhow::anyhow!("{branch:?} at {p} Pa, {t} K: {e:?}"))?;
                let want = [t.ln(), (rho * m).ln(), x.h / m, x.w.ln()];
                for (a, b) in node.iter().zip(want) {
                    anyhow::ensure!(
                        (a - b).abs() < 1e-9 * b.abs().max(1.0),
                        "{branch:?} {p} Pa {t} K: {node:?} vs {want:?}"
                    );
                }
            }
        }
        // No vapour has a liquid's entropy at low pressure.
        let liquid = molar(fluid, liquid_density(fluid, 5e6, 80.0)?, 80.0);
        anyhow::ensure!(
            Branch::Vapour.node(fluid, 1e4, liquid.s, None).is_err(),
            "vapour"
        );
        anyhow::ensure!(INVALID.iter().all(|v| v.is_nan()), "invalid node");
        Ok(())
    }
}
