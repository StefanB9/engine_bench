//! Building one fluid's tables over an envelope of pressure and temperature.

use anyhow::Context;

use super::{
    Tables,
    branch::{Branch, End, INVALID, Node, liquid_density},
    grid::{Axis, Grid2},
    saturation::SaturationTable,
};
use crate::fluid::{
    Fluid, PureFluid,
    frozen::vapour_density,
    helmholtz::molar,
    saturation::{dome_top, saturation_si},
};

/// Node spacing in ln p.
const LN_P_STEP: f64 = 0.1;
/// Node spacing in s, as a fraction of the ideal-gas cp at 300 K: ln T then
/// changes by about this much from node to node.
const S_STEP: f64 = 0.02;
/// Nodes beyond the envelope on every side, so that its edges have full
/// stencils.
const MARGIN: u32 = 2;
/// Nodes dropped where a branch ends at its spinodal: the fields turn
/// singular there (w → 0), so states that close come from the exact equation.
const SPINODAL_GUARD: usize = 3;

/// The states a table covers: pressures in Pa and temperatures in K.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Envelope {
    pub(super) p_min: f64,
    pub(super) p_max: f64,
    pub(super) t_min: f64,
    pub(super) t_max: f64,
}

impl Envelope {
    /// Cold-gas propulsion: vacuum nozzle exits (and a valve outlet draining
    /// to vacuum) to 500 bar tanks, from the
    /// lowest liquid or frozen vapour up to 900 K.
    pub(super) const DEFAULT: Self = Self {
        p_min: 0.1,
        p_max: 5e7,
        t_min: 0.0,
        t_max: 900.0,
    };
}

/// Builds `fluid`'s tables over `envelope`, rows in parallel.
pub(super) fn tables(fluid: PureFluid, envelope: &Envelope) -> anyhow::Result<Tables> {
    let data = fluid.data();
    let m = data.data.molar_mass;
    let p_max = envelope.p_max.min(data.data.p_max);
    let t_low = envelope.t_min.max(data.data.t_triple);
    let s_min = molar(data, liquid_density(data, p_max, t_low)?, t_low).s / m;
    let t_high = envelope.t_max.min(data.data.t_max);
    let s_max = molar(data, vapour_density(data, envelope.p_min, t_high)?, t_high).s / m;
    let cp = molar(data, 1e-6, 300.0).cp / m;
    let x = axis(envelope.p_min.ln(), p_max.ln(), LN_P_STEP);
    let y = axis(s_min, s_max, S_STEP * cp);
    let (p_triple, _, _) = saturation_si(data, data.data.t_triple)?;
    let (p_top, _, _) = saturation_si(data, dome_top(data))?;

    let workers = std::thread::available_parallelism().map_or(1, std::num::NonZero::get);
    let mut rows: Vec<(Row, Row)> = vec![(Row::empty(y), Row::empty(y)); x.len()];
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|k| {
                scope.spawn(move || {
                    (k..x.len())
                        .step_by(workers)
                        .map(|i| {
                            let p = x.node(i).exp();
                            let liquid_wanted = p >= p_triple * (-3.0 * LN_P_STEP).exp();
                            (
                                i,
                                row(data, p, y, Branch::Vapour),
                                if liquid_wanted {
                                    row(data, p, y, Branch::Liquid)
                                } else {
                                    Row::empty(y)
                                },
                            )
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        for handle in handles {
            // A worker that panicked leaves its rows invalid: the exact
            // equation then answers there.
            if let Ok(done) = handle.join() {
                for (i, vapour, liquid) in done {
                    if let Some(slot) = rows.get_mut(i) {
                        *slot = (vapour, liquid);
                    }
                }
            }
        }
    });
    let (vapour, liquid): (Vec<_>, Vec<_>) =
        rows.into_iter().map(|(v, l)| (v.nodes, l.nodes)).unzip();
    let critical = molar(data, data.data.rho_crit, data.data.t_crit);
    let gas = molar(data, vapour_density(data, 1e5, 300.0)?, 300.0);
    Ok(Tables {
        fluid,
        vapour: Grid2::from_rows(x, y, vapour),
        liquid: Grid2::from_rows(x, y, liquid),
        saturation: SaturationTable::build(data).context("saturation table")?,
        dome: (p_triple, p_top),
        s_crit: critical.s / m,
        rho_crit: data.data.rho_crit * m,
        t_triple: data.data.t_triple,
        t_crit: data.data.t_crit,
        reference: (gas.s / m, gas.cp / m, data.data.gas_constant / m),
    })
}

/// An axis covering `start ..= end` in steps of about `step`, with the margin
/// beyond both ends.
fn axis(start: f64, end: f64, step: f64) -> Axis {
    let cells = ((end - start) / step).ceil().max(1.0);
    let margin = f64::from(MARGIN) * step;
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a positive, modest node count"
    )]
    let len = cells as usize + 1 + 2 * MARGIN as usize;
    Axis::new(start - margin, start + cells * step + margin, len)
}

/// One row of a branch table.
#[derive(Debug, Clone, PartialEq)]
struct Row {
    nodes: Vec<Node>,
}

impl Row {
    fn empty(y: Axis) -> Self {
        Self {
            nodes: vec![INVALID; y.len()],
        }
    }
}

/// One row of `branch` at `p` Pa: swept from the branch's far end (hot vapour,
/// cold liquid) towards the dome, each node warm-started from the last, until
/// the branch ends. Before a spinodal the last nodes are dropped.
fn row(fluid: &Fluid, p: f64, y: Axis, branch: Branch) -> Row {
    let mut out = Row::empty(y);
    if p.is_nan() || p > fluid.data.p_max {
        return out;
    }
    let m = fluid.data.molar_mass;
    // Entropies beyond the branch's temperature range at this pressure.
    let t_edge = match branch {
        Branch::Vapour => fluid.data.t_max,
        Branch::Liquid => fluid.data.t_triple,
    };
    let edge_ok = |s: f64, edge: f64| match branch {
        Branch::Vapour => s <= edge,
        Branch::Liquid => s >= edge,
    };
    let edge = branch
        .density(fluid, p, t_edge)
        .map_or(f64::NAN, |rho| molar(fluid, rho, t_edge).s / m);
    let order: Vec<usize> = match branch {
        Branch::Vapour => (0..y.len()).rev().collect(),
        Branch::Liquid => (0..y.len()).collect(),
    };
    let mut guess = None;
    let mut filled: Vec<usize> = Vec::with_capacity(y.len());
    for j in order {
        let s = y.node(j);
        if edge.is_finite() && !edge_ok(s, edge) {
            continue;
        }
        match branch.node(fluid, p, s * m, guess) {
            Ok(node) => {
                guess = Some(node[0].exp());
                filled.push(j);
                if let Some(slot) = out.nodes.get_mut(j) {
                    *slot = node;
                }
            }
            Err(end) if !filled.is_empty() => {
                if end == End::Spinodal {
                    for &j in filled.iter().rev().take(SPINODAL_GUARD) {
                        if let Some(slot) = out.nodes.get_mut(j) {
                            *slot = INVALID;
                        }
                    }
                }
                break;
            }
            Err(_) => {}
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_axis_covers_its_range_with_margins() {
        let a = axis(1.0, 2.0, 0.1);
        let margin = MARGIN as usize;
        assert_eq!(a.len(), 11 + 2 * margin);
        assert!((a.node(margin) - 1.0).abs() < 1e-12);
        assert!((a.node(margin + 10) - 2.0).abs() < 1e-12);
    }
}
