//! Plot panel: selected channels, one plot per quantity kind, in display units.

use eframe::egui::Ui;
use egui_plot::{Line, Plot, PlotPoints, VLine};
use engine_bench_core::spec::{ChannelRole, ChannelSpec, QuantityKind};

use crate::units::UnitSelection;

/// Which channels are plotted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PlotSelection {
    selected: Vec<bool>,
}

impl PlotSelection {
    /// The first two observable channels of `channels`.
    pub(crate) fn new(channels: &[ChannelSpec]) -> Self {
        let mut remaining = 2;
        let selected = channels
            .iter()
            .map(|c| {
                let pick = remaining > 0 && c.role == ChannelRole::Observable;
                remaining -= usize::from(pick);
                pick
            })
            .collect();
        Self { selected }
    }

    /// Whether channel `index` is plotted.
    pub(crate) fn is_selected(&self, index: usize) -> bool {
        self.selected.get(index).copied().unwrap_or(false)
    }

    /// Checkboxes for the observable channels.
    pub(crate) fn picker(&mut self, ui: &mut Ui, channels: &[ChannelSpec]) {
        ui.horizontal_wrapped(|ui| {
            for (selected, channel) in self.selected.iter_mut().zip(channels) {
                if channel.role == ChannelRole::Observable {
                    ui.checkbox(selected, channel.label);
                }
            }
        });
    }
}

/// Draws the selected channels; `points(index)` gives `[time, SI value]` pairs
/// and `cursor` marks a time (replay).
pub(crate) fn plots(
    ui: &mut Ui,
    channels: &[ChannelSpec],
    selection: &PlotSelection,
    units: &UnitSelection,
    points: impl Fn(usize) -> Vec<[f64; 2]>,
    cursor: Option<f64>,
) {
    let mut kinds: Vec<QuantityKind> = Vec::with_capacity(channels.len());
    for (index, channel) in channels.iter().enumerate() {
        if selection.is_selected(index) && !kinds.contains(&channel.quantity) {
            kinds.push(channel.quantity);
        }
    }
    if kinds.is_empty() {
        ui.weak("Select channels to plot.");
        return;
    }
    ui.columns(kinds.len(), |columns| {
        for (column, kind) in columns.iter_mut().zip(&kinds) {
            let unit = units.unit(*kind);
            let title = format!("{} [{}]", crate::units::quantity_name(*kind), unit.symbol);
            column.label(&title);
            Plot::new(title.as_str())
                .height(column.available_height().max(120.0))
                .show(column, |plot| {
                    for (index, channel) in channels.iter().enumerate() {
                        if selection.is_selected(index) && channel.quantity == *kind {
                            let converted: Vec<[f64; 2]> = points(index)
                                .into_iter()
                                .map(|[t, v]| [t, unit.to_display(v)])
                                .collect();
                            plot.line(Line::new(channel.label, PlotPoints::from(converted)));
                        }
                    }
                    if let Some(time) = cursor {
                        plot.vline(VLine::new("cursor", time));
                    }
                });
        }
    });
}

#[cfg(test)]
mod tests {
    use engine_bench_core::engine::EngineId;

    use super::*;

    #[test]
    fn test_plot_selection_starts_with_two_observables() -> anyhow::Result<()> {
        let channels = EngineId::TwoTankDemo.channels();
        let selection = PlotSelection::new(channels);
        let picked: Vec<usize> = (0..channels.len())
            .filter(|i| selection.is_selected(*i))
            .collect();
        anyhow::ensure!(picked == [0, 1], "{picked:?}");
        anyhow::ensure!(!selection.is_selected(99), "out of range");
        Ok(())
    }
}
