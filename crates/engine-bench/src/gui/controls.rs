//! Side-panel widgets. They return [`Command`]s instead of sending them, so
//! tests can drive them without a simulation.

use eframe::egui::{self, Slider, Ui};
use engine_bench_core::spec::{ControlId, ControlKind, ControlSpec, QuantityKind};

use super::simulator::{Command, Snapshot};
use crate::units::{Preset, UnitSelection, quantity_name, units};

/// Sliders, toggles and buttons for `specs` at `values` (SI), in display
/// `units`.
pub(crate) fn engine_controls(
    ui: &mut Ui,
    specs: &[ControlSpec],
    values: &[f64],
    units: &UnitSelection,
    commands: &mut Vec<Command>,
) {
    for (index, spec) in specs.iter().enumerate() {
        let id = ControlId::new(index);
        let value = values
            .get(index)
            .copied()
            .unwrap_or_else(|| spec.default_value());
        match spec.kind {
            ControlKind::Continuous {
                quantity, min, max, ..
            } => {
                let unit = units.unit(quantity);
                let mut shown = unit.to_display(value);
                let range = unit.to_display(min)..=unit.to_display(max);
                let suffix = if unit.symbol.is_empty() {
                    String::new()
                } else {
                    format!(" {}", unit.symbol)
                };
                let slider = Slider::new(&mut shown, range)
                    .text(spec.label)
                    .suffix(suffix);
                if ui.add(slider).changed() {
                    // Converting back can land a hair outside the SI range; the
                    // engine would reject that.
                    let si = unit.to_si(shown).clamp(min, max);
                    commands.push(Command::SetControl(id, si));
                }
            }
            ControlKind::Toggle { .. } => {
                let mut on = value >= 0.5;
                if ui.checkbox(&mut on, spec.label).changed() {
                    commands.push(Command::SetControl(id, f64::from(u8::from(on))));
                }
            }
            ControlKind::Trigger => {
                if ui.button(spec.label).clicked() {
                    commands.push(Command::SetControl(id, 1.0));
                }
            }
        }
    }
}

/// Pause / resume, single step, reset and time scale.
pub(crate) fn transport(ui: &mut Ui, snapshot: &Snapshot, commands: &mut Vec<Command>) {
    ui.horizontal(|ui| {
        if snapshot.paused {
            if ui.button("Resume").clicked() {
                commands.push(Command::Resume);
            }
            if ui.button("Step").clicked() {
                commands.push(Command::Step);
            }
        } else if ui.button("Pause").clicked() {
            commands.push(Command::Pause);
        }
        if ui.button("Restart").clicked() {
            commands.push(Command::Reset);
        }
    });
    let mut scale = snapshot.time_scale;
    let slider = Slider::new(&mut scale, 0.001..=100.0)
        .logarithmic(true)
        .text("Time scale")
        .suffix("×");
    if ui.add(slider).changed() {
        commands.push(Command::SetTimeScale(scale));
    }
}

fn preset_name(preset: Preset) -> &'static str {
    match preset {
        Preset::Si => "SI",
        Preset::Engineering => "Engineering",
        Preset::Imperial => "Imperial",
    }
}

/// Preset menu and per-quantity unit overrides.
pub(crate) fn unit_picker(ui: &mut Ui, selection: &mut UnitSelection) {
    egui::ComboBox::from_label("Preset")
        .selected_text(preset_name(selection.preset()))
        .show_ui(ui, |ui| {
            for preset in [Preset::Si, Preset::Engineering, Preset::Imperial] {
                if ui
                    .selectable_label(selection.preset() == preset, preset_name(preset))
                    .clicked()
                {
                    selection.set_preset(preset);
                }
            }
        });
    ui.collapsing("Per quantity", |ui| {
        for kind in QuantityKind::ALL {
            let current = selection.unit(kind).symbol;
            egui::ComboBox::from_label(quantity_name(kind))
                .selected_text(if current.is_empty() { "–" } else { current })
                .show_ui(ui, |ui| {
                    for unit in units(kind) {
                        let label = if unit.symbol.is_empty() {
                            "–"
                        } else {
                            unit.symbol
                        };
                        if ui.selectable_label(unit.symbol == current, label).clicked() {
                            // Every listed unit belongs to `kind`, so this
                            // cannot fail.
                            let _ = selection.set_override(kind, unit.symbol);
                        }
                    }
                });
        }
    });
}

#[cfg(test)]
mod tests {
    use egui_kittest::{
        Harness,
        kittest::{By, Queryable},
    };
    use engine_bench_core::engine::EngineId;

    use super::*;

    struct Panel {
        specs: Vec<ControlSpec>,
        values: Vec<f64>,
        units: UnitSelection,
        snapshot: Snapshot,
        commands: Vec<Command>,
    }

    fn panel() -> Panel {
        let engine = EngineId::TwoTankDemo;
        Panel {
            specs: engine.controls().to_vec(),
            values: engine
                .controls()
                .iter()
                .map(ControlSpec::default_value)
                .collect(),
            units: UnitSelection::new(Preset::Engineering),
            snapshot: Snapshot::empty(engine),
            commands: Vec::new(),
        }
    }

    fn harness(panel: Panel) -> Harness<'static, Panel> {
        Harness::new_ui_state(
            |ui, p: &mut Panel| {
                engine_controls(ui, &p.specs, &p.values, &p.units, &mut p.commands);
                transport(ui, &p.snapshot, &mut p.commands);
                unit_picker(ui, &mut p.units);
            },
            panel,
        )
    }

    #[test]
    fn test_controls_show_every_control_in_display_units() -> anyhow::Result<()> {
        let mut h = harness(panel());
        h.run();
        for label in [
            "Valve opening",
            "Initial pressure A (on reset)",
            "Heater power, tank B",
            "Reset",
        ] {
            anyhow::ensure!(
                h.query_all_by_label_contains(label).next().is_some(),
                "missing {label}"
            );
        }
        // 8 bar default shown in the engineering preset.
        anyhow::ensure!(
            h.query_all(
                By::new()
                    .predicate(|node| { node.value().is_some_and(|value| value.contains("bar")) })
            )
            .next()
            .is_some(),
            "pressure not in bar"
        );
        anyhow::ensure!(h.state().commands.is_empty(), "commands without input");
        Ok(())
    }

    #[test]
    fn test_controls_trigger_button_sends_command() -> anyhow::Result<()> {
        let mut h = harness(panel());
        h.run();
        h.get_by_label("Reset").click();
        h.run();
        let reset = ControlId::new(4);
        anyhow::ensure!(
            h.state()
                .commands
                .iter()
                .any(|c| matches!(c, Command::SetControl(id, _) if *id == reset)),
            "{:?}",
            h.state().commands
        );
        Ok(())
    }

    #[test]
    fn test_transport_pauses_steps_and_resumes() -> anyhow::Result<()> {
        let mut h = harness(panel());
        h.run();
        h.get_by_label("Pause").click();
        h.run();
        anyhow::ensure!(
            h.state().commands.contains(&Command::Pause),
            "{:?}",
            h.state().commands
        );
        h.state_mut().snapshot.paused = true;
        h.state_mut().commands.clear();
        h.run();
        h.get_by_label("Step").click();
        h.run();
        h.get_by_label("Resume").click();
        h.run();
        anyhow::ensure!(
            h.state().commands.contains(&Command::Step)
                && h.state().commands.contains(&Command::Resume),
            "{:?}",
            h.state().commands
        );
        h.get_by_label("Restart").click();
        h.run();
        anyhow::ensure!(
            h.state().commands.contains(&Command::Reset),
            "{:?}",
            h.state().commands
        );
        Ok(())
    }
}
