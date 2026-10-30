//! The engine abstraction: every simulated engine implements [`Engine`];
//! [`EngineId`] names them, holds their static descriptions and draws their
//! scenes; [`AnyEngine`] dispatches statically over all of them.

use anyhow::Context;
use engine_bench_physics::units::Time;

use crate::{
    CoreError,
    engines::{TwoTankDemo, two_tank},
    scene::Scene,
    spec::{ChannelSpec, ControlId, ControlSpec, StepContext},
};

/// Static description of an engine model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EngineDescriptor {
    /// Which engine.
    pub id: EngineId,
    /// Display name.
    pub name: &'static str,
    /// Position in the 38-engine catalogue, `None` for demos.
    pub catalogue_number: Option<u8>,
    /// One-sentence description.
    pub summary: &'static str,
    /// Incremented whenever the channel list changes, so old recordings are
    /// recognised.
    pub schema_version: u32,
}

/// A simulated engine.
///
/// Controls and channels are fixed per engine type. The scene is drawn by
/// [`EngineId::draw`] from a channel sample, so live view and replay share one
/// code path.
pub trait Engine {
    /// Static description.
    fn descriptor(&self) -> &'static EngineDescriptor;
    /// Controls, indexed by [`ControlId`].
    fn controls(&self) -> &'static [ControlSpec];
    /// Current value of every control, in control order.
    fn control_values(&self) -> &[f64];
    /// Reported channels, in sample order.
    fn channels(&self) -> &'static [ChannelSpec];
    /// Sets a control after validating it against its spec.
    fn set_control(&mut self, id: ControlId, value: f64) -> anyhow::Result<()>;
    /// Advances the physics by `dt` (> 0) at the accuracy in `context`.
    fn advance(&mut self, dt: Time, context: &StepContext) -> anyhow::Result<()>;
    /// Simulation time since the last reset.
    fn time(&self) -> Time;
    /// Writes every channel's current value (SI) into `out`, whose length must
    /// equal the channel count.
    fn sample(&self, out: &mut [f64]) -> anyhow::Result<()>;
    /// Returns to the initial state given by the current controls.
    fn reset(&mut self) -> anyhow::Result<()>;
}

/// Every engine model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EngineId {
    /// Two tanks and a valve; exercises the runtime before the catalogue
    /// engines exist.
    TwoTankDemo,
}

impl EngineId {
    /// Every engine, in menu order.
    pub const ALL: [Self; 1] = [Self::TwoTankDemo];

    /// Stable command-line and recording key, e.g. `two-tank-demo`.
    pub fn key(self) -> &'static str {
        match self {
            Self::TwoTankDemo => "two-tank-demo",
        }
    }

    /// The engine with command-line `key`.
    pub fn from_key(key: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|id| id.key() == key)
    }

    /// Static description.
    pub fn descriptor(self) -> &'static EngineDescriptor {
        match self {
            Self::TwoTankDemo => TwoTankDemo::descriptor_static(),
        }
    }

    /// Channels in sample order, without needing an instance (for recordings).
    pub fn channels(self) -> &'static [ChannelSpec] {
        match self {
            Self::TwoTankDemo => &two_tank::CHANNELS,
        }
    }

    /// Draws the scene for one channel `sample` (live or recorded) into
    /// `scene`, replacing its contents.
    pub fn draw(self, sample: &[f64], scene: &mut Scene) -> anyhow::Result<()> {
        let expected = self.channels().len();
        if sample.len() != expected {
            return Err(CoreError::SampleLength {
                expected,
                actual: sample.len(),
            })
            .with_context(|| format!("drawing {}", self.key()));
        }
        match self {
            Self::TwoTankDemo => two_tank::draw(sample, scene),
        }
    }
}

/// Any engine, statically dispatched.
#[derive(Debug, Clone)]
pub enum AnyEngine {
    /// See [`EngineId::TwoTankDemo`].
    TwoTankDemo(TwoTankDemo),
}

impl AnyEngine {
    /// A new instance of engine `id` with default controls.
    pub fn new(id: EngineId) -> anyhow::Result<Self> {
        Ok(match id {
            EngineId::TwoTankDemo => Self::TwoTankDemo(TwoTankDemo::new()?),
        })
    }
}

macro_rules! dispatch {
    ($self:ident, $engine:ident => $body:expr) => {
        match $self {
            AnyEngine::TwoTankDemo($engine) => $body,
        }
    };
}

impl Engine for AnyEngine {
    fn descriptor(&self) -> &'static EngineDescriptor {
        dispatch!(self, e => e.descriptor())
    }

    fn controls(&self) -> &'static [ControlSpec] {
        dispatch!(self, e => e.controls())
    }

    fn control_values(&self) -> &[f64] {
        dispatch!(self, e => e.control_values())
    }

    fn channels(&self) -> &'static [ChannelSpec] {
        dispatch!(self, e => e.channels())
    }

    fn set_control(&mut self, id: ControlId, value: f64) -> anyhow::Result<()> {
        dispatch!(self, e => e.set_control(id, value))
    }

    fn advance(&mut self, dt: Time, context: &StepContext) -> anyhow::Result<()> {
        dispatch!(self, e => e.advance(dt, context))
    }

    fn time(&self) -> Time {
        dispatch!(self, e => e.time())
    }

    fn sample(&self, out: &mut [f64]) -> anyhow::Result<()> {
        dispatch!(self, e => e.sample(out))
    }

    fn reset(&mut self) -> anyhow::Result<()> {
        dispatch!(self, e => e.reset())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_engine_id_keys_round_trip() -> anyhow::Result<()> {
        for id in EngineId::ALL {
            anyhow::ensure!(!id.key().is_empty(), "{id:?} has no key");
            anyhow::ensure!(EngineId::from_key(id.key()) == Some(id), "{id:?}");
            anyhow::ensure!(id.descriptor().id == id, "{id:?} descriptor");
        }
        anyhow::ensure!(
            EngineId::from_key("no-such-engine").is_none(),
            "unknown key"
        );
        Ok(())
    }

    #[test]
    fn test_any_engine_matches_static_description() -> anyhow::Result<()> {
        for id in EngineId::ALL {
            let engine = AnyEngine::new(id)?;
            anyhow::ensure!(engine.descriptor().id == id, "{id:?}");
            anyhow::ensure!(engine.channels() == id.channels(), "{id:?} channels");
            anyhow::ensure!(
                engine.control_values().len() == engine.controls().len(),
                "{id:?} control values"
            );
            let mut sample = vec![0.0; engine.channels().len()];
            engine.sample(&mut sample)?;
            let mut scene = Scene::new();
            id.draw(&sample, &mut scene)?;
            anyhow::ensure!(!scene.shapes().is_empty(), "{id:?} drew nothing");
        }
        Ok(())
    }

    #[test]
    fn test_engine_id_draw_rejects_wrong_sample_length() -> anyhow::Result<()> {
        let mut scene = Scene::new();
        let Err(e) = EngineId::TwoTankDemo.draw(&[1.0, 2.0], &mut scene) else {
            anyhow::bail!("accepted a 2-value sample");
        };
        anyhow::ensure!(
            matches!(
                e.downcast_ref::<CoreError>(),
                Some(CoreError::SampleLength { actual: 2, .. })
            ),
            "wrong error: {e:#}"
        );
        Ok(())
    }
}
