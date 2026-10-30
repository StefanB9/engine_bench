//! What an engine exposes: controls it accepts, channels it reports, the
//! quantity each one measures, and how accurately to step it.

mod channel;
mod control;
mod quantity;
mod step;

pub use channel::{ChannelRole, ChannelSpec, find_channel};
pub use control::{ControlId, ControlKind, ControlSpec, find_control};
pub use quantity::QuantityKind;
pub use step::{Accuracy, StepContext};
