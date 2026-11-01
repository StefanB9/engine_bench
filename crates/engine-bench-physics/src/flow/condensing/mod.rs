//! Condensation at a finite rate in nozzle flows: vapour expanding past
//! saturation nucleates droplets that grow, releasing latent heat into the
//! flow (plan 0005).

mod contour;
mod growth;
mod march;
mod model;
mod nucleation;
mod properties;
mod vapour;

pub use contour::Contour;
pub use march::{CondensingFlow, CondensingNozzle, FlowSample};
pub use model::CondensationModel;
pub use properties::{Condensable, CondensationProperties, VapourTransport};
pub use vapour::{CondensingFluid, VapourPoint};
pub(crate) use vapour::{kernel, saturation_temperature};
