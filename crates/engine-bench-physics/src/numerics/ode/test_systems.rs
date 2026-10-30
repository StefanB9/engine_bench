//! Reference systems with known solutions, shared by the integrator tests.

use super::{Direction, EventSet, OdeSystem};
use crate::PhysicsError;

pub(super) struct Decay {
    pub(super) rate: f64,
}

impl OdeSystem<1> for Decay {
    fn derivatives(&mut self, _t: f64, y: &[f64; 1], dydt: &mut [f64; 1]) {
        dydt[0] = -self.rate * y[0];
    }
}

/// y = [x, v], unit angular frequency: x = cos t.
pub(super) struct Oscillator;

impl OdeSystem<2> for Oscillator {
    fn derivatives(&mut self, _t: f64, y: &[f64; 2], dydt: &mut [f64; 2]) {
        dydt[0] = y[1];
        dydt[1] = -y[0];
    }
}

/// y = [height, vertical velocity] under gravity.
pub(super) struct Ball;

pub(super) const G: f64 = 9.81;

impl OdeSystem<2> for Ball {
    fn derivatives(&mut self, _t: f64, y: &[f64; 2], dydt: &mut [f64; 2]) {
        dydt[0] = y[1];
        dydt[1] = -G;
    }
}

pub(super) struct Ground(pub(super) Direction);

impl EventSet<2, 1> for Ground {
    fn directions(&self) -> [Direction; 1] {
        [self.0]
    }

    fn values(&mut self, _t: f64, y: &[f64; 2], g: &mut [f64; 1]) {
        g[0] = y[0];
    }
}

/// Fires at apogee (velocity falling through zero) and on ground contact.
pub(super) struct ApogeeAndGround;

impl EventSet<2, 2> for ApogeeAndGround {
    fn directions(&self) -> [Direction; 2] {
        [Direction::Falling, Direction::Falling]
    }

    fn values(&mut self, _t: f64, y: &[f64; 2], g: &mut [f64; 2]) {
        *g = [y[1], y[0]];
    }
}

/// Finite only for t < 1.
pub(super) struct BlowsUp;

impl OdeSystem<1> for BlowsUp {
    fn derivatives(&mut self, t: f64, _y: &[f64; 1], dydt: &mut [f64; 1]) {
        dydt[0] = (1.0 - t).sqrt().recip();
    }
}

pub(super) struct AlwaysNan;

impl OdeSystem<1> for AlwaysNan {
    fn derivatives(&mut self, _t: f64, _y: &[f64; 1], dydt: &mut [f64; 1]) {
        dydt[0] = f64::NAN;
    }
}

pub(super) fn physics_error(e: &anyhow::Error) -> Option<PhysicsError> {
    e.downcast_ref::<PhysicsError>().copied()
}
