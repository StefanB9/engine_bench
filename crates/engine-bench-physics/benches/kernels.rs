//! Hot-kernel benchmarks: a baseline for the real-time loop, not a pass/fail
//! gate.

use criterion::{Criterion, criterion_group, criterion_main};

fn benches(_c: &mut Criterion) {}

criterion_group!(kernels, benches);
criterion_main!(kernels);
