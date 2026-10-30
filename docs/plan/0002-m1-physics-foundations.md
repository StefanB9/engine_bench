# Plan 0002 — M1 Physics Foundations

Parent plan: `0001-architecture-and-roadmap.md`, milestone M1. All code lives in
`engine-bench-physics`. No engine models, no GUI.

## Decisions (from planning Q&A)

| Topic | Decision |
|---|---|
| Thermo data | NASA Glenn 9-coefficient (`thermo.inp`, NASA/TP-2002-211556), as used by CEA. Supersedes "7-coefficient" in plan 0001 |
| Species representation | Compile-time species sets; mixtures are `[f64; N]`, no heap |
| Hot-path errors | Validate at the edges: constructors and solver steps return `Result`; property/flux kernels are infallible on validated inputs (CLAUDE.md gets this exception) |
| Nozzle scope | Full quasi-1D regimes including normal shock in the nozzle; flow separation flagged |
| Control-volume state | Conserved variables: species masses + internal energy; T recovered by safeguarded Newton |
| Nozzle/orifice gas model | Constant γ frozen at upstream state; thermally perfect expansion added in M5 |
| Benchmarks | criterion benches for hot kernels from M1 |
| Delivery | 4 sequential PRs |

## Scope

In: species data, mixture properties, ideal-gas control volume, orifice flow, isentropic and
normal-shock relations, converging–diverging nozzle, ODE integrators, root finders.

Out (later milestones): chemical equilibrium (M5), steam/IAPWS-IF97 (M3), 1D/2D solvers,
heat-transfer correlations, real-gas effects.

## Modules

```
engine-bench-physics/
├── build.rs                     Parses data/thermo.inp → const tables in OUT_DIR
├── data/thermo.inp              Vendored subset (12 species) + source/attribution header
└── src/
    ├── units.rs                 (exists)
    ├── error.rs                 PhysicsError (thiserror)
    ├── thermo/
    │   ├── mod.rs               re-exports
    │   ├── nasa9.rs             cp/R, h/RT, s/R for one species (infallible)
    │   ├── species.rs           Species (name, molar mass, elements, intervals); generated tables
    │   ├── species_set.rs       SpeciesSet<const N> trait; Air3, Major7 sets
    │   └── mixture.rs           Mixture<S, N>: cp, cv, h, u, s, R, γ, M; T from u / h
    ├── flow/
    │   ├── isentropic.rs        p/p0, T/T0, A/A*(M), M(A/A*) both branches, Prandtl–Meyer
    │   ├── shock.rs             Normal-shock relations
    │   ├── orifice.rs           Compressible orifice (choked/unchoked, Cd, bidirectional)
    │   └── nozzle.rs            CD nozzle operating point: regime, ṁ, exit state, thrust, Isp
    ├── gas/
    │   └── control_volume.rs    ControlVolume<S, N>: masses + U + V → p, T; derivative assembly
    └── numerics/
        ├── roots.rs             Brent, safeguarded Newton
        └── ode.rs               OdeSystem<const N>, Rk4, DormandPrince45 (adaptive, PI control, events)
```

Species vendored now: N2, O2, Ar, H2O, CO2, CO, H2 (used in M1) plus OH, H, O, NO, N
(data only, for M5). Sets: `Air3` (N2, O2, Ar), `Major7` (all seven majors).

## Key design points

- **Species sets on stable Rust:** `trait SpeciesSet<const N: usize> { const SPECIES: [&'static Species; N]; }`
  and `Mixture<S: SpeciesSet<N>, const N: usize>`, with aliases (`type AirMixture = Mixture<Air3, 3>`).
  Avoids the incomplete `generic_const_exprs` feature.
- **Generated tables:** `build.rs` parses the vendored `thermo.inp` (fixed-width Fortran
  format, `D` exponents) into `const` arrays. No runtime parsing, exact match with CEA
  data, and a new species only needs a data entry. The build script returns `Result` (lints apply).
- **Out-of-range temperature:** below 200 K, extrapolate with constant cp (linear h). Above
  the top interval, use the last interval. Kernels never fail; `ControlVolume` checks T at
  each step against the set's valid range and returns `PhysicsError::TemperatureOutOfRange`.
- **T inversion:** Newton on u(T) or h(T), bracketed to [T_min, T_max], warm-started from the
  previous T (one or two iterations per step in practice).
- **Orifice:** isentropic flow with frozen γ at upstream conditions, choked above the
  critical pressure ratio. Below a small Δp it switches to a C¹ linearisation, so the ODE
  Jacobian stays bounded (no infinite slope at Δp = 0).
- **Nozzle:** given p0, T0, gas, A*, Ae, pb → regime enum (Subsonic, ShockInNozzle{x/A},
  ShockAtExit, OverExpanded, Ideal, UnderExpanded) plus ṁ, exit p/T/v/M, thrust, Isp.
  Separation flagged with the Summerfield criterion (pe < 0.4·pb). Calorically perfect
  gas with γ frozen at chamber state; the thermally perfect variant arrives in M5.
- **ODE:** fixed-size state `[f64; N]`, `derivatives(&mut self, t, &y, &mut dydt)` is
  infallible. `step()` returns `Result`, failing on non-finite state or step-size underflow.
  Dormand–Prince with dense output and event location (root of g(t, y)), needed in M3 for
  "water exhausted" and "apogee".
- **Errors:** `PhysicsError` (thiserror): InvalidComposition, TemperatureOutOfRange,
  NonFiniteState, StepSizeUnderflow, NoConvergence, InvalidGeometry. Display tested per variant.

## Validation (TDD; golden values from CEA / JANAF / textbook)

- **Species:** cp, h, s of N2, O2, H2O, CO2 at 298.15, 1000, 3000 K vs JANAF tables
  (tolerance about 0.1 %). Interval continuity at 1000/6000 K (proptest).
- **Mixture:** air γ ≈ 1.400 at 300 K, R = 287.0 J/kg/K. Proptest: h = u + RT,
  dh/dT = cp (finite difference), T(u(T)) round-trip over 200–6000 K, 1000 cases.
- **Isentropic/shock:** A/A* = 1.6875 at M = 2 (γ = 1.4); M2 = 0.5774, p2/p1 = 4.5 at M1 = 2.
  Round-trip M(A/A*(M)), proptest.
- **Nozzle:** regime boundaries and shock location vs Anderson, *Modern Compressible Flow*,
  worked examples. Thrust continuity across regime boundaries (proptest).
- **ODE:** decay and harmonic oscillator vs analytic. Order-of-convergence check (RK4 error
  ratio ≈ 16 per step halving, DP45 ≈ 32). Event time accuracy.
- **Control volume:** proptest conservation of total mass and energy for two volumes joined
  by an orifice (relative 1e-12). Adiabatic compression conserves entropy. Constant-γ tank
  blowdown vs analytic solution.

## Dependencies

- `thiserror` (physics): domain error enum. Required by CLAUDE.md; no alternative considered.
- `proptest` (dev): property tests, required by CLAUDE.md.
- `criterion` (dev): kernel benchmarks (`benches/`), features `cargo_bench_support` only.
- Versions: thiserror 2.0.21, proptest 1.12.0 (feature `std`), criterion 0.8.2; all
  `default-features = false`.

## Benchmarks

`crates/engine-bench-physics/benches/kernels.rs`: species/mixture property evaluation,
T-from-u inversion, orifice mass flow, nozzle operating point, RK4/DP45 step on a
control-volume system. Run with `cargo bench -p engine-bench-physics`. They are a baseline,
not a gate.

## CLAUDE.md changes

- Error handling: "Hot-path kernels (property evaluation, fluxes, derivatives) are
  infallible on validated inputs; validation happens in constructors and solver steps,
  which return `Result`."
- Project structure: new physics modules.

## Delivery

Four sequential PRs into `dev`, each merged before the next one starts:

1. `feature/m1-numerics`: this plan, CLAUDE.md hot-path rule, `PhysicsError`, roots, ODE
   integrators, benchmark harness.
2. `feature/m1-thermo`: build.rs + vendored data, species, sets, mixture.
3. `feature/m1-flow`: isentropic, shock, orifice, nozzle.
4. `feature/m1-control-volume`: control volume, conservation tests.

Each PR updates the CLAUDE.md project structure for the modules it adds.
