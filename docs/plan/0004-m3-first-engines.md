# Plan 0004 — M3 First Engines: Aeolipile, Cold-Gas Thruster, Water Rocket

Parent plan: `0001-architecture-and-roadmap.md`, milestone M3 (catalogue engines 1–3).

## Decisions (from planning Q&A)

| Topic | Decision |
|---|---|
| Steam/water | **Full IAPWS-IF97**, all five regions, with verification against the IF97 tables |
| Cold, dense gas | **Reference multiparameter Helmholtz EOS** for the four pure propellants, **Peng–Robinson** for mixtures; both used in volumes **and through nozzles and orifices** (revised at the start of PR 4) |
| Cold-gas propellants | **N₂, Ar, He, CO₂**. CO₂ is stored as liquid + vapour |
| Below the NASA data range | **Extended on physical grounds**: the cold-gas set is valid down to 60 K through the constant-cp ideal-gas extrapolation (exact for monatomic He and Ar, ~0.1 % for N₂; CO₂ is bounded by its 216.6 K triple point) |
| Water rocket | **2D flight with launch tube and parachute**, as a **3-DOF planar rigid body** (decided at the start of PR 1) |
| Camera | **Scene sets its view**: an optional explicit view rectangle, still drawn from channels |
| Real-gas shocks | **Full real-gas Rankine–Hugoniot** in the nozzle's shock regime |
| Pressure regulator | **Optional** in the cold-gas thruster (setpoint 0 = pure blowdown) |
| Delivery | **Physics first, then engines** (8 PRs below) |

## Data sources (network: IAPWS and NIST are unreachable; GitHub and PyPI are reachable)

- **IF97 coefficients**: transcribed from `CoolProp/IF97` (`IF97.h`, MIT, © Ian H. Bell) into
  `engine-bench-physics/src/steam/coefficients/`, with attribution in each file. Checked against
  the verification values of IAPWS R7-97(2012) as reproduced in that repository.
- **Reference values for tests**: generated offline with `CoolProp` 8.0 (Python, MIT) from its
  reference equations of state (Span–Wagner for N₂ and CO₂; Tegeler for Ar; Ortiz-Vega for He) and its
  IF97 backend, then pasted into tests as golden values with the generating script committed
  under `tools/reference/` for reproducibility. CoolProp is not a build dependency.
- **Helium**: its NASA Glenn record is vendored into `data/thermo.inp` like the other 12 species.

## Physics additions (`engine-bench-physics`)

1. **`steam`: IAPWS-IF97.**
   - Forward equations for regions 1 and 2 (Gibbs), 3 (Helmholtz), 4 (saturation line) and
     5 (Gibbs). Region selection from (p, T), and backward equations T(p, h), T(p, s).
   - A state solver from (u, v), which is what a volume conserves: single-phase via Newton on
     (T, ρ), two-phase via T_sat with quality from v.
   - Hot kernels stay infallible; validation at the edges.
2. **`gas::eos`: equation-of-state trait.**
   - `ControlVolume<S, N, E: EquationOfState>`. `IdealGas` is today's behaviour, zero-cost;
     existing code keeps `IdealGas` as the default type parameter.
   - `PengRobinson` for mixtures, with critical constants and acentric factors (N₂, O₂, Ar, He,
     CO₂, H₂O, H₂, CO), van der Waals one-fluid mixing (k_ij = 0), and u, h, s departures added
     to the NASA ideal-gas parts. Single-phase: a state inside the spinodal is reported as
     invalid.
   - **Revised at the start of PR 4.** Measured against CoolProp's reference equations, PR misses
     Z by 2–4 % at 300 K and up to 10 % for N₂ and Ar at 150 K and 300 bar, so the pure
     propellants get **reference multiparameter Helmholtz equations** instead: Span et al. (N₂),
     Tegeler et al. (Ar), Ortiz-Vega et al. (He) and Span–Wagner (CO₂), coefficients from
     CoolProp's fluid files (MIT). Vapour–liquid equilibrium (equal Gibbs energy on the
     reference equation) makes a CO₂ volume two-phase; PR's own pure-species VLE is no longer
     needed.
3. **`flow::real`: real-fluid isentropic flow** for any `RealFluid` (IF97 water, PR mixtures).
   - Isentropic states from (p, s₀). The homogeneous equilibrium model (HEM) applies inside the
     two-phase dome.
   - Choked mass flux = max over p of ρ·√(2(h₀ − h)). This gives a real-fluid orifice and
     nozzle (subsonic, choked, under-, ideally and over-expanded).
   - Normal shocks are solved with the real EOS: conservation of mass, momentum and energy across
     the shock (Rankine–Hugoniot), a 2D Newton per shock. The shock is located in the nozzle so
     that the exit state matches the back pressure, as the perfect-gas nozzle does today.
   - Same API shape as `flow::Orifice` / `flow::Nozzle`.
4. **`flight`**, used by the water rocket:
   - International Standard Atmosphere, −500 m to 20 km (troposphere and the isothermal lower
     stratosphere).
   - Planar 3-DOF rigid body (x, y, pitch) with mass, inertia and centre of gravity supplied per
     evaluation, so they can change as propellant leaves.
   - Aerodynamics against the air-relative velocity (wind included): axial drag, normal force at
     the centre of pressure, pitch damping, and jet damping from the exhaust mass flow.
   - Barrowman normal-force slope and centre of pressure from nose, body and fin geometry, so
     stability emerges from the geometry: a rocket with its centre of pressure ahead of its
     centre of gravity tumbles.
5. **Thermo data**: helium record and a `ColdGas` species set (N₂, Ar, He, CO₂).

## Core additions (`engine-bench-core`)

- `Scene::set_view(min, max)` / `Scene::view()`. The painter fits the explicit view when it is
  set, otherwise the scene bounds as now.

## Engines

### #2 Cold-gas thruster
- **Hardware**: propellant selector (N₂, Ar, He, CO₂), PR tank (CO₂ liquid + vapour), optional
  pressure regulator (setpoint; 0 = unregulated blowdown) feeding a small plenum, real-fluid nozzle
  into an adjustable ambient (vacuum to 1 atm), fire valve, wall heat transfer as in the demo.
- **Channels**: tank p, T and liquid fraction; thrust; Isp; ṁ; total impulse; exit Mach;
  nozzle regime; separation flag.
- **Scene**: tank cutaway with liquid level, valve and regulator, nozzle, and a plume whose
  length scales with thrust. Real plume structure (shock diamonds) waits for 2D CFD in M4.

### #3 Water rocket
- **Configuration**: 2 L PET bottle; water fill; pressure; nozzle diameter; launch angle; launch-tube
  length; drag coefficient; parachute (deploys at apogee plus a delay, with its own area and C_d).
- **Phases, as DP45 events**:
  1. Launch tube: air pushes the rocket off the tube; no water flows.
  2. Water expulsion: unsteady Bernoulli. Air expands as an ideal-gas `ControlVolume` of growing
     volume. Water density comes from IF97.
  3. Air blowdown: real-fluid orifice. Air at a few bar is near-ideal, so `IdealGas` is used.
  4. Coast.
  5. Parachute descent.
  6. Ground contact.
- **Flight**: 3-DOF planar rigid body. Thrust acts along the body axis; attitude follows from
  the aerodynamic moments rather than being assumed. Geometry controls cover nose, body, fins
  (count, root and tip chord, span, sweep) and ballast, so an unstable configuration tumbles.
  CG and inertia shift as the water leaves.
- **Channels**: x, y, velocity, attitude, water mass, air p and T, thrust, phase.
- **Scene**: the camera follows the rocket through `Scene::set_view`, with ground line, launch
  tube, altitude grid and trajectory trace, a bottle cutaway showing the water level, and the
  parachute when deployed.

### #1 Aeolipile
- **Physics**:
  - Sphere boiler as a two-phase IF97 volume, heated by a burner (power control) with heat loss
    to ambient.
  - Saturated vapour leaves through two tangential real-fluid nozzles.
  - Rotor torque τ = ṁ·r·(v_jet − ω·r) (Euler turbine equation), against bearing friction
    (Coulomb + viscous) and air drag on the arms.
  - Rotor angle is a scene-state channel.
- **Channels**: p, T, water mass, quality, steam flow, jet velocity, ω (rpm), torque, efficiency
  (shaft power / burner power).
- **Scene**: rotating sphere and arms, water level inside, steam jets, flame.

## Performance budget

Real-fluid nozzle evaluations involve root solves. Targets, measured by criterion benches:
- PR state from (u, v): under 2 µs.
- Real-fluid orifice / nozzle point: under 20 µs.

If an engine's real-time step exceeds about 30 % of wall time, the plan is to tabulate isentropes
per propellant at reset, without changing the API.

Measured after PR 5 (criterion, `bench-unwind`, one point each):

| Kernel | Time |
|--------|------|
| Reference N2 state at (p, s) | ~125 µs |
| Real orifice, reference N2 at 20 MPa, unchoked / choked | ~1.6 ms / ~1.3 ms |
| Real orifice, wet steam (IF97), choked | ~120 µs |
| Real nozzle, reference N2, under-expanded | ~5.7 ms |
| Real nozzle, reference N2, shock in the nozzle | ~70 ms |

All points miss the 20 µs target by 6× (IF97) to more than 1000×. The cost is the
reference EOS: each (p, s) state nests a temperature Newton around a density solve (~10 µs,
~12 Helmholtz evaluations), and a choke search makes ~10 of them (a nozzle shock search ~400).
Tuning the solves cannot close this gap, so the isentrope tabulation above is expected in PR 6, as
soon as the cold-gas thruster's real-time step is measured.

## Validation

- **IF97**: official verification values for every region and the saturation line. Round trips
  (p, T) ↔ (p, h) ↔ (p, s) ↔ (u, v). Continuity across region boundaries (proptest).
- **Peng–Robinson**: the analytic critical compressibility Z_c = 0.3074; Z against CoolProp with
  the measured PR accuracy stated per species; thermodynamic consistency (Maxwell relations,
  c_v = ∂u/∂T) and the ideal-gas limit at low density (proptests).
- **Reference equations**: the published check values of each equation, and CoolProp goldens
  for p, u, h, s, c_p, w and saturation (CO₂) to ~1e-9.
- **Real-fluid flow**: equals the perfect-gas orifice/nozzle in the ideal limit; choked flux
  matches the analytic perfect-gas value; real-gas shock jumps conserve mass, momentum and energy
  and reduce to the perfect-gas relations at low density; thrust stays continuous across regime
  boundaries; steam choked flow against CoolProp-derived isentropes;
  mass and energy conserved through coupled volumes.
- **Engines**:
  - Cold gas: ideal-limit thrust and Isp against `Nozzle`; total impulse ≈ ∫F dt; blowdown
    conservation; isothermal vs adiabatic limits.
  - Water rocket: analytic vacuum-free-fall coast; zero-drag 45° launch range against
    projectile motion; a stable configuration weathervanes into the wind and an unstable one
    tumbles; water phase against the incompressible Bernoulli + adiabatic-air ODE.
  - Aeolipile: steady ω with no friction equals v_jet / r; boiler energy balance.

## Delivery: 8 PRs, physics first

1. `feature/m3-scene-view-flight`: `Scene` view, ISA atmosphere, 3-DOF rigid body,
   aerodynamics, Barrowman.
2. `feature/m3-if97-basic`: IF97 regions 1, 2, 4, 5 + verification.
3. `feature/m3-if97-region3`: region 3, (u, v) solver, two-phase water/steam volume.
4. `feature/m3-peng-robinson`: EOS trait, `ControlVolume` generic over the EOS, PR for mixtures,
   He, the cold-gas species set and its low-temperature range.
4b. `feature/m3-reference-eos`: reference Helmholtz equations for N₂, Ar, He, CO₂; two-phase
   pure-fluid volumes (CO₂).
5. `feature/m3-real-flow`: real-fluid isentropic flow, orifice, nozzle, HEM.
6. `feature/m3-cold-gas-thruster`: engine #2.
7. `feature/m3-water-rocket`: engine #3.
8. `feature/m3-aeolipile`: engine #1.
