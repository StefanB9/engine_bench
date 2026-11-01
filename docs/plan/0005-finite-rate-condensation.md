# Plan 0005 — Finite-Rate Condensation in Nozzles

Parent plans: `0004-m3-first-engines.md` (frozen nozzle vapour, cold-gas thruster). Frozen vapour
fails where real vapour condenses: argon through the thruster's default nozzle (area ratio 49)
would reach ~13 K, below where the frozen model trusts the equation of state, and CO₂ reaches its
spinodal or the solid's region within a few bar. This plan replaces the frozen limit with
condensation at a finite rate.

## Decisions (from planning Q&A)

| Topic | Decision |
|---|---|
| Nozzle shape | **Conical**: converging cone, throat radius of curvature, diverging cone with a half-angle control (default 15°) |
| Droplet population | **Droplet classes**: each marching step's nuclei form a class whose radius grows downstream (polydisperse, exact in steady 1D flow); similar classes merged |
| Condensate | **Supercooled liquid** (Ostwald's step rule): liquid properties extrapolated below the triple point |
| Validation | **Measured data**: N₂ condensation onset (Yale nozzle B) and, if its data become reachable, steam in Moore's nozzle B; physics limits throughout |

## Physics (`engine-bench-physics`, `flow::condensing`)

- **Geometry**: `ConicalNozzle` gives A(x) and dA/dx: converging cone (half-angle fixed at 30°)
  into a circular-arc throat (radius of curvature a multiple of the throat radius), diverging
  cone to the exit. Planar nozzles (constant depth, as the validation nozzle) as a second shape
  behind the same contour trait.
- **Properties** per fluid: surface tension σ(T) from the CoolProp correlations (generated
  offline by a `tools/reference` script, like the equations of state), liquid density, latent
  heat and vapour pressure from the reference equations (N₂, Ar, He, CO₂) or IF97 with the
  IAPWS surface tension (water); below the triple point extrapolated as supercooled liquid.
- **Nucleation**: classical theory (Volmer–Frenkel–Zeldovich) with the Kantrowitz
  non-isothermal correction; critical radius r* = 2σ / (ρ_l·R·T·ln S).
- **Growth**: Knudsen-bridged heat-transfer-limited law (Gyarmathy / Young), droplet temperature
  with the capillary (Kelvin) correction.
- **Flow**: steady quasi-1D, inviscid: mass, momentum and energy of the vapour and the droplet
  classes, marched in x by Dormand–Prince from just downstream of the throat, where the frozen
  solution still holds. The vapour on the reference equation (metastable vapour branch).
- **Shocks**: a normal shock in the nozzle passes the droplets through unchanged; behind it the
  same growth law evaporates them in the subsonic flow.
- **Outputs**: exit state, thrust, liquid mass fraction, droplet spectrum (mean radius), Wilson
  point (where condensation sets in).

## Engine (`engine-bench-core`)

- Cold-gas thruster: mass flow from the throat as now (condensation starts downstream); thrust
  and exit state from the condensing march. New control: divergence half-angle. New channels:
  exit liquid fraction, mean droplet radius, Wilson point position. Argon works; CO₂ joins the
  propellants.
- Real-time cost measured; if a march per step is too slow, thrust is evaluated less often and
  interpolated (the fidelity rule: cut only where measured real time fails).

## Validation

- **N₂ onset, Yale nozzle B** (Wegener et al., *Study of experiments on condensation of
  nitrogen…*, final report to NASA Langley, grant NSG-1612, May 1980, NASA CR-163217;
  transcribed from the scanned Tables 1, 2 and 4):
  - Nozzle B, planar: throat height h* = 0.40 cm, depth b = 0.70 cm, A* = 0.28 cm², throat
    radius of curvature R* = 2 cm, divergence α = 6.8° between the walls (3.4° each: the only
    reading consistent with the length), length l = 8.0 cm, exit height h = 1.33 cm,
    A/A* = 3.33 (inviscid exit M = 2.75 for γ = 1.4).
  - Onset of condensation observed (light scattering) at the nozzle exit for these supply
    states (p₀ in torr, T₀ in K): #2 463/113, #6 244/107, #9 244/107, #11 338/110, #12 405/111,
    #13 409/111, #20 671/122, #22 468/116, #23 502/116, #27 710/118, #28 973/124, #29 744/121,
    #30 893/122, #31 871/121, #32 1150/129, #37 1209/132, #40 1357/134. Exit states bracketed
    between the inviscid limit (p_e/p₀ = 0.0398) and a room-temperature calibration
    (p_e/p₀ ≈ 0.054–0.066).
  - Mean onset states (Table 4): about 48–53 K at 16–52 torr, Mach 2.5–2.7, ~10 K below
    saturation, at cooling rates of ~0.1 °C/µs.
  - Test: for each supply state the predicted Wilson point lies near the exit. The tolerance is
    set from the measured agreement and reported (classical nucleation theory is known to be
    off by orders of magnitude in rate, much less in onset temperature).
- **Steam, Moore et al. (1973) nozzle B** (p₀ = 25 kPa, T₀ = 358 K; pressure along the nozzle and
  mean droplet radius): only if the measured data become available (the International Wet
  Steam Modelling Project report is not reachable from this environment).
- **Physics limits**: no nucleation reduces to the frozen solution; fast growth approaches
  equilibrium (HEM); mass and energy conserved along the march; droplet classes merged without
  changing the liquid mass.

## Delivery

1. `feature/condensation-properties`: surface tension (generated), supercooled-liquid
   properties, `ConicalNozzle` (and the planar contour).
2. `feature/condensing-nozzle`: nucleation, growth, the marching solver with droplet classes,
   shocks, validation against the Yale data.
3. `feature/thruster-condensation`: the thruster on the condensing nozzle, half-angle control,
   new channels, CO₂.
