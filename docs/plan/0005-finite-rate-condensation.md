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
- **Transport**: the growth law needs the vapour's viscosity (mean free path, Knudsen number)
  and thermal conductivity, in the dilute-gas limit (condensing vapour is far below its
  critical density): Lemmon and Jacobsen (2004) for N₂ and Ar, Arp et al. (1998) and Hands
  (1981) for He, Laesecke and Muzny (2017) and Huber et al. (2016) for CO₂, IAPWS (2008, 2011)
  for water; coefficients from CoolProp's sources, checked against CoolProp goldens.
- **Model parameters**: Young's growth law with his α = 9 and β = 0 (Young 1982; the project
  found α ≈ 11 best for the low-pressure cases), condensation coefficient 1, Kantrowitz
  correction on, flat-film surface tension. Fitted to steam; for the other fluids no data
  fixes them, so the same values.
- **Formulation**: conservation held algebraically (mass flow, stagnation enthalpy, impulse
  function F = pA + ṁu with dF/dx = p·dA/dx), so the march conserves to round-off; the ODE
  carries F and the droplet classes. The droplets nucleated within one integration step
  form one class (number and liquid mass), closed when the step is accepted; with at most
  64 classes the two most alike are merged.
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
- **Steam, Moses and Stein nozzle, cases 252 and 257** (Starzmann et al., *Results of the
  International Wet Steam Modeling Project*, Proc. IMechE Part A 232(5), 2018,
  doi:10.1177/0957650918758779; measured points read from its Figures 4 and 5, the markers
  located on the page rendered at 300 dpi, to about 0.5 % of each axis; the paper gives no
  tables). Planar, throat 10 × 10 mm. The measured droplet radius is ±20 %.
  - Case 252: p₀ = 40.05 kPa, T₀ = 374.3 K. p/p₀ at x (cm): 2.054 0.3969, 2.259 0.3829,
    2.469 0.3687, 2.682 0.3544, 2.902 0.3402, 3.141 0.3259, 3.371 0.3119, 3.559 0.3004,
    3.732 0.2970, 3.911 0.2971, 4.085 0.2989, 4.259 0.2985, 4.424 0.2953, 4.560 0.2896,
    4.770 0.2798, 5.068 0.2661, 5.412 0.2508, 5.691 0.2381, 5.993 0.2241. Sauter radius (nm)
    at x (cm): 3.465 10.57, 3.665 11.90, 3.853 12.92, 4.023 13.75, 4.209 14.53, 4.411 15.24,
    4.604 15.81, 4.791 16.35, 4.979 16.78, 5.165 17.13, 5.377 17.44, 5.575 17.75, 5.773 18.03,
    5.970 18.19.
  - Case 257: p₀ = 67.66 kPa, T₀ = 376.7 K. p/p₀: 1.957 0.4000, 2.212 0.3897, 2.840 0.3814,
    3.446 0.3841, 3.920 0.3557, 4.480 0.3236, 5.111 0.2870. Sauter radius (nm): 2.495 16.83,
    3.895 29.24, 6.003 34.34.
  - Geometry: the wall arcs (5.3 cm into 68.6 cm) are given without the point where they meet,
    and the small throat makes boundary-layer blockage significant. As the paper recommends
    for 1D validation, the effective area comes from the dry (pre-condensation) pressures: a
    throat arc whose radius of curvature is fitted to the dry points of each case (A/A* − 1 ∝
    x² there).
  - Moore nozzle B is not used: its inlet is not fully specified and its curvature
    discontinuity sends waves across the flow that no 1D model follows.
- **Nucleation of the cryogenic fluids** (decided after the first results): classical theory
  matches steam (the project's Fig. 14 to 0.1 decade; Moses and Stein's pressures within 0.004
  to 0.021 of p₀) but needs 20–25 K of subcooling for nitrogen where the Yale nozzle shows onset
  at 7–10 K. Kalikmanov's mean-field kinetic nucleation theory (J. Chem. Phys. 124, 124505,
  2006) raises nitrogen's rate by ~11 orders but fits steam worse (pressures off by up to
  0.03). Fitting a rate factor to the Yale onsets would take ~10¹⁷ — Kalikmanov finds his
  theory 4–7 orders low for nitrogen (Iland's pulse chamber) and Sinha et al. (J. Chem. Phys.
  132, 064304, 2010) within one order for argon in a cryogenic nozzle. The 1980 onsets are
  therefore taken as not homogeneous (impurities or clusters in the supply). **Decision**:
  classical theory with Kantrowitz for water; MKNT as published (with Kantrowitz, the vapour
  being pure) for N₂, Ar and CO₂ (`CondensationModel::cryogenic`); the Yale gap kept as a test.
  A rate factor exists for calibration but is 1 by default. Sinha et al.'s argon data, if
  obtained, become the cryogenic validation.
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
