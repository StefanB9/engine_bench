//! Reference equations of state transcribed from `CoolProp`'s fluid files by
//! `tools/reference/gen_fluids.py` (MIT, © Ian H. Bell and other `CoolProp`
//! developers; licence in each file).

pub(crate) mod argon;
pub(crate) mod carbon_dioxide;
pub(crate) mod helium;
pub(crate) mod nitrogen;

/// A reference Helmholtz equation of state, molar basis:
/// α = α⁰(δ, τ) + αʳ(δ, τ), δ = `ρ/ρ_c`, τ = `T_c/T`.
#[derive(Debug)]
pub(crate) struct FluidData {
    /// NASA species name of the fluid.
    pub(crate) species: &'static str,
    /// The equation's own molar gas constant, J/(mol·K).
    pub(crate) gas_constant: f64,
    /// kg/mol.
    pub(crate) molar_mass: f64,
    /// Critical (reducing) temperature, K.
    pub(crate) t_crit: f64,
    /// Critical (reducing) density, mol/m³.
    pub(crate) rho_crit: f64,
    /// K.
    pub(crate) t_triple: f64,
    /// Upper temperature of the equation, K.
    pub(crate) t_max: f64,
    /// Upper pressure of the equation, Pa.
    pub(crate) p_max: f64,
    /// α⁰ term a·ln τ.
    pub(crate) log_tau: f64,
    /// α⁰ terms n·τ^t.
    pub(crate) ideal_power: &'static [(f64, f64)],
    /// α⁰ terms n·ln(1 − exp(−θτ)), as (n, θ).
    pub(crate) planck_einstein: &'static [(f64, f64)],
    pub(crate) power: &'static [Power],
    pub(crate) gaussian: &'static [Gaussian],
    pub(crate) non_analytic: &'static [NonAnalytic],
    /// Saturated liquid density ancillary (a starting point for the
    /// equilibrium, not the equation's own value).
    pub(crate) liquid_density: Ancillary,
    pub(crate) vapour_density: Ancillary,
}

/// αʳ term n δ^d τ^t exp(−δ^l) (no exponential when l = 0).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Power {
    pub(crate) n: f64,
    pub(crate) d: i32,
    pub(crate) t: f64,
    pub(crate) l: i32,
}

/// αʳ term n δ^d τ^t exp(−η(δ − ε)² − β(τ − γ)²).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Gaussian {
    pub(crate) n: f64,
    pub(crate) d: i32,
    pub(crate) t: f64,
    pub(crate) eta: f64,
    pub(crate) epsilon: f64,
    pub(crate) beta: f64,
    pub(crate) gamma: f64,
}

/// Non-analytic critical-region term of Span and Wagner (and IAPWS-95):
/// n Δ^b δ ψ, Δ = θ² + B((δ−1)²)^a, θ = (1 − τ) + A((δ−1)²)^(1/(2β)),
/// ψ = exp(−C(δ−1)² − D(τ−1)²).
#[derive(Debug, Clone, Copy)]
pub(crate) struct NonAnalytic {
    pub(crate) n: f64,
    pub(crate) a: f64,
    pub(crate) b: f64,
    pub(crate) beta: f64,
    pub(crate) big_a: f64,
    pub(crate) big_b: f64,
    pub(crate) big_c: f64,
    pub(crate) big_d: f64,
}

/// A saturation ancillary in θ = 1 − `T/T_r`: Σ nᵢ θ^tᵢ, optionally times
/// `T_r/T`, mapped to a pressure or density by `form`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Ancillary {
    pub(crate) form: AncillaryForm,
    pub(crate) using_tau: bool,
    pub(crate) t_reducing: f64,
    pub(crate) reducing_value: f64,
    pub(crate) terms: &'static [(f64, f64)],
}

/// How an ancillary's sum maps to its value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AncillaryForm {
    /// value = reducing · (1 + sum).
    LiquidDensity,
    /// value = reducing · exp(sum).
    VapourDensity,
}
