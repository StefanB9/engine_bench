//! IAPWS-IF97 coefficient tables, generated from CoolProp/IF97 by
//! `tools/reference/gen_if97.py` (MIT, © Ian H. Bell; licence in each file).

pub(crate) mod backward;
pub(crate) mod gibbs;
pub(crate) mod region3;
pub(crate) mod saturation;

/// One term nᵢ · x^Iᵢ · y^Jᵢ of an IF97 polynomial. `I` is `i32`, except in the
/// region 2a T(p, s) backward equation, whose pressure exponents are
/// fractional.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Term<I> {
    /// Exponent of the first (pressure-like) variable.
    pub(crate) i: I,
    /// Exponent of the second (temperature-, enthalpy- or entropy-like)
    /// variable.
    pub(crate) j: i32,
    /// Coefficient.
    pub(crate) n: f64,
}

impl<I> Term<I> {
    pub(crate) const fn new(i: I, j: i32, n: f64) -> Self {
        Self { i, j, n }
    }
}

/// One subregion of the region 3 backward equation v(p, T) (IAPWS
/// SR5-05(2016)): v / v* = (Σ nᵢ ((π − a)^c)^Iᵢ ((θ − b)^d)^Jᵢ)^e, or exp(Σ nᵢ
/// (π − a)^Iᵢ (θ − b)^Jᵢ) when `exponential` (subregion 3n), with π = p / p*, θ
/// = T / T*.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Subregion {
    pub(crate) v_star: f64,
    /// Reducing pressure, Pa.
    pub(crate) p_star: f64,
    pub(crate) t_star: f64,
    pub(crate) a: f64,
    pub(crate) b: f64,
    pub(crate) c: f64,
    pub(crate) d: f64,
    pub(crate) e: f64,
    pub(crate) exponential: bool,
    pub(crate) terms: &'static [Term<i32>],
}
