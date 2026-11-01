// IAPWS-IF97 coefficients, transcribed by `tools/reference/gen_if97.py` from
// CoolProp/IF97 `IF97.h` v2.4.0 (https://github.com/CoolProp/IF97).
// Copyright (C) 2015 Ian H. Bell, MIT licence: permission is hereby granted,
// free of charge, to any person obtaining a copy of this software and
// associated documentation files (the "Software"), to deal in the Software
// without restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies of the
// Software, and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions: The above copyright notice and this
// permission notice shall be included in all copies or substantial portions of
// the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED.

//! Gibbs free energy coefficients of regions 1, 2 and 5 (IAPWS R7-97(2012)
//! Tables 2, 10, 11, 37, 38).

use super::Term;

/// Region 1: γ = Σ nᵢ (7.1 − π)^Iᵢ (τ − 1.222)^Jᵢ (Table 2).
pub(crate) const REGION1: [Term<i32>; 34] = [
    Term::new(0, -2, 0.14632971213167),
    Term::new(0, -1, -0.84548187169114),
    Term::new(0, 0, -3.756360367204),
    Term::new(0, 1, 3.3855169168385),
    Term::new(0, 2, -0.95791963387872),
    Term::new(0, 3, 0.15772038513228),
    Term::new(0, 4, -0.016616417199501),
    Term::new(0, 5, 0.00081214629983568),
    Term::new(1, -9, 0.00028319080123804),
    Term::new(1, -7, -0.00060706301565874),
    Term::new(1, -1, -0.018990068218419),
    Term::new(1, 0, -0.032529748770505),
    Term::new(1, 1, -0.021841717175414),
    Term::new(1, 3, -5.283835796993e-5),
    Term::new(2, -3, -0.00047184321073267),
    Term::new(2, 0, -0.00030001780793026),
    Term::new(2, 1, 4.7661393906987e-5),
    Term::new(2, 3, -4.4141845330846e-6),
    Term::new(2, 17, -7.2694996297594e-16),
    Term::new(3, -4, -3.1679644845054e-5),
    Term::new(3, 0, -2.8270797985312e-6),
    Term::new(3, 6, -8.5205128120103e-10),
    Term::new(4, -5, -2.2425281908e-6),
    Term::new(4, -2, -6.5171222895601e-7),
    Term::new(4, 10, -1.4341729937924e-13),
    Term::new(5, -8, -4.0516996860117e-7),
    Term::new(8, -11, -1.2734301741641e-9),
    Term::new(8, -6, -1.7424871230634e-10),
    Term::new(21, -29, -6.8762131295531e-19),
    Term::new(23, -31, 1.4478307828521e-20),
    Term::new(29, -38, 2.6335781662795e-23),
    Term::new(30, -39, -1.1947622640071e-23),
    Term::new(31, -40, 1.8228094581404e-24),
    Term::new(32, -41, -9.3537087292458e-26),
];

/// Region 2 residual part: γʳ = Σ nᵢ π^Iᵢ (τ − 0.5)^Jᵢ (Table 11).
pub(crate) const REGION2_RESIDUAL: [Term<i32>; 43] = [
    Term::new(1, 0, -0.0017731742473213),
    Term::new(1, 1, -0.017834862292358),
    Term::new(1, 2, -0.045996013696365),
    Term::new(1, 3, -0.057581259083432),
    Term::new(1, 6, -0.05032527872793),
    Term::new(2, 1, -3.3032641670203e-5),
    Term::new(2, 2, -0.00018948987516315),
    Term::new(2, 4, -0.0039392777243355),
    Term::new(2, 7, -0.043797295650573),
    Term::new(2, 36, -2.6674547914087e-5),
    Term::new(3, 0, 2.0481737692309e-8),
    Term::new(3, 1, 4.3870667284435e-7),
    Term::new(3, 3, -3.227767723857e-5),
    Term::new(3, 6, -0.0015033924542148),
    Term::new(3, 35, -0.040668253562649),
    Term::new(4, 1, -7.8847309559367e-10),
    Term::new(4, 2, 1.2790717852285e-8),
    Term::new(4, 3, 4.8225372718507e-7),
    Term::new(5, 7, 2.2922076337661e-6),
    Term::new(6, 3, -1.6714766451061e-11),
    Term::new(6, 16, -0.0021171472321355),
    Term::new(6, 35, -23.895741934104),
    Term::new(7, 0, -5.905956432427e-18),
    Term::new(7, 11, -1.2621808899101e-6),
    Term::new(7, 25, -0.038946842435739),
    Term::new(8, 8, 1.1256211360459e-11),
    Term::new(8, 36, -8.2311340897998),
    Term::new(9, 13, 1.9809712802088e-8),
    Term::new(10, 4, 1.0406965210174e-19),
    Term::new(10, 10, -1.0234747095929e-13),
    Term::new(10, 14, -1.0018179379511e-9),
    Term::new(16, 29, -8.0882908646985e-11),
    Term::new(16, 50, 0.10693031879409),
    Term::new(18, 57, -0.33662250574171),
    Term::new(20, 20, 8.9185845355421e-25),
    Term::new(20, 35, 3.0629316876232e-13),
    Term::new(20, 48, -4.2002467698208e-6),
    Term::new(21, 21, -5.9056029685639e-26),
    Term::new(22, 53, 3.7826947613457e-6),
    Term::new(23, 39, -1.2768608934681e-15),
    Term::new(24, 26, 7.3087610595061e-29),
    Term::new(24, 40, 5.5414715350778e-17),
    Term::new(24, 58, -9.436970724121e-7),
];

/// Region 2 ideal-gas part: γ° = ln π + Σ n°ᵢ τ^J°ᵢ (Table 10).
pub(crate) const REGION2_IDEAL: [(i32, f64); 9] = [
    (0, -9.6927686500217),
    (1, 10.086655968018),
    (-5, -0.005608791128302),
    (-4, 0.071452738081455),
    (-3, -0.40710498223928),
    (-2, 1.4240819171444),
    (-1, -4.383951131945),
    (2, -0.28408632460772),
    (3, 0.021268463753307),
];

/// Ideal-gas part of the metastable-vapour equation, Eq. (18): Table 10
/// with n°₁ and n°₂ readjusted (its footnote).
pub(crate) const METASTABLE_IDEAL: [(i32, f64); 9] = [
    (0, -9.6937268393049),
    (1, 10.087275970006),
    (-5, -0.005608791128302),
    (-4, 0.071452738081455),
    (-3, -0.40710498223928),
    (-2, 1.4240819171444),
    (-1, -4.383951131945),
    (2, -0.28408632460772),
    (3, 0.021268463753307),
];

/// Residual part of the metastable-vapour equation: γʳ = Σ nᵢ π^Iᵢ
/// (τ − 0.5)^Jᵢ (Table 16).
pub(crate) const METASTABLE_RESIDUAL: [Term<i32>; 13] = [
    Term {
        i: 1,
        j: 0,
        n: -0.73362260186506e-2,
    },
    Term {
        i: 1,
        j: 2,
        n: -0.88223831943146e-1,
    },
    Term {
        i: 1,
        j: 5,
        n: -0.72334555213245e-1,
    },
    Term {
        i: 1,
        j: 11,
        n: -0.40813178534455e-2,
    },
    Term {
        i: 2,
        j: 1,
        n: 0.20097803380207e-2,
    },
    Term {
        i: 2,
        j: 7,
        n: -0.53045921898642e-1,
    },
    Term {
        i: 2,
        j: 16,
        n: -0.76190409086970e-2,
    },
    Term {
        i: 3,
        j: 4,
        n: -0.63498037657313e-2,
    },
    Term {
        i: 3,
        j: 16,
        n: -0.86043093028588e-1,
    },
    Term {
        i: 4,
        j: 7,
        n: 0.75321581522770e-2,
    },
    Term {
        i: 4,
        j: 10,
        n: -0.79238375446139e-2,
    },
    Term {
        i: 5,
        j: 9,
        n: -0.22888160778447e-3,
    },
    Term {
        i: 5,
        j: 10,
        n: -0.26456501482810e-2,
    },
];

/// Region 5 residual part: γʳ = Σ nᵢ π^Iᵢ τ^Jᵢ (Table 38).
pub(crate) const REGION5_RESIDUAL: [Term<i32>; 6] = [
    Term::new(1, 1, 0.0015736404855259),
    Term::new(1, 2, 0.00090153761673944),
    Term::new(1, 3, -0.0050270077677648),
    Term::new(2, 3, 2.2440037409485e-6),
    Term::new(2, 9, -4.1163275453471e-6),
    Term::new(3, 7, 3.7919454822955e-8),
];

/// Region 5 ideal-gas part: γ° = ln π + Σ n°ᵢ τ^J°ᵢ (Table 37).
pub(crate) const REGION5_IDEAL: [(i32, f64); 6] = [
    (0, -13.179983674201),
    (1, 6.8540841634434),
    (-3, -0.024805148933466),
    (-2, 0.36901534980333),
    (-1, -3.1161318213925),
    (2, -0.32961626538917),
];
