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

//! Saturation line (region 4) and the B23 and B2bc boundaries (IAPWS
//! R7-97(2012) Tables 1, 19, 34).

/// Region 4 saturation equation coefficients n₁ … n₁₀ (Table 34).
pub(crate) const SATURATION: [f64; 10] = [
    1167.0521452767,
    -724213.16703206,
    -17.073846940092,
    12020.82470247,
    -3232555.0322333,
    14.91510861353,
    -4823.2657361591,
    405113.40542057,
    -0.23855557567849,
    650.17534844798,
];

/// B23 boundary coefficients n₁ … n₅, π = n₁ + n₂θ + n₃θ² (Table 1).
pub(crate) const B23: [f64; 5] = [
    348.05185628969,
    -1.1671859879975,
    0.0010192970039326,
    572.54459862746,
    13.91883977887,
];

/// B2bc boundary coefficients n₁ … n₅, π = n₁ + n₂η + n₃η² (Table 19).
pub(crate) const B2BC: [f64; 5] = [
    905.84278514723,
    -0.67955786399241,
    0.00012809002730136,
    2652.6571908428,
    4.5257578905948,
];
