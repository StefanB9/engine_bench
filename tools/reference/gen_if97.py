"""Transcribes IAPWS-IF97 coefficient tables from CoolProp/IF97 `IF97.h` into Rust consts."""
import re, sys

src = open(sys.argv[1]).read()
out_dir = sys.argv[2]

def table(name):
    m = re.search(r"static const \w+ " + name + r"\[\] = \{(.*?)\n\s*\};", src, re.S)
    body = re.sub(r"//[^\n]*", "", m.group(1))
    return [tuple(x.strip() for x in t.split(",")) for t in re.findall(r"\{([^{}]*)\}", body)]

def flat(name):
    m = re.search(r"static const double " + name + r"\[\] = \{(.*?)\};", src, re.S)
    return [x.strip() for x in re.sub(r"//[^\n]*", "", m.group(1)).split(",") if x.strip()]

def num(s):
    v = float(s)
    r = repr(v)
    if "e" in r or "E" in r:
        mant, exp = r.split("e")
        if "." not in mant: mant += ".0"
        return f"{mant}e{int(exp)}"
    return r if "." in r else r + ".0"

HEADER = """// IAPWS-IF97 coefficients, transcribed by `tools/reference/gen_if97.py` from
// CoolProp/IF97 `IF97.h` v2.4.0 (https://github.com/CoolProp/IF97).
// Copyright (C) 2015 Ian H. Bell, MIT licence: permission is hereby granted, free of charge, to
// any person obtaining a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
// Software, and to permit persons to whom the Software is furnished to do so, subject to the
// following conditions: The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
// WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.
"""

def residual(rust, doc, rows, ifloat=False, sign=lambda i: 1):
    ity = "f64" if ifloat else "i32"
    lines = [f"/// {doc}", f"pub(crate) const {rust}: [Term<{ity}>; {len(rows)}] = ["]
    for i, j, n in rows:
        ii = num(i) if ifloat else str(int(float(i)))
        nn = float(n) * sign(int(float(i)) if not ifloat else 0)
        lines.append(f"    Term::new({ii}, {int(float(j))}, {num(repr(nn))}),")
    lines.append("];")
    return "\n".join(lines)

def ideal(rust, doc, rows):
    lines = [f"/// {doc}", f"pub(crate) const {rust}: [(i32, f64); {len(rows)}] = ["]
    for j, n in rows:
        lines.append(f"    ({int(j)}, {num(n)}),")
    lines.append("];")
    return "\n".join(lines)

def consts(rust, doc, vals):
    return "\n".join([f"/// {doc}", f"pub(crate) const {rust}: [f64; {len(vals)}] = ["]
                     + [f"    {num(v)}," for v in vals] + ["];"])

gibbs = [HEADER, "//! Gibbs free energy coefficients of regions 1, 2 and 5 (IAPWS R7-97(2012) Tables 2, 10, 11, 37, 38).",
  "", "use super::Term;", "",
  # CoolProp folds (-1)^I into n to write (π − 7.1)^I; undo it to match IAPWS Table 2 term for term.
  residual("REGION1", "Region 1: γ = Σ nᵢ (7.1 − π)^Iᵢ (τ − 1.222)^Jᵢ (Table 2).", table("Region1residdata"), sign=lambda i: (-1) ** i),
  "", residual("REGION2_RESIDUAL", "Region 2 residual part: γʳ = Σ nᵢ π^Iᵢ (τ − 0.5)^Jᵢ (Table 11).", table("Region2residdata")),
  "", ideal("REGION2_IDEAL", "Region 2 ideal-gas part: γ° = ln π + Σ n°ᵢ τ^J°ᵢ (Table 10).", table("Region2idealdata")),
  "", residual("REGION5_RESIDUAL", "Region 5 residual part: γʳ = Σ nᵢ π^Iᵢ τ^Jᵢ (Table 38).", table("Region5residdata")),
  "", ideal("REGION5_IDEAL", "Region 5 ideal-gas part: γ° = ln π + Σ n°ᵢ τ^J°ᵢ (Table 37).", table("Region5idealdata")),
  ""]
open(f"{out_dir}/gibbs.rs", "w").write("\n".join(gibbs))

sat = [HEADER, "//! Saturation line (region 4) and the B23 and B2bc boundaries (IAPWS R7-97(2012) Tables 1, 19, 34).",
  "", consts("SATURATION", "Region 4 saturation equation coefficients n₁ … n₁₀ (Table 34).", [n for _, n in table("sat")]),
  "", consts("B23", "B23 boundary coefficients n₁ … n₅, π = n₁ + n₂θ + n₃θ² (Table 1).", flat("Region23data")),
  "", consts("B2BC", "B2bc boundary coefficients n₁ … n₅, π = n₁ + n₂η + n₃η² (Table 19).", flat("Region2b2cdata")),
  ""]
open(f"{out_dir}/saturation.rs", "w").write("\n".join(sat))

bk = [HEADER, "//! Backward equations T(p, h) and T(p, s) for regions 1 and 2 (IAPWS R7-97(2012) Tables 6, 8, 20–23, 25–27).",
  "", "use super::Term;", "",
  residual("REGION1_PH", "Region 1 T(p, h): θ = Σ nᵢ π^Iᵢ (η + 1)^Jᵢ (Table 6).", table("Coeff1H")),
  "", residual("REGION1_PS", "Region 1 T(p, s): θ = Σ nᵢ π^Iᵢ (σ + 2)^Jᵢ (Table 8).", table("Coeff1S")),
  "", residual("REGION2A_PH", "Region 2a T(p, h): θ = Σ nᵢ π^Iᵢ (η − 2.1)^Jᵢ (Table 20).", table("Coeff2aH")),
  "", residual("REGION2B_PH", "Region 2b T(p, h): θ = Σ nᵢ (π − 2)^Iᵢ (η − 2.6)^Jᵢ (Table 21).", table("Coeff2bH")),
  "", residual("REGION2C_PH", "Region 2c T(p, h): θ = Σ nᵢ (π + 25)^Iᵢ (η − 1.8)^Jᵢ (Table 22).", table("Coeff2cH")),
  "", residual("REGION2A_PS", "Region 2a T(p, s): θ = Σ nᵢ π^Iᵢ (σ − 2)^Jᵢ, fractional Iᵢ (Table 25).", table("Coeff2aS"), ifloat=True),
  "", residual("REGION2B_PS", "Region 2b T(p, s): θ = Σ nᵢ π^Iᵢ (10 − σ)^Jᵢ (Table 26).", table("Coeff2bS")),
  "", residual("REGION2C_PS", "Region 2c T(p, s): θ = Σ nᵢ π^Iᵢ (2 − σ)^Jᵢ (Table 27).", table("Coeff2cS")),
  ""]
open(f"{out_dir}/backward.rs", "w").write("\n".join(bk))
print("ok")

# --- Region 3: Helmholtz terms, SR5-05 v(p, T) subregions and their boundary curves ---
helm = table("Region3residdata")
n1 = helm[0][2]
r3 = [HEADER, "//! Region 3: Helmholtz free energy (IAPWS R7-97(2012) Table 30) and the backward",
      "//! equations v(p, T) of IAPWS SR5-05(2016) with their subregion boundaries (Tables 2, 3, 4, 6–9).",
      "", "use super::{Subregion, Term};", "",
      f"/// n₁ of φ = n₁ ln δ + Σ nᵢ δ^Iᵢ τ^Jᵢ (Table 30).\npub(crate) const HELMHOLTZ_LN: f64 = {num(n1)};", "",
      residual("HELMHOLTZ", "Terms 2 … 40 of φ = n₁ ln δ + Σ nᵢ δ^Iᵢ τ^Jᵢ (Table 30).", helm[1:]), ""]
params = {}
for m in re.finditer(r"Region3(\w)\(\) : Region3BackwardsRegion\(Region3\wdata, (\d+)\)\{\s*(.*?)\};", src, re.S):
    letter, body = m.group(1), m.group(3)
    kv = dict(re.findall(r"(\w+) = ([\d.]+)(?:\*p_fact)?", body))
    params[letter.upper()] = kv
subs = []
for letter in "abcdefghijklmnopqrstuvwxyz":
    L = letter.upper()
    rows = table(f"Region3{L}data")
    kv = params[L]
    name = f"TERMS_3{L}"
    r3 += [residual(name, f"Subregion 3{letter} (SR5-05 Tables 2 and 4–9).", rows), ""]
    exp_form = L == "N"
    c, d, e = (kv.get("c", "1"), kv.get("d", "1"), kv.get("e", "1"))
    subs.append(f"    // 3{letter}\n    Subregion {{ v_star: {num(kv['v_star'])}, p_star: {num(float(kv['p_star'])*1e6)}, "
                f"t_star: {num(kv['T_star'])}, a: {num(kv['a'])}, b: {num(kv['b'])}, c: {num(c)}, d: {num(d)}, "
                f"e: {num(e)}, exponential: {'true' if exp_form else 'false'}, terms: &{name} }},")
r3 += ["/// The 26 subregions 3a … 3z, in order.", f"pub(crate) const SUBREGIONS: [Subregion; 26] = ["] + subs + ["];", ""]
for line in ["AB", "CD", "GH", "IJ", "JK", "MN", "OP", "QU", "RX", "UV", "WX"]:
    rows = table(f"{line}data")
    log_form = line in ("AB", "OP", "WX")
    var = "ln π" if log_form else "π"
    r3 += [ideal(f"LINE_{line}", f"Boundary {line[0].lower()}{line[1].lower()}: T = Σ nᵢ ({var})^Iᵢ K (SR5-05 Table 3).", rows), ""]
open(f"{out_dir}/region3.rs", "w").write("\n".join(r3))
print("region 3 ok", len(helm), sum(len(table(f'Region3{L}data')) for L in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
