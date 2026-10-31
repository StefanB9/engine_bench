//! Parses the vendored NASA Glenn `data/thermo.inp` into const species tables
//! (`$OUT_DIR/gases.rs`).
//!
//! Record layout (McBride et al., `NASA/TP-2002-211556`, appendix A), columns
//! 0-based:
//! - line 1: species name (first token).
//! - line 2: interval count `[0..2]`, formula as 5 × (symbol A2, count F6.2)
//!   `[10..50]`, phase `[50..52]`, molar mass in g/mol `[52..65]`, heat of
//!   formation in J/mol `[65..80]`.
//! - per interval, line 3: `T_min [0..11]`, `T_max [11..22]`, coefficient count
//!   `[22]`, T exponents 8 × F5.1 `[23..63]`; lines 4-5: a1..a5, then a6, a7,
//!   (blank), b1, b2 as 16-character Fortran `D` floats.

use std::{env, error::Error, fmt::Write as _, fs, io::Write as _, path::PathBuf};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

const SOURCE: &str = "data/thermo.inp";
const EXPECTED_EXPONENTS: [f64; 8] = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0];
const ELEMENTS: [(&str, &str); 6] = [
    ("H", "H"),
    ("C", "C"),
    ("N", "N"),
    ("O", "O"),
    ("AR", "Ar"),
    ("HE", "He"),
];

struct Interval {
    t_min: f64,
    t_max: f64,
    a: [f64; 7],
    b: [f64; 2],
}

struct Record {
    name: String,
    elements: Vec<(&'static str, f64)>,
    molar_mass_g: f64,
    heat_of_formation: f64,
    intervals: Vec<Interval>,
}

fn main() -> Result<()> {
    writeln!(std::io::stdout(), "cargo::rerun-if-changed={SOURCE}")?;
    let text = fs::read_to_string(SOURCE).map_err(|e| format!("reading {SOURCE}: {e}"))?;
    let records = parse(&text)?;
    let out = PathBuf::from(env::var("OUT_DIR")?).join("gases.rs");
    fs::write(out, generate(&records)?)?;
    Ok(())
}

fn parse(text: &str) -> Result<Vec<Record>> {
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('!'))
        .skip_while(|l| !l.starts_with("thermo"))
        .skip(2)
        .peekable();
    let mut records = Vec::new();
    while let Some(first) = lines.next() {
        if first.starts_with("END") {
            return Ok(records);
        }
        let name = first
            .split_whitespace()
            .next()
            .ok_or("empty species line")?;
        let record = parse_record(name, &mut lines).map_err(|e| format!("species {name}: {e}"))?;
        records.push(record);
    }
    Err("missing END line".into())
}

fn parse_record<'a>(name: &str, lines: &mut impl Iterator<Item = &'a str>) -> Result<Record> {
    let header = lines.next().ok_or("missing header line")?;
    let count: usize = field(header, 0, 2)?.parse()?;
    if count == 0 {
        return Err("no temperature intervals".into());
    }
    let mut elements = Vec::new();
    for slot in 0..5 {
        let start = 10 + 8 * slot;
        let symbol = field(header, start, start + 2)?;
        let amount: f64 = field(header, start + 2, start + 8)?.parse()?;
        if amount != 0.0 {
            let variant = ELEMENTS
                .iter()
                .find(|(s, _)| *s == symbol)
                .map(|(_, v)| *v)
                .ok_or_else(|| format!("unsupported element {symbol:?}"))?;
            elements.push((variant, amount));
        }
    }
    if field(header, 50, 52)? != "0" {
        return Err("only gas-phase species are supported".into());
    }
    let molar_mass_g = field(header, 52, 65)?.parse()?;
    let heat_of_formation = field(header, 65, 80)?.parse()?;

    let mut intervals = Vec::with_capacity(count);
    for _ in 0..count {
        let range = lines.next().ok_or("missing range line")?;
        let t_min = field(range, 0, 11)?.parse()?;
        let t_max = field(range, 11, 22)?.parse()?;
        if field(range, 22, 23)? != "7" {
            return Err("expected 7 coefficients".into());
        }
        let exponents = field(range, 23, 63)?
            .split_whitespace()
            .map(str::parse)
            .collect::<std::result::Result<Vec<f64>, _>>()?;
        if exponents != EXPECTED_EXPONENTS {
            return Err(format!("unexpected exponents {exponents:?}").into());
        }
        let first = lines.next().ok_or("missing coefficient line 1")?;
        let second = lines.next().ok_or("missing coefficient line 2")?;
        let mut a = [0.0; 7];
        for (i, slot) in a.iter_mut().take(5).enumerate() {
            *slot = fortran(field(first, 16 * i, 16 * i + 16)?)?;
        }
        a[5] = fortran(field(second, 0, 16)?)?;
        a[6] = fortran(field(second, 16, 32)?)?;
        let b = [
            fortran(field(second, 48, 64)?)?,
            fortran(field(second, 64, 80)?)?,
        ];
        intervals.push(Interval { t_min, t_max, a, b });
    }
    if intervals
        .windows(2)
        .any(|w| (w[0].t_max - w[1].t_min).abs() > 1e-9)
    {
        return Err("intervals are not contiguous".into());
    }
    Ok(Record {
        name: name.to_owned(),
        elements,
        molar_mass_g,
        heat_of_formation,
        intervals,
    })
}

fn field(line: &str, start: usize, end: usize) -> Result<&str> {
    let end = end.min(line.len());
    Ok(line
        .get(start..end)
        .ok_or_else(|| format!("line too short: {line:?}"))?
        .trim())
}

fn fortran(text: &str) -> Result<f64> {
    Ok(text.replace(['D', 'd'], "E").parse()?)
}

fn const_name(species: &str) -> String {
    species.to_ascii_uppercase()
}

fn generate(records: &[Record]) -> Result<String> {
    let mut out = String::with_capacity(512 * records.len());
    writeln!(
        out,
        "// Generated by build.rs from {SOURCE}. Do not edit.\n"
    )?;
    for r in records {
        writeln!(
            out,
            "/// {} (NASA Glenn data, see `data/thermo.inp`).",
            r.name
        )?;
        writeln!(
            out,
            "pub const {}: Species = Species {{",
            const_name(&r.name)
        )?;
        writeln!(out, "    name: {:?},", r.name)?;
        writeln!(out, "    molar_mass: {:?},", r.molar_mass_g * 1e-3)?;
        let elements: Vec<String> = r
            .elements
            .iter()
            .map(|(e, n)| format!("(Element::{e}, {n:?})"))
            .collect();
        writeln!(out, "    elements: &[{}],", elements.join(", "))?;
        writeln!(out, "    enthalpy_of_formation: {:?},", r.heat_of_formation)?;
        writeln!(out, "    intervals: &[")?;
        for i in &r.intervals {
            writeln!(
                out,
                "        Interval {{ t_min: {:?}, t_max: {:?}, a: {:?}, b: {:?} }},",
                i.t_min, i.t_max, i.a, i.b
            )?;
        }
        writeln!(out, "    ],\n}};\n")?;
    }
    let all: Vec<String> = records
        .iter()
        .map(|r| format!("&{}", const_name(&r.name)))
        .collect();
    writeln!(
        out,
        "/// Every species in the vendored data, in file order."
    )?;
    writeln!(
        out,
        "pub const ALL: [&Species; {}] = [{}];",
        records.len(),
        all.join(", ")
    )?;
    Ok(out)
}
