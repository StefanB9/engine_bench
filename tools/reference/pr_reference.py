"""Critical constants and golden compressibilities for the Peng-Robinson EOS
(`crates/engine-bench-physics/src/gas/eos/peng_robinson`), from CoolProp's reference
equations of state.

    venv/bin/python tools/reference/pr_reference.py

CoolProp (MIT) is used offline only; it is not a build dependency.
"""
import CoolProp
import CoolProp.CoolProp as CP

SPECIES = [  # NASA name, CoolProp fluid
    ("N2", "Nitrogen"), ("O2", "Oxygen"), ("Ar", "Argon"), ("He", "Helium"),
    ("CO2", "CarbonDioxide"), ("H2O", "Water"), ("H2", "Hydrogen"), ("CO", "CarbonMonoxide"),
]


def main():
    print(f"// Critical constants from CoolProp {CoolProp.__version__} reference EOS.")
    for name, fluid in SPECIES:
        tc, pc, w = (CP.PropsSI(k, fluid) for k in ("Tcrit", "pcrit", "acentric"))
        print(f'    ("{name}", Critical {{ temperature: {tc!r}, pressure: {pc!r}, acentric_factor: {w!r} }}),')
    print()
    print("// Reference Z: species, T K, p Pa, Z.")
    for name, fluid in SPECIES[:5]:
        for t, p in ((300.0, 1e7), (300.0, 3e7), (200.0, 1e7), (150.0, 3e7)):
            if fluid == "CarbonDioxide":
                t, p = {(300.0, 1e7): (300.0, 5e6), (300.0, 3e7): (350.0, 2e7),
                        (200.0, 1e7): (400.0, 1e7), (150.0, 3e7): (350.0, 5e6)}[(t, p)]
            z = CP.PropsSI("Z", "T", t, "P", p, fluid)
            print(f'    ("{name}", {t!r}, {p!r}, {z!r}),')


if __name__ == "__main__":
    main()
