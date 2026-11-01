"""Golden dilute-gas viscosity and thermal conductivity, from CoolProp's transport correlations
at a vanishing density (1 µmol/m³, where every density-dependent term is below 1e-12).

    pip install CoolProp==8.0.0
    python3 tools/reference/transport_golden.py
"""
import CoolProp.CoolProp as CP

TEMPERATURES = {
    "Nitrogen": [70.0, 100.0, 200.0, 300.0],
    "Argon": [90.0, 150.0, 300.0],
    "Helium": [3.0, 5.0, 20.0, 100.0],
    "CarbonDioxide": [220.0, 300.0, 500.0],
    "Water": [280.0, 373.15, 500.0],
}

for fluid, temperatures in TEMPERATURES.items():
    print(f"// {fluid}: (T K, viscosity Pa s, thermal conductivity W/(m K))")
    for t in temperatures:
        mu = CP.PropsSI("V", "T", t, "Dmolar", 1e-6, fluid)
        k = CP.PropsSI("L", "T", t, "Dmolar", 1e-6, fluid)
        print(f"({t!r}, {mu!r}, {k!r}),")
