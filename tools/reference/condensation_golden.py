"""Golden values for condensation properties: surface tension (Mulero et al. 2012, as CoolProp
uses it), latent heat, saturated liquid density and vapour pressure from CoolProp's reference
equations, between the triple point and 95 % of the critical temperature.

    pip install CoolProp==8.0.0
    python3 tools/reference/condensation_golden.py
"""
import CoolProp.CoolProp as CP

FLUIDS = ["Nitrogen", "Argon", "Helium", "CarbonDioxide"]

for fluid in FLUIDS:
    t_triple = CP.PropsSI("Ttriple", fluid)
    t_crit = CP.PropsSI("Tcrit", fluid)
    print(f"// {fluid}: (T K, sigma N/m, latent heat J/kg, liquid density kg/m3, p_sat Pa)")
    for fraction in (0.05, 0.4, 0.8):
        t = t_triple + fraction * (0.95 * t_crit - t_triple)
        sigma = CP.PropsSI("I", "T", t, "Q", 0, fluid)
        h_l = CP.PropsSI("H", "T", t, "Q", 0, fluid)
        h_v = CP.PropsSI("H", "T", t, "Q", 1, fluid)
        rho_l = CP.PropsSI("D", "T", t, "Q", 0, fluid)
        p = CP.PropsSI("P", "T", t, "Q", 0, fluid)
        print(f"({t!r}, {sigma!r}, {h_v - h_l!r}, {rho_l!r}, {p!r}),")
