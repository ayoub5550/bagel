#!/usr/bin/env python3
"""Techno-economic model (TEA) for acid-stable C-phycocyanin from heterotrophic
Galdieria sulphuraria, fed with sugar from low-grade Algerian dates.

Every parameter carries a reference key from docs/references.md, "computed",
or UNVERIFIED (engineering assumption to be replaced by a quote/measurement).

Outputs (all regenerated on each run, never edit by hand):
  data/01_parameters.tsv        every distribution used, with its source key
  data/02_route_comparison.tsv  Monte Carlo cost per kg PC for 6 process routes x 4 scales
  data/03_base_breakdown.tsv    deterministic cost breakdown (mode values), recommended route
  data/04_dcf.tsv               NPV / IRR / P(NPV>0) at 3 price anchors, recommended route
  data/05_tornado.tsv           one-at-a-time sensitivity, recommended route, 100 t/yr
  data/06_lab_targets.tsv       unit cost over a grid of PC content x biomass productivity
  data/07_scenarios.tsv         deterministic scenarios + break-even price of low-grade dates vs glucose
  data/fig_*.png                figures

Run:  python projects/galdieria-blue/scripts/01_tea_model.py   (numpy, pandas, matplotlib)
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)
SEED = 20261007
N = 100_000                      # Monte Carlo draws per route x scale
SCALES_T = [25, 50, 100, 200]    # t of pure C-phycocyanin (PC) sold per year
FX_DZD_PER_USD = 134.40          # 2026-10-07 (repo FX note in HANDOFF)
PRODUCT_PC_FRACTION = 0.34       # typical PC in galdieria extract blue (fr2025_galdieria)
E18_PC_FRACTION = 0.24           # E18 grade ~24 % C-PC on carrier (tint2026)
CRF_YEARS = 15                   # plant life (ruiz2022 uses 15 y)
TAX_RATE = 0.19                  # IBS, production of goods (dgi_ibs)

# --------------------------------------------------------------------------------------
# Price anchors, $/kg of pure PC (computed from spirulina-blue market quotes)
# --------------------------------------------------------------------------------------
PRICE_ANCHORS = {
    "low_china_E18_bulk": 32.0 / E18_PC_FRACTION,            # $32/kg E18, 25 kg MOQ (botanicalcube2026)
    "mid_china_E25": 115.0 / (0.24 + (25 - 18) / (40 - 18) * (0.60 - 0.24)),  # $115/kg E25 (botanicalcube2026), PC% interpolated E18 24% -> E40 60% (tint2026)
    "high_eu_E18_250kg": 0.69e-3 * 180_000 / E18_PC_FRACTION,  # $0.69 per 1000 CU, E18 = 180 CU/g (tint2026)
}

# --------------------------------------------------------------------------------------
# Parameters: (low, mode, high) triangular, unit, source.  Route overrides below.
# --------------------------------------------------------------------------------------
P = {
    # biology
    "Px":        ((15, 25, 40), "g biomass/L/d at scale", "graverholt2007 17.5 fed-batch; burns2020 31.8-42.8 (3.5 L); scale-down UNVERIFIED"),
    "c_pc":      ((15, 25, 40), "mg PC/g biomass", "graverholt2007 26.7 (glucose, fed-batch); burns2020 31.4 glucose fed-batch, 29.9 glycerol 42C, max 45"),
    "Y_glc":     ((0.41, 0.48, 0.55), "g biomass/g sugar", "schmidt2005 0.48-0.50; graverholt2007 0.41-0.43"),
    "f_date_Y":  ((0.80, 0.95, 1.05), "x of glucose yield on date syrup", "UNVERIFIED: stage-1 gate is >=0.80 (no study on dates)"),
    "f_date_c":  ((0.70, 0.95, 1.10), "x of glucose PC content on date syrup", "UNVERIFIED: fructose/impurities may repress pigment (sloth2006, perezsaura2022)"),
    "uptime":    ((0.85, 0.90, 0.93), "fraction of 365 d", "ruiz2022 330 d/y = 0.90"),
    "R_dsp":     ((0.50, 0.65, 0.80), "PC recovered / PC in broth", "sorensen2013 (purification works); overall loss UNVERIFIED"),
    # feedstock
    "date_dzd":  ((15, 30, 50), "DZD/kg low-grade dates (rebuts)", "UNVERIFIED: below processed date feed 68-72 DZD/kg (tsa2026); common dates 50-250 DZD/kg (ouargla_market)"),
    "s_frac":    ((0.45, 0.55, 0.64), "kg sugar/kg low-grade dates", "scirp2021_hchef 45.3%; ifstj2021 51-63.6%"),
    "eta_ext":   ((0.85, 0.90, 0.95), "sugar extraction yield", "UNVERIFIED"),
    "proc_date": ((0.02, 0.05, 0.10), "$/kg sugar to extract+clarify", "UNVERIFIED"),
    "p_glc":     ((0.60, 0.70, 0.90), "$/kg glucose landed Algeria", "procurement2026 China FOB $476/t + 30% duty (dz_tarif_1702) + freight UNVERIFIED"),
    # nutrients & utilities
    "n_frac":    ((0.075, 0.09, 0.12), "g N/g biomass", "burns2020 Yx/N 11.7; abiusi2024 7.8-12% N"),
    "p_nh3":     ((0.45, 0.55, 0.75), "$/kg NH3 delivered", "gtaic2026 ~516 $/t import proxy"),
    "nut_other": ((0.02, 0.04, 0.08), "$/kg biomass (P,K,Mg,trace,acid,antifoam)", "UNVERIFIED"),
    "e_O2":      ((0.6, 1.0, 1.8), "kWh/kg O2 transferred (incl. VPSA)", "vpsa 0.3-0.6 kWh/Nm3 (doer2026); transfer UNVERIFIED"),
    "e_misc":    ((0.3, 0.5, 0.9), "kWh/kg biomass (pumps, DSP, cooling fans)", "UNVERIFIED"),
    "p_el":      ((0.035, 0.045, 0.08), "$/kWh", "gpp2025 business 4.68 DZD = $0.035; upside = subsidy removal"),
    "p_gas":     ((0.003, 0.01, 0.03), "$/kWh thermal", "gpp2025 gas $0.003; 2025 large-user repricing (dz_gas2025)"),
    # downstream
    "dsp_cons":  ((8, 15, 30), "$/kg PC (enzyme, membranes, filter aid, QC, pack)", "UNVERIFIED"),
    "p_carrier": ((0.6, 1.2, 2.0), "$/kg carrier", "trehalose2026 $0.65-1.8/kg"),
    "logistics": ((1, 2, 4), "$/kg product to EU/US", "UNVERIFIED"),
    "v_bio":     ((0.10, 0.20, 0.30), "$/kg spent biomass, net of drying", "tridge2026 soybean meal DZ import $0.29-0.34/kg as ceiling"),
    # capex
    "C200":      ((0.6, 1.0, 1.7), "M$ per 200 m3-working vessel", "ruiz2022 EUR397k/250.7 m3 (2020) x cepci; lever2023 $1.7M"),
    "f_vessel":  ((0.7, 0.9, 1.0), "x vessel cost (acid, non-autoclaved, lined)", "UNVERIFIED; burns2020 plastics, no autoclave"),
    "vpsa_cost": ((600, 1200, 2500), "$ per Nm3/h O2", "doer2026 / ecer quotes 300-10,000 Nm3/h"),
    "f_O2":      ((0.3, 0.5, 0.8), "share of O2 from VPSA", "UNVERIFIED"),
    "dsp_m3":    ((4000, 7000, 10000), "$ DSP equipment per m3 working volume (at 800 m3)", "cxbio2025 3,000-10,000 $/m3"),
    "date_plant":((1.0, 2.0, 3.5), "M$ date-syrup plant at 20 kt sugar/y", "UNVERIFIED"),
    "lang":      ((2.8, 3.3, 4.5), "installed/equipment", "ruiz2022 factors = 3.24; lever2023 3.5 (cites Humbird 4.5); cxbio2025 2-4"),
    "foak":      ((1.0, 1.2, 1.6), "first-of-a-kind overrun", "UNVERIFIED; cxbio2025 real plants 74k-734k $/m3"),
    "maint":     ((0.015, 0.03, 0.04), "fraction of FCI per year", "ruiz2022 4% of equipment"),
    # labour & fixed
    "fte_cost":  ((9000, 14000, 22000), "$/FTE/y loaded", "almohandiss2026 45-110k DZD/month; CNAS employer 26% (cnas2025)"),
    "expat":     ((0.3, 0.6, 1.0), "M$/y foreign technical staff", "UNVERIFIED"),
    "fixed_qa":  ((0.5, 1.0, 2.0), "M$/y QA, regulatory, sales", "UNVERIFIED"),
    "rate":      ((0.10, 0.12, 0.15), "discount rate", "UNVERIFIED (Algeria risk premium)"),
}

# Process routes: overrides of the base distributions + design flags.
ROUTES = {
    "R0_eu_sterile_str_glucose": dict(
        label="Reference: EU, sterile stirred tank, glucose (incumbent-like)",
        sugar="glucose", o2_enriched=False, eu=True,
        over={"p_glc": (0.45, 0.50, 0.60), "p_el": (0.10, 0.12, 0.16), "p_gas": (0.03, 0.05, 0.08),
              "f_vessel": (1.07, 1.5, 2.0), "fte_cost": (45000, 55000, 70000), "expat": (0, 0, 0.01),
              "rate": (0.07, 0.08, 0.10), "logistics": (0.5, 1, 2), "e_O2": (1.0, 1.6, 2.5),
              "Px": (12, 18, 30), "c_pc": (15, 25, 35)},
        sources="ruiz2022 (EUR 0.44/kg glucose, EUR 0.122/kWh, salaries); f_vessel: agitator EUR576k (ruiz2022)"),
    "R1_naive_batch_dates": dict(
        label="Naive: batch on date syrup, sugar in excess, air only",
        sugar="date", o2_enriched=False, eu=False,
        over={"Px": (8, 15, 25), "c_pc": (3, 5, 11), "e_O2": (1.0, 1.6, 2.5)},
        sources="schmidt2005 3-4 mg/g; sloth2006 4-5 mg/g on sugars; molasses 11.2 (burns2020)"),
    "R2_climited_glucose": dict(
        label="Recommended design on imported glucose",
        sugar="glucose", o2_enriched=True, eu=False, over={}, sources="see base"),
    "R3_climited_dates_fedbatch": dict(
        label="RECOMMENDED: carbon-limited fed-batch on date syrup, high DO (O2-enriched), 42C, pH 2",
        sugar="date", o2_enriched=True, eu=False, over={}, sources="see base"),
    "R4_climited_dates_continuous": dict(
        label="Same, semi-continuous / chemostat cycles",
        sugar="date", o2_enriched=True, eu=False,
        over={"Px": (20, 32, 45), "c_pc": (15, 25, 35), "uptime": (0.80, 0.87, 0.92)},
        sources="burns2020 continuous 42.8 g/L/d, 29.9 mg/g (glycerol); revertants lower PC; pleissner2025_gald 40 d non-sterile"),
    "R5_dates_with_biomass_credit": dict(
        label="R3 + spent biomass sold as feed protein",
        sugar="date", o2_enriched=True, eu=False, over={}, credit=True, sources="v_bio"),
}
RECOMMENDED = "R3_climited_dates_fedbatch"


def tri(rng, lo, mode, hi, n):
    if hi <= lo:
        return np.full(n, float(mode))
    return rng.triangular(lo, mode, hi, n)


def draw(route: dict, rng, n: int, mode_only: bool = False, fix: dict | None = None) -> dict:
    """Draw every parameter for a route; mode_only returns the modes (deterministic base)."""
    out = {}
    for k, (rngdef, _u, _s) in P.items():
        lo, mo, hi = route["over"].get(k, rngdef)
        out[k] = np.full(n, float(mo)) if mode_only else tri(rng, lo, mo, hi, n)
    for k, v in (fix or {}).items():
        out[k] = np.full(n, float(v))
    return out


def o2_per_biomass(Y_g):
    """kg O2 per kg biomass from degree-of-reduction balance (computed).
    Sugar CH2O (30 g/Cmol, gamma 4); biomass CH1.8O0.5N0.2 (24.6 g/Cmol, gamma 4.2, NH3 as N)."""
    Y_cmol = Y_g * 30.0 / 24.6
    o2_mol_per_cmol_s = (4.0 - Y_cmol * 4.2) / 4.0
    return o2_mol_per_cmol_s * 32.0 / 30.0 / Y_g


def crf(r, n=CRF_YEARS):
    return r * (1 + r) ** n / ((1 + r) ** n - 1)


def unit_cost(route: dict, p: dict, P_final_kg: float) -> dict:
    """Annual mass balance, capex and opex at full capacity -> $/kg PC (vectorised)."""
    date = route["sugar"] == "date"
    c = p["c_pc"] / 1000.0 * (p["f_date_c"] if date else 1.0)
    Y = p["Y_glc"] * (p["f_date_Y"] if date else 1.0)
    days = 365.0 * p["uptime"]
    pc_broth = P_final_kg / p["R_dsp"]
    X = pc_broth / c                                   # kg biomass / y
    Vw = X / (p["Px"] * days)                          # m3 working volume (Px g/L/d = kg/m3/d)
    sugar = X / Y                                      # kg sugar / y
    o2pb = o2_per_biomass(Y)
    O2 = X * o2pb                                      # kg O2 / y
    otr = p["Px"] * o2pb * 1000 / 32 / 24              # mmol O2 / L / h
    heat_MW = O2 * 14.4e6 / (days * 86400) / 1e6       # ~14.4 MJ/kg O2 oxycaloric equivalent (cooney1969 correlation; coefficient UNVERIFIED)
    evap_m3 = O2 * 14.4 / 2.4 / 1000                   # cooling-tower evaporation, m3/y (latent ~2.4 MJ/kg)
    product = P_final_kg / PRODUCT_PC_FRACTION

    # sugar price
    if date:
        p_sugar = (p["date_dzd"] / FX_DZD_PER_USD) / (p["s_frac"] * p["eta_ext"]) + p["proc_date"]
        dates_t = sugar / (p["s_frac"] * p["eta_ext"]) / 1000
    else:
        p_sugar = p["p_glc"]
        dates_t = np.zeros_like(sugar)

    # ---- capex (bottom-up, M$) ----
    v_size = np.clip(Vw / 2.0, 20.0, 200.0)
    n_vessels = np.ceil(Vw / v_size)
    ferm = n_vessels * p["C200"] * (v_size / 200.0) ** 0.6 * p["f_vessel"]
    seed = 0.15 * ferm                                 # lever2023: seed train 1.37 / (6 x 1.7) = 0.13
    compressor = 0.23 * n_vessels * (v_size / 200.0) ** 0.6   # lever2023 $1.375M for 6 x 200 m3
    if route["o2_enriched"]:
        nm3h = O2 / (days * 24) / 1.429 / 0.6 * p["f_O2"]       # 60 % O2 utilisation UNVERIFIED
        vpsa = nm3h * p["vpsa_cost"] / 1e6
    else:
        vpsa = np.zeros_like(ferm)
    dsp = p["dsp_m3"] * 800.0 * (Vw / 800.0) ** 0.7 / 1e6
    feed_plant = (p["date_plant"] if date else 0.3) * (sugar / 2.0e7) ** 0.6
    mec = ferm + seed + compressor + vpsa + dsp + feed_plant
    fci = mec * p["lang"] * p["foak"]                  # M$

    # ---- opex ($/y) ----
    c_sugar = sugar * p_sugar
    c_nut = X * (p["n_frac"] * 17.0 / 14.0 * p["p_nh3"] + p["nut_other"])
    kwh = O2 * p["e_O2"] + X * p["e_misc"]
    c_el = kwh * p["p_el"]
    th_kwh = (product * 4.0 * 4.5 / 3.6                  # spray drying ~4 kg water/kg product, 4.5 MJ/kg
              + (sugar * 3.0 * 2.6 / 2.5 / 3.6 if date else 0)   # syrup evaporation, 3-effect
              + (X * (1 - c) * 3.0 * 3.0 / 3.6 if route.get("credit") else 0))
    c_th = th_kwh * p["p_gas"]
    c_dsp = P_final_kg * p["dsp_cons"] + product * 0.40 * p["p_carrier"] + product * p["logistics"]
    fte = 30 + 30 * (Vw / 1000.0) ** 0.6
    c_lab = fte * p["fte_cost"] + p["expat"] * 1e6
    c_maint = p["maint"] * fci * 1e6
    c_ins = 0.012 * fci * 1e6
    c_ovh = 0.55 * (c_lab + c_maint)
    c_fixed = p["fixed_qa"] * 1e6
    c_cap = crf(p["rate"]) * fci * 1e6 * 1.10          # +10 % start-up & working capital
    credit = X * (1 - c) * p["v_bio"] if route.get("credit") else np.zeros_like(X)

    parts = {"sugar": c_sugar, "nutrients": c_nut, "electricity": c_el, "thermal": c_th,
             "dsp_carrier_logistics": c_dsp, "labour": c_lab, "maintenance": c_maint,
             "insurance": c_ins, "overheads": c_ovh, "qa_sales_fixed": c_fixed,
             "capital_charge": c_cap, "biomass_credit": -credit}
    total = sum(parts.values())
    cash_opex = total - c_cap
    return dict(unit=total / P_final_kg, parts={k: v / P_final_kg for k, v in parts.items()},
                fci=fci, Vw=Vw, n_vessels=n_vessels, sugar_t=sugar / 1000, dates_t=dates_t, X_t=X / 1000,
                otr=otr, heat_MW=heat_MW, evap_m3=evap_m3, p_sugar=p_sugar, cash_opex=cash_opex,
                pc_content=c * 1000, P_pc=p["Px"] * c)


def pct(a, q):
    return float(np.percentile(a, q))


def route_comparison(rng):
    rows = []
    for rk, route in ROUTES.items():
        for s in SCALES_T:
            p = draw(route, rng, N)
            r = unit_cost(route, p, s * 1000.0)
            rows.append(dict(route=rk, label=route["label"], scale_t_pc=s,
                             cost_pc_P10=pct(r["unit"], 10), cost_pc_P50=pct(r["unit"], 50), cost_pc_P90=pct(r["unit"], 90),
                             cost_product34_P50=pct(r["unit"], 50) * PRODUCT_PC_FRACTION,
                             cost_E18eq_P50=pct(r["unit"], 50) * E18_PC_FRACTION,
                             p_below_low=float(np.mean(r["unit"] < PRICE_ANCHORS["low_china_E18_bulk"])),
                             p_below_mid=float(np.mean(r["unit"] < PRICE_ANCHORS["mid_china_E25"])),
                             p_below_high=float(np.mean(r["unit"] < PRICE_ANCHORS["high_eu_E18_250kg"])),
                             capex_M_P50=pct(r["fci"], 50), capex_M_P10=pct(r["fci"], 10), capex_M_P90=pct(r["fci"], 90),
                             Vw_m3_P50=pct(r["Vw"], 50), vessels_P50=pct(r["n_vessels"], 50),
                             sugar_t_P50=pct(r["sugar_t"], 50), dates_t_P50=pct(r["dates_t"], 50),
                             biomass_t_P50=pct(r["X_t"], 50), OTR_mmol_L_h_P50=pct(r["otr"], 50),
                             OTR_P90=pct(r["otr"], 90), cooling_MW_P50=pct(r["heat_MW"], 50),
                             evap_m3_y_P50=pct(r["evap_m3"], 50), sugar_usd_kg_P50=pct(r["p_sugar"], 50)))
    return pd.DataFrame(rows)


def base_breakdown():
    rows = []
    for s in SCALES_T:
        for rk in ("R2_climited_glucose", RECOMMENDED):
            route = ROUTES[rk]
            p = draw(route, None, 1, mode_only=True)
            r = unit_cost(route, p, s * 1000.0)
            row = dict(route=rk, scale_t_pc=s, total_usd_per_kg_pc=float(r["unit"][0]),
                       capex_M=float(r["fci"][0]), Vw_m3=float(r["Vw"][0]))
            row.update({k: float(v[0]) for k, v in r["parts"].items()})
            rows.append(row)
    return pd.DataFrame(rows)


def irr_vec(cf):
    """IRR by bisection for many cash-flow rows (NaN if no sign change in [-0.5, 1.0])."""
    lo = np.full(cf.shape[0], -0.5); hi = np.full(cf.shape[0], 1.0)
    t = np.arange(cf.shape[1])
    f = lambda r: (cf / (1 + r[:, None]) ** t).sum(1)
    flo, fhi = f(lo), f(hi)
    ok = np.sign(flo) != np.sign(fhi)
    for _ in range(60):
        mid = (lo + hi) / 2; fm = f(mid)
        left = np.sign(fm) == np.sign(flo)
        lo = np.where(left, mid, lo); flo = np.where(left, fm, flo); hi = np.where(left, hi, mid)
    out = (lo + hi) / 2
    out[~ok] = np.nan
    return out


def dcf(rng, n=20_000):
    rows = []
    route = ROUTES[RECOMMENDED]
    for s in SCALES_T:
        p = draw(route, rng, n)
        r = unit_cost(route, p, s * 1000.0)
        fci = r["fci"] * 1e6
        for pk, price in PRICE_ANCHORS.items():
            years = 2 + CRF_YEARS
            cf = np.zeros((n, years))
            cf[:, 0] = -0.4 * fci; cf[:, 1] = -0.6 * fci
            wc = 0.25 * r["cash_opex"]
            cf[:, 1] -= wc
            loss_cf = np.zeros(n)
            for y in range(CRF_YEARS):
                ramp = (0.5, 0.8)[y] if y < 2 else 1.0
                rev = price * s * 1000.0 * ramp
                # variable part of cash opex scales with ramp; ~35 % is fixed (labour, maintenance, QA) UNVERIFIED
                opex = r["cash_opex"] * (0.35 + 0.65 * ramp)
                dep = fci / 10.0 if y < 10 else 0.0
                ebt = rev - opex - dep
                taxable = ebt - loss_cf
                loss_cf = np.where(taxable < 0, -taxable, 0.0)
                tax = np.where(taxable > 0, taxable * TAX_RATE, 0.0)
                cf[:, 2 + y] = rev - opex - tax
            cf[:, -1] += wc
            disc = (1 + p["rate"][:, None]) ** np.arange(years)
            npv = (cf / disc).sum(1)
            irr = irr_vec(cf)
            rows.append(dict(scale_t_pc=s, price_anchor=pk, price_usd_per_kg_pc=round(price, 1),
                             revenue_M=price * s / 1000.0, capex_M_P50=pct(r["fci"], 50),
                             NPV_M_P10=pct(npv, 10) / 1e6, NPV_M_P50=pct(npv, 50) / 1e6, NPV_M_P90=pct(npv, 90) / 1e6,
                             P_NPV_pos=float(np.mean(npv > 0)), IRR_P50=float(np.nanmedian(irr)) if np.isfinite(irr).any() else np.nan))
    return pd.DataFrame(rows)


def tornado(scale_t=100):
    route = ROUTES[RECOMMENDED]
    base = draw(route, None, 1, mode_only=True)
    b = float(unit_cost(route, base, scale_t * 1000.0)["unit"][0])
    rows = []
    for k, (rngdef, unit, src) in P.items():
        lo, mo, hi = route["over"].get(k, rngdef)
        if lo == hi:
            continue
        vals = []
        for v in (lo, hi):
            p = draw(route, None, 1, mode_only=True, fix={k: v})
            vals.append(float(unit_cost(route, p, scale_t * 1000.0)["unit"][0]))
        rows.append(dict(param=k, unit=unit, low=lo, high=hi, cost_at_low=vals[0], cost_at_high=vals[1],
                         swing=abs(vals[1] - vals[0]), base=b, source=src))
    return pd.DataFrame(rows).sort_values("swing", ascending=False)


def lab_targets(scale_t=100):
    route = ROUTES[RECOMMENDED]
    rows = []
    for c_pc in (5, 10, 15, 20, 25, 30, 35, 40, 50):
        for px in (10, 15, 20, 25, 30, 40):
            p = draw(route, None, 1, mode_only=True, fix={"c_pc": c_pc, "Px": px, "f_date_c": 1.0})
            rows.append(dict(c_pc_mg_g=c_pc, Px_g_L_d=px, P_pc_g_L_d=c_pc * px / 1000,
                             cost_usd_per_kg_pc=float(unit_cost(route, p, scale_t * 1000.0)["unit"][0])))
    return pd.DataFrame(rows)


def scenarios(scale_t=100):
    """Deterministic what-ifs (all other parameters at mode) and the date-price break-even."""
    rows = []
    cases = {
        "base_mode_dates": (RECOMMENDED, {}),
        "base_mode_glucose": ("R2_climited_glucose", {}),
        "dates_no_penalty": (RECOMMENDED, {"f_date_c": 1.0, "f_date_Y": 1.0}),
        "breakthrough_dates": (RECOMMENDED, {"c_pc": 40, "Px": 30, "R_dsp": 0.75, "f_date_c": 1.0, "f_date_Y": 1.0}),
        "breakthrough_dates_cheap_rebuts": (RECOMMENDED, {"c_pc": 40, "Px": 30, "R_dsp": 0.75, "f_date_c": 1.0,
                                                          "f_date_Y": 1.0, "date_dzd": 15}),
        "pessimistic_dates": (RECOMMENDED, {"c_pc": 15, "Px": 15, "R_dsp": 0.5, "f_date_c": 0.8}),
    }
    for name, (rk, fix) in cases.items():
        p = draw(ROUTES[rk], None, 1, mode_only=True, fix=fix)
        r = unit_cost(ROUTES[rk], p, scale_t * 1000.0)
        rows.append(dict(case=name, route=rk, fixed=str(fix), usd_per_kg_pc=float(r["unit"][0]),
                         usd_per_kg_product34=float(r["unit"][0]) * PRODUCT_PC_FRACTION,
                         usd_per_kg_E18eq=float(r["unit"][0]) * E18_PC_FRACTION, capex_M=float(r["fci"][0])))
    # date price at which the date route equals the glucose route (no biological penalty)
    g = float(unit_cost(ROUTES["R2_climited_glucose"], draw(ROUTES["R2_climited_glucose"], None, 1, mode_only=True),
                        scale_t * 1000.0)["unit"][0])
    lo, hi = 0.0, 150.0
    for _ in range(50):
        mid = (lo + hi) / 2
        p = draw(ROUTES[RECOMMENDED], None, 1, mode_only=True, fix={"f_date_c": 1.0, "f_date_Y": 1.0, "date_dzd": mid})
        d = float(unit_cost(ROUTES[RECOMMENDED], p, scale_t * 1000.0)["unit"][0])
        lo, hi = (mid, hi) if d < g else (lo, mid)
    rows.append(dict(case="breakeven_date_price_vs_glucose", route=RECOMMENDED, fixed="f_date=1",
                     usd_per_kg_pc=g, usd_per_kg_product34=np.nan, usd_per_kg_E18eq=np.nan,
                     capex_M=np.nan, date_dzd_breakeven=(lo + hi) / 2))
    return pd.DataFrame(rows)


def figures(cmp_df, tor, lab):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = cmp_df[cmp_df.scale_t_pc == 100].copy()
    fig, ax = plt.subplots(figsize=(9, 4.8))
    y = np.arange(len(d))
    ax.barh(y, d.cost_pc_P50, xerr=[d.cost_pc_P50 - d.cost_pc_P10, d.cost_pc_P90 - d.cost_pc_P50],
            color=["#999" if "R3" not in r else "#1f5fa8" for r in d.route], capsize=3)
    ax.set_yticks(y); ax.set_yticklabels(d.route, fontsize=8); ax.set_xscale("log")
    for name, v in PRICE_ANCHORS.items():
        ax.axvline(v, ls="--", lw=1, color="#c0392b")
        ax.text(v * 1.02, -0.45, f"{name} ${v:.0f}", fontsize=7, color="#c0392b", rotation=90, va="bottom")
    ax.set_xlabel("$ per kg pure C-phycocyanin (P10-P50-P90), 100 t PC/yr")
    ax.set_title("Unit cost by process route vs spirulina-blue price anchors")
    fig.tight_layout(); fig.savefig(DATA / "fig_cost_by_route.png", dpi=130); plt.close(fig)

    t = tor.head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 5))
    base = t.base.iloc[0]
    ax.barh(t.param, t.cost_at_low - base, left=base, color="#1f5fa8", label="param at low")
    ax.barh(t.param, t.cost_at_high - base, left=base, color="#e67e22", label="param at high")
    ax.axvline(base, color="k", lw=1)
    ax.set_xlabel("$ per kg PC (recommended route, 100 t/yr, others at mode)")
    ax.set_title("What moves the cost (tornado)"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(DATA / "fig_tornado.png", dpi=130); plt.close(fig)

    piv = lab.pivot(index="c_pc_mg_g", columns="Px_g_L_d", values="cost_usd_per_kg_pc")
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(piv.values, origin="lower", aspect="auto", cmap="RdYlGn_r", vmin=50, vmax=600)
    ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns)
    ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, f"{piv.values[i, j]:.0f}", ha="center", va="center", fontsize=7)
    ax.set_xlabel("biomass productivity at scale, g/L/day"); ax.set_ylabel("PC content, mg/g biomass")
    ax.set_title("$/kg PC on date syrup, 100 t/yr: what the lab must prove")
    fig.colorbar(im, ax=ax, label="$/kg PC")
    fig.tight_layout(); fig.savefig(DATA / "fig_lab_targets.png", dpi=130); plt.close(fig)


def main():
    rng = np.random.default_rng(SEED)
    prm = [dict(param=k, low=v[0][0], mode=v[0][1], high=v[0][2], unit=v[1], source=v[2]) for k, v in P.items()]
    for rk, r in ROUTES.items():
        for k, v in r["over"].items():
            prm.append(dict(param=f"{rk}:{k}", low=v[0], mode=v[1], high=v[2], unit=P[k][1], source=r["sources"]))
    prm += [dict(param=f"price:{k}", low=v, mode=v, high=v, unit="$/kg PC", source="computed (botanicalcube2026, tint2026)")
            for k, v in PRICE_ANCHORS.items()]
    pd.DataFrame(prm).to_csv(DATA / "01_parameters.tsv", sep="\t", index=False)

    cmp_df = route_comparison(rng); cmp_df.round(3).to_csv(DATA / "02_route_comparison.tsv", sep="\t", index=False)
    bb = base_breakdown(); bb.round(2).to_csv(DATA / "03_base_breakdown.tsv", sep="\t", index=False)
    dc = dcf(rng); dc.round(3).to_csv(DATA / "04_dcf.tsv", sep="\t", index=False)
    tor = tornado(); tor.round(2).to_csv(DATA / "05_tornado.tsv", sep="\t", index=False)
    lab = lab_targets(); lab.round(1).to_csv(DATA / "06_lab_targets.tsv", sep="\t", index=False)
    sc = scenarios(); sc.round(2).to_csv(DATA / "07_scenarios.tsv", sep="\t", index=False)
    figures(cmp_df, tor, lab)

    pd.set_option("display.width", 220); pd.set_option("display.max_columns", 30)
    print("price anchors $/kg PC:", {k: round(v, 1) for k, v in PRICE_ANCHORS.items()})
    print(cmp_df[["route", "scale_t_pc", "cost_pc_P10", "cost_pc_P50", "cost_pc_P90", "cost_product34_P50",
                  "p_below_low", "p_below_mid", "p_below_high", "capex_M_P50", "Vw_m3_P50", "dates_t_P50",
                  "OTR_mmol_L_h_P50", "sugar_usd_kg_P50"]].round(2).to_string(index=False))
    print(bb.round(1).to_string(index=False))
    print(dc.round(2).to_string(index=False))
    print(sc.round(1).to_string(index=False))
    print(tor[["param", "low", "high", "cost_at_low", "cost_at_high", "swing"]].head(15).round(1).to_string(index=False))


if __name__ == "__main__":
    main()
