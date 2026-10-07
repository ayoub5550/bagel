#!/usr/bin/env python3
"""05 — A targeted HRMS/MS2 list for the saponins of *Anabasis articulata*.

Why this is the cheapest decisive experiment
--------------------------------------------
Everything published on this plant's saponins comes from four isolation papers; no one
has run high-resolution MS over the crude saponin fraction. A single LC-HRMS/MS run
(one injection, ~50-100 mg of dried aerial parts) either
  (a) shows only the already-published compositions  -> the "undescribed saponins"
      hypothesis is dead, and the project should stop spending on it, or
  (b) shows compositions that are not in the literature -> there is something new,
      and the next spend (preparative isolation + NMR) is justified.

This script builds the list needed to read that run: the exact masses of every
plausible glycoform of the aglycones this plant is *known* to make, the adduct ions to
watch in both polarities, and the MS2 Y-ion ladders that identify the sugar chain.

Inputs : data/01b_triterpenoids.csv  (aglycone formulas, from constructed structures)
Outputs: data/05_ms_targets.csv            full enumerated space
         data/05_ms_inclusion_list.csv     the short list to put on the instrument
         data/05_ms_method_notes.md        how to run and read it

Nothing here is a measurement. Masses are arithmetic from atomic constants; retention
times are unknown and deliberately left empty.
"""

from __future__ import annotations

import csv
import itertools
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

ROOT = Path(__file__).resolve().parents[1]
IN_CSV = ROOT / "data" / "01b_triterpenoids.csv"
OUT_ALL = ROOT / "data" / "05_ms_targets.csv"
OUT_LIST = ROOT / "data" / "05_ms_inclusion_list.csv"
OUT_MD = ROOT / "data" / "05_ms_method_notes.md"

# CODATA/NIST monoisotopic masses (u)
M = {"C": 12.0, "H": 1.00782503207, "O": 15.9949146196, "N": 14.0030740048,
     "Na": 22.9897692809, "K": 38.96370668, "Cl": 34.96885268}
PROTON = 1.007276466879
ELECTRON = 0.000548579909

# Glycosyl residues: mass added to the aglycone per sugar (free sugar - H2O)
RESIDUES = {
    "Hex": {"C": 6, "H": 10, "O": 5, "aka": "Glc / Gal"},
    "Pen": {"C": 5, "H": 8, "O": 4, "aka": "Xyl / Ara"},
    "GlcA": {"C": 6, "H": 8, "O": 6, "aka": "glucuronic acid"},
    "dHex": {"C": 6, "H": 10, "O": 4, "aka": "Rha (not yet reported here)"},
}

# Compositions already in the literature, as (aglycone_id, {residue: n}) -> reference
PUBLISHED = {
    ("AGL-1", ("Hex", 1)): "metwally2012 (oleanolic acid 3-O-Glc)",
    ("AGL-1", ("Hex", 1), ("Pen", 1)): "metwally2012 (3-O-Glc-28-O-Xyl)",
    ("AGL-2", ("Hex", 1)): "salaheldine2019 cpd 1 (28-O-Glc ester)",
    ("AGL-2", ("Hex", 2)): "salaheldine2019 cpd 2 (3-O-Gal, 28-O-Glc)",
    ("AGL-3", ("Hex", 1), ("Pen", 1)): "salaheldine2019 cpd 3 (3-O-Xyl, 28-O-Glc)",
    ("AGL-6", ("Hex", 1)): "gamal2022 cpd 5 (23-O-Glc ester)",
    ("AGL-6", ("GlcA", 1), ("Hex", 1)): "gamal2022 cpd 6 (3-O-GlcA, 28-O-Glc)",
    ("AGL-5", ("GlcA", 1), ("Hex", 1)): "gamal2022 cpd 7 (3-O-GlcA, 28-O-Glc)",
}

MAX_RESIDUES = 4
LIMITS = {"Hex": 3, "Pen": 2, "GlcA": 1, "dHex": 1}


def formula_mass(counts: dict[str, int]) -> float:
    return sum(M[el] * n for el, n in counts.items())


def parse_formula(formula: str) -> dict[str, int]:
    out: dict[str, int] = {}
    token, num = "", ""
    for ch in formula + "|":
        if ch.isupper() or ch == "|":
            if token:
                out[token] = out.get(token, 0) + (int(num) if num else 1)
            token, num = ("" if ch == "|" else ch), ""
        elif ch.islower():
            token += ch
        elif ch.isdigit():
            num += ch
    return out


def load_aglycones() -> list[dict[str, object]]:
    keep = {"AGL-1", "AGL-2", "AGL-3", "AGL-4", "AGL-5", "AGL-6"}
    rows = []
    with IN_CSV.open() as fh:
        for row in csv.DictReader(fh):
            if row["id"] not in keep:
                continue
            mol = Chem.MolFromSmiles(row["smiles"])
            counts = parse_formula(CalcMolFormula(mol))
            rows.append({
                "id": row["id"], "name": row["literature_name"],
                "skeleton": row["skeleton"], "source_key": row["source_key"],
                "counts": counts, "mass": formula_mass(counts),
                "formula": CalcMolFormula(mol),
            })
    return rows


def enumerate_glycoforms():
    """All residue combinations from 0 to MAX_RESIDUES sugars."""
    names = list(RESIDUES)
    for combo in itertools.product(*[range(LIMITS[n] + 1) for n in names]):
        total = sum(combo)
        if total == 0 or total > MAX_RESIDUES:
            continue
        yield {n: c for n, c in zip(names, combo) if c}


def composition_key(agl_id: str, comp: dict[str, int]):
    return (agl_id, *sorted(comp.items()))


def adducts(neutral: float) -> dict[str, float]:
    return {
        "[M-H]-": neutral - PROTON,
        "[M+HCOO]-": neutral + formula_mass({"C": 1, "H": 2, "O": 2}) - PROTON,
        "[M+Cl]-": neutral + M["Cl"] + ELECTRON,
        "[M+Na-2H]-": neutral + M["Na"] - 2 * PROTON - ELECTRON,
        "[M+H]+": neutral + PROTON,
        "[M+Na]+": neutral + M["Na"] - ELECTRON,
        "[M+K]+": neutral + M["K"] - ELECTRON,
        "[M+NH4]+": neutral + formula_mass({"N": 1, "H": 4}) + ELECTRON * 0,
    }


def y_ions(agl_mass: float, comp: dict[str, int]) -> list[str]:
    """Y-ion ladder in negative mode: sequential neutral losses of the sugars."""
    out = []
    remaining = dict(comp)
    mass = agl_mass + sum(formula_mass({k: v for k, v in RESIDUES[r].items()
                                        if k in ("C", "H", "O")}) * n
                          for r, n in comp.items())
    order = [r for r in ("Hex", "dHex", "Pen", "GlcA") for _ in range(remaining.get(r, 0))]
    for res in order:
        mass -= formula_mass({k: v for k, v in RESIDUES[res].items() if k in ("C", "H", "O")})
        out.append(f"-{res}: {mass - PROTON:.4f}")
    return out


def main() -> int:
    aglycones = load_aglycones()
    if not aglycones:
        print("[05] no aglycones found — run scripts/01b_build_triterpenoids.py first")
        return 1

    rows: list[dict[str, object]] = []
    for agl in aglycones:
        # the free aglycone itself
        rows.append(make_row(agl, {}, "aglycone (free)"))
        for comp in enumerate_glycoforms():
            key = composition_key(str(agl["id"]), comp)
            rows.append(make_row(agl, comp, PUBLISHED.get(key, "")))

    with OUT_ALL.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # ---- the practical inclusion list ---------------------------------------
    short = [r for r in rows
             if int(r["n_sugars"]) <= 3 and str(r["residues"]) != "dHex1"
             and "dHex" not in str(r["residues"])]
    short.sort(key=lambda r: (r["aglycone_id"], float(r["neutral_mass"])))
    list_cols = ["target_id", "aglycone_id", "aglycone_name", "composition", "formula",
                 "neutral_mass", "mz_M-H", "mz_M+HCOO", "mz_M+Na", "status",
                 "published_as", "ms2_y_ions_neg", "rt_min"]
    with OUT_LIST.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list_cols)
        writer.writeheader()
        for i, r in enumerate(short, 1):
            writer.writerow({
                "target_id": f"T{i:03d}",
                "aglycone_id": r["aglycone_id"],
                "aglycone_name": r["aglycone_name"],
                "composition": r["residues"] or "aglycone",
                "formula": r["formula"],
                "neutral_mass": r["neutral_mass"],
                "mz_M-H": r["[M-H]-"],
                "mz_M+HCOO": r["[M+HCOO]-"],
                "mz_M+Na": r["[M+Na]+"],
                "status": r["status"],
                "published_as": r["published_as"],
                "ms2_y_ions_neg": r["ms2_y_ions_neg"],
                "rt_min": "",
            })

    n_pub = sum(1 for r in rows if r["status"] == "published")
    n_pred = sum(1 for r in rows if r["status"] == "predicted")
    write_notes(aglycones, len(rows), len(short), n_pub, n_pred)

    print(f"[05] aglycone scaffolds           : {len(aglycones)}")
    print(f"[05] enumerated glycoforms        : {len(rows)} "
          f"({n_pub} already published, {n_pred} predicted/unreported)")
    print(f"[05] inclusion list (<=3 sugars)  : {len(short)} targets")
    print(f"[05] wrote {OUT_ALL.relative_to(ROOT)}, {OUT_LIST.relative_to(ROOT)}, "
          f"{OUT_MD.relative_to(ROOT)}")
    for r in rows:
        if r["status"] == "published":
            print(f"       published  {r['aglycone_id']} {str(r['residues']):16s} "
                  f"{r['formula']:12s} [M-H]- {r['[M-H]-']}  {r['published_as']}")
    return 0


def make_row(agl: dict[str, object], comp: dict[str, int], published: str) -> dict[str, object]:
    counts = dict(agl["counts"])  # type: ignore[arg-type]
    for res, n in comp.items():
        for el in ("C", "H", "O"):
            counts[el] = counts.get(el, 0) + RESIDUES[res][el] * n  # type: ignore[index]
    neutral = formula_mass(counts)
    ions = adducts(neutral)
    formula = "".join(f"{el}{counts[el]}" for el in ("C", "H", "O") if counts.get(el))
    row: dict[str, object] = {
        "aglycone_id": agl["id"],
        "aglycone_name": agl["name"],
        "skeleton": agl["skeleton"],
        "aglycone_source_key": agl["source_key"],
        "residues": "".join(f"{r}{n}" for r, n in sorted(comp.items())),
        "n_sugars": sum(comp.values()),
        "formula": formula,
        "neutral_mass": round(neutral, 4),
        "status": "published" if published and comp else
                  ("aglycone" if not comp else "predicted"),
        "published_as": published if comp else "",
        "ms2_y_ions_neg": " | ".join(y_ions(float(agl["mass"]), comp)) if comp else "",
    }
    for name, mz in ions.items():
        row[name] = round(mz, 4)
    return row


def write_notes(aglycones, n_all, n_short, n_pub, n_pred) -> None:
    lines = [
        "# 05 — targeted HRMS/MS2 list: method notes",
        "",
        "Generated by `scripts/05_ms_target_list.py`. Masses are arithmetic from NIST "
        "monoisotopic atomic masses; **no retention time is given because none is known**.",
        "",
        "## What goes on the instrument",
        "",
        f"- `data/05_ms_inclusion_list.csv` — {n_short} targets (<= 3 sugars, no "
        "deoxyhexose). Negative-mode ESI is the primary channel for saponins: watch "
        "`[M-H]-` and `[M+HCOO]-` (formate adducts dominate when the mobile phase "
        "contains ammonium formate).",
        f"- Full enumerated space (up to 4 sugars, incl. rhamnose): "
        f"`data/05_ms_targets.csv`, {n_all} rows — {n_pub} compositions already "
        f"published, {n_pred} not reported for this species.",
        "",
        "## Aglycone scaffolds used (all with NMR-level literature support)",
        "",
        "| id | aglycone | formula | monoisotopic | source |",
        "|---|---|---|---|---|",
    ]
    for a in aglycones:
        lines.append(f"| {a['id']} | {a['name']} | {a['formula']} | "
                     f"{a['mass']:.4f} | {a['source_key']} |")
    lines += [
        "",
        "## Minimum acceptable run",
        "",
        "1. Dried, milled aerial parts, 50-100 mg; 70% MeOH or 80% EtOH extract, "
        "filtered. A defatted n-BuOH partition enriches the saponins, which is what the "
        "isolation papers used.",
        "2. C18 column, water/acetonitrile + 0.1% formic acid or 5 mM ammonium formate, "
        "20-40 min gradient.",
        "3. HRMS (Q-TOF or Orbitrap), **both polarities**, mass accuracy < 5 ppm, with "
        "data-dependent MS2 at 2-3 collision energies (20/35/50 eV).",
        "4. Report measured m/z, ppm error against this list, isotope pattern, and the "
        "MS2 Y-ion ladder.",
        "",
        "## How to read the result (decision rule, set before the run)",
        "",
        "| observation | conclusion | next step |",
        "|---|---|---|",
        "| only published compositions detected | the saponin pool is already described |"
        " stop spending on saponin discovery; pivot to quantification (content %) |",
        "| new composition, < 5 ppm, with a clean Y-ion ladder and the right isotope "
        "pattern | a glycoform not reported for this species | preparative isolation of "
        "that single peak, then NMR; only NMR proves a structure |",
        "| a mass matches but MS2 shows no sugar losses | probably an isobaric "
        "non-saponin | discard |",
        "| nothing matches | the literature saponin set may be wrong, or the material "
        "differs (season, population, chemotype) | re-check the voucher and the "
        "extraction before blaming the plant |",
        "",
        "## Limits of this list",
        "",
        "- A mass match is **not** an identification. Sugar identity (Glc vs Gal, Xyl vs "
        "Ara), linkage position and anomeric configuration are invisible to MS: "
        "hexose/pentose counts are all this list can decide.",
        "- Isomers collapse: `Hex1` on C-3 and `Hex1` on C-28 have the same mass.",
        "- Aglycone coverage is limited to the six scaffolds with literature support. A "
        "genuinely new aglycone would appear as an unexplained mass and would need "
        "untargeted processing (MZmine/GNPS) rather than this list.",
        "- Boussingoside E was also reported from this plant (salaheldine2019). Its "
        "structure was not reconstructed here, so it is absent from the list: UNVERIFIED.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
