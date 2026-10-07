#!/usr/bin/env python3
"""02 — Drug-likeness descriptors and structural alerts for the whole inventory.

Inputs : data/01_resolved.csv  (PubChem-resolved inventory)
         data/01b_triterpenoids.csv  (constructed triterpenoids)
Output : data/02_descriptors.csv
         data/02_descriptors_summary.md  (what the numbers mean, per compound class)

All numbers come from RDKit (version printed at runtime). Nothing here is a measurement:
descriptors are *filters*, they say which molecules could plausibly act as oral
small-molecule drugs — not whether anything in this plant works.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
IN_RESOLVED = ROOT / "data" / "01_resolved.csv"
IN_TRITERP = ROOT / "data" / "01b_triterpenoids.csv"
OUT_CSV = ROOT / "data" / "02_descriptors.csv"
OUT_MD = ROOT / "data" / "02_descriptors_summary.md"

SUGAR_SMARTS = Chem.MolFromSmarts("[OX2][CX4H1]1[OX2][CX4][CX4][CX4][CX4]1")  # pyranoside


def alert_catalog() -> tuple[FilterCatalog, FilterCatalog]:
    pains = FilterCatalogParams()
    pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    brenk = FilterCatalogParams()
    brenk.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog(pains), FilterCatalog(brenk)


def load_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with IN_RESOLVED.open() as fh:
        for row in csv.DictReader(fh):
            rows.append({
                "id": row["cid"] and f"CID{row['cid']}" or "-",
                "name": row["name"],
                "class": row["chem_class"],
                "evidence": row["evidence"],
                "source_key": row["source_key"],
                "smiles": row["smiles"],
                "smiles_source": row["smiles_source"],
                "origin": "01_resolved",
            })
    if IN_TRITERP.exists():
        with IN_TRITERP.open() as fh:
            for row in csv.DictReader(fh):
                rows.append({
                    "id": row["id"],
                    "name": row["literature_name"],
                    "class": f"{row['skeleton']} / {row['kind']}",
                    "evidence": row["evidence"],
                    "source_key": row["source_key"],
                    "smiles": row["smiles"],
                    "smiles_source": "constructed (01b)",
                    "origin": "01b_triterpenoids",
                })
    return [r for r in rows if r["smiles"]]


def describe(mol: Chem.Mol, pains: FilterCatalog, brenk: FilterCatalog) -> dict[str, object]:
    mw = Descriptors.MolWt(mol)
    clogp = Crippen.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    lipinski = sum([mw > 500, clogp > 5, hbd > 5, hba > 10])
    veber_ok = rotb <= 10 and tpsa <= 140
    pains_hits = [e.GetDescription() for e in pains.GetMatches(mol)]
    brenk_hits = [e.GetDescription() for e in brenk.GetMatches(mol)]
    return {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(mw, 2),
        "exact_mass": round(Descriptors.ExactMolWt(mol), 4),
        "clogp": round(clogp, 2),
        "tpsa": round(tpsa, 1),
        "hbd": hbd,
        "hba": hba,
        "rot_bonds": rotb,
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "arom_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "qed": round(QED.qed(mol), 3),
        "lipinski_violations": lipinski,
        "veber_oral_ok": "yes" if veber_ok else "no",
        "n_pyranose_units": len(mol.GetSubstructMatches(SUGAR_SMARTS)),
        "pains_alerts": "; ".join(pains_hits) or "-",
        "brenk_alerts": "; ".join(brenk_hits) or "-",
    }


def main() -> int:
    pains, brenk = alert_catalog()
    rows = load_rows()
    out: list[dict[str, object]] = []
    failed: list[str] = []

    for row in rows:
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            failed.append(row["name"])
            continue
        rec = dict(row)
        rec.pop("smiles")
        rec.update(describe(mol, pains, brenk))
        rec["smiles"] = row["smiles"]
        out.append(rec)

    cols = ["id", "name", "class", "evidence", "source_key", "origin", "smiles_source",
            "formula", "mw", "exact_mass", "clogp", "tpsa", "hbd", "hba", "rot_bonds",
            "rings", "arom_rings", "fraction_csp3", "heavy_atoms", "qed",
            "lipinski_violations", "veber_oral_ok", "n_pyranose_units",
            "pains_alerts", "brenk_alerts", "smiles"]
    with OUT_CSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        writer.writerows(out)

    # ---- summary ------------------------------------------------------------
    drug_like = [r for r in out if r["lipinski_violations"] == 0 and r["veber_oral_ok"] == "yes"]
    saponins = [r for r in out if int(r["n_pyranose_units"]) >= 1 and float(r["mw"]) > 550]
    pains_flagged = [r for r in out if r["pains_alerts"] != "-"]
    artefacts = [r for r in out if r["evidence"] == "artifact_suspect"]

    lines = [
        "# 02 — descriptor summary (generated by scripts/02_descriptors.py)",
        "",
        f"RDKit {rdkit.__version__}. Structures scored: **{len(out)}** "
        f"(failed to parse: {len(failed)}).",
        "",
        "| group | n | meaning |",
        "|---|---|---|",
        f"| passes Lipinski *and* Veber | {len(drug_like)} | could plausibly be an oral "
        "small molecule; says nothing about activity |",
        f"| glycosides with MW > 550 (saponins) | {len(saponins)} | far outside oral "
        "small-molecule space: expect low passive absorption, gut hydrolysis, and "
        "detergent-like membrane effects rather than clean receptor pharmacology |",
        f"| PAINS alerts | {len(pains_flagged)} | assay-interference motifs (catechols, "
        "polyphenols); these are exactly the compounds that light up DPPH-type assays |",
        f"| flagged as suspected artefacts in data/compounds.tsv | {len(artefacts)} | "
        "BHT/BHA/phthalates/nitrosamine — plasticiser and additive contamination, not "
        "plant chemistry |",
        "",
        "## Compounds that pass both oral filters",
        "",
        "| name | MW | cLogP | QED | PAINS |",
        "|---|---|---|---|---|",
    ]
    for r in sorted(drug_like, key=lambda r: -float(r["qed"])):
        lines.append(f"| {r['name']} | {r['mw']} | {r['clogp']} | {r['qed']} | "
                     f"{'yes' if r['pains_alerts'] != '-' else 'no'} |")
    lines += [
        "",
        "## The saponins, measured against oral drug space",
        "",
        "| id | name | MW | cLogP | TPSA | HBD | rot. bonds | sugars | Lipinski viol. |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in sorted(saponins, key=lambda r: float(r["mw"])):
        lines.append(f"| {r['id']} | {str(r['name'])[:60]} | {r['mw']} | {r['clogp']} | "
                     f"{r['tpsa']} | {r['hbd']} | {r['rot_bonds']} | "
                     f"{r['n_pyranose_units']} | {r['lipinski_violations']} |")
    lines += [
        "",
        "Read this table as a warning, not a result: every saponin breaks two to four "
        "Lipinski rules and the Veber limits. Any claim that one of them is an oral drug "
        "candidate has to explain the absorption problem first.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))

    print(f"[02] RDKit {rdkit.__version__}; scored {len(out)} structures, {len(failed)} failures")
    print(f"[02] oral-filter pass: {len(drug_like)} | saponins: {len(saponins)} | "
          f"PAINS: {len(pains_flagged)}")
    for name in failed:
        print(f"       PARSE-FAIL {name}")
    print(f"[02] wrote {OUT_CSV.relative_to(ROOT)} and {OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
