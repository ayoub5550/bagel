#!/usr/bin/env python3
"""03 — What has actually been *measured* for these molecules (ChEMBL).

Two questions, answered with measured data only:

  1. Which inventory compounds exist in ChEMBL, and what are their measured activities
     with a pChEMBL value (i.e. real dose-response numbers, not screening noise)?
  2. For the triterpenoids this plant makes (which are absent from ChEMBL), what are the
     nearest measured analogues, and on which targets? -> hypotheses to test, nothing more.

Inputs : data/01_resolved.csv, data/01b_triterpenoids.csv
Outputs: data/03_chembl_molecules.csv   one row per inventory compound
         data/03_chembl_activities.csv  one row per (compound, target)
         data/03_similar_actives.csv    nearest ChEMBL analogues of the triterpenoids
         data/03_chembl_summary.md
Cache  : data/cache/chembl/
"""

from __future__ import annotations

import csv
import json
import statistics
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN_RESOLVED = ROOT / "data" / "01_resolved.csv"
IN_TRITERP = ROOT / "data" / "01b_triterpenoids.csv"
OUT_MOL = ROOT / "data" / "03_chembl_molecules.csv"
OUT_ACT = ROOT / "data" / "03_chembl_activities.csv"
OUT_SIM = ROOT / "data" / "03_similar_actives.csv"
OUT_MD = ROOT / "data" / "03_chembl_summary.md"
CACHE = ROOT / "data" / "cache" / "chembl"

API = "https://www.ebi.ac.uk/chembl/api/data"
UA = "bagel-dry-lab/1.0 (open research; https://github.com/ayoub5550/bagel)"
MAX_ACT_PAGES = 5          # 1000 activities per compound is plenty for a ranking
SIMILARITY = 70            # % Tanimoto for the analogue search


def fetch(url: str, key: str) -> dict | None:
    CACHE.mkdir(parents=True, exist_ok=True)
    cached = CACHE / f"{key}.json"
    if cached.exists():
        payload = json.loads(cached.read_text())
        return None if payload.get("__status__") == 404 else payload
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                payload = json.loads(resp.read().decode())
            cached.write_text(json.dumps(payload))
            time.sleep(0.2)
            return payload
        except urllib.error.HTTPError as err:
            if err.code == 404:
                cached.write_text(json.dumps({"__status__": 404}))
                return None
            if attempt < 3:
                time.sleep(2 ** attempt + 1)
                continue
            print(f"[03] HTTP {err.code} for {url}")
            return None
        except Exception as err:  # noqa: BLE001 - network flakiness
            if attempt < 3:
                time.sleep(2 ** attempt + 1)
                continue
            print(f"[03] failed {url}: {err}")
            return None
    return None


def molecule_by_inchikey(inchikey: str) -> dict | None:
    url = f"{API}/molecule.json?molecule_structures__standard_inchi_key={inchikey}"
    payload = fetch(url, f"mol_{inchikey}")
    if not payload:
        return None
    mols = payload.get("molecules", [])
    return mols[0] if mols else None


def activities(chembl_id: str) -> list[dict]:
    out: list[dict] = []
    for page in range(MAX_ACT_PAGES):
        url = (f"{API}/activity.json?molecule_chembl_id={chembl_id}"
               f"&pchembl_value__isnull=false&limit=200&offset={page * 200}")
        payload = fetch(url, f"act_{chembl_id}_{page}")
        if not payload:
            break
        batch = payload.get("activities", [])
        out.extend(batch)
        if len(batch) < 200:
            break
    return out


def similar_molecules(smiles: str, tag: str) -> list[dict]:
    url = (f"{API}/similarity/{urllib.parse.quote(smiles, safe='')}/{SIMILARITY}.json"
           f"?limit=25")
    payload = fetch(url, f"sim_{tag}")
    return payload.get("molecules", []) if payload else []


def target_label(act: dict) -> tuple[str, str, str]:
    return (act.get("target_chembl_id") or "-",
            act.get("target_pref_name") or "-",
            act.get("target_organism") or "-")


def target_type(target_chembl_id: str) -> str:
    """SINGLE PROTEIN / CELL-LINE / ORGANISM ... — a cell-line IC50 is not a target."""
    if not target_chembl_id or target_chembl_id == "-":
        return "-"
    payload = fetch(f"{API}/target/{target_chembl_id}.json", f"tgt_{target_chembl_id}")
    return (payload or {}).get("target_type", "-") or "-"


def main() -> int:
    # ---------------------------------------------------------------- 1. inventory
    with IN_RESOLVED.open() as fh:
        inventory = [r for r in csv.DictReader(fh) if r["inchikey"]]

    mol_rows: list[dict[str, object]] = []
    act_rows: list[dict[str, object]] = []

    for row in inventory:
        mol = molecule_by_inchikey(row["inchikey"])
        rec = {
            "name": row["name"], "inchikey": row["inchikey"],
            "evidence": row["evidence"], "source_key": row["source_key"],
            "chembl_id": "", "chembl_pref_name": "", "max_phase": "",
            "n_activities_pchembl": 0, "n_targets": 0, "best_pchembl": "",
            "best_target": "",
        }
        if mol:
            rec["chembl_id"] = mol.get("molecule_chembl_id", "")
            rec["chembl_pref_name"] = mol.get("pref_name") or ""
            rec["max_phase"] = mol.get("max_phase") if mol.get("max_phase") is not None else ""
            acts = activities(str(rec["chembl_id"]))
            per_target: dict[tuple[str, str, str], list[float]] = {}
            for act in acts:
                try:
                    pch = float(act.get("pchembl_value"))
                except (TypeError, ValueError):
                    continue
                per_target.setdefault(target_label(act), []).append(pch)
            rec["n_activities_pchembl"] = sum(len(v) for v in per_target.values())
            rec["n_targets"] = len(per_target)
            best = None
            for (tid, tname, torg), values in sorted(per_target.items(),
                                                     key=lambda kv: -max(kv[1])):
                ttype = target_type(tid)
                act_rows.append({
                    "compound": row["name"], "chembl_id": rec["chembl_id"],
                    "target_chembl_id": tid, "target": tname, "target_type": ttype,
                    "organism": torg, "n": len(values),
                    "median_pchembl": round(statistics.median(values), 2),
                    "max_pchembl": round(max(values), 2),
                })
                if best is None and ttype == "SINGLE PROTEIN":
                    best = (f"{tname} ({torg})", max(values))
            if best:
                rec["best_target"] = best[0]
                rec["best_pchembl"] = round(best[1], 2)
        mol_rows.append(rec)
        print(f"[03] {row['name'][:45]:45s} {rec['chembl_id'] or '-':15s} "
              f"targets={rec['n_targets']:4d} best_pChEMBL={rec['best_pchembl'] or '-'}")

    with OUT_MOL.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(mol_rows[0].keys()))
        writer.writeheader()
        writer.writerows(mol_rows)
    with OUT_ACT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["compound", "chembl_id",
                                                "target_chembl_id", "target",
                                                "target_type", "organism", "n",
                                                "median_pchembl", "max_pchembl"])
        writer.writeheader()
        writer.writerows(sorted(act_rows, key=lambda r: -float(r["max_pchembl"])))

    # ------------------------------------------------- 2. analogues of the triterpenoids
    sim_rows: list[dict[str, object]] = []
    if IN_TRITERP.exists():
        with IN_TRITERP.open() as fh:
            triterpenoids = list(csv.DictReader(fh))
        for row in triterpenoids:
            hits = similar_molecules(row["smiles"], row["id"])
            print(f"[03] similarity {row['id']:7s} -> {len(hits):2d} ChEMBL molecules "
                  f"at >={SIMILARITY}%")
            for hit in hits:
                cid = hit.get("molecule_chembl_id", "")
                acts = activities(cid)
                per_target: dict[tuple[str, str, str], list[float]] = {}
                for act in acts:
                    try:
                        pch = float(act.get("pchembl_value"))
                    except (TypeError, ValueError):
                        continue
                    per_target.setdefault(target_label(act), []).append(pch)
                if not per_target:
                    continue
                proteins = {k: v for k, v in per_target.items()
                            if target_type(k[0]) == "SINGLE PROTEIN"}
                if not proteins:
                    continue
                (tid, tname, torg), values = max(proteins.items(),
                                                 key=lambda kv: max(kv[1]))
                sim_rows.append({
                    "query_id": row["id"], "query_name": row["literature_name"],
                    "similarity_pct": round(float(hit.get("similarity") or 0), 1),
                    "analogue_chembl_id": cid,
                    "analogue_name": hit.get("pref_name") or "",
                    "analogue_max_phase": hit.get("max_phase"),
                    "best_target": tname, "target_chembl_id": tid, "organism": torg,
                    "max_pchembl": round(max(values), 2),
                    "n_targets_with_pchembl": len(per_target),
                })
    if sim_rows:
        with OUT_SIM.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(sim_rows[0].keys()))
            writer.writeheader()
            writer.writerows(sorted(sim_rows, key=lambda r: (str(r["query_id"]),
                                                             -float(r["max_pchembl"]))))

    # ---------------------------------------------------------------- 3. summary
    in_chembl = [r for r in mol_rows if r["chembl_id"]]
    with_data = [r for r in mol_rows if int(r["n_activities_pchembl"]) > 0]
    lines = [
        "# 03 — measured activity (ChEMBL) summary",
        "",
        f"Inventory compounds queried: **{len(mol_rows)}** | present in ChEMBL: "
        f"**{len(in_chembl)}** | with at least one pChEMBL value: **{len(with_data)}**.",
        "",
        "A pChEMBL value is a real dose-response number (-log10 of IC50/EC50/Ki/Kd). "
        "No pChEMBL means nobody has measured a proper curve for that molecule, not that "
        "it is inactive.",
        "",
        "## Inventory compounds with measured data",
        "",
        "Only `SINGLE PROTEIN` targets are used for the \"best target\" column: an IC50 "
        "on a cell line is a toxicity reading, not a target.",
        "",
        "| compound | ChEMBL | max phase | targets | best pChEMBL (protein) | best protein target |",
        "|---|---|---|---|---|---|",
    ]
    for r in sorted(with_data, key=lambda r: -float(r["best_pchembl"] or 0)):
        lines.append(f"| {r['name']} | {r['chembl_id']} | {r['max_phase']} | "
                     f"{r['n_targets']} | {r['best_pchembl']} | {r['best_target']} |")
    lines += [
        "",
        "## Nearest measured analogues of the plant's own triterpenoids",
        "",
        f"The saponins and nor-triterpenoids of this plant are **not in ChEMBL**. The "
        f"table below is a >= {SIMILARITY}% Tanimoto search: these are analogues, and "
        "their targets are hypotheses for the plant's compounds, not evidence about them.",
        "",
        "| query | analogue | ChEMBL | similarity | best target | max pChEMBL |",
        "|---|---|---|---|---|---|",
    ]
    for r in sim_rows[:40]:
        lines.append(f"| {r['query_id']} | {str(r['analogue_name'])[:40] or '-'} | "
                     f"{r['analogue_chembl_id']} | {r['similarity_pct']} | "
                     f"{r['best_target']} | {r['max_pchembl']} |")
    lines.append("")
    OUT_MD.write_text("\n".join(lines))

    print(f"[03] in ChEMBL: {len(in_chembl)}/{len(mol_rows)}; with pChEMBL data: "
          f"{len(with_data)}; analogue rows: {len(sim_rows)}")
    print(f"[03] wrote {OUT_MOL.relative_to(ROOT)}, {OUT_ACT.relative_to(ROOT)}, "
          f"{OUT_SIM.relative_to(ROOT) if sim_rows else '(no analogues)'}, "
          f"{OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
