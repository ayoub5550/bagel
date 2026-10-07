#!/usr/bin/env python3
"""01 — Resolve the manual compound inventory to real chemical structures.

Input : data/compounds.tsv           (manually curated, one row per compound)
Output: data/01_resolved.csv         (one row per compound + structure)
Cache : data/cache/pubchem/*.json    (raw PubChem answers, so reruns are free)

Rules (AGENTS.md):
  * every structure carries its provenance: `smiles_source` = pubchem | manual
  * a name PubChem does not know is NOT silently dropped; it is reported and,
    where the literature gives an unambiguous structure, taken from
    data/manual_smiles.tsv with smiles_source=manual.

The script is idempotent: delete data/cache/pubchem to force a refetch.
"""

from __future__ import annotations

import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IN_TSV = ROOT / "data" / "compounds.tsv"
MANUAL_TSV = ROOT / "data" / "manual_smiles.tsv"
OUT_CSV = ROOT / "data" / "01_resolved.csv"
CACHE = ROOT / "data" / "cache" / "pubchem"

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PROPS = "SMILES,ConnectivitySMILES,InChIKey,MolecularFormula,MolecularWeight"
UA = "bagel-dry-lab/1.0 (open research; https://github.com/ayoub5550/bagel)"
SLEEP = 0.25  # PubChem asks for <= 5 requests/second


def slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")[:120]


def fetch_json(url: str, cache_key: str) -> dict | None:
    """GET a PubChem URL with an on-disk cache. None = resource not found (404)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    cached = CACHE / f"{cache_key}.json"
    if cached.exists():
        raw = json.loads(cached.read_text())
        return None if raw.get("__status__") == 404 else raw

    req = urllib.request.Request(url, headers={"User-Agent": UA})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=40) as resp:
                payload = json.loads(resp.read().decode())
            cached.write_text(json.dumps(payload))
            time.sleep(SLEEP)
            return payload
        except urllib.error.HTTPError as err:
            if err.code == 404:
                cached.write_text(json.dumps({"__status__": 404}))
                time.sleep(SLEEP)
                return None
            if err.code in (429, 503) and attempt < 3:
                time.sleep(2 ** attempt + 1)
                continue
            raise
        except (urllib.error.URLError, TimeoutError):
            if attempt < 3:
                time.sleep(2 ** attempt + 1)
                continue
            raise
    return None


def by_name(name: str) -> dict | None:
    url = f"{PUG}/compound/name/{urllib.parse.quote(name)}/property/{PROPS}/JSON"
    payload = fetch_json(url, f"name_{slug(name)}")
    if not payload:
        return None
    props = payload.get("PropertyTable", {}).get("Properties", [])
    return props[0] if props else None


def load_manual() -> dict[str, dict[str, str]]:
    if not MANUAL_TSV.exists():
        return {}
    with MANUAL_TSV.open() as fh:
        return {row["name"].strip(): row for row in csv.DictReader(fh, delimiter="\t")}


def main() -> int:
    with IN_TSV.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    manual = load_manual()

    out: list[dict[str, object]] = []
    unresolved: list[str] = []
    delegated: list[str] = []

    for row in rows:
        name = row["name"].strip()
        built_elsewhere = (row.get("structure_source") or "").strip()
        aliases = [a.strip() for a in (row.get("aliases") or "").split(";") if a.strip()]
        hit, via = None, ""
        if built_elsewhere:
            delegated.append(f"{name} -> {built_elsewhere}")
        for candidate in ([] if built_elsewhere else [name, *aliases]):
            hit = by_name(candidate)
            if hit:
                via = "name" if candidate == name else f"alias:{candidate}"
                break

        record = {
            "name": name,
            "chem_class": row["chem_class"],
            "plant_part": row["plant_part"],
            "evidence": row["evidence"],
            "source_key": row["source_key"],
            "cid": "",
            "smiles": "",
            "connectivity_smiles": "",
            "inchikey": "",
            "formula": "",
            "mw_pubchem": "",
            "resolved_via": "",
            "smiles_source": "",
            "structure_source": built_elsewhere,
            "notes": row.get("notes", ""),
        }

        if hit:
            record.update(
                cid=hit.get("CID", ""),
                smiles=hit.get("SMILES", "") or hit.get("ConnectivitySMILES", ""),
                connectivity_smiles=hit.get("ConnectivitySMILES", ""),
                inchikey=hit.get("InChIKey", ""),
                formula=hit.get("MolecularFormula", ""),
                mw_pubchem=hit.get("MolecularWeight", ""),
                resolved_via=via,
                smiles_source="pubchem",
            )
        elif built_elsewhere:
            record["resolved_via"] = f"delegated:{built_elsewhere}"
            record["smiles_source"] = "see_structure_source"
        elif name in manual:
            record.update(
                smiles=manual[name]["smiles"],
                connectivity_smiles=manual[name]["smiles"],
                resolved_via=f"manual:{manual[name].get('basis', 'literature')}",
                smiles_source="manual",
            )
        else:
            unresolved.append(name)
            record["resolved_via"] = "UNRESOLVED"
            record["smiles_source"] = "none"

        out.append(record)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        writer.writeheader()
        writer.writerows(out)

    pubchem = sum(1 for r in out if r["smiles_source"] == "pubchem")
    manual_n = sum(1 for r in out if r["smiles_source"] == "manual")
    print(f"[01] input compounds        : {len(out)}")
    print(f"[01] resolved via PubChem   : {pubchem}")
    print(f"[01] structures from manual : {manual_n}")
    print(f"[01] built by another script: {len(delegated)}")
    for item in delegated:
        print(f"       DELEGATED   {item}")
    print(f"[01] unresolved             : {len(unresolved)}")
    for name in unresolved:
        print(f"       UNRESOLVED  {name}")
    print(f"[01] wrote {OUT_CSV.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
