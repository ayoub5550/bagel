#!/usr/bin/env python3
"""Which high-value plant molecules have a producing species in the Algerian flora?

Pipeline (public, key-less APIs):
  1. name -> PubChem CID + InChIKey
  2. InChIKey -> Wikidata/LOTUS "found in taxon" (P703) -> producing taxa
  3. Algerian vascular-plant checklist from GBIF occurrences (country=DZ, Tracheophyta):
     every species and genus with >=1 record, with record counts
  4. intersect: producing species present in DZ (with GBIF record count) and producing genera present in DZ

Caveats written into the output:
  - LOTUS pairs are literature claims, some are artefacts (e.g. huperzine A "in" Citrullus colocynthis).
  - GBIF record count measures how often a species was recorded, not its abundance; cultivated plants count too.

Run:   python3 projects/scouting/scripts/02_high_value_dz_sources.py           (~10 min; caches in work/)
Output: projects/scouting/data/04_high_value_dz_sources.tsv
"""
import csv, json, os, time, urllib.parse, concurrent.futures as cf
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "..", "work")
OUT = os.path.join(HERE, "..", "data", "04_high_value_dz_sources.tsv")
S = requests.Session()
S.headers["User-Agent"] = "bagel-research/1.0 (dry-lab scouting)"

# High-value or supply-limited plant molecules (drugs, API precursors, rare reagents, fragrance materials).
NAMES = [
    "thapsigargin", "trilobolide", "nortrilobolide", "paclitaxel", "10-deacetylbaccatin III", "vinblastine",
    "vincristine", "camptothecin", "podophyllotoxin", "etoposide", "colchicine", "galantamine", "artemisinin",
    "digoxin", "homoharringtonine", "harringtonine", "triptolide", "ingenol mebutate", "pilocarpine",
    "scopolamine", "hyoscyamine", "anisodamine", "littorine", "emetine", "physostigmine", "huperzine A",
    "quinine", "reserpine", "sanguinarine", "tubocurarine", "noscapine", "withaferin A", "celastrol",
    "forskolin", "ginkgolide B", "cytisine", "sparteine", "lobeline", "parthenolide", "shikonin", "alkannin",
    "hypericin", "berberine", "maytansine", "ouabain", "calotropin", "cucurbitacin E", "cucurbitacin I",
    "rotenone", "deguelin", "silibinin", "psoralen", "xanthotoxin", "bergapten", "khellin", "visnagin",
    "peganine", "harmine", "solasodine", "diosgenin", "hecogenin", "glycyrrhizin", "mogroside V",
    "rebaudioside A", "crocin", "safranal", "picrocrocin", "bruceantin", "betulinic acid", "ursolic acid",
    "salvinorin A", "yohimbine", "vinorelbine", "anisodine", "cocaine", "ephedrine", "morphine", "codeine",
    "thebaine", "papaverine", "capsaicin", "cannabidiol", "tetrahydrocannabinol", "conessine", "usnic acid",
    "sclareol", "nootkatone", "cis-jasmone", "methyl jasmonate", "rosmarinic acid", "carnosic acid",
    "verbascoside", "echinacoside",
]
CONTROLLED = {"cocaine", "morphine", "codeine", "thebaine", "tetrahydrocannabinol", "ephedrine", "salvinorin A"}


def cached(name, fn):
    os.makedirs(WORK, exist_ok=True)
    p = os.path.join(WORK, name)
    if os.path.exists(p):
        return json.load(open(p))
    data = fn()
    json.dump(data, open(p, "w"))
    return data


def pubchem(name):
    u = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + urllib.parse.quote(name) + "/property/InChIKey/JSON"
    for a in range(3):
        try:
            r = S.get(u, timeout=30)
            if r.status_code == 200:
                p = r.json()["PropertyTable"]["Properties"][0]
                return p["CID"], p["InChIKey"]
            if r.status_code == 404:
                return None, None
        except Exception:
            pass
        time.sleep(1 + a)
    return None, None


def lotus_taxa(inchikey):
    q = ('SELECT ?name WHERE { ?c wdt:P235 "%s" . ?c p:P703 ?st . ?st ps:P703 ?t . ?t wdt:P225 ?name . }' % inchikey)
    for a in range(4):
        try:
            r = S.get("https://query.wikidata.org/sparql", params={"query": q, "format": "json"}, timeout=60)
            if r.status_code == 200:
                return sorted({b["name"]["value"] for b in r.json()["results"]["bindings"]})
        except Exception:
            pass
        time.sleep(2 + 2 * a)
    return None


def fetch_lotus():
    out = {}
    for n in NAMES:
        cid, ik = pubchem(n)
        out[n] = dict(cid=cid, inchikey=ik, taxa=lotus_taxa(ik) if ik else None)
        time.sleep(0.5)
    return out


def fetch_flora():
    def facet(field):
        r = S.get("https://api.gbif.org/v1/occurrence/search", timeout=300, params={
            "country": "DZ", "phylumKey": 7707728, "limit": 0, "facet": field, "facetLimit": 20000}).json()
        return [(int(c["name"]), c["count"]) for c in r["facets"][0]["counts"]]

    def name(k):
        for a in range(4):
            try:
                j = S.get(f"https://api.gbif.org/v1/species/{k}", timeout=60).json()
                return k, j.get("canonicalName") or j.get("scientificName"), j.get("family")
            except Exception:
                time.sleep(1 + a)
        return k, None, None

    out = {}
    for field in ("genusKey", "speciesKey"):
        f = facet(field)
        with cf.ThreadPoolExecutor(16) as ex:
            names = {k: (n, fam) for k, n, fam in ex.map(name, [k for k, _ in f])}
        out[field] = [dict(key=k, count=c, name=names[k][0], family=names[k][1]) for k, c in f]
    return out


def main():
    lotus = cached("lotus_taxa.json", fetch_lotus)
    flora = cached("dz_flora.json", fetch_flora)
    sp = {d["name"]: d for d in flora["speciesKey"] if d["name"]}
    ge = {d["name"] for d in flora["genusKey"] if d["name"]}
    rows = []
    for comp in NAMES:
        d = lotus.get(comp) or {}
        taxa = d.get("taxa") or []
        hits = sorted(((t, sp[t]["count"], sp[t]["family"]) for t in taxa if t in sp), key=lambda x: -x[1])
        genera = sorted({t.split()[0] for t in taxa if t.split()[0] in ge})
        rows.append(dict(
            compound=comp, pubchem_cid=d.get("cid"), inchikey=d.get("inchikey"), lotus_taxa=len(taxa),
            dz_species_n=len(hits), dz_species="; ".join(f"{t} ({c} GBIF rec, {f})" for t, c, f in hits[:15]),
            dz_genera="; ".join(genera), controlled_substance="yes" if comp in CONTROLLED else ""))
    with open(OUT, "w", newline="") as f:
        f.write("# generated by scripts/02_high_value_dz_sources.py. LOTUS pairs are literature claims (some artefacts); "
                "GBIF record counts measure recording effort, not abundance; cultivated species included.\n")
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} compounds; {sum(1 for r in rows if r['dz_species_n'])} with a producing species recorded in DZ "
          f"(flora: {len(sp)} species, {len(ge)} genera) -> {OUT}")


if __name__ == "__main__":
    main()
