#!/usr/bin/env python3
"""Who has already studied Ghars dates? Literature + public sequence-data audit.

Outputs (re-run freely; only public, key-less APIs are used):
  data/01_literature_ghars.tsv      papers whose title/abstract names the cultivar, with a layer label
  data/01_public_sequence_data.tsv  every NCBI BioSample of Phoenix dactylifera whose metadata names
                                    Ghars under any spelling (Ghars / Rhars / Gharss / Ghers), plus
                                    the per-country count of all date-palm BioSamples

Spelling matters: the only Ghars genome in NCBI is registered as "Rhars" (French transliteration).
Europe PMC fuzzy-matches "Ghars" to unrelated words, so every hit is re-checked with a word-boundary regex.
"""
import csv, json, re, time, collections, xml.etree.ElementTree as ET
import requests

S = requests.Session()
S.headers["User-Agent"] = "ghars-research/1.0 (mailto:research@example.org)"
NAME = re.compile(r"\b(ghars|rhars|gharss|ghers)\b", re.I)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def openalex(query, pages=3):
    out, cur = [], "*"
    for _ in range(pages):
        r = S.get("https://api.openalex.org/works", params={"search": query, "per-page": 100, "cursor": cur}, timeout=60).json()
        for w in r.get("results", []):
            ai = w.get("abstract_inverted_index") or {}
            ab = " ".join(t for p, t in sorted((p, t) for t, ps in ai.items() for p in ps))
            src = ((w.get("primary_location") or {}).get("source") or {})
            out.append(dict(year=w.get("publication_year"), title=w.get("title") or "", venue=src.get("display_name") or "",
                            doi=(w.get("doi") or "").replace("https://doi.org/", "").lower(), abstract=ab, found_via="openalex"))
        cur = r.get("meta", {}).get("next_cursor")
        if not cur:
            break
        time.sleep(0.2)
    return out


def europepmc(query):
    r = S.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
              params={"query": query, "format": "json", "pageSize": 200, "resultType": "core"}, timeout=60).json()
    return [dict(year=x.get("pubYear"), title=x.get("title", ""), venue=x.get("journalTitle", "") or "",
                 doi=(x.get("doi") or "").lower(), abstract=re.sub("<[^>]+>", "", x.get("abstractText", "")), found_via="europepmc")
            for x in r.get("resultList", {}).get("result", [])]


# ordered: first match wins
LAYERS = [
    ("genetics", r"\brapd\b|\bissr\b|\bssr\b|microsatell|genome|resequenc|chloroplast genome|mitochondrial genome|\bsnp"),
    ("in_silico", r"docking|in silico"),
    ("human", r"healthy (human )?subjects|post.?prandial|glyc(a)?emic|blood glucose"),
    ("animal", r"\brats?\b|\bmice\b|in vivo|hepatotox"),
    ("pests", r"pyral|ectomyelois|parlatoria|cochineal|cochenille|date moth|insect|entomolog|bird depredation|arthropod"),
    ("materials_energy", r"adsorp|activated carbon|corrosion|pyrolysis|graphene|desalination|\bdye\b|biomass fuel|bioethanol|alkylpolyglucoside|tensioactif"),
    ("microbial_biotech", r"citric acid|bacteriocin|enterococcus|aspergillus|fermentable sugars"),
    ("food_tech", r"syrup|sirop|bread|valoris|fructose-rich|browning|brunissement|polyphenoloxidase|functional ingredient"),
    ("chemistry_bioassay", r"phenol|flavono|antioxid|lc-|hplc|gc-?ms|volatil|tocopherol|fatty acid|seed oil|phytochem|chemical compos|biochemical"),
    ("agronomy", r"pollen|pollinat|farm|maturit|sun exposure|temperature|salinity|biodiversity|typology|inventaire|vegetative|bagging|ensachage|تكييس|palmiers mâles|éclaircissage"),
]
NOT_DATES = re.compile(r"ni.?ma|ne.?mah|udayan|khaja|hadith|ḥadīth|encyclopaedia|children|methods and protocols", re.I)
# "Ghars" is also an author surname (M.A. Ghars, plant stress papers) and a classical Arabic writer;
# those records are dropped by NOT_DATES or by the date-palm context filter.


def literature():
    recs = {}
    for q in ["Ghars date", "Ghars Phoenix dactylifera", "Ghars cultivar", '"Ghars"', "Rhars date palm",
              "Ghars dates paste", "Ghars transcriptome", "Ghars genome", "Ghars RNA-seq", "Ghars metabolomic",
              "Ghars proteomic", "Ghars microsatellite", "Ghars sucrose"]:
        for r in openalex(q):
            recs.setdefault(r["doi"] or r["title"].lower()[:80], r)
    for q in ['"Ghars"', '"Ghars" AND dactylifera', '"Rhars"']:
        for r in europepmc(q):
            recs.setdefault(r["doi"] or r["title"].lower()[:80], r)
    rows = []
    for r in recs.values():
        blob = f"{r['title']} {r['abstract']}"
        if r["title"].strip().lower() == "ghars" or not NAME.search(blob) or NOT_DATES.search(r["title"]) or not re.search(r"dactylifera|date|datte|palm|ghars|tamr|phoenicicole", blob, re.I):
            continue
        b = blob.lower()
        layer = next((name for name, rx in LAYERS if re.search(rx, b)), "other")
        rows.append(dict(year=r["year"], layer=layer, title=re.sub(r"\s+", " ", r["title"]), venue=r["venue"], doi=r["doi"], found_via=r["found_via"]))
    rows.sort(key=lambda x: (x["layer"], -int(x["year"] or 0)))
    with open("data/01_literature_ghars.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print("papers naming the cultivar:", len(rows), dict(collections.Counter(r["layer"] for r in rows)))


def biosamples():
    r = S.get(EUTILS + "esearch.fcgi", params={"db": "biosample", "term": "Phoenix dactylifera[Organism]", "retmax": 5000,
                                              "usehistory": "y", "retmode": "json"}, timeout=60).json()["esearchresult"]
    xml = ""
    for start in range(0, int(r["count"]), 500):
        time.sleep(0.5)
        xml += S.get(EUTILS + "efetch.fcgi", params={"db": "biosample", "query_key": r["querykey"], "WebEnv": r["webenv"],
                                                     "retstart": start, "retmax": 500}, timeout=120).text
    xml = "<BioSampleSet>" + re.sub(r"<\?xml[^>]*\?>|</?BioSampleSet>", "", xml) + "</BioSampleSet>"
    root = ET.fromstring(xml)
    hits, countries = [], collections.Counter()
    for b in root.findall("BioSample"):
        attrs = {a.get("attribute_name"): (a.text or "") for a in b.findall("Attributes/Attribute")}
        geo = (attrs.get("geo_loc_name") or attrs.get("country") or "").split(":")[0].strip()
        countries[geo or "unknown"] += 1
        blob = " ".join([b.findtext("Description/Title") or ""] + list(attrs.values()))
        if NAME.search(blob):
            acc = b.get("accession")
            time.sleep(0.4)
            runs = S.get("https://www.ebi.ac.uk/ena/portal/api/filereport", params={
                "accession": acc, "result": "read_run", "format": "tsv",
                "fields": "run_accession,study_accession,library_strategy,instrument_model,read_count,base_count"}, timeout=60).text.splitlines()[1:]
            for line in runs or [""]:
                p = (line.split("\t") + [""] * 7)[:7]
                hits.append(dict(biosample=acc, name_in_record=attrs.get("cultivar", ""), country=geo,
                                 submitter=b.findtext("Owner/Name") or "", released=b.get("publication_date", "")[:10],
                                 run=p[0], bioproject=p[1], strategy=p[2], instrument=p[3], read_pairs=p[4], bases=p[5]))
    with open("data/01_public_sequence_data.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(hits[0].keys()) if hits else ["none"], delimiter="\t")
        w.writeheader(); w.writerows(hits)
        f.write("\n# all Phoenix dactylifera BioSamples by country (" + time.strftime("%Y-%m-%d") + "): "
                + json.dumps(dict(countries.most_common()), ensure_ascii=False) + "\n")
    print("BioSamples scanned:", sum(countries.values()), "| naming Ghars:", len({h['biosample'] for h in hits}))
    for h in hits:
        print("  ", h)


if __name__ == "__main__":
    literature()
    biosamples()
