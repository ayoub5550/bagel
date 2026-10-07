#!/usr/bin/env python3
"""Which Saharan / arid Algerian plant is still open? Literature + public sequence-data gap scan.

For each candidate (all accepted names and common synonyms are OR-ed) count:
  epmc        Europe PMC hits (title, abstract AND open full text — reviews that merely mention a name count too)
  epmc_dz     same, AND (Algeria OR Algerian OR Algérie)
  epmc_omics  same, AND (transcriptome OR RNA-seq OR genome sequence/assembly OR biosynthe*)
  pubmed      PubMed title/abstract hits
  sra         NCBI SRA runs ([Organism])
  assembly    NCBI Assembly records for the species; genus_assembly for the whole genus
  nuccore     NCBI Nucleotide records

Counts are a map of attention, not of value: the value column in projects/scouting/README.md is
judgement, with sources. Re-run freely (key-less public APIs, ~3 min):
  python3 projects/scouting/scripts/01_desert_gap_scan.py
Output: projects/scouting/data/01_desert_gap_scan.tsv
"""
import csv, os, time, threading, concurrent.futures as cf, datetime
import requests

SPECIES = [  # (key, Arabic/local name, names searched)
    ('Cistanche phelypaea', 'الذنون', ['Cistanche phelypaea', 'Cistanche tinctoria', 'Cistanche lutea']),
    ('Cistanche violacea', 'الذنون البنفسجي', ['Cistanche violacea']),
    ('Cistanche tubulosa', 'الذنون (الصيني/الشرق)', ['Cistanche tubulosa']),
    ('Thapsia garganica', 'الدرياس/بونافع', ['Thapsia garganica']),
    ('Ammodaucus leucotrichus', 'المسوفة/كمون الصحراء', ['Ammodaucus leucotrichus']),
    ('Cymbopogon schoenanthus', 'الإذخر', ['Cymbopogon schoenanthus']),
    ('Hyoscyamus muticus', 'السكران/فلزلز', ['Hyoscyamus muticus', 'Hyoscyamus falezlez']),
    ('Balanites aegyptiaca', 'الهجليج/تيبوراق', ['Balanites aegyptiaca']),
    ('Citrullus colocynthis', 'الحنظل', ['Citrullus colocynthis']),
    ('Peganum harmala', 'الحرمل', ['Peganum harmala']),
    ('Calotropis procera', 'الكرنكة/العشار', ['Calotropis procera']),
    ('Pergularia tomentosa', 'الغلقة', ['Pergularia tomentosa']),
    ('Zygophyllum album', 'العقاية', ['Zygophyllum album', 'Tetraena alba']),
    ('Haloxylon scoparium', 'الرمث', ['Haloxylon scoparium', 'Hammada scoparia']),
    ('Retama raetam', 'الرتم', ['Retama raetam']),
    ('Artemisia herba-alba', 'الشيح', ['Artemisia herba-alba']),
    ('Artemisia judaica', 'البعيثران', ['Artemisia judaica']),
    ('Ziziphus lotus', 'السدرة/النبق', ['Ziziphus lotus']),
    ('Pistacia atlantica', 'البطم', ['Pistacia atlantica']),
    ('Argania spinosa', 'الأرقان (تندوف)', ['Argania spinosa']),
    ('Olea europaea subsp. laperrinei', 'زيتون لابرين (الهقار)', ['Olea europaea subsp. laperrinei', 'Olea laperrinei', 'laperrinei']),
    ('Cupressus dupreziana', 'سرو الطاسيلي (تارووت)', ['Cupressus dupreziana']),
    ('Myrtus nivellei', 'ريحان الطوارق', ['Myrtus nivellei']),
    ('Lawsonia inermis', 'الحنة', ['Lawsonia inermis']),
    ('Salvadora persica', 'الأراك', ['Salvadora persica']),
    ('Anvillea radiata', 'النقد', ['Anvillea radiata', 'Anvillea garcinii']),
    ('Deverra scoparia', 'القزاح', ['Deverra scoparia', 'Pituranthos scoparius']),
    ('Calligonum comosum', 'العرطة', ['Calligonum comosum', 'Calligonum polygonoides']),
    ('Atriplex halimus', 'القطف', ['Atriplex halimus']),
    ('Solenostemma argel', 'الحرجل', ['Solenostemma argel', 'Solenostemma arghel']),
    ('Senna italica', 'السنا', ['Senna italica', 'Cassia italica']),
    ('Nitraria retusa', 'الغرقد', ['Nitraria retusa']),
    ('Euphorbia guyoniana', 'لبينة', ['Euphorbia guyoniana']),
    ('Ferula vesceritensis', 'الكلخ الصحراوي', ['Ferula vesceritensis']),
    ('Matricaria pubescens', 'الوزيعة', ['Matricaria pubescens']),
    ('Saccocalyx satureioides', 'الزعتر الصحراوي', ['Saccocalyx satureioides']),
    ('Withania adpressa', '—', ['Withania adpressa']),
    ('Cotula cinerea', 'الشويحية/القرطوفة', ['Cotula cinerea', 'Brocchia cinerea']),
    ('Bubonium graveolens', 'التافس', ['Bubonium graveolens', 'Nauplius graveolens', 'Asteriscus graveolens']),
    ('Fredolia aretioides', 'الدّقف', ['Fredolia aretioides']),
    ('Oudneya africana', 'الحنّة الجمل', ['Oudneya africana']),
    ('Limoniastrum guyonianum', 'الزيتة', ['Limoniastrum guyonianum']),
    ('Rhanterium adpressum', 'العرفج', ['Rhanterium adpressum']),
    ('Tamarix aphylla', 'الأثل', ['Tamarix aphylla']),
    ('Ephedra alata', 'العلندة', ['Ephedra alata']),
    ('Capparis spinosa', 'الكبار', ['Capparis spinosa']),
    ('Glycyrrhiza glabra', 'عرق السوس', ['Glycyrrhiza glabra']),
    ('Apteranthes europaea', '—', ['Caralluma europaea', 'Apteranthes europaea']),
    ('Vachellia tortilis raddiana', 'الطلح', ['Acacia raddiana', 'Vachellia tortilis', 'Acacia tortilis']),
    ('Terfezia (fungus, not plant)', 'الترفاس', ['Terfezia', 'Tirmania']),
]

S = requests.Session()
S.headers["User-Agent"] = "bagel-research/1.0"
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_lock, _last = threading.Lock(), [0.0]


def ncbi(db, term):
    for attempt in range(6):
        with _lock:  # NCBI allows 3 requests/s without a key
            wait = 0.36 - (time.time() - _last[0])
            if wait > 0:
                time.sleep(wait)
            _last[0] = time.time()
        try:
            r = S.get(EU, params=dict(db=db, term=term, retmode="json", retmax=0), timeout=60)
            if r.status_code == 200:
                return int(r.json()["esearchresult"]["count"])
        except Exception:
            pass
        time.sleep(1 + attempt)
    return -1


def epmc(query):
    for attempt in range(4):
        try:
            r = S.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                      params=dict(query=query, format="json", pageSize=1), timeout=60)
            return int(r.json()["hitCount"])
        except Exception:
            time.sleep(1 + attempt)
    return -1


def scan(entry):
    key, ar, names = entry
    q = " OR ".join(f'"{n}"' for n in names)
    org = " OR ".join(f'"{n}"[Organism]' for n in names)
    tiab = " OR ".join(f'"{n}"[tiab]' for n in names)
    genus = names[0].split()[0]
    return dict(
        species=key, local_name=ar, names_searched="; ".join(names),
        epmc=epmc(f"({q})"),
        epmc_dz=epmc(f"({q}) AND (Algeria OR Algerian OR Algérie)"),
        epmc_omics=epmc(f'({q}) AND (transcriptome OR RNA-seq OR "genome sequence" OR "genome assembly" OR biosynthe*)'),
        pubmed=ncbi("pubmed", tiab), sra=ncbi("sra", org), assembly=ncbi("assembly", org),
        genus_assembly=ncbi("assembly", f'"{genus}"[Organism]'), nuccore=ncbi("nuccore", org))


def main():
    with cf.ThreadPoolExecutor(6) as ex:
        rows = list(ex.map(scan, SPECIES))
    out = os.path.join(os.path.dirname(__file__), "..", "data", "01_desert_gap_scan.tsv")
    with open(out, "w", newline="") as f:
        f.write(f"# generated by scripts/01_desert_gap_scan.py on {datetime.date.today()}; -1 = query failed\n")
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in sorted(rows, key=lambda r: r["epmc"]):
        print(f"{r['species'][:34]:34s} epmc={r['epmc']:5d} dz={r['epmc_dz']:4d} omics={r['epmc_omics']:4d} sra={r['sra']:4d} asm={r['assembly']}")


if __name__ == "__main__":
    main()
