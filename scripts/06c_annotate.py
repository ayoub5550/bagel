"""Step 5: annotate HMM hits -> gene-level candidate table.

Inputs : mine/hmm_hits.tsv, mine/hits.faa, mine/hits.sprot.tsv (DIAMOND vs Swiss-Prot, k=25),
         mine/hits.dmnd.tsv (DIAMOND vs refdb incl. quinoa/spinach/beet proteomes, k=5), quant/q/quant.sf
Outputs: mine/orfs_annot.tsv (one row per ORF), mine/genes.tsv (best ORF per rnaSPAdes gene x family)
"""
import csv, re, collections, statistics

FAMS = [("SQHop_cyclase_N", "OSC"), ("SQHop_cyclase_C", "OSC"), ("p450", "CYP"), ("UDPGT", "UGT"),
        ("Cellulose_synt", "CSL"), ("Transferase", "BAHD"), ("Terpene_synth", "TPS"), ("Terpene_synth_C", "TPS")]
FAMMAP = dict(FAMS)

def family(hmmset):
    for h, f in FAMS:
        if h in hmmset:
            return f
    return "?"

# --- ORF headers
orf = {}
for l in open("mine/hits.faa"):
    if not l.startswith(">"):
        continue
    p = l[1:].split()
    oid = p[0]
    tx = oid.rsplit("_ORF.", 1)[0]
    m = re.search(r"_g(\d+)_i(\d+)$", tx)
    coords = re.match(r"\[(\d+)-(\d+)\]\((.)\)", p[1])
    typ = p[2].split(":")[1]
    nt = int(p[3].split(":")[1])
    orf[oid] = dict(orf=oid, transcript=tx, gene=f"g{m.group(1)}", iso=int(m.group(2)),
                    tx_len=int(re.search(r"_length_(\d+)_", tx).group(1)),
                    strand=coords.group(3), orf_type=typ, length_aa=nt // 3 - (1 if typ in ("complete", "5-prime-partial") else 0))

hmm = collections.defaultdict(dict)
for r in csv.DictReader(open("mine/hmm_hits.tsv"), delimiter="\t"):
    hmm[r["orf"]][r["hmm"]] = r

# --- Swiss-Prot hits: best overall + best "named" hit (CYP/UGT nomenclature or curated function)
NAME_RE = [
    ("CYP", re.compile(r"\bGN=(CYP\d+[A-Z]*\d*)")),
    ("CYP", re.compile(r"Cytochrome P450 (\d+[A-Z]+\d+)")),
    ("UGT", re.compile(r"\bGN=(UGT\d+[A-Z]+\d*)")),
    ("UGT", re.compile(r"(?:glycosyltransferase|glucosyl transferase|glucosyltransferase) (\d+[A-Z]+\d+)", re.I)),
]

def parse_title(t):
    os_ = re.search(r"OS=(.+?) OX=", t)
    desc = t.split(" OS=")[0].split(" ", 1)[1] if " " in t else t
    return desc, (os_.group(1) if os_ else "")

def named(t):
    for fam, rx in NAME_RE:
        m = rx.search(t)
        if m:
            n = m.group(1)
            if fam == "CYP" and not n.startswith("CYP"):
                n = "CYP" + n
            if fam == "UGT" and not n.startswith("UGT"):
                n = "UGT" + n
            return n
    return None

sp = collections.defaultdict(list)
for l in open("mine/hits.sprot.tsv"):
    f = l.rstrip("\n").split("\t")
    sp[f[0]].append(f)
rd = collections.defaultdict(list)
for l in open("mine/hits.dmnd.tsv"):
    f = l.rstrip("\n").split("\t")
    rd[f[0]].append(f)

# --- expression
tpm_tx = {}
for r in csv.DictReader(open("quant/q/quant.sf"), delimiter="\t"):
    tpm_tx[r["Name"]] = float(r["TPM"])
gene_tpm = collections.defaultdict(float)
for tx, v in tpm_tx.items():
    g = "g" + re.search(r"_g(\d+)_i\d+$", tx).group(1)
    gene_tpm[g] += v
expressed = sorted(v for v in gene_tpm.values() if v > 0)

def pct_rank(v):
    import bisect
    return round(100.0 * bisect.bisect_left(expressed, v) / len(expressed), 1)

rows = []
for oid, o in orf.items():
    hs = hmm.get(oid, {})
    if not hs:
        continue
    fam = family(set(hs))
    best_sp = sp[oid][0] if sp[oid] else None
    best_named = next((h for h in sp[oid] if named(h[10])), None)
    best_rd = rd[oid][0] if rd[oid] else None
    r = dict(o)
    r.update(family=fam, pfam=";".join(sorted(f"{k}:{float(v['dom_ievalue']):.1e}" for k, v in hs.items())))
    if best_sp:
        d, os_ = parse_title(best_sp[10])
        r.update(sp_acc=best_sp[1].split("|")[1], sp_desc=d, sp_org=os_, sp_pid=float(best_sp[2]),
                 sp_qcov=float(best_sp[6]), sp_scov=float(best_sp[7]), sp_evalue=best_sp[8], sp_slen=int(best_sp[5]))
    if best_named:
        r.update(named_hit=named(best_named[10]), named_acc=best_named[1].split("|")[1],
                 named_org=parse_title(best_named[10])[1], named_pid=float(best_named[2]),
                 named_qcov=float(best_named[6]))
    if best_rd:
        d, os_ = parse_title(best_rd[10])
        r.update(ref_acc=best_rd[1].split("|")[1], ref_desc=d, ref_org=os_, ref_pid=float(best_rd[2]),
                 ref_qcov=float(best_rd[6]), ref_slen=int(best_rd[5]))
    r.update(tx_tpm=round(tpm_tx.get(o["transcript"], 0.0), 2), gene_tpm=round(gene_tpm[o["gene"]], 2),
             gene_tpm_pct=pct_rank(gene_tpm[o["gene"]]))
    rows.append(r)

cols = ["orf", "transcript", "gene", "iso", "tx_len", "strand", "orf_type", "length_aa", "family", "pfam",
        "sp_acc", "sp_desc", "sp_org", "sp_pid", "sp_qcov", "sp_scov", "sp_slen", "sp_evalue",
        "named_hit", "named_acc", "named_org", "named_pid", "named_qcov",
        "ref_acc", "ref_desc", "ref_org", "ref_pid", "ref_qcov", "ref_slen",
        "tx_tpm", "gene_tpm", "gene_tpm_pct"]
with open("mine/orfs_annot.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in sorted(rows, key=lambda r: (r["family"], -r["length_aa"])):
        w.writerow(r)

# gene-level: per (gene, family) keep longest ORF on + strand (stranded RF library => sense = '+'), else longest
best = {}
for r in rows:
    k = (r["gene"], r["family"])
    key = (r["strand"] == "+", r["length_aa"], r["tx_tpm"])
    if k not in best or key > best[k][0]:
        best[k] = (key, r)
genes = [v[1] for v in best.values()]
niso = collections.Counter((r["gene"], r["family"]) for r in rows)
with open("mine/genes.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols + ["n_orfs"], delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in sorted(genes, key=lambda r: (r["family"], -r["length_aa"])):
        r["n_orfs"] = niso[(r["gene"], r["family"])]
        w.writerow(r)

print("expressed genes:", len(expressed), "median TPM:", round(statistics.median(expressed), 2))
print("ORF rows:", len(rows), "gene x family rows:", len(genes))
print(collections.Counter(r["family"] for r in genes))
print("minus-strand ORFs:", sum(r["strand"] == "-" for r in rows))
