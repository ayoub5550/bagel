"""Step 6: function-labelled phylogenies for OSC, CYP716, CYP72 and CSL candidates.

References are characterized enzymes (Swiss-Prot PE=1/2 or the papers named in the label);
labels carry the *published* activity, never the UniProt automatic name.
Anabasis sequences: best ORF per rnaSPAdes gene (>= MINLEN aa), plus the in-silico joined bAS-A model.
Alignment: pyfamsa; tree: VeryFastTree (-lg -gamma, SH-like local supports).
"""
import gzip, os, re, subprocess, urllib.request, csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyfamsa

os.makedirs("trees", exist_ok=True)
VFT = ".venv/lib/python3.13/site-packages/veryfasttree/bin/VeryFastTree-avx2"

# ---------- references (accession -> label)
REF = {
 "OSC": {
  # beta-amyrin synthases
  "A0AAW1L0L7": "bAS|Saponaria_SobAS1", "O82140": "bAS|Panax_PNY1", "O82146": "bAS|Panax_PNY2",
  "Q8W3Z1": "bAS|Betula_OSCBPY", "A8CDT2": "bAS|Bruguiera_BgbAS", "Q9MB42": "bAS|Glycyrrhiza_GgbAS1",
  "B6EXY6": "bAS|Arabidopsis_AtBAS", "A0A0S2IHL6": "bAS|Kalopanax_KsBAS1", "A0A0U2U4F3": "bAS|Barbarea_BAS1",
  "E7DN63": "bAS|Solanum_TTS1", "A0A067FI21": "bAS|Citrus_OSC2", "Q9LRH8": "bAS|Pisum_OSCPSY",
  "Q9LRH7": "mixed-amyrin|Pisum_OSCPSM", "E7DN64": "delta-amyrin|Solanum_TTS2",
  "A0A5B8NBN0": "tirucalladienol|Azadirachta_OSC1",
  # lupeol synthases
  "Q2XPU7": "LUP|Ricinus_LUS", "Q8W3Z2": "LUP|Betula_OSCBPW", "A8CDT3": "LUP|Bruguiera_LUS",
  "E2IUA9": "LUP|Kalanchoe_LUS", "Q764T8": "LUP|Glycyrrhiza_LUS1", "P0DXI1": "LUP|Ailanthus_OSC3",
  "A0A067DDU9": "LUP|Citrus_OSC3", "Q9C5M3": "LUP-multi|Arabidopsis_LUP1", "Q8RWT0": "amyrin-multi|Arabidopsis_LUP2",
  "P0C8Y0": "camelliol|Arabidopsis_CAMS1",
  # sterol cyclases (outgroup)
  "P38605": "CAS|Arabidopsis_CAS1", "O82139": "CAS|Panax_CAS", "Q8W3Z4": "CAS|Betula_CAS1",
  "Q9SXV6": "CAS|Glycyrrhiza_CAS1", "Q6BE25": "CAS|Cucurbita_CPX", "Q2XPU6": "CAS|Ricinus_CAS",
  "Q6Z2X6": "CAS|Oryza_CAS", "Q1G1A4": "LAS|Arabidopsis_LAS1",
 },
 "CYP716": {
  "A0AAW1NEA3": "C28ox|Saponaria_CYP716A378", "A0AAW1JA93": "C28ox+C16a|Saponaria_CYP716A379",
  "Q2MJ20": "C28ox|Medicago_CYP716A12", "F6H9N6": "C28ox|Vitis_CYP716A15", "A5BFI4": "C28ox|Vitis_CYP716A17",
  "I7C6E8": "C28ox|Panax_CYP716A52v2", "A0A3Q7HBJ5": "C28ox|Solanum_CYP716A44", "K4CEE8": "C28ox|Solanum_CYP716A46",
  "A0A0S2II38": "C28ox|Kalopanax_CYP716A94", "A0A0U2U8U5": "C28ox|Barbarea_CYP716A81",
  "A0A1I9Q5Z0": "C16b|Platycodon_CYP716A141", "H2DH16": "C12|Panax_CYP716A47", "I7CT85": "C6|Panax_CYP716A53v2",
  "A0A3Q7HS74": "C6b|Solanum_CYP716E26", "A0A0B4L1W8": "Maesa_CYP716A75",
  "Q50EK1": "Picea_CYP716B1", "Q93Z79": "outgroup|Arabidopsis_CYP714A1",
 },
 "CYP72": {
  "A0AAW1J8D7": "C23ox|Saponaria_CYP72A984", "Q2MJ19": "C23ox|Medicago_CYP72A68",
  "A0A0S2IHL2": "C23ox|Kalopanax_CYP72A397", "H1A981": "C30ox|Medicago_CYP72A63",
  "H1A988": "C30ox|Glycyrrhiza_CYP72A154", "Q2MJ21": "Medicago_CYP72A67", "A0A481NR20": "Barbarea_CYP72A552",
  "A0A517FNC4": "Paris_CYP72A616", "A0A517FNC9": "Trigonella_CYP72A613", "W8JWW3": "Catharanthus_CYP72A225",
  "Q9LUC5": "Arabidopsis_CYP72A15", "Q9LUC9": "Arabidopsis_CYP72A11",
  "Q9SHG5": "outgroup|Arabidopsis_CYP72C1", "Q9ASR3": "outgroup|Arabidopsis_CYP709B1",
 },
 "CSL": {
  "A0AAW1HA02": "C3-sugar|Saponaria_SoCSL1", "A0A7I6P8N6": "C3-GlcA|Glycyrrhiza_GuCSyGT",
  "A0A7I6PD81": "C3-GlcA|Glycine_GmCSyGT1", "A0A7I6P8A1": "C3-GlcA|Lotus_LjCSyGT",
  "A0A7I6PDX0": "C3-GlcA|Glycine_GmCSLM1", "A0A7I6PAZ0": "C3-GlcA|Glycine_GmCSLM2",
  "Q570S7": "CslG|Arabidopsis_CSLG1", "Q8VYR4": "CslG|Arabidopsis_CSLG2", "Q0WVN5": "CslG|Arabidopsis_CSLG3",
  "Q8VZK9": "CslE|Arabidopsis_CSLE1", "O80898": "CslB|Arabidopsis_CSLB1",
  "Q9M9M4": "outgroup|Arabidopsis_CSLD3",
 },
}

def load_sprot(accs):
    out = {}
    with gzip.open("ref/sprot.fa.gz", "rt") as fh:
        for r in SeqIO.parse(fh, "fasta"):
            a = r.id.split("|")[1]
            if a in accs:
                out[a] = str(r.seq)
    return out

def fetch_uniprot(acc):
    t = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=60).read().decode()
    return "".join(t.split("\n")[1:])

allacc = {a for d in REF.values() for a in d}
seqs = load_sprot(allacc)
for a in sorted(allacc - set(seqs)):
    seqs[a] = fetch_uniprot(a)
missing = allacc - set(seqs)
assert not missing, missing

# ---------- Anabasis candidates
genes = list(csv.DictReader(open("mine/genes.tsv"), delimiter="\t"))
orfseq = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse("mine/hits.faa", "fasta")}

def pick(fam, pred, minlen):
    out = []
    for g in genes:
        if g["family"] == fam and int(g["length_aa"]) >= minlen and pred(g):
            out.append((f"Aart_{g['gene']}_{g['length_aa']}aa_TPM{float(g['gene_tpm']):.0f}", orfseq[g["orf"]]))
    return out

# joined bAS-A model: g7405_i1 ORF (aa 1-555) + g2508_i2 ORF (internal Met..stop), 238-aa overlap at 100% identity
a = next(v for k, v in orfseq.items() if "_g7405_i1_ORF" in k)
b = next(v for k, v in orfseq.items() if "_g2508_i2_ORF" in k)
ov = b[:238]
i = a.find(ov[:30])
assert i > 0 and a[i:] == b[:len(a) - i], "overlap mismatch"
basA = a[:i] + b
SeqIO.write([SeqRecord(Seq(basA), id="Aart_bAS-A_joined_g7405+g2508", description=f"len={len(basA)} joined in silico")],
            "trees/bAS-A_joined.faa", "fasta")
print("bAS-A joined length:", len(basA), "junction at aa", i)

sets = {
 "OSC": pick("OSC", lambda g: g["gene"] not in ("g7405", "g2508"), 450) + [("Aart_bAS-A_joined_g7405+g2508_TPM~900", basA)],
 "CYP716": pick("CYP", lambda g: (g.get("named_hit") or "").startswith("CYP716"), 140),
 "CYP72": pick("CYP", lambda g: (g.get("named_hit") or "").startswith("CYP72A") and float(g["named_qcov"] or 0) > 80, 350),
 "CSL": pick("CSL", lambda g: not (g.get("sp_desc") or "").startswith(("Cellulose synthase A", "Probable cellulose synthase A")), 600),
}

for name, cands in sets.items():
    recs = [(f"REF|{lab}|{acc}", seqs[acc]) for acc, lab in REF[name].items()] + cands
    seqs_fa = [pyfamsa.Sequence(n.encode(), s.encode()) for n, s in recs]
    aln = pyfamsa.Aligner(guide_tree="upgma", threads=16).align(seqs_fa)
    with open(f"trees/{name}.aln.faa", "w") as fh:
        for s in aln:
            fh.write(f">{s.id.decode()}\n{s.sequence.decode()}\n")
    env = {k: v for k, v in os.environ.items() if k not in ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT")}
    env["OMP_NUM_THREADS"] = "16"
    subprocess.run([VFT, "-lg", "-gamma", "-threads", "16", "-out", f"trees/{name}.nwk", f"trees/{name}.aln.faa"],
                   check=True, env=env, stdout=subprocess.DEVNULL, stderr=open(f"trees/{name}.log", "w"))
    print(name, "n_ref", len(REF[name]), "n_anabasis", len(cands))

# ---------- OSC product-determining residue (Kushiro et al. 2000: Leu->Trp turns a lupeol synthase into a bAS)
# column = SobAS1 Trp257 (motif MWCYCR at 256-261); reported for every sequence in the OSC alignment.
recs = list(SeqIO.parse("trees/OSC.aln.faa", "fasta"))
ref = next(r for r in recs if "SobAS1" in r.id)
pos = 0
for col, ch in enumerate(str(ref.seq)):
    if ch != "-":
        pos += 1
        if pos == 257:
            wcol = col
            break
with open("trees/OSC_motif.tsv", "w") as fh:
    fh.write("sequence\tmotif_window\tresidue_at_SobAS1_W257\n")
    for r in recs:
        fh.write(f"{r.id}\t{str(r.seq)[wcol-4:wcol+6]}\t{str(r.seq)[wcol]}\n")
print("wrote trees/OSC_motif.tsv")
