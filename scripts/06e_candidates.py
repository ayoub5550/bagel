"""Step 7: curated saponin-pathway candidate table for Anabasis articulata (SRR6435311 assembly).

Everything numeric is read from the pipeline outputs (mine/genes.tsv, quant/q/quant.sf, trees/*.nwk);
only the candidate list, the proposed step and the notes are curated here (and documented in docs/03 §6).

Usage (from the transcriptome work dir):  .venv/bin/python scripts/06e_candidates.py /path/to/repo
"""
import csv, re, sys, collections, os, shutil
from Bio import Phylo, SeqIO

REPO = sys.argv[1] if len(sys.argv) > 1 else "."
genes = {(r["gene"], r["family"]): r for r in csv.DictReader(open("mine/genes.tsv"), delimiter="\t")}
orfseq = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse("mine/hits.faa", "fasta")}
joined = {r.id: str(r.seq) for r in SeqIO.parse("trees/bAS-A_joined.faa", "fasta")}

gt = collections.defaultdict(float)
for r in csv.DictReader(open("quant/q/quant.sf"), delimiter="\t"):
    gt["g" + re.search(r"_g(\d+)_i\d+$", r["Name"]).group(1)] += float(r["TPM"])
ranked = sorted(gt, key=lambda g: -gt[g])
rank = {g: i + 1 for i, g in enumerate(ranked)}
n_expr = sum(1 for v in gt.values() if v > 0)

trees = {}
def tree(name):
    if name not in trees:
        t = Phylo.read(f"trees/{name}.nwk", "newick")
        outs = [c for c in t.get_terminals() if "outgroup" in c.name or c.name.startswith(("REF|CAS", "REF|LAS"))]
        t.root_with_outgroup(outs if len(outs) > 1 else outs[0])
        trees[name] = t
    return trees[name]

def support(name, leaf_sub, ref_sub):
    t = tree(name)
    a = next(c for c in t.get_terminals() if leaf_sub in c.name)
    b = next(c for c in t.get_terminals() if ref_sub in c.name)
    ca = t.common_ancestor(a, b)
    n = len(ca.get_terminals())
    return f"with {b.name.split('|')[2]} ({b.name.split('|')[1]}) in a clade of {n} leaves, SH-like support {ca.confidence}"

# id, genes(list; first = representative ORF), family, tree, ref-in-tree, proposed step, confidence, notes
C = [
 ("AaOSC-bAS-A", ["g7405", "g2508"], "OSC", "OSC", "SobAS1", "2,3-oxidosqualene -> beta-amyrin (oleanane branch)", "medium-high",
  "ORF joined in silico from two rnaSPAdes genes (g7405 aa 1-555 + g2508 internal-Met ORF) over a 238-aa overlap with 100% identity; bAS-type Trp at the Kushiro-2000 position (IWCYCR). Junction must be confirmed by RT-PCR."),
 ("AaOSC-bAS-B", ["g14739", "g10797"], "OSC", "OSC", "SobAS1", "beta-amyrin (or mixed) synthase", "medium",
  "3'-partial (aa 1-493). g14739/g10797 99.4% identical (alleles or assembly variants); 95-96% identical to bAS-A over the overlap = second paralog. Phe (IFCYCR) instead of the bAS Trp at the Kushiro-2000 position -> product not predictable."),
 ("AaOSC-LUP-1", ["g7093"], "OSC", "OSC", "Betula_OSCBPW", "2,3-oxidosqualene -> lupeol (lupane branch)", "medium",
  "Complete. Clusters with characterised lupeol synthases of Betula and Glycyrrhiza; carries the lupeol-type Leu (TLCYCR). Fits the lupane saponins of gamal2022."),
 ("AaOSC-LUP-2", ["g6626"], "OSC", "OSC", "Betula_OSCBPW", "lupeol (or other pentacyclic) synthase", "low-medium",
  "Complete. Sister to AaOSC-LUP-1 inside the LUP clade, but has Trp (NFWCYCR) at the Kushiro-2000 position -> product uncertain."),
 ("AaOSC-CAS", ["g3483"], "OSC", "OSC", "Ricinus_CAS", "cycloartenol (sterols; housekeeping control)", "high",
  "Complete; sterol cyclase, not saponin pathway. Included as an internal control (WCHCR motif like all CAS)."),
 ("AaCYP716A", ["g8081", "g28922", "g33474"], "CYP", "CYP716", "CYP716A378", "C-28 oxidation (beta-amyrin -> oleanolic acid; lupeol -> betulinic acid is a hypothesis)", "high",
  "Complete 483 aa; 81.5% identical to soapwort CYP716A378 (characterised C-28 oxidase, jo2025). g28922 and g33474 are 100%-identical fragments of the same gene; TPM below is the sum of the three."),
 ("AaCYP72A-1", ["g31"], "CYP", "CYP72", "CYP72A984", "C-23 oxidation (-> 23-CHO / 23-COOH)", "medium-high",
  "Complete 532 aa; 60.9% to soapwort CYP72A984 (C-23 methyl -> carboxyl, jo2025). Subfamily-level identity, but CYP72A functions vary (C-2, C-23, C-30). Matches the 23-aldehyde/23-acid lupane saponins of gamal2022."),
 ("AaCYP72A-2..5", ["g6754", "g8816", "g1975", "g8362"], "CYP", "CYP72", "CYP72A984", "paralogs of AaCYP72A-1 (function open)", "low",
  "Complete 521-532 aa, 56-60% to CYP72A984; expanded Anabasis-specific CYP72A cluster in the tree. Lower expression."),
 ("AaCSL-1", ["g936"], "CSL", "CSL", "SoCSL1", "C-3 first sugar (glucuronic acid by analogy)", "high",
  "Complete 704 aa; 68.4% to soapwort CSL1 (starts the C-3 sugar chain, jo2025). Cellulose-synthase-derived enzymes add the C-3 GlcA in legumes (chung2020) and quinoa (zhang2025). gamal2022 compounds 6-7 carry 3-O-GlcA."),
 ("AaCSL-M1..3", ["g6759", "g7475", "g7869"], "CSL", "CSL", "GmCSLM1", "alternative C-3 glucuronosyltransferases (CslM-type)", "low-medium",
  "Complete 748-762 aa; sister to soybean CSLM1/2 (beta-amyrin glucuronosyltransferases, jozwiak2020). Low expression."),
 ("AaUGT74-1", ["g8469"], "UGT", None, None, "C-28 ester first sugar (Glc)", "low-medium",
  "Complete 487 aa; 51.9% to soapwort UGT74CD1 (starts the C-28 chain, jo2025). Family-level identity only. Quinoa C-28 ester glucosyltransferase is also UGT74 (CqUGT74BB3, zhang2025)."),
 ("AaUGT73-1", ["g7054"], "UGT", None, None, "sugar-chain elongation (open)", "low",
  "Complete 500 aa; 41.5% to soapwort UGT73CC6 (branches the C-3 chain with Xyl). High expression, weak homology."),
 ("AaUGT73-2", ["g13108"], "UGT", None, None, "sugar-chain elongation (open)", "low",
  "Complete 510 aa; 50.1% to soapwort UGT73DL1 (adds Gal to the C-3 GlcA)."),
 ("AaUGT73-3", ["g12002"], "UGT", None, None, "sugar-chain elongation (open)", "low",
  "Complete 485 aa; 51.3% to soapwort UGT73M2 (C-28 chain Xyl)."),
 ("AaUGT79-1", ["g16827"], "UGT", None, None, "sugar-chain elongation (open)", "low",
  "3'-partial 418 aa; 55.2% to soapwort UGT79T1 (C-28 chain deoxyhexose)."),
]

def seq_for(cid, g, fam):
    if cid == "AaOSC-bAS-A":
        return next(iter(joined.values()))
    return orfseq[genes[(g, fam)]["orf"]]

import subprocess, tempfile
REFSEQ = {}
for n in ("OSC", "CYP716", "CYP72", "CSL"):
    for rec in SeqIO.parse(f"trees/{n}.aln.faa", "fasta"):
        if rec.id.startswith("REF|"):
            REFSEQ[rec.id] = str(rec.seq).replace("-", "")

def pid_to_ref(seq, ref_sub):
    rid = next(k for k in REFSEQ if ref_sub in k)
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/q.fa", "w").write(f">q\n{seq}\n"); open(f"{d}/r.fa", "w").write(f">r\n{REFSEQ[rid]}\n")
        subprocess.run(["tools/diamond", "makedb", "--in", f"{d}/r.fa", "-d", f"{d}/r", "--quiet"], check=True)
        out = subprocess.run(["tools/diamond", "blastp", "-d", f"{d}/r", "-q", f"{d}/q.fa", "--quiet", "--more-sensitive",
                              "-k", "1", "--max-hsps", "1", "--outfmt", "6", "pident", "qcovhsp", "scovhsp"],
                             check=True, capture_output=True, text=True).stdout.split()
    return rid.split("|")[2], float(out[0]), float(out[1]), float(out[2])

rows, faa = [], []
for cid, gl, fam, tname, ref, step, conf, notes in C:
    r = genes[(gl[0], fam)]
    s = seq_for(cid, gl[0], fam)
    tpm = sum(gt[g] for g in gl)
    if cid == "AaOSC-bAS-A":
        length, comp = len(s), "joined_in_silico"
    else:
        length = int(r["length_aa"]); comp = "complete" if r["orf_type"] == "complete" else "PARTIAL(" + r["orf_type"] + ")"
    leaf = "bAS-A_joined" if cid == "AaOSC-bAS-A" else f"Aart_{gl[0]}_"
    placement = support(tname, leaf, ref) if tname else "no tree (UGT: family-level only)"
    if tname:
        rname, pid, qcov, scov = pid_to_ref(s, ref)
    else:
        rname, pid, qcov, scov = (r.get("named_hit") or r["sp_desc"]), float(r.get("named_pid") or r["sp_pid"]), float(r.get("named_qcov") or r["sp_qcov"]), None
    rows.append(dict(candidate_id=cid, rnaspades_genes=";".join(gl), transcript_id=r["transcript"], length_aa=length,
                     complete_or_partial=comp, family=fam, pfam=r["pfam"],
                     reference=rname, reference_acc=(rname.split("|")[-1] if False else (next(k for k in REFSEQ if ref in k).split("|")[3] if tname else (r.get("named_acc") or r["sp_acc"]))),
                     pct_identity=pid, query_cov=qcov,
                     best_amaranthaceae_hit=f"{r['ref_acc']} {r['ref_org']} {r['ref_pid']}%" if "Spinacia" in r["ref_org"] or "Chenopodium" in r["ref_org"] or "Beta vulgaris" in r["ref_org"] else "",
                     tree_placement=placement, tpm_sum=round(tpm, 1),
                     expression_rank=f"{min(rank[g] for g in gl)} / {n_expr}", proposed_step=step, confidence=conf, notes=notes))
    if cid == "AaOSC-bAS-A":
        rows[-1]["transcript_id"] = "NODE_20735(g7405_i1) + NODE_5413(g2508_i2), joined"
        rows[-1]["pfam"] = "SQHop_cyclase_N + SQHop_cyclase_C (across the two parts)"
    faa.append(f">{cid} {';'.join(gl)} len={len(s)} {comp}\n{s}\n")

os.makedirs(f"{REPO}/data/06_trees", exist_ok=True)
with open(f"{REPO}/data/06_pathway_candidates.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
open(f"{REPO}/data/06_pathway_candidates.faa", "w").write("".join(faa))
for n in ("OSC", "CYP716", "CYP72", "CSL"):
    shutil.copy(f"trees/{n}.nwk", f"{REPO}/data/06_trees/{n}.nwk")
    shutil.copy(f"trees/{n}.aln.faa", f"{REPO}/data/06_trees/{n}.aln.faa")
for r in rows:
    print(r["candidate_id"], r["length_aa"], r["complete_or_partial"], r["pct_identity"], r["tpm_sum"], r["expression_rank"], "|", r["tree_placement"])
