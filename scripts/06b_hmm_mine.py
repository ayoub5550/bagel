"""Step 4: HMM mining of saponin-pathway families in the Anabasis articulata ORFs.

Input : orf/orfs.faa (orfipy, both strands, >=100 aa, partials allowed)
Output: mine/hmm_hits.tsv (per ORF x HMM best domain), mine/hits.faa
"""
import os, re, csv, collections
import pyhmmer
from pyhmmer.easel import SequenceFile, Alphabet
from pyhmmer.plan7 import HMMFile

ALPHA = Alphabet.amino()
S = lambda x: x.decode() if isinstance(x, bytes) else x
IEVAL = 1e-5
os.makedirs("mine", exist_ok=True)

hmms = []
for f in sorted(os.listdir("hmm")):
    with HMMFile(f"hmm/{f}") as hf:
        for h in hf:
            hmms.append(h)
print("HMMs:", [(S(h.accession), S(h.name), h.M) for h in hmms])

with SequenceFile("orf/orfs.faa", digital=True, alphabet=ALPHA) as sf:
    seqs = sf.read_block()
print("ORFs:", len(seqs))

rows = []
for top in pyhmmer.hmmsearch(hmms, seqs, cpus=16, E=1e-3, domE=1e-3):
    qname = S(top.query.name)
    for hit in top:
        best = None
        for d in hit.domains:
            if d.i_evalue < IEVAL and (best is None or d.score > best.score):
                best = d
        if best is None:
            continue
        ali = best.alignment
        rows.append(dict(orf=S(hit.name), hmm=qname, evalue=hit.evalue, score=hit.score,
                         dom_ievalue=best.i_evalue, hmm_from=ali.hmm_from, hmm_to=ali.hmm_to,
                         hmm_len=top.query.M,
                         env_from=best.env_from, env_to=best.env_to))
print("hits:", len(rows))
cnt = collections.Counter(r["hmm"] for r in rows)
print(cnt)
with open("mine/hmm_hits.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

keep = {r["orf"] for r in rows}
n = 0
with open("orf/orfs.faa") as fi, open("mine/hits.faa", "w") as fo:
    write = False
    for line in fi:
        if line.startswith(">"):
            write = line[1:].split()[0] in keep
            n += write
        if write:
            fo.write(line)
print("wrote", n, "seqs")
