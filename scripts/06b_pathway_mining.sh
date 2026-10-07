#!/usr/bin/env bash
# 06b_pathway_mining.sh — steps 4-7 on the rnaSPAdes assembly made by 06_transcriptome_assembly.sh
#
# Exact sequence run on 2026-10-07 (05:00-05:10 UTC, 16 threads, ~10 min wall time in total):
#   4. ORFs (orfipy, both strands, >=100 aa, partials kept)
#   5. HMM mining (pyhmmer, domain i-Evalue < 1e-5) + DIAMOND vs Swiss-Prot and vs refdb (quinoa/spinach/beet)
#   6. salmon quantification (library type auto-detected; expected ISR for this RF library)
#   7. annotation table, function-labelled trees (pyfamsa + VeryFastTree), curated candidate table
# Outputs committed to the repo: data/06_pathway_candidates.{tsv,faa}, data/06_trees/*
# Usage: WORK=work/tx bash scripts/06b_pathway_mining.sh      (run from the repo root)
set -euo pipefail
REPO=$(pwd)
WORK=${WORK:-work/tx}
THREADS=${THREADS:-16}
cd "$WORK"
mkdir -p orf mine trees quant logs
NOOMP=(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT OMP_NUM_THREADS="$THREADS")   # see AGENTS.md §5

# 4. ORFs
[ -s orf/orfs.faa ] || "${NOOMP[@]}" .venv/bin/orfipy asm/transcripts.fasta --pep orfs.faa --min 300 \
    --strand b --partial-3 --partial-5 --procs "$THREADS" --outdir orf > logs/orfipy.log 2>&1
# extra Pfam model for cellulose-synthase-like glycosyltransferases (C-3 glucuronidation; chung2020, jo2025)
[ -s hmm/PF03552.hmm ] || curl -sSL "https://www.ebi.ac.uk/interpro/api/entry/pfam/PF03552?annotation=hmm" | gunzip > hmm/PF03552.hmm

# 5. HMM + DIAMOND
"${NOOMP[@]}" .venv/bin/python "$REPO/scripts/06b_hmm_mine.py"
[ -s ref/sprot.dmnd ] || tools/diamond makedb --in ref/sprot.fa.gz -d ref/sprot -p "$THREADS" > /dev/null 2>&1
OF="6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore stitle"
tools/diamond blastp -d ref/sprot -q mine/hits.faa --more-sensitive -e 1e-10 -k 25 -p "$THREADS" -o mine/hits.sprot.tsv --outfmt $OF 2>/dev/null
tools/diamond blastp -d ref/refdb -q mine/hits.faa --more-sensitive -e 1e-10 -k 5 -p "$THREADS" -o mine/hits.dmnd.tsv --outfmt $OF 2>/dev/null

# 6. salmon (mapping rate 94.99%, detected library type ISR on 2026-10-07)
S=$(ls -d tools/salmon*/bin/salmon | head -1)
[ -d quant/idx ] || "${NOOMP[@]}" "$S" index -t asm/transcripts.fasta -i quant/idx -p "$THREADS" > logs/salmon_index.log 2>&1
[ -s quant/q/quant.sf ] || "${NOOMP[@]}" "$S" quant -i quant/idx -l A -1 qc/R1.trim.fq.gz -2 qc/R2.trim.fq.gz \
    -p "$THREADS" --validateMappings -o quant/q > logs/salmon.log 2>&1

# 7. annotation, trees, candidates
.venv/bin/python "$REPO/scripts/06c_annotate.py"
"${NOOMP[@]}" .venv/bin/python "$REPO/scripts/06d_trees.py"
.venv/bin/python "$REPO/scripts/06e_candidates.py" "$REPO"
