#!/usr/bin/env bash
# 06_transcriptome_assembly.sh — first de-novo transcriptome of Anabasis articulata (SRR6435311)
#
# Steps 1-3 below are the EXACT commands that were run on 2026-10-07 (download, QC, assembly),
# plus the reference/HMM downloads used for gene mining. Steps 4-7 (ORFs, HMM + DIAMOND, tree,
# salmon) are NOT automated yet — their plan is in HANDOFF.md §2.
#
# Idempotent: every step is skipped if its output already exists. Large files go to $WORK,
# which is git-ignored (work/). Needs: curl, python3 (>=3.8), uv, ~30 GB disk, >=64 GB RAM.
#
# Usage:   bash scripts/06_transcriptome_assembly.sh
#          WORK=/path/to/scratch THREADS=16 MEM_GB=300 bash scripts/06_transcriptome_assembly.sh
# Resume a killed assembly:   (cd $WORK && tools/SPAdes-4.0.0-Linux/bin/rnaspades.py --continue -o asm)
set -euo pipefail

WORK=${WORK:-work/tx}
THREADS=${THREADS:-$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)}
THREADS=$(( THREADS > 16 ? 16 : THREADS ))
MEM_GB=${MEM_GB:-300}
SS=${SS:-rf}   # strandedness passed to rnaSPAdes. UNVERIFIED for this library — see HANDOFF.md §2 caveat 1.
               # Use SS=none to assemble as unstranded.

mkdir -p "$WORK"/{raw,qc,tools,logs,hmm,ref}
cd "$WORK"

# --- tools (static binaries; versions used on 2026-10-07) -------------------------------------
cd tools
[ -x fastp ] || { curl -sSL -o fastp http://opengene.org/fastp/fastp && chmod +x fastp; }          # 1.4.0
[ -d SPAdes-4.0.0-Linux ] || curl -sSL https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz | tar xz
[ -x seqkit ] || curl -sSL https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz | tar xz
[ -x diamond ] || curl -sSL https://github.com/bbuchfink/diamond/releases/download/v2.1.10/diamond-linux64.tar.gz | tar xz
ls -d salmon*/ >/dev/null 2>&1 || curl -sSL https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz | tar xz
cd ..
[ -d .venv ] || { uv venv -q .venv && uv pip install -q -p .venv pyhmmer pyfamsa biopython orfipy veryfasttree; }

# --- 1. download from ENA (no SRA toolkit needed) ---------------------------------------------
# ENA filereport: https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR6435311&result=read_run&fields=fastq_ftp,fastq_md5
B=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR643/001/SRR6435311
[ -s raw/R1.fastq.gz ] || curl -sS -o raw/R1.fastq.gz $B/SRR6435311_1.fastq.gz
[ -s raw/R2.fastq.gz ] || curl -sS -o raw/R2.fastq.gz $B/SRR6435311_2.fastq.gz
md5sum raw/R1.fastq.gz raw/R2.fastq.gz | tee raw/md5.txt
# expected (2026-10-07): b66db361ec6541d9414d4c97641720c8 R1 (1,483,376,308 B)
#                        bb06e1ed5f2166f0374368ad3ebc341c R2 (1,491,351,188 B)

# --- 2. QC + trimming ----------------------------------------------------------------------------
# NOTE: base qualities in this run are binned (every base = '?' = Q30), so -q 20 filters nothing;
# only adapters, N content and length actually filter. Expected: 64,188,264 -> 63,972,036 reads.
[ -s qc/R2.trim.fq.gz ] || tools/fastp -i raw/R1.fastq.gz -I raw/R2.fastq.gz \
    -o qc/R1.trim.fq.gz -O qc/R2.trim.fq.gz --detect_adapter_for_pe -q 20 -l 36 -w "$THREADS" \
    -j fastp.json -h fastp.html > logs/fastp.log 2>&1

# --- 3. assembly (rnaSPAdes 4.0.0; k = 41,61 chosen automatically from read length) -------------
# Slow in a gVisor sandbox: the K41 graph alone took > 1 h with 16 threads. Run it in the
# background (nohup) — a tool/session timeout that kills the shell kills the assembler too.
SS_ARG=(--ss "$SS"); [ "$SS" = none ] && SS_ARG=()
if [ ! -s asm/transcripts.fasta ]; then
  tools/SPAdes-4.0.0-Linux/bin/rnaspades.py "${SS_ARG[@]}" -1 qc/R1.trim.fq.gz -2 qc/R2.trim.fq.gz \
      -t "$THREADS" -m "$MEM_GB" -o asm > logs/rnaspades.stdout 2>&1
fi
tools/seqkit stats -a -T asm/transcripts.fasta | tee assembly_stats.txt

# --- references for gene mining (used in steps 5-6, see HANDOFF.md) ---------------------------
for p in PF13243 PF13249 PF00067 PF00201 PF01397 PF03936 PF02458; do
  [ -s hmm/$p.hmm ] || curl -sSL "https://www.ebi.ac.uk/interpro/api/entry/pfam/$p?annotation=hmm" | gunzip > hmm/$p.hmm
done
cd ref
[ -s sprot.fa.gz ] || curl -sS -o sprot.fa.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
# proteomes: quinoa UP000596660, spinach UP000813463, beet UP000035740 (same family, Amaranthaceae)
for u in UP000596660 UP000813463 UP000035740; do
  [ -s $u.fa.gz ] || curl -sS -o $u.fa.gz "https://rest.uniprot.org/uniprotkb/stream?query=proteome:$u&format=fasta&compressed=true"
done
[ -s refdb.dmnd ] || { zcat sprot.fa.gz UP*.fa.gz > refdb.fa && ../tools/diamond makedb --in refdb.fa -d refdb -p "$THREADS"; }
# Curated function sets (bAS.fa, CAS.fa, LUP.fa, CYP716.fa, CYP72.fa, CYPother.fa, UGT.fa) were pulled
# from the UniProt REST search API. WARNING: the query 'protein_name:"beta-amyrin synthase"' also
# returns lupeol/delta-amyrin synthases — label every reference by its Swiss-Prot function before
# building a tree, never by the search term that fetched it.
echo "steps 1-3 done; continue with HANDOFF.md §2 (steps 4-7)"
