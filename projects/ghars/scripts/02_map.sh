#!/usr/bin/env bash
# Map public WGS reads of Ghars + controls to the Barhee BC4 genome and keep chr14 hits.
# Then run 02_invertase_cnv.py to genotype the sucrose-locus invertase deletions.
#
#   WORK=work bash scripts/02_map.sh          (≈ 40 min on 16 threads + downloads, ~17 GB disk)
#
# Reads: first 20 M reads per mate from ENA (~5x depth) — enough to tell 0 / 0.5 / 1 copies over a
# 35-kb window (thousands of reads expected per window). Barhee has only 12.4 M pairs in total.
set -euo pipefail
WORK=${WORK:-work}
T=${THREADS:-16}
mkdir -p "$WORK"/{tools,ref,reads,paf,logs}
MM2="$WORK/tools/minimap2-2.28_x64-linux/minimap2"

[ -x "$MM2" ] || curl -sL https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -xj -C "$WORK/tools"

R=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/389/715/GCF_009389715.1_palm_55x_up_171113_PBpolish2nd_filt_p/GCF_009389715.1_palm_55x_up_171113_PBpolish2nd_filt_p
[ -s "$WORK/ref/barhee.fna.gz" ] || curl -s -o "$WORK/ref/barhee.fna.gz" "${R}_genomic.fna.gz"
[ -s "$WORK/ref/barhee.gff.gz" ] || curl -s -o "$WORK/ref/barhee.gff.gz" "${R}_genomic.gff.gz"
[ -s "$WORK/ref/barhee.sr.mmi" ] || env -u OMP_NUM_THREADS "$MM2" -x sr -t "$T" -d "$WORK/ref/barhee.sr.mmi" "$WORK/ref/barhee.fna.gz" 2> "$WORK/logs/index.log"

# run  ENA-path        cultivar (role)
SAMPLES="SRR2511369:SRR251/009:Rhars=Ghars(test)
SRR2577995:SRR257/005:Deglet_Noor_NYUAD(pre-registered_positive_control;_FAILED_->_label_check)
SRR2559387:SRR255/007:Medjool(negative_control)
SRR2577516:SRR257/006:Barhee(negative_control)"

for line in $SAMPLES; do
  run=${line%%:*}; rest=${line#*:}; dir=${rest%%:*}
  for m in 1 2; do
    f="$WORK/reads/${run}_${m}.fq.gz"
    [ -s "$f" ] || curl -s "https://ftp.sra.ebi.ac.uk/vol1/fastq/${dir}/${run}/${run}_${m}.fastq.gz" | gzip -dc | head -n 80000000 | gzip -1 > "$f"
  done
  out="$WORK/paf/${run}.chr14.paf"
  if [ ! -s "$out" ]; then
    env -u OMP_NUM_THREADS "$MM2" -x sr -t "$T" --secondary=no "$WORK/ref/barhee.sr.mmi" \
      "$WORK/reads/${run}_1.fq.gz" "$WORK/reads/${run}_2.fq.gz" 2> "$WORK/logs/${run}.map.log" \
      | awk -F'\t' -v o="$out" '{n++} $6=="NC_052405.1"{print > o} END{print n > "/dev/stderr"}' 2> "$WORK/logs/${run}.n_alignments"
  fi
  echo "$run mapped: $(wc -l < "$out") chr14 alignments of $(cat "$WORK/logs/${run}.n_alignments") total"
done

# Independent Deglet Noor positive controls, added AFTER the pre-registered one (SRR2577995) did not
# show the deletion (see docs/02-sugar-locus.md §3). Declared post hoc.
#   SRR6439416 = dnPdF, Weill Cornell Qatar (PRJNA427409), 151-bp PE: first 12 M reads per mate
#   SRR121605  = Deglet Noor female PI 4611, WCMQ 2011 (PRJNA40349), 84-bp single-end: whole file
for m in 1 2; do
  f="$WORK/reads/SRR6439416_${m}.fq.gz"
  [ -s "$f" ] || curl -s "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR643/006/SRR6439416/SRR6439416_${m}.fastq.gz" | gzip -dc | head -n 48000000 | gzip -1 > "$f"
done
[ -s "$WORK/reads/SRR121605.fq.gz" ] || curl -s -o "$WORK/reads/SRR121605.fq.gz" "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR121/SRR121605/SRR121605.fastq.gz"
map_extra() { run=$1; shift; out="$WORK/paf/${run}.chr14.paf"
  [ -s "$out" ] && return 0
  env -u OMP_NUM_THREADS "$MM2" -x sr -t "$T" --secondary=no "$WORK/ref/barhee.sr.mmi" "$@" 2> "$WORK/logs/${run}.map.log" \
    | awk -F'\t' -v o="$out" '{n++} $6=="NC_052405.1"{print > o} END{print n > "/dev/stderr"}' 2> "$WORK/logs/${run}.n_alignments"; }
map_extra SRR6439416 "$WORK/reads/SRR6439416_1.fq.gz" "$WORK/reads/SRR6439416_2.fq.gz"
map_extra SRR121605 "$WORK/reads/SRR121605.fq.gz"

python3 scripts/02_invertase_cnv.py "$WORK"
python3 scripts/02b_plot_invertase.py   # needs matplotlib
