#!/usr/bin/env python3
"""Genotype the sucrose-locus invertase deletions (Hazzouri 2019; Malek 2020) from public
whole-genome resequencing reads by normalised read depth.

Design (decided before looking at the data):
  * reference  : Barhee BC4 (GCF_009389715.1), the assembly used by Hazzouri 2019.
  * loci       : DEL1 ~40 kb at chr14 (NC_052405.1) 2.467-2.507 Mb, removes the 5' half of
                 CWINV1 (LOC103698975); DEL2 ~5 kb at 3.088-3.093 Mb, A/N-INV1 promoter.
  * controls   : Deglet Noor  = sucrose-type, published homozygous deletion (Malek 2020)  -> expect ~0
                 Barhee, Medjool = reducing-sugar types                                      -> expect >=0.5
  * test       : Rhars (= Ghars, Algeria) from the same 2015 sequencing batch (PRJNA296800).
  * RESULT OF THE FIRST RUN (2026-10-07): the pre-registered positive control SRR2577995 (NYUAD
    "Deglet noor", Algeria) did NOT show the deletion (DEL1 0.88, like the negative controls).
    Two independent Deglet Noor genomes were then added as positive controls (post hoc, declared):
      SRR6439416 = dnPdF (Weill Cornell Qatar, PRJNA427409; the sample family Malek 2020 used),
      SRR121605  = Deglet Noor female PI 4611 (Weill Cornell Qatar 2011, PRJNA40349; 84-bp SE).
    The run is valid only if these call del/del and the negative controls do not.
  * call rule  : mean normalised depth in the core window < 0.15 -> homozygous deletion (del/del);
                 0.30-0.70 -> heterozygous; > 0.80 -> no deletion; otherwise ambiguous.
  * normalisation: median depth of 1-kb bins on chr14 outside 2.0-3.3 Mb (MAPQ >= 20 only).

Usage: 02_invertase_cnv.py WORKDIR   (expects WORKDIR/paf/<sample>.chr14.paf from 02_map.sh)
"""
import sys, os, statistics, csv, json

CHR = "NC_052405.1"
BIN = 1000
CHR_LEN = None
SAMPLES = [  # run, label, published sugar type (source), role
    ("SRR2511369", "Rhars (Ghars)", "reported reducing-sugar-rich; fresh-fruit sucrose not cleanly measured (mimouni2014, mimouni2026 vs derkaoui2021)", "test"),
    ("SRR6439416", "Deglet Noor dnPdF (WCMQ)", "sucrose-type, homozygous deletion (malek2020)", "positive"),
    ("SRR121605", "Deglet Noor PI 4611 (WCMQ 2011)", "sucrose-type, homozygous deletion (malek2020)", "positive"),
    ("SRR2577995", "'Deglet noor' NYUAD 2015", "label says sucrose-type; failed as positive control", "label_check"),
    ("SRR2559387", "Medjool", "reducing-type, high invert sugars (terada2026)", "negative"),
    ("SRR2577516", "Barhee", "reducing-type (malek2020)", "negative"),
]
WINDOWS = {
    "DEL1_core 2.470-2.505 Mb": (2_470_000, 2_505_000),
    "CWINV1 gene (LOC103698975)": (2_466_060, 2_469_686),
    "CWINV3 gene (LOC103713368)": (2_508_976, 2_513_070),
    "DEL2 3.088-3.093 Mb": (3_088_000, 3_093_000),
    "A/N-INV1 gene (LOC103706133)": (3_079_319, 3_087_117),
    "flank_left 2.30-2.40 Mb": (2_300_000, 2_400_000),
    "flank_right 2.60-2.70 Mb": (2_600_000, 2_700_000),
}

def depth_bins(paf, minq=20):
    bins = {}
    n = 0
    with open(paf) as f:
        for line in f:
            p = line.split("\t")
            if len(p) < 12:  # truncated line
                continue
            if p[5] != CHR or int(p[11]) < minq or "tp:A:P" not in line:
                continue
            s, e = int(p[7]), int(p[8])
            n += 1
            b = s // BIN
            while b * BIN < e:
                lo, hi = max(s, b * BIN), min(e, (b + 1) * BIN)
                bins[b] = bins.get(b, 0) + (hi - lo)
                b += 1
    return {b: v / BIN for b, v in bins.items()}, n

def call(x):
    if x < 0.15: return "del/del"
    if 0.30 <= x <= 0.70: return "del/+"
    if x > 0.80: return "+/+"
    return "ambiguous"

def main(work):
    out_rows, profiles = [], {}
    for run, label, sugar, role in SAMPLES:
        paf = os.path.join(work, "paf", f"{run}.chr14.paf")
        if not os.path.exists(paf):
            print("missing", paf); continue
        d, n = depth_bins(paf)
        norm_bins = [v for b, v in d.items() if not (2_000 <= b < 3_300) and v > 0]
        med = statistics.median(norm_bins)
        maxb = max(d)
        prof = [d.get(b, 0.0) / med for b in range(maxb + 1)]
        profiles[run] = prof
        row = {"run": run, "cultivar": label, "role": role, "published_sugar_type": sugar,
               "chr14_alignments_q20": n, "median_depth_x": round(med, 2)}
        for w, (s, e) in WINDOWS.items():
            vals = prof[s // BIN: e // BIN]
            row[w] = round(sum(vals) / len(vals), 3)
        # descriptive uncertainty (added after the first run): mean +/- 1.96*SE over the 35 1-kb
        # bins of the DEL1 window. Bins are not fully independent, so read this as approximate.
        s1, e1 = WINDOWS["DEL1_core 2.470-2.505 Mb"]
        v1 = prof[s1 // BIN: e1 // BIN]
        se = statistics.stdev(v1) / len(v1) ** 0.5
        row["DEL1_ci95"] = f"{max(0, sum(v1)/len(v1) - 1.96*se):.2f}-{sum(v1)/len(v1) + 1.96*se:.2f}"
        row["call_DEL1"] = call(row["DEL1_core 2.470-2.505 Mb"])
        row["call_DEL2"] = call(row["DEL2 3.088-3.093 Mb"])
        # diagnostic: same windows with MAPQ >= 0. If depth comes back at MAPQ 0, the window is a
        # repeat / duplicated sequence (reads present but ambiguous), not a deletion.
        d0, _ = depth_bins(paf, minq=0)
        med0 = statistics.median([v for b, v in d0.items() if not (2_000 <= b < 3_300) and v > 0])
        for w in ("DEL1_core 2.470-2.505 Mb", "DEL2 3.088-3.093 Mb", "CWINV1 gene (LOC103698975)"):
            s0, e0 = WINDOWS[w]
            vals = [d0.get(b, 0.0) / med0 for b in range(s0 // BIN, e0 // BIN)]
            row[w.split()[0] + "_anyMAPQ"] = round(sum(vals) / len(vals), 3)
        out_rows.append(row)
        print(json.dumps(row, ensure_ascii=False))
    # validity gate (decided before reading the test sample)
    pos = [r for r in out_rows if r["role"] == "positive"]
    neg = [r for r in out_rows if r["role"] == "negative"]
    ok = bool(pos) and all(r["call_DEL1"] == "del/del" for r in pos) and all(r["call_DEL1"] != "del/del" for r in neg)
    print("VALIDITY:", "PASS" if ok else "FAIL",
          {r["cultivar"]: r["call_DEL1"] for r in pos + neg})
    os.makedirs("data", exist_ok=True)
    with open("data/02_invertase_genotypes.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(out_rows)
    # 1-kb normalised depth profile 2.30-3.15 Mb for plotting / re-use
    with open("data/02_invertase_depth_profile.tsv", "w") as f:
        runs = [r["run"] for r in out_rows]
        f.write("bin_start\t" + "\t".join(runs) + "\n")
        for b in range(2_300, 3_150):
            f.write(f"{b*BIN}\t" + "\t".join(f"{profiles[r][b]:.3f}" for r in runs) + "\n")
    return out_rows

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "work")
