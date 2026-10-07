#!/usr/bin/env python3
"""Plot normalised read depth across the chr14 sucrose locus (output of 02_invertase_cnv.py).

Left column : zoom on DEL1 (2.44-2.53 Mb), raw 1-kb bins.
Right column: whole region 2.30-3.15 Mb (DEL1 + DEL2), 5-kb running mean.
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = list(csv.DictReader(open("data/02_invertase_depth_profile.tsv"), delimiter="\t"))
geno = {r["run"]: r for r in csv.DictReader(open("data/02_invertase_genotypes.tsv"), delimiter="\t")}
runs = [k for k in rows[0] if k != "bin_start"]
x = [int(r["bin_start"]) / 1e6 for r in rows]
ROLE = {"test": "TEST", "positive": "positive control", "negative": "negative control",
        "label_check": "label check"}
COLOR = {"test": "#d62728", "positive": "#2ca02c", "negative": "#4c78a8", "label_check": "#9467bd"}
genes = [("CWINV1", 2.46606, 2.46969), ("CWINV3", 2.50898, 2.51307), ("A/N-INV1", 3.07932, 3.08712)]


def smooth(v, k=5):
    out = []
    for i in range(len(v)):
        w = v[max(0, i - k // 2): i + k // 2 + 1]
        out.append(sum(w) / len(w))
    return out


fig, axes = plt.subplots(len(runs), 2, figsize=(13, 1.9 * len(runs)), sharex="col",
                         gridspec_kw={"width_ratios": [1, 1.6]})
for (axz, axw), run in zip(axes, runs):
    g = geno[run]
    col = COLOR.get(g["role"], "#4c78a8")
    y = [float(r[run]) for r in rows]
    for ax, yy in ((axz, y), (axw, smooth(y))):
        ax.fill_between(x, 0, yy, step="mid", color=col, alpha=0.85, linewidth=0)
        ax.axhline(1.0, color="grey", lw=0.6, ls="--")
        ax.axhline(0.5, color="grey", lw=0.4, ls=":")
        ax.axvspan(2.470, 2.505, color="#e45756", alpha=0.13)   # DEL1 decision window
        ax.axvspan(3.088, 3.093, color="#e45756", alpha=0.25)   # DEL2
        for _, s, e in genes:
            ax.axvspan(s, e, color="#f58518", alpha=0.35)
        ax.set_ylim(0, 2.2)
    axz.set_xlim(2.44, 2.53)
    axw.set_xlim(2.30, 3.15)
    axz.set_ylabel("depth / median", fontsize=7)
    axz.set_title(f"{g['cultivar']} — {ROLE.get(g['role'], g['role'])}", fontsize=8.5, loc="left")
    axw.set_title(f"{run} · DEL1 window {float(g['DEL1_core 2.470-2.505 Mb']):.2f} → {g['call_DEL1']}"
                  f"   ·   DEL2 {g['call_DEL2']} (unreliable: low mappability)", fontsize=8, loc="left")
    for ax in (axz, axw):
        ax.tick_params(labelsize=7)
axes[-1][0].set_xlabel("chr14 Mb (zoom DEL1)", fontsize=8)
axes[-1][1].set_xlabel("Barhee BC4 chr14 (NC_052405.1), Mb · light red = DEL1 decision window 2.470–2.505 & DEL2 "
                       "(Hazzouri 2019) · orange = invertase genes", fontsize=7.5)
plt.tight_layout()
plt.savefig("data/02_invertase_depth.png", dpi=130)
print("wrote data/02_invertase_depth.png")
