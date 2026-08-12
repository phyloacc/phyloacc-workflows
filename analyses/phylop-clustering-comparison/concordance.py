#!/usr/bin/env python3
"""
concordance.py - score the clustering methods against a phastCons reference and tune each
method's parameters to best reproduce phastCons on the same 2 Mb chr1 slice.

Reference: phastCons --most-conserved elements, run on the real 241-species alignment of
chr1:26-28Mb (data/phastcons_reference.0based.bed, coordinates shifted to the 0-based
mini-chromosome to match data/chr1_26-28Mb.conserved-sites.bed).

Objective: bp-level F1 between a method's conserved regions and the phastCons elements
  precision = |method ∩ phastCons| / |method|
  recall    = |method ∩ phastCons| / |phastCons|
  F1        = 2·P·R / (P+R)

For each method we grid-search its knobs, keep the max-F1 setting, and report default vs
tuned. The HMM stays the simple 2-state model - we only pick better values for its existing
5 probabilities. Run in phyloacc-workflows (numpy, matplotlib, hdbscan).

Outputs: figures/concordance_*.png, concordance_results.csv
"""

import csv
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import lib.clustering as C
import lib.intervals as I

HERE = os.path.dirname(os.path.abspath(__file__))
LENGTH = 2_000_000
SLICE = os.path.join(HERE, "data", "chr1_26-28Mb.conserved-sites.bed")
REF = os.path.join(HERE, "data", "phastcons_reference.0based.bed")


def load_sites():
    out = []
    with open(SLICE) as f:
        for ln in f:
            _, s, e = ln.split("\t")
            out.append((int(s), int(e)))
    return out


def load_regions(path):
    """Tolerant BED loader: first 3 cols -> (start, end), clipped to [0, LENGTH)."""
    regs = []
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith(("#", "track")):
                continue
            p = ln.split()
            s, e = int(p[1]), int(p[2])
            s, e = max(0, s), min(LENGTH, e)
            if e > s:
                regs.append((s, e))
    return regs


def mask(regions):
    m = np.zeros(LENGTH, dtype=bool)
    for r in regions:
        m[r[0]:r[1]] = True
    return m


def prf(regions, ref_mask):
    """bp-level precision/recall/F1 of `regions` against the reference mask, plus the raw
    true-positive / false-positive / false-negative base-pair counts.
      TP = called bp that ARE in a phastCons element
      FP = called bp that are NOT (over-calls / false positives)
      FN = phastCons element bp the method MISSED (false negatives)
    """
    ref = int(ref_mask.sum())
    if not regions:
        return dict(precision=0.0, recall=0.0, f1=0.0, coverage_pct=0.0, n_regions=0,
                    tp_bp=0, fp_bp=0, fn_bp=ref)
    m = mask(regions)
    inter = int(np.logical_and(m, ref_mask).sum())
    called = int(m.sum())
    p = inter / called if called else 0.0
    r = inter / ref if ref else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return dict(precision=p, recall=r, f1=f1, coverage_pct=100.0 * called / LENGTH,
                n_regions=len(regions), tp_bp=inter, fp_bp=called - inter, fn_bp=ref - inter)


# method -> (default callable, param grid as list of (label, callable)); each callable: sites -> regions
def make_methods(sites):
    def gap(a, b, c):
        return [(s, e) for (_, s, e, _, _) in I.cluster_sites([("s", x, y) for x, y in sites], a, b, c)[0]]

    windowed_grid = [
        (f"w={w},min={m}", (lambda s, w=w, m=m: [(a, b) for a, b, _ in C.windowed(s, LENGTH, w, m)]))
        for w in (20, 50, 100, 200) for m in (2, 3, 5, 10)
    ]
    hmm_grid = [
        (f"t11={t11},e11={e11},min={ml}",
         (lambda s, t11=t11, e11=e11, ml=ml:
          [(a, b) for a, b, _ in C.hmm(s, LENGTH, t0_0=0.9, t1_1=t11, e0_0=0.8, e1_1=e11, s0=0.9,
                                       min_len=ml, max_len=100000)]))
        for t11 in (0.8, 0.9, 0.95, 0.99) for e11 in (0.3, 0.5, 0.7) for ml in (1, 20)
    ]
    # hdbscan given a fair, thorough tune: min_cluster_size x min_samples x selection method
    # ('leaf' = finer/smaller clusters, less gap-bridging, than default 'eom').
    hdb_grid = [
        (f"mcs={mcs},ms={ms},{csm}",
         (lambda s, mcs=mcs, ms=ms, csm=csm: [(a, b) for a, b, _ in C.hdbscan_cluster(
             s, min_cluster_size=mcs, min_samples=ms, cluster_selection_method=csm)]))
        for mcs in (5, 15, 50) for ms in (None, 10) for csm in ("eom", "leaf")
    ]
    return {
        "windowed": dict(default=lambda s: [(a, b) for a, b, _ in C.windowed(s, LENGTH, 50, 5)], grid=windowed_grid),
        "hmm": dict(default=lambda s: [(a, b) for a, b, _ in C.hmm(s, LENGTH)], grid=hmm_grid),
        "hdbscan": dict(default=lambda s: [(a, b) for a, b, _ in C.hdbscan_cluster(s, min_cluster_size=10)], grid=hdb_grid),
        "gap-merge (current)": dict(default=lambda s: gap(20, 5, 20),
                                    grid=[(f"gap={g},min={mn},len={ln}", (lambda s, g=g, mn=mn, ln=ln: gap(g, mn, ln)))
                                          for g in (10, 20, 50) for mn in (3, 5, 10) for ln in (1, 20)]),
    }


def main():
    sites = load_sites()
    ref = load_regions(REF)
    ref_mask = mask(ref)
    print(f"phastCons reference: {len(ref)} elements, {100*ref_mask.sum()/LENGTH:.1f}% of slice\n")

    methods = make_methods(sites)
    rows = []
    tuned_regions = {}
    for name, spec in methods.items():
        dflt = prf(spec["default"](sites), ref_mask)
        best_label, best, best_regs = None, dict(f1=-1), []
        for label, fn in spec["grid"]:
            regs = fn(sites)
            sc = prf(regs, ref_mask)
            if sc["f1"] > best["f1"]:
                best, best_label, best_regs = sc, label, regs
        tuned_regions[name] = best_regs
        med = int(np.median([e - s for s, e in best_regs])) if best_regs else 0
        rows.append(dict(method=name,
                         default_f1=round(dflt["f1"], 3), default_cov=round(dflt["coverage_pct"], 1),
                         tuned_f1=round(best["f1"], 3), tuned_precision=round(best["precision"], 3),
                         tuned_recall=round(best["recall"], 3), tuned_cov=round(best["coverage_pct"], 1),
                         tuned_fp_kb=round(best["fp_bp"] / 1000), tuned_fn_kb=round(best["fn_bp"] / 1000),
                         tuned_n_regions=len(best_regs), tuned_median_bp=med, tuned_params=best_label))
        print(f"{name:<20} default F1={dflt['f1']:.3f} (cov {dflt['coverage_pct']:.0f}%)  ->  "
              f"tuned F1={best['f1']:.3f} (P={best['precision']:.2f} R={best['recall']:.2f} "
              f"FP={best['fp_bp']//1000}kb FN={best['fn_bp']//1000}kb "
              f"cov {best['coverage_pct']:.0f}%, {len(best_regs)} regions, median {med}bp)  [{best_label}]")

    with open(os.path.join(HERE, "concordance_results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    names = [r["method"] for r in rows]
    colors = {"windowed": "#4C72B0", "hmm": "#C44E52", "hdbscan": "#55A868", "gap-merge (current)": "#888888"}

    # Figure: default vs tuned F1
    fig, ax = plt.subplots(figsize=(8.5, 5))
    x = np.arange(len(names))
    ax.bar(x - 0.2, [r["default_f1"] for r in rows], 0.4, label="default params", color="#cccccc")
    ax.bar(x + 0.2, [r["tuned_f1"] for r in rows], 0.4, label="tuned (max F1)",
           color=[colors[n] for n in names])
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=25, ha="right")
    ax.set_ylabel("F1 vs phastCons elements"); ax.set_ylim(0, 1)
    ax.set_title("Concordance with phastCons: default vs tuned parameters")
    ax.legend(); ax.grid(True, axis="y", alpha=0.15)
    for xi, r in zip(x, rows):
        ax.text(xi + 0.2, r["tuned_f1"] + 0.02, f"{r['tuned_f1']:.2f}", ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "figures", "concordance_f1.png"), dpi=150)
    plt.close(fig)

    # ---- tuned descriptive figures, with phastCons as a reference series ----
    PC = "#B8860B"
    series = [("phastCons (reference)", ref, PC)] + [(n, tuned_regions[n], colors[n]) for n in names]
    starts = np.array([s for s, _ in sites])

    # Track view: 60 kb window - raw sites, then phastCons + each tuned method
    w_lo, w_hi = 1_000_000, 1_060_000
    fig, axes = plt.subplots(len(series) + 1, 1, figsize=(11, 6.2), sharex=True,
                             gridspec_kw={"height_ratios": [1.5] + [1] * len(series)})
    axes[0].hist(starts[(starts >= w_lo) & (starts < w_hi)],
                 bins=np.arange(w_lo, w_hi + 1000, 1000), color="#555555")
    axes[0].set_ylabel("sites\n/1 kb", fontsize=8)
    axes[0].set_title(f"Tuned methods vs phastCons — {w_lo:,}-{w_hi:,} (0-based). Raw sites, then each element set.")
    for ax, (label, regs, color) in zip(axes[1:], series):
        for s, e, *_ in regs:
            if e > w_lo and s < w_hi:
                ax.axvspan(max(s, w_lo), min(e, w_hi), color=color, alpha=0.85)
        ax.set_ylim(0, 1); ax.set_yticks([])
        ax.set_ylabel(label.replace(" (reference)", "\n(ref)").replace(" (current)", ""),
                      fontsize=8, rotation=0, ha="right", va="center")
    axes[-1].set_xlabel("position (bp, 0-based within slice)")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "tuned_track.png"), dpi=150); plt.close(fig)

    # Region-size distributions (tuned + phastCons)
    fig, ax = plt.subplots(figsize=(8.5, 5))
    bins = np.logspace(0, 4, 40)
    for label, regs, color in series:
        if regs:
            ax.hist([e - s for s, e, *_ in regs], bins=bins, histtype="step", lw=2, color=color,
                    label=label, ls="--" if "phastCons" in label else "-")
    ax.set_xscale("log"); ax.set_xlabel("element length (bp)"); ax.set_ylabel("number of elements")
    ax.set_title("Element-length distributions: tuned methods vs phastCons")
    ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.15)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "tuned_sizes.png"), dpi=150); plt.close(fig)

    # Jaccard (bp) among tuned methods + phastCons
    labs = [s[0].replace(" (reference)", "").replace(" (current)", "") for s in series]
    masks = [mask(s[1]) for s in series]
    k = len(series); J = np.ones((k, k))
    for a in range(k):
        for b in range(k):
            inter = np.logical_and(masks[a], masks[b]).sum(); union = np.logical_or(masks[a], masks[b]).sum()
            J[a, b] = inter / union if union else 0.0
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    im = ax.imshow(J, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(k)); ax.set_xticklabels(labs, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(k)); ax.set_yticklabels(labs, fontsize=8)
    for a in range(k):
        for b in range(k):
            ax.text(b, a, f"{J[a,b]:.2f}", ha="center", va="center",
                    color="white" if J[a, b] < 0.6 else "black", fontsize=8)
    ax.set_title("Agreement among tuned methods + phastCons (Jaccard of covered bp)")
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "tuned_jaccard.png"), dpi=150); plt.close(fig)

    print("\nwrote concordance_results.csv + figures/{concordance_f1,tuned_track,tuned_sizes,tuned_jaccard}.png")


if __name__ == "__main__":
    main()
