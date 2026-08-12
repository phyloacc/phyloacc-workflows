#!/usr/bin/env python3
"""
cluster_compare.py - run the conserved-region clustering methods on a real phyloP slice
and compare their results.

Input: a 2 Mb slice of the echolocation-v3 (Zoonomia 241) chr1 FDR-significant conserved
sites (data/chr1_26-28Mb.conserved-sites.bed; 0-based, shifted to a mini-chromosome).
Methods (lib/clustering.py, plus the current gap-merge baseline in lib/intervals.py):
  windowed  - fixed 50 bp windows, >= 5 sites/window, merge adjacent
  hmm       - 2-state online HMM (legacy defaults)
  hdbscan   - 1-D density clustering, min_cluster_size=10
  gap-merge - current naive method (max_gap 20, min 5 sites, min_len 20)  [baseline]

Outputs: figures/*.png and cluster_metrics.csv. Run in the phyloacc-workflows env
(needs numpy, matplotlib, hdbscan).
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import lib.clustering as C
import lib.intervals as I

HERE = os.path.dirname(os.path.abspath(__file__))
SLICE = os.path.join(HERE, "data", "chr1_26-28Mb.conserved-sites.bed")
LENGTH = 2_000_000
REGION_LABEL = "chr1:26,000,000-28,000,000 (human GRCh38; Zoonomia 241 phyloP), 0-based"

# method -> (callable(sites, length) -> [(start,end,n_sites)], color, params-string)
METHODS = [
    ("windowed", lambda s: C.windowed(s, LENGTH, window_bp=50, min_sites_per_window=5),
     "#4C72B0", "50 bp windows, >=5 sites/window"),
    ("hmm", lambda s: C.hmm(s, LENGTH),
     "#C44E52", "2-state HMM (t00=.9,t11=.8,e00=.8,e11=.5,s0=.9; 20<len<10000)"),
    ("hdbscan", lambda s: C.hdbscan_cluster(s, min_cluster_size=10),
     "#55A868", "min_cluster_size=10"),
    ("gap-merge (current)", lambda s: [(st, en, n) for (_, st, en, n, _)
                                        in I.cluster_sites([("s", a, b) for a, b in s], 20, 5, 20)[0]],
     "#888888", "max_gap 20, >=5 sites, min_len 20  [baseline]"),
]


def load_sites():
    sites = []
    with open(SLICE) as f:
        for line in f:
            _, s, e = line.rstrip("\n").split("\t")
            sites.append((int(s), int(e)))
    return sites


def coverage_mask(regions):
    m = np.zeros(LENGTH, dtype=bool)
    for s, e, _ in regions:
        m[s:e] = True
    return m


def metrics(regions, n_total_sites):
    if not regions:
        return dict(n_regions=0, coverage_bp=0, coverage_pct=0.0, median_len=0,
                    mean_len=0.0, max_len=0, sites_in_regions=0, pct_sites=0.0)
    lens = np.array([e - s for s, e, _ in regions])
    sites_in = sum(n for _, _, n in regions)
    cov = int(coverage_mask(regions).sum())
    return dict(
        n_regions=len(regions), coverage_bp=cov, coverage_pct=100.0 * cov / LENGTH,
        median_len=int(np.median(lens)), mean_len=float(lens.mean()), max_len=int(lens.max()),
        sites_in_regions=int(sites_in), pct_sites=100.0 * sites_in / n_total_sites,
    )


def main():
    sites = load_sites()
    starts = np.array([s for s, _ in sites])
    n_total = len(sites)
    print(f"{n_total} sites over {LENGTH} bp ({REGION_LABEL})\n")

    results, masks, rows = {}, {}, []
    for name, fn, color, params in METHODS:
        regs = fn(sites)
        results[name] = (regs, color, params)
        masks[name] = coverage_mask(regs)
        m = metrics(regs, n_total)
        rows.append({"method": name, "params": params, **m})
        print(f"{name:<20} {m['n_regions']:>6} regions | cov {m['coverage_pct']:5.1f}% | "
              f"median {m['median_len']:>4} bp | sites captured {m['pct_sites']:5.1f}%")

    # ---- metrics CSV ----
    with open(os.path.join(HERE, "cluster_metrics.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    names = [n for n, *_ in METHODS]
    colors = {n: c for n, _, c, _ in METHODS}

    # ---- Figure 1: region-length distributions ----
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    bins = np.logspace(0, 4, 40)
    for name in names:
        regs = results[name][0]
        if regs:
            lens = [e - s for s, e, _ in regs]
            ax.hist(lens, bins=bins, histtype="step", lw=2, color=colors[name], label=name)
    ax.set_xscale("log")
    ax.set_xlabel("conserved-region length (bp)")
    ax.set_ylabel("number of regions")
    ax.set_title("Region-length distributions by clustering method")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.15)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "figures", "region_length_dist.png"), dpi=150)
    plt.close(fig)

    # ---- Figure 2: track view of a 60 kb sub-window ----
    w_lo, w_hi = 1_000_000, 1_060_000
    fig, axes = plt.subplots(len(names) + 1, 1, figsize=(11, 5.6), sharex=True,
                             gridspec_kw={"height_ratios": [1.6] + [1] * len(names)})
    sub = starts[(starts >= w_lo) & (starts < w_hi)]
    axes[0].hist(sub, bins=np.arange(w_lo, w_hi + 1000, 1000), color="#555555")
    axes[0].set_ylabel("sites\n/1 kb", fontsize=8)
    axes[0].set_title(f"Zoomed track view: {w_lo:,}-{w_hi:,} (0-based) — raw significant sites, then each method's regions")
    for k, name in enumerate(names, start=1):
        ax = axes[k]
        for s, e, _ in results[name][0]:
            if e > w_lo and s < w_hi:
                ax.axvspan(max(s, w_lo), min(e, w_hi), color=colors[name], alpha=0.85)
        ax.set_ylim(0, 1); ax.set_yticks([])
        ax.set_ylabel(name, fontsize=8, rotation=0, ha="right", va="center")
    axes[-1].set_xlabel("position (bp, 0-based within slice)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "figures", "track_view.png"), dpi=150)
    plt.close(fig)

    # ---- Figure 3: pairwise Jaccard (bp overlap) ----
    n = len(names)
    J = np.ones((n, n))
    for a in range(n):
        for b in range(n):
            ma, mb = masks[names[a]], masks[names[b]]
            inter = np.logical_and(ma, mb).sum()
            union = np.logical_or(ma, mb).sum()
            J[a, b] = inter / union if union else 0.0
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    im = ax.imshow(J, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(names, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(n)); ax.set_yticklabels(names, fontsize=8)
    for a in range(n):
        for b in range(n):
            ax.text(b, a, f"{J[a,b]:.2f}", ha="center", va="center",
                    color="white" if J[a, b] < 0.6 else "black", fontsize=8)
    ax.set_title("Pairwise agreement (Jaccard of covered bp)")
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "figures", "jaccard_heatmap.png"), dpi=150)
    plt.close(fig)

    # ---- Figure 4: coverage % and region count ----
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.2))
    cov = [next(r for r in rows if r["method"] == n)["coverage_pct"] for n in names]
    nreg = [next(r for r in rows if r["method"] == n)["n_regions"] for n in names]
    a1.bar(range(n), cov, color=[colors[x] for x in names])
    a1.set_xticks(range(n)); a1.set_xticklabels(names, rotation=35, ha="right", fontsize=8)
    a1.set_ylabel("% of slice covered"); a1.set_title("Coverage")
    a2.bar(range(n), nreg, color=[colors[x] for x in names])
    a2.set_xticks(range(n)); a2.set_xticklabels(names, rotation=35, ha="right", fontsize=8)
    a2.set_ylabel("number of regions"); a2.set_title("Region count")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "figures", "summary_bars.png"), dpi=150)
    plt.close(fig)

    print(f"\nwrote figures/ and cluster_metrics.csv")


if __name__ == "__main__":
    main()
