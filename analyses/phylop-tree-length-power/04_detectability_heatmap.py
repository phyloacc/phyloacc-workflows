#!/usr/bin/env python3
"""
04_detectability_heatmap.py - the two-variable (tree length x number of sites) view.

Whether ANY conserved site can survive genome-wide BH-FDR is a race between two
quantities:

  * the CEILING on the conserved score, set by total tree length T:
        ceiling(T) = -log10( 2*Phi(-sqrt(LRT_max(T))) ),  LRT_max ~ kappa*T
  * the BAR a site must clear, set by the number of tested sites N (and alpha):
        bar(N) = -log10(alpha / N)          (BH rank-1 / genome-wide threshold)

A conserved site is detectable only where ceiling(T) >= bar(N). This script maps the
"detectability margin" ceiling(T) - bar(N) (in orders of magnitude of p) over the
(T, N) plane, draws the zero contour (the detectability boundary), and marks where the
four real dataset families fall.

ceiling(T) is computed from a real fitted model (the 241-taxon mammal neutral model,
uniformly rescaled) via phylop_power - no simulated data. kappa is mildly
model-dependent (~2-3), which shifts the boundary by a fraction of a tree-length unit;
noted on the figure.

Outputs: figures/detectability_heatmap.png
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import phylop_power as pp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
MODEL = os.path.join(REPO, "data/mammals-v3/echolocation-v3/02-neutral-model/phylofit/autosomes/chr1-corrected.mod")

ALPHA = 0.05

# real families: (label, tree_length, representative per-chromosome N, color)
# turtle N is the measured site count of a real scaffold (6.32e6); others use a
# representative large-chromosome value - the boundary is steep in T, so the exact
# N (within 1e6-1e8) does not change which side of it a family lands on.
FAMILIES = [
    ("hamster", 0.93, 9.5e7, "#4C72B0"),
    ("birds", 1.09, 2.0e7, "#55A868"),
    ("turtle", 2.81, 6.32e6, "#C44E52"),
    ("mammal", 16.07, 1.9e8, "#8172B3"),
]


def ceiling_curve(T_grid):
    md = pp.parse_mod(MODEL)
    Tbase = pp.tree_length(pp.parse_newick(md["tree_newick"]))
    return np.array([pp.ceiling(md, rho_scale=t / Tbase)["neglog10p"] for t in T_grid])


def main():
    T = np.geomspace(0.3, 30, 160)
    N = np.geomspace(1e1, 1e9, 160)
    ceil = ceiling_curve(T)                        # shape (len(T),)
    bar = -np.log10(ALPHA / N)                     # shape (len(N),)

    # margin[i, j] = ceiling(T_j) - bar(N_i)  (orders of magnitude of p)
    margin = ceil[None, :] - bar[:, None]

    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    vlim = 12
    pcm = ax.pcolormesh(T, N, margin, cmap="RdBu", vmin=-vlim, vmax=vlim, shading="auto")
    cbar = fig.colorbar(pcm, ax=ax, extend="both")
    cbar.set_label("detectability margin  ceiling(T) − bar(N)\n(orders of magnitude of p; >0 = detectable)")

    # detectability boundary, drawn explicitly (ceiling is monotonic in T, so invert
    # T = ceiling^{-1}(bar(N)) by interpolation) - a complete, clean curve.
    T_boundary = np.interp(bar, ceil, T, left=np.nan, right=np.nan)
    ax.plot(T_boundary, N, color="black", lw=2.4, zorder=5)
    ax.text(6.6, 1.2e5, "detectability\nboundary", fontsize=8.5, rotation=68,
            color="black", ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("total neutral tree length  T  (expected substitutions / site)")
    ax.set_ylabel("number of sites tested  N  (genome-wide FDR)")
    ax.set_title("Can any conserved site survive genome-wide FDR?\n"
                 "A race between the tree-length ceiling and the multiple-testing bar (α=0.05)")

    # region labels
    ax.text(0.62, 5e6, "NO conserved site\ncan pass\n(ceiling below bar)",
            fontsize=11, color="#7a1f1f", ha="center", va="center", weight="bold",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.65))
    ax.text(20, 8e3, "conserved sites\ndetectable",
            fontsize=11, color="#1f3a7a", ha="center", va="center", weight="bold",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.65))

    # real families
    for label, Tf, Nf, color in FAMILIES:
        ax.scatter([Tf], [Nf], s=90, color=color, edgecolor="black", linewidth=0.8, zorder=6)
        ax.annotate(f"{label}\nT={Tf}", (Tf, Nf), textcoords="offset points",
                    xytext=(8, 8), fontsize=8.5, color=color, weight="bold")

    ax.text(0.32, 2e1, "ceiling(T) from real mammal model rescaled; κ≈2–3 is mildly\n"
                       "model-dependent, shifting the boundary by <~1 tree-length unit",
            fontsize=6.8, color="gray")

    fig.tight_layout()
    out = os.path.join(HERE, "figures", "detectability_heatmap.png")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")

    # print the boundary N for a few T, and the family verdicts
    print("\ndetectability boundary (min detectable ceiling vs T):")
    for t in [1, 2, 5, 10, 16, 24]:
        c = float(ceiling_curve(np.array([t]))[0])
        Nmax = ALPHA * 10 ** c   # bar(N)=ceiling -> N = alpha*10^ceiling
        print(f"  T={t:>2}: ceiling={c:5.2f} -> largest N still detectable = {Nmax:.2e} sites")
    print("\nfamily verdicts (at their representative N):")
    for label, Tf, Nf, _ in FAMILIES:
        c = float(ceiling_curve(np.array([Tf]))[0])
        bar_f = -np.log10(ALPHA / Nf)
        verdict = "DETECTABLE" if c >= bar_f else "undetectable"
        print(f"  {label:<8} T={Tf:<5} N={Nf:.1e}: ceiling={c:5.2f}  bar={bar_f:5.2f}  -> {verdict} (margin {c-bar_f:+.1f})")


if __name__ == "__main__":
    main()
