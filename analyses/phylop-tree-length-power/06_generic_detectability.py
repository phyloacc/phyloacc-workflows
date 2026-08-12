#!/usr/bin/env python3
"""
06_generic_detectability.py - a GENERIC detectability nomogram (no dataset markers).

Companion to 04_detectability_heatmap.py, built to be reusable in a tutorial /
walkthrough: a reader locates their own (tree depth T, test burden M) and reads off
whether per-site phyloP can make a call. A conserved site is detectable only where the
tree-length CEILING clears the multiple-testing BAR:

    ceiling(T, model)  >=  log10(M / alpha)

The ceiling rises with total tree length T *and* with taxon sampling (more taxa at fixed
T = more independent "no substitution" evidence). It also depends on the substitution
model: a synthetic JC/star model badly *under*-estimates real ceilings (a 250-taxon
JC-star tree scores like a ~15-taxon real model), so this figure uses REAL fitted neutral
models, uniformly rescaled in T, for the band edges - sparse (~15-taxon) and dense
(~241-taxon) vertebrate models. The band between them is the honest taxon-sampling spread.

Approximation: rescaling a fitted model far from its own depth extrapolates its shape;
fine for an order-of-magnitude nomogram. Reads two .mod files from this repo's data/
(same convention as 04); to reuse elsewhere, point MODELS at your own neutral models.

Output: figures/generic_detectability.png
"""
import os
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import phylop_power as pp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
ALPHA = 0.05

# Band edges: real fitted neutral models of contrasting taxon count, rescaled in T.
# (Labeled by taxon count, not organism, to keep the figure generic.)
SPARSE_MOD = os.path.join(REPO, "data/hamsters/workflow-tests/test-full/"
                          "02-neutral-model/phylofit/autosomes/CM001006.3-corrected.mod")
DENSE_MOD = os.path.join(REPO, "data/mammals-v3/echolocation-v3/"
                         "02-neutral-model/phylofit/autosomes/chr1-corrected.mod")


def ceiling_vs_T(mod_path, T_grid):
    md = pp.parse_mod(mod_path)
    Tbase = pp.tree_length(pp.parse_newick(md["tree_newick"]))
    n_tips = pp.count_tips(pp.parse_newick(md["tree_newick"]))
    ceil = np.array([pp.ceiling(md, rho_scale=t / Tbase)["neglog10p"] for t in T_grid])
    return ceil, n_tips


def boundary_T(M_grid, ceil, T_grid):
    """Invert the monotone ceiling(T): minimum T whose ceiling clears bar(M)."""
    bar = np.log10(M_grid / ALPHA)
    return np.interp(bar, ceil, T_grid, left=np.nan, right=np.nan)


def main():
    T_grid = np.geomspace(0.1, 40, 90)
    M_grid = np.geomspace(1e1, 1e9, 200)
    ceil_sparse, n_sparse = ceiling_vs_T(SPARSE_MOD, T_grid)
    ceil_dense, n_dense = ceiling_vs_T(DENSE_MOD, T_grid)
    Tb_sparse = boundary_T(M_grid, ceil_sparse, T_grid)  # right edge (fewer taxa -> harder)
    Tb_dense = boundary_T(M_grid, ceil_dense, T_grid)     # left edge  (more taxa  -> easier)
    # For the zone fills, clamp a boundary that runs off the right edge (NaN) to T=40 so
    # the band/regions tile the whole plane; the boundary *lines* keep NaN (stop cleanly).
    Tb_sparse_f = np.where(np.isnan(Tb_sparse), 40.0, Tb_sparse)
    Tb_dense_f = np.where(np.isnan(Tb_dense), 40.0, Tb_dense)

    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    ax.fill_betweenx(M_grid, Tb_sparse_f, 40, color="#c7dcef", alpha=0.7, zorder=0)
    ax.fill_betweenx(M_grid, 0.1, Tb_dense_f, color="#f2c9c4", alpha=0.7, zorder=0)
    ax.fill_betweenx(M_grid, Tb_dense_f, Tb_sparse_f, color="#ece3b0", alpha=0.75, zorder=0)
    ax.plot(Tb_dense, M_grid, color="#8a6d1a", lw=1.6, zorder=3)
    ax.plot(Tb_sparse, M_grid, color="#8a6d1a", lw=1.6, zorder=3)

    ax.text(22, 3e2, "single-site calls\nDETECTABLE", fontsize=12, weight="bold",
            color="#1f3a7a", ha="center", va="center")
    ax.text(0.35, 5e6, "single-site power\nlost\n(no site can pass,\nhowever conserved)",
            fontsize=10.5, weight="bold", color="#7a1f1f", ha="center", va="center")
    ax.text(15.5, 1.2e7, "depends on\ntaxon sampling", fontsize=9, style="italic",
            color="#6b5310", ha="center", va="center", rotation=64)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("total neutral tree depth  T  (sum of branch lengths, subs/site)")
    ax.set_ylabel("number of sites tested  M  (multiple-testing burden)")
    ax.set_title("Where can per-site phyloP make a call?\n"
                 "detectable iff  ceiling(T, model) ≥ log10(M / α),  α=0.05")
    ax.set_xlim(0.1, 40)
    ax.set_ylim(1e1, 1e9)
    ax.legend(handles=[Patch(fc="#f2c9c4", label="power-limited (any taxon sampling)"),
                       Patch(fc="#ece3b0", label=f"boundary band: ~{n_sparse}→~{n_dense} taxa"),
                       Patch(fc="#c7dcef", label="detectable (any taxon sampling)")],
              loc="lower right", fontsize=8.5, framealpha=0.9)
    ax.text(0.11, 1.4e1,
            f"band edges = real fitted neutral models (~{n_sparse} and ~{n_dense} taxa) rescaled in T;\n"
            "real substitution structure sets the ceiling (a synthetic JC model underestimates it).",
            fontsize=6.6, color="gray")

    fig.tight_layout()
    out = os.path.join(HERE, "figures", "generic_detectability.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}  (sparse={n_sparse} taxa, dense={n_dense} taxa)")

    print("\nminimum tree depth T to enter the detectable zone:")
    for m in (1e2, 1e6, 1e7, 1e8, 1e9):
        i = int(np.argmin(np.abs(M_grid - m)))
        print(f"  M={m:.0e}:  dense (~{n_dense} taxa) T≈{Tb_dense[i]:5.1f}   sparse (~{n_sparse} taxa) T≈{Tb_sparse[i]:5.1f}")


if __name__ == "__main__":
    main()
