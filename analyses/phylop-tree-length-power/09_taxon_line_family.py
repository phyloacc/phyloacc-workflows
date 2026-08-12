#!/usr/bin/env python3
"""
09_taxon_line_family.py - detectability boundary for EVERY taxon count n = 1..20.

Same axes as 08 (tree depth T x sites tested M), but no shading: just one dashed
detectability boundary per taxon count n, colored by n. Each line is the
detectable/undetectable boundary for that n (detectable to its right / below its wall):

    detectable iff  ceiling(n, T) >= log10(M / alpha),   kappa = 1 (JC).

Reading it: lines march rightward as n falls (fewer taxa -> deeper tree needed), and each
one flattens into its own horizontal taxon wall at M_wall(n) = alpha * 10^(2(n-1)ln4 cap)
above which no depth detects. Small n never clear the bar within the plot; the walls fan
upward as n grows.

Output: figures/taxon_line_family.png
"""
import os
import math
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

import phylop_power as pp

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = 0.05
XMAX = 40.0
NMAX = 20


def ceiling_jc_star(n, T):
    x = math.exp(-4.0 * T / (3.0 * n))
    p_same, p_diff = 0.25 + 0.75 * x, 0.25 - 0.25 * x
    return pp.neglog10_from_lrt(-2.0 * math.log(p_same ** n + 3.0 * p_diff ** n))


def boundary_T(M_grid, n, T_grid):
    ceil = np.array([ceiling_jc_star(n, t) for t in T_grid])
    return np.interp(np.log10(M_grid / ALPHA), ceil, T_grid, left=np.nan, right=np.nan)


def main():
    T_grid = np.geomspace(0.1, 600, 500)   # wide so each curve reaches its wall
    M_grid = np.geomspace(1e1, 1e9, 400)

    fig, ax = plt.subplots(figsize=(9.6, 6.7))
    norm = Normalize(vmin=1, vmax=NMAX)
    cmap = plt.cm.viridis

    for n in range(1, NMAX + 1):
        Tb = boundary_T(M_grid, n, T_grid)
        if np.all(np.isnan(Tb)):
            continue  # n so small it never clears the bar anywhere in range (n=1, etc.)
        ax.plot(Tb, M_grid, ls="--", lw=1.3, color=cmap(norm(n)), zorder=3)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("total neutral tree depth  T  (sum of branch lengths, subs/site)")
    ax.set_ylabel("number of sites tested  M  (multiple-testing burden)")
    ax.set_title("Detectability boundary for every taxon count n = 1..20  (JC, κ=1, α=0.05)\n"
                 "each dashed line: detectable to its right; flattens into its taxon wall above")
    ax.set_xlim(0.1, XMAX)
    ax.set_ylim(1e1, 1e9)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.015)
    cbar.set_label("number of taxa  n")
    cbar.set_ticks(range(2, NMAX + 1, 2))

    fig.tight_layout()
    out = os.path.join(HERE, "figures", "taxon_line_family.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
    print("\ntaxon wall M_wall(n) = alpha * 10^cap(n):")
    for n in range(1, NMAX + 1):
        cap = pp.neglog10_from_lrt(2.0 * (n - 1) * math.log(4.0))
        print(f"  n={n:2d}:  cap ≈ {cap:6.2f}   M_wall ≈ {ALPHA*10**cap:.1e}")


if __name__ == "__main__":
    main()
