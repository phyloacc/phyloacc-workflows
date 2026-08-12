#!/usr/bin/env python3
"""
10_rate_effect.py - how the substitution rate max(-Q_bb) shifts the detectability map.

Companion to 08_taxon_walls.py (which fixes the rate at max(-Q_bb)=1). The depth-limited
part of the ceiling is  LRT ~= 2 * max(-Q_bb) * T,  so the rate enters ONLY as a horizontal
rescaling of the T axis: a faster model slides every depth-limited boundary LEFT (needs less
depth) by the factor 1/rate. The taxon wall  2(n-1)ln4  is rate-INDEPENDENT and does not move.

Shown for a single taxon count (n=12): the rising boundary fans left with higher rate, but
all curves meet the SAME horizontal wall - a faster model lets you reach the wall with less
depth, yet cannot raise it (rate can't fix too-few-species).

Output: figures/rate_effect.png
"""
import os
import math
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import phylop_power as pp  # neglog10_from_lrt: LRT -> -log10 p

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = 0.05
XMAX = 40.0
N = 12
# Okabe-Ito colorblind-safe trio, plus distinct linestyles so color isn't the only cue.
RATES = [(1.0, "#0072B2", "-"), (1.5, "#E69F00", "--"), (2.0, "#009E73", "-.")]


def ceiling(n, T, rate):
    lrt = min(2.0 * rate * T, 2.0 * (n - 1) * math.log(4.0))   # depth-limited vs taxon wall
    return pp.neglog10_from_lrt(lrt)


def boundary_T(M_grid, n, rate, T_grid):
    ceil = np.array([ceiling(n, t, rate) for t in T_grid])
    return np.interp(np.log10(M_grid / ALPHA), ceil, T_grid, left=np.nan, right=np.nan)


def main():
    T_grid = np.geomspace(0.1, 400, 800)
    M_grid = np.geomspace(1e1, 1e9, 300)
    wall = ALPHA * 10 ** ceiling(N, 1e9, 1.0)   # M where n=12 saturates (rate-independent)

    fig, ax = plt.subplots(figsize=(9.2, 6.4))

    # shade the detectable region for the slowest model (rate=1) as a faint reference
    Tb1 = boundary_T(M_grid, N, 1.0, T_grid)
    ax.fill_betweenx(M_grid, np.where(np.isnan(Tb1), XMAX, Tb1), XMAX,
                     color="#0072B2", alpha=0.06, zorder=0)

    for rate, col, ls in RATES:
        Tb = boundary_T(M_grid, N, rate, T_grid)
        ax.plot(Tb, M_grid, color=col, ls=ls, lw=2.1, zorder=4,
                label=f"max(−Q_bb) = {rate:g}")

    # the taxon wall (same for every rate)
    ax.axhline(wall, color="#444", lw=1.0, ls=":", zorder=3)
    ax.text(XMAX * 0.97, wall * 1.28,
            f"taxon wall for n={N}:  2(n−1)ln4   (M≈{wall:.0e})\nrate cannot raise this",
            color="#444", fontsize=7.6, ha="right", va="bottom")

    # arrow: higher rate -> boundary slides left
    ax.annotate("", xy=(0.30, 0.30), xytext=(0.52, 0.30), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#333", lw=1.6))
    ax.text(0.41, 0.335, "faster model →\nless depth needed", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=9, style="italic", color="#333")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.1, XMAX); ax.set_ylim(1e1, 1e9)
    ax.set_xlabel("total neutral tree depth  T  (sum of branch lengths, subs/site)")
    ax.set_ylabel("number of sites tested  M  (multiple-testing burden)")
    ax.set_title(f"What the substitution rate does to the boundary   (held at n={N} species)\n"
                 "a faster model slides the depth-limited boundary LEFT (needs less depth); the taxon wall does not move\n"
                 "LRT ≈ min(2·T·max(−Q_bb), 2(n−1)ln4);  equal base frequencies;  α=0.05")
    ax.legend(loc="lower right", fontsize=9, title="max exit rate  (=1 is the equal-rate floor)",
              title_fontsize=8.5, framealpha=0.95)

    fig.tight_layout()
    out = os.path.join(HERE, "figures", "rate_effect.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
    for rate, _, _ in RATES:
        # depth needed to be detectable at M=1e5, as a function of rate
        T_at = boundary_T(np.array([1e5]), N, rate, T_grid)[0]
        print(f"  max(−Q_bb)={rate:g}:  T needed at M=1e5 ≈ {T_at:.2f} subs/site")


if __name__ == "__main__":
    main()
