#!/usr/bin/env python3
"""
08_taxon_walls.py - detectability nomogram as a "how many taxa do you need" map.

Same axes as 07 (tree depth T x number of sites tested M), but here EACH contour is the
detectable/undetectable boundary for a DIFFERENT taxon count n. The detectable region is
not fixed: with fewer taxa the boundary moves RIGHT (you need a deeper tree), because

    ceiling(n, T) = -log10 erfc(sqrt(LRT/2)),   LRT ~= min( 2*kappa*T , 2*(n-1)*ln4 )

is increasing in BOTH T and n. So the plane splits into nested zones by the MINIMUM taxon
count that can make a call there:

    left of the n->inf curve : impossible for ANY n (depth too small for the test burden)
    then, moving right        : need >20 taxa ... 12-20 ... 8-12 ... <=8 taxa (easy)

Two limits are visible at once:
  * kappa (=1 here, JC) only stretches the T axis - fixed at 1 so the contours isolate
    the TAXON effect. (06/07 show the kappa/model band.)
  * n sets a HARD CEILING 2(n-1)ln4 on the score. Once log10(M/alpha) exceeds it, that n's
    boundary flattens into a horizontal "taxon wall" at M_wall(n)=alpha*10^cap(n) and its
    zone pinches shut - no tree depth rescues it. Fewer taxa -> lower wall.

Output: figures/taxon_walls.png
"""
import os
import math
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import phylop_power as pp  # neglog10_from_lrt: LRT -> -log10 p

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = 0.05
XMAX = 40.0


def ceiling_many_taxon(T, rate=1.0):
    return pp.neglog10_from_lrt(2.0 * rate * T)           # n -> inf limit (LRT -> 2*rate*T), rate = max exit rate


def ceiling_jc_star(n, T):
    x = math.exp(-4.0 * T / (3.0 * n))                    # JC star tree, n tips, depth T
    p_same, p_diff = 0.25 + 0.75 * x, 0.25 - 0.25 * x
    return pp.neglog10_from_lrt(-2.0 * math.log(p_same ** n + 3.0 * p_diff ** n))


def taxon_cap_ceiling(n):
    return pp.neglog10_from_lrt(2.0 * (n - 1) * math.log(4.0))   # saturated: LRT -> 2(n-1)ln4


def boundary_T(M_grid, ceil_of_T, T_grid):
    ceil = np.array([ceil_of_T(t) for t in T_grid])
    return np.interp(np.log10(M_grid / ALPHA), ceil, T_grid, left=np.nan, right=np.nan)


def main():
    T_grid = np.geomspace(0.1, 400, 400)   # wide, so finite-n curves reach their walls
    M_grid = np.geomspace(1e1, 1e9, 300)

    # boundaries, easiest (most taxa) -> hardest (fewest taxa); nested left -> right
    Tb_inf = boundary_T(M_grid, ceiling_many_taxon, T_grid)
    Tb = {n: boundary_T(M_grid, lambda t, n=n: ceiling_jc_star(n, t), T_grid) for n in (20, 12, 8)}
    walls = {n: ALPHA * 10 ** taxon_cap_ceiling(n) for n in (20, 12, 8)}
    clamp = lambda a: np.where(np.isnan(a), XMAX, a)
    Bi, B20, B12, B8 = clamp(Tb_inf), clamp(Tb[20]), clamp(Tb[12]), clamp(Tb[8])

    fig, ax = plt.subplots(figsize=(9.8, 6.7))
    # zones by MINIMUM taxa needed: pale (needs many) near the red edge -> strong (few taxa) at right
    ax.fill_betweenx(M_grid, 0.1, Bi, color="#f2c9c4", alpha=0.8, zorder=0)   # impossible, any n
    ax.fill_betweenx(M_grid, Bi, B20, color="#d9e7f3", alpha=0.9, zorder=0)  # >20
    ax.fill_betweenx(M_grid, B20, B12, color="#9fc3e0", alpha=0.9, zorder=0)  # 12-20
    ax.fill_betweenx(M_grid, B12, B8, color="#5f9bd1", alpha=0.9, zorder=0)   # 8-12
    ax.fill_betweenx(M_grid, B8, XMAX, color="#2f6fb0", alpha=0.85, zorder=0) # <=8 (easy)

    # faint reference lines: the detectability boundary for EVERY n = 1..20 (thin, grey)
    for n in range(1, 21):
        if n in (8, 12, 20):
            continue  # drawn thick below
        Tn = boundary_T(M_grid, lambda t, n=n: ceiling_jc_star(n, t), T_grid)
        if not np.all(np.isnan(Tn)):
            ax.plot(Tn, M_grid, ls=":", lw=0.6, color="#3a3a3a", alpha=0.45, zorder=2)

    # the contour lines = the moving boundary, one per taxon count
    ax.plot(Tb_inf, M_grid, color="#20456e", lw=1.7, zorder=4)
    for n, col in ((20, "#20456e"), (12, "#20456e"), (8, "#20456e")):
        ax.plot(Tb[n], M_grid, color=col, lw=1.7, zorder=4)
        w = walls[n]
        if M_grid[0] < w < M_grid[-1]:
            ax.axhline(w, color="#20456e", lw=1.0, ls=":", zorder=4)
            ax.text(XMAX * 0.97, w * 1.35, f"taxon wall, n={n}  (M≈{w:.0e}): no depth helps above",
                    color="#20456e", fontsize=7.3, ha="right", va="bottom", zorder=6)

    # label each boundary with the taxon count it belongs to
    logM = np.log10(M_grid)
    T_at = lambda Tarr, M: float(np.interp(np.log10(M), logM, Tarr))
    ax.text(T_at(Tb_inf, 4e8), 4e8, "n→∞", color="#20456e", fontsize=8.5, weight="bold",
            ha="right", va="center", rotation=90, zorder=6)
    ax.text(T_at(Tb[20], 4e8), 4e8, "n=20", color="#20456e", fontsize=8.5, weight="bold",
            ha="right", va="center", rotation=90, zorder=6)
    ax.text(T_at(Tb[12], 3e5), 3e5, "n=12", color="#123", fontsize=8.5, weight="bold",
            ha="left", va="bottom", zorder=6)
    ax.text(T_at(Tb[8], 1.2e3), 1.2e3, "n=8", color="#fff", fontsize=8.5, weight="bold",
            ha="left", va="bottom", zorder=6)

    # --- rate effect on the n→∞ edge (dashed red) -------------------------------------
    # A larger max exit rate scales the depth term 2*rate*T, sliding the rising part of
    # EVERY contour left by the factor 1/rate. Only the n→∞ edge is drawn here (it is all
    # rising part, no wall, so it slides cleanly and stays clear of the taxon contours);
    # the same left-shift applies to the rising part of the finite-n curves, and taxon
    # walls are rate-independent.
    for rate_r, ls_r in ((1.5, (0, (6, 3))), (2.0, (0, (3, 3)))):
        Tb_r = clamp(boundary_T(M_grid, lambda t, r=rate_r: ceiling_many_taxon(t, r), T_grid))
        ax.plot(Tb_r, M_grid, color="#20456e", lw=1.4, ls=ls_r, zorder=5)
        ax.text(T_at(Tb_r, 3e7), 3e7, f"×{rate_r:g}", color="#20456e", fontsize=8, weight="bold",
                ha="right", va="center", rotation=90, zorder=6)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("total neutral tree depth  T  (sum of branch lengths, subs/site)")
    ax.set_ylabel("number of sites tested  M  (multiple-testing burden)")
    ax.set_title("How many taxa do you need?  each contour is the detectability boundary for one taxon count\n"
                 "detectable iff  −log10 erfc(√(LRT/2)) ≥ log10(M/α);   LRT ≈ min(2·T·max(−Q_bb), 2(n−1)ln4)\n"
                 "contours & zones at max(−Q_bb)=1;   dashed = n→∞ edge at max(−Q_bb) 1.5 & 2;   α=0.05")
    ax.set_xlim(0.1, XMAX); ax.set_ylim(1e1, 1e9)

    ax.legend(handles=[
        Patch(fc="#f2c9c4", label="impossible for any n (depth too small)"),
        Patch(fc="#d9e7f3", label="detectable only with >20 taxa"),
        Patch(fc="#9fc3e0", label="need 12–20 taxa"),
        Patch(fc="#5f9bd1", label="need 8–12 taxa"),
        Patch(fc="#2f6fb0", label="detectable with ≤8 taxa (easy)"),
        Line2D([0], [0], color="#20456e", lw=1.0, ls=":", label="taxon wall (score cap 2(n−1)ln4)"),
        Line2D([0], [0], color="#20456e", lw=1.4, ls="--", label="n→∞ edge at max(−Q_bb)=1.5, 2 (faster → left)"),
    ], loc="upper left", fontsize=7.8, framealpha=0.95, title="minimum taxa needed", title_fontsize=8)

    ax.annotate("fewer taxa → boundary moves right\n(need a deeper tree, or you can't call it at all)",
                xy=(0.42, 0.42), xycoords="axes fraction", ha="center", va="center",
                fontsize=8.5, style="italic", color="#5a2a2a")
    ax.text(0.11, 1.25e1,
            "purely theoretical, equal rates and equal base frequencies. Each contour is a full "
            "detectable/undetectable boundary; with fewer taxa\nit shifts right and develops a horizontal wall at "
            "2(n−1)ln4 above which no depth detects. A faster top substitution rate only stretches the T axis.",
            fontsize=6.6, color="gray", zorder=6)

    fig.tight_layout()
    out = os.path.join(HERE, "figures", "taxon_walls.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
    for n in (8, 12, 20, 50):
        cap = taxon_cap_ceiling(n)
        print(f"  n={n:3d}:  cap ceiling ≈ {cap:5.1f}   M_wall ≈ {ALPHA*10**cap:.1e}")


if __name__ == "__main__":
    main()
