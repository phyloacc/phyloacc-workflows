#!/usr/bin/env python3
"""
07_theoretical_detectability.py - a PURELY THEORETICAL detectability nomogram.

Same axes as 06 (tree depth T x number of sites tested M), but the boundary is drawn
from first principles - NO fitted models, no dataset benchmarks. A site is detectable iff

    ceiling(T)  >=  log10(M / alpha)

In the many-taxon limit the score ceiling has a closed form with a single model knob:

    ceiling(T) = -log10 erfc( sqrt(kappa * T) )

  * kappa = max per-base substitution rate (rate matrix normalized to mean rate 1), so
    kappa >= 1, with kappa = 1 exactly for a flat / Jukes-Cantor model. kappa just
    stretches the T axis by 1/kappa; it does not change the shape. For real nucleotide
    models kappa is a narrow ~1.3-1.5.
The band spans kappa = 1 (flat model, conservative) to kappa = 1.5 (typical nucleotide
structure). Taxon count enters as a *saturating* finite-sample correction (matters only
below ~50 taxa); a dashed line shows a sparse 20-taxon Jukes-Cantor tree for reference.

Output: figures/theoretical_detectability.png
"""
import os
import math
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import phylop_power as pp  # only for neglog10_from_lrt (the LRT -> -log10 p map)

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = 0.05


def ceiling_many_taxon(T, kappa):
    # many-taxon limit: LRT_max -> 2*kappa*T ; ceiling = -log10 erfc(sqrt(kappa*T))
    return pp.neglog10_from_lrt(2.0 * kappa * T)


def ceiling_jc_star(n, T):
    # closed-form Jukes-Cantor on a star tree of n tips, total depth T (kappa = 1)
    x = math.exp(-4.0 * T / (3.0 * n))
    p_same = 0.25 + 0.75 * x
    p_diff = 0.25 - 0.25 * x
    lrt = -2.0 * math.log(p_same ** n + 3.0 * p_diff ** n)
    return pp.neglog10_from_lrt(lrt)


def boundary_T(M_grid, ceil_of_T, T_grid):
    ceil = np.array([ceil_of_T(t) for t in T_grid])
    bar = np.log10(M_grid / ALPHA)
    return np.interp(bar, ceil, T_grid, left=np.nan, right=np.nan)


def main():
    T_grid = np.geomspace(0.1, 40, 200)
    M_grid = np.geomspace(1e1, 1e9, 200)

    Tb_flat = boundary_T(M_grid, lambda t: ceiling_many_taxon(t, 1.0), T_grid)   # kappa=1 (right edge)
    Tb_struct = boundary_T(M_grid, lambda t: ceiling_many_taxon(t, 1.5), T_grid)  # kappa=1.5 (left edge)
    Tb_sparse = boundary_T(M_grid, lambda t: ceiling_jc_star(20, t), T_grid)     # 20-taxon JC (dashed)
    clamp = lambda a: np.where(np.isnan(a), 40.0, a)
    Tb_flat_f, Tb_struct_f = clamp(Tb_flat), clamp(Tb_struct)

    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    ax.fill_betweenx(M_grid, Tb_flat_f, 40, color="#c7dcef", alpha=0.7, zorder=0)     # detectable (any kappa)
    ax.fill_betweenx(M_grid, 0.1, Tb_struct_f, color="#f2c9c4", alpha=0.7, zorder=0)  # power-limited (any kappa)
    ax.fill_betweenx(M_grid, Tb_struct_f, Tb_flat_f, color="#ece3b0", alpha=0.75, zorder=0)
    ax.plot(Tb_struct, M_grid, color="#8a6d1a", lw=1.6, zorder=3)
    ax.plot(Tb_flat, M_grid, color="#8a6d1a", lw=1.6, zorder=3)
    ax.plot(Tb_sparse, M_grid, color="#555555", lw=1.5, ls="--", zorder=4)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("total neutral tree depth  T  (sum of branch lengths, subs/site)")
    ax.set_ylabel("number of sites tested  M  (multiple-testing burden)")
    ax.set_title("Where can per-site phyloP make a call?  (first-principles)\n"
                 "detectable iff  −log10 erfc(√(κT)) ≥ log10(M/α),  α=0.05")
    ax.set_xlim(0.1, 40)
    ax.set_ylim(1e1, 1e9)
    ax.legend(handles=[Patch(fc="#f2c9c4", label="power-limited (any model)"),
                       Patch(fc="#ece3b0", label="band: κ = 1 (flat) → 1.5 (structured)"),
                       Patch(fc="#c7dcef", label="detectable (any model)"),
                       Line2D([0], [0], color="#555555", ls="--", label="20-taxon JC (finite-taxa shift)")],
              loc="upper left", fontsize=8, framealpha=0.92)

    # horizontal region labels (no rotation)
    ax.text(24, 1.5e2, "single-site calls\nDETECTABLE", fontsize=12, weight="bold",
            color="#1f3a7a", ha="center", va="center")
    ax.text(0.55, 3e4, "single-site power lost\n(no site can pass,\nhowever conserved)",
            fontsize=10.5, weight="bold", color="#7a1f1f", ha="center", va="center")
    ax.text(0.11, 1.25e1,
            "purely theoretical: ceiling = −log10 erfc(√(κT)), many-taxon limit. κ = max per-base rate "
            "(=1 for a flat/JC model);\nreal nucleotide models κ≈1.3–1.5. Taxon count saturates by ~50 (dashed = 20).",
            fontsize=6.6, color="gray")

    # Finalize layout/transforms, then align the two curve-following labels to the
    # actual display slope of the curve they annotate (so they sit ON the curve).
    fig.tight_layout()
    fig.canvas.draw()
    logM = np.log10(M_grid)

    def T_at(Tarr, M):
        return float(np.interp(np.log10(M), logM, Tarr))

    def disp_angle(Tarr, M):
        p0 = ax.transData.transform((T_at(Tarr, M / 1.8), M / 1.8))
        p1 = ax.transData.transform((T_at(Tarr, M * 1.8), M * 1.8))
        return math.degrees(math.atan2(p1[1] - p0[1], p1[0] - p0[0]))

    Tmid = (Tb_struct_f + Tb_flat_f) / 2.0
    Tsp = clamp(Tb_sparse)
    Mb = 4e6
    ax.text(T_at(Tmid, Mb), Mb, "depends on substitution structure (κ)", fontsize=8.5, style="italic",
            color="#6b5310", ha="center", va="center", rotation=disp_angle(Tmid, Mb), rotation_mode="anchor")
    Ms = 4e6
    ax.text(T_at(Tsp, Ms) * 1.18, Ms, "20-taxon JC", fontsize=8, color="#3a3a3a",
            ha="center", va="center", rotation=disp_angle(Tsp, Ms), rotation_mode="anchor")

    out = os.path.join(HERE, "figures", "theoretical_detectability.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
    print("\nboundary tree depth T at a few M (kappa=1 flat / kappa=1.5 structured / 20-taxon JC):")
    for m in (1e2, 1e6, 1e8):
        i = int(np.argmin(np.abs(M_grid - m)))
        print(f"  M={m:.0e}:  flat T≈{Tb_flat[i]:5.1f}   structured T≈{Tb_struct[i]:5.1f}   20-taxon T≈{Tb_sparse[i]:5.1f}")


if __name__ == "__main__":
    main()
