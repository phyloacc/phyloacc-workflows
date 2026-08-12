#!/usr/bin/env python3
"""
05_element_size.py - predict the minimum detectable ELEMENT size vs tree length.

Extends the single-base ceiling to multi-base elements. This is exactly the framework
the 29-mammals paper (Lindblad-Toh et al. 2011) used with its P1 / P12 quantities.

Key idea - evidence adds across sites:
  P1(T)  = probability a NEUTRAL single base is invariant (identical in all species)
           = the "single-base" surprise; computed from the real model (phylop_power).
  P_L    = probability a neutral L-mer is entirely invariant ~ P1(T)**L  (independent sites).
An element of L perfectly-conserved bases is detectable when a neutral L-mer being that
conserved is rarer than the genome-wide bar: P1(T)**L < alpha/N. Solving for L:

  L_min(T, N) = ln(N / alpha) / ln(1 / P1(T)).

So minimum element size scales as 1 / ln(1/P1(T)) (shrinks fast as the tree deepens) and
only ~ln(N) in the number of tests. This is the perfect-conservation floor; real elements
(not perfectly invariant) need to be somewhat larger.

Cross-check vs the paper: at T=4.5 (29 mammals) our P1 ~ 0.023 (paper: P1 < 0.02) and
P1**12 ~ 3e-20 (paper: P12 < 1e-25) - same regime; 12-mers hugely detectable, single
bases not. Reproduced from real fitted models, no simulated data.

Outputs: figures/min_element_size.png
"""

import math
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

# real families (own models) + the two mammalian landmark datasets (by tree length)
FAMILIES = [
    ("hamster", 0.93, "data/hamsters/workflow-tests/test-full/02-neutral-model/phylofit/autosomes/CM000995.3-corrected.mod", "#4C72B0"),
    ("birds", 1.09, "data/birds/test1/03-phylofit/group1/chr10-corrected.mod", "#55A868"),
    ("turtle", 2.81, "data/turtles/small-test2/03-phylofit/shortest-scaffolds/NC_050095.1-corrected.mod", "#C44E52"),
    ("mammal / Zoonomia 241", 16.07, "data/mammals-v3/echolocation-v3/02-neutral-model/phylofit/autosomes/chr1-corrected.mod", "#8172B3"),
]
LANDMARKS = [("29 mammals\n(Lindblad-Toh 2011)", 4.5), ("Zoonomia 241\n(Christmas 2023)", 16.0)]


def L_min_exact(P1, N):
    """Min element size using the EXACT invariant probability P1 (P_L = P1^L).
    This is the 29-mammals paper's criterion - the strictest / most conservative."""
    if P1 >= 1.0:
        return float("inf")
    return math.log(N / ALPHA) / math.log(1.0 / P1)


def L_min_chi2(lrt_per_site, N):
    """Min element size using phyloP's actual asymptotic (chi^2) score: an element of
    L perfectly-conserved sites has LRT ~ L*lrt_per_site and p = 2*Phi(-sqrt(LRT)).
    This is what the pipeline AND Zoonomia actually use; it is less conservative than
    the exact-P1 criterion (phyloP's tail is anticonservative). Solve for L such that
    -log10(p_element) = -log10(alpha/N)."""
    target = math.log10(N / ALPHA)
    lo, hi = 1e-3, 1e5
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if pp.neglog10_from_lrt(mid * lrt_per_site) < target:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def main():
    md = pp.parse_mod(MODEL)
    Tbase = pp.tree_length(pp.parse_newick(md["tree_newick"]))
    N = 1e8

    T = np.geomspace(0.5, 30, 200)
    P1_T = np.array([pp.invariant_probability(md, t / Tbase) for t in T])
    LRT_T = np.array([pp.ceiling(md, t / Tbase)["LRT_max"] for t in T])  # per-site ceiling LRT
    L_exact = np.array([L_min_exact(p, N) for p in P1_T])
    L_chi2 = np.array([L_min_chi2(s, N) for s in LRT_T])

    fig, ax = plt.subplots(figsize=(9.4, 6.2))

    # two p-value conventions (two curves), with a light band between them
    ax.fill_between(T, L_chi2, L_exact, color="#B0B7C6", alpha=0.18, zorder=1)
    ax.plot(T, L_exact, "-", color="#333333", lw=1.8, zorder=3,
            label="exact invariant-probability criterion P₁  (stricter; 29-mammals paper)")
    ax.plot(T, L_chi2, "-", color="#C44E52", lw=2.4, zorder=4,
            label="phyloP's actual scoring (χ²)  — what phyloP achieves")

    ax.axhline(1, color="crimson", ls="--", lw=1.3)
    ax.text(0.52, 1.05, "single-base resolution (L = 1)", color="crimson", fontsize=8.5)

    xcross = float(np.interp(1.0, L_chi2[::-1], T[::-1]))
    ax.plot([xcross], [1.0], "o", color="crimson", ms=7, zorder=8)
    ax.annotate(f"phyloP reaches single base at T ≈ {xcross:.0f}\n→ Zoonomia (T=16) clears it",
                (xcross, 1.0), textcoords="offset points", xytext=(-6, 20), fontsize=8,
                color="crimson", ha="center")

    for lab, Tl in LANDMARKS:
        ax.axvline(Tl, color="black", ls="-", lw=0.8, alpha=0.30)
        ax.text(Tl * 1.02, 50, lab, fontsize=8, rotation=90, va="top")

    # dataset points = what phyloP can ACTUALLY resolve → placed on the phyloP (χ²) curve
    ax.scatter([], [], s=85, color="none", edgecolor="black", label="dataset (on the phyloP curve = actually resolvable)")
    print(f"{'dataset':<22} {'T':>5} {'L_min χ²':>9} {'L_min P₁':>9}   (min conserved element, bp; N=1e8)")
    for name, Tf, mod, color in FAMILIES:
        m = pp.parse_mod(os.path.join(REPO, mod))
        p1 = pp.invariant_probability(m)
        s = pp.ceiling(m)["LRT_max"]
        le, lc = L_min_exact(p1, N), L_min_chi2(s, N)
        short = name.split(" /")[0]
        ax.scatter([Tf], [lc], s=90, color=color, edgecolor="black", linewidth=0.8, zorder=6)
        dx = 9 if Tf < 12 else -11
        ha = "left" if Tf < 12 else "right"
        lab = f"{short}: ~{lc:.0f} bp" if lc >= 1 else f"{short}: single base"
        ax.annotate(lab, (Tf, lc), textcoords="offset points", xytext=(dx, 4),
                    fontsize=8.5, color=color, weight="bold", ha=ha)
        print(f"{name:<22} {Tf:>5.2f} {lc:>8.1f} {le:>9.1f}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("total neutral tree length  T  (expected substitutions / site)")
    ax.set_ylabel("minimum detectable element size  L_min  (bp)")
    ax.set_title("Minimum detectable conserved-element size vs tree length  (N=10⁸)")
    ax.set_ylim(0.55, 55)
    ax.set_xlim(0.5, 30)
    ax.legend(loc="upper right", fontsize=8.0)
    ax.grid(True, which="both", alpha=0.15)
    ax.text(0.53, 1.55, "dots = each dataset on the phyloP (red) curve = the smallest element phyloP can\n"
                        "actually resolve there. Both curves are perfect-conservation floors (lower bounds):\n"
                        "nothing smaller is detectable at any level; real imperfect elements need to be larger.",
            fontsize=7.0, color="gray")

    fig.tight_layout()
    out = os.path.join(HERE, "figures", "min_element_size.png")
    fig.savefig(out, dpi=150)
    print(f"\nwrote {out}")
    print(f"chi2 (real-method) single-base threshold: T ~ {xcross:.1f} subs/site  (Zoonomia T=16 clears it)")
    print(f"cross-check P1 at T=4.5: {pp.invariant_probability(md, 4.5/Tbase):.4f} "
          f"(paper P1<0.02); P1^12 = {pp.invariant_probability(md, 4.5/Tbase)**12:.2e} (paper P12<1e-25)")


if __name__ == "__main__":
    main()
