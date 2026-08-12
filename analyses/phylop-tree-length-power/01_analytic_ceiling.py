#!/usr/bin/env python3
"""
01_analytic_ceiling.py - the centerpiece analytic figure.

Sweeps total tree length T (by uniformly rescaling each family's real neutral
tree) and plots the maximum achievable per-site conserved score,
-log10(p_ceiling), against T. Overlays:

  * one curve per dataset family (hamster / birds / turtle / mammal), each from
    its own real fitted GTR model - they nearly coincide because the ceiling is
    governed by T, with mild model-dependent (kappa) scatter;
  * a marker at each family's actual T;
  * the genome-wide BH-FDR detectability bar -log10(alpha/N) for a range of
    per-chromosome site counts N - the most-conserved site can only clear FDR
    when its ceiling rises above this bar.

No simulated data: every point evaluates a real model's own likelihood for a
fully-conserved column (see phylop_power.py).

Outputs:
  figures/analytic_ceiling_vs_treelen.png
  analytic_ceiling_table.csv
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import phylop_power as pp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))

# One representative corrected neutral model per family (the 4 the user named).
FAMILIES = [
    ("hamster (15 sp)",  "data/hamsters/workflow-tests/test-full/02-neutral-model/phylofit/autosomes/CM000995.3-corrected.mod", "#4C72B0"),
    ("birds (44 sp)",    "data/birds/test1/03-phylofit/group1/chr10-corrected.mod",                                            "#55A868"),
    ("turtle (22 sp)",   "data/turtles/small-test2/03-phylofit/shortest-scaffolds/NC_050095.1-corrected.mod",                 "#C44E52"),
    ("mammal (241 sp)",  "data/mammals-v3/echolocation-v3/02-neutral-model/phylofit/autosomes/chr1-corrected.mod",            "#8172B3"),
]

# BH bar: a conserved site can only survive genome-wide BH-FDR if its p-value
# clears ~alpha/N, i.e. -log10(p) must exceed -log10(alpha/N).
ALPHA = 0.05
N_VALUES = [1e6, 1e7, 1e8]

T_GRID = np.geomspace(0.3, 32, 60)


def sweep_family(mod_path):
    md = pp.parse_mod(os.path.join(REPO, mod_path))
    root = pp.parse_newick(md["tree_newick"])
    T_base = pp.tree_length(root)
    Ts, ceils, kappas = [], [], []
    for T in T_GRID:
        s = T / T_base
        c = pp.ceiling(md, rho_scale=s)
        Ts.append(c["tree_length"])
        ceils.append(c["neglog10p"])
        kappas.append(c["kappa"])
    # the family's actual (unscaled) point
    c0 = pp.ceiling(md, rho_scale=1.0)
    return np.array(Ts), np.array(ceils), np.array(kappas), c0


def main():
    fig, ax = plt.subplots(figsize=(8.2, 5.6))

    rows = []
    for label, mod, color in FAMILIES:
        Ts, ceils, kappas, c0 = sweep_family(mod)
        ax.plot(Ts, ceils, color=color, lw=1.8, alpha=0.9, label=f"{label}")
        ax.scatter([c0["tree_length"]], [c0["neglog10p"]], color=color,
                   s=70, zorder=5, edgecolor="black", linewidth=0.6)
        rows.append({
            "family": label, "mod": mod, "tips": c0["tips"],
            "tree_length": round(c0["tree_length"], 3),
            "LRT_max": round(c0["LRT_max"], 3),
            "kappa": round(c0["kappa"], 3),
            "ceiling_neglog10p": round(c0["neglog10p"], 3),
            "ceiling_p": c0_p(c0["neglog10p"]),
        })

    # BH detectability bars
    for N in N_VALUES:
        bar = -np.log10(ALPHA / N)
        ax.axhline(bar, color="gray", ls="--", lw=1.0, alpha=0.7)
        ax.text(0.32, bar + 0.15, f"BH bar, N=10^{int(np.log10(N))} sites  (need −log10 p > {bar:.1f})",
                fontsize=7.5, color="gray")

    # shade the "no conserved site can ever pass" region (below the gentlest bar)
    gentlest = -np.log10(ALPHA / min(N_VALUES))
    ax.axhspan(0, gentlest, color="gray", alpha=0.06)

    ax.set_xscale("log")
    ax.set_xlabel("total neutral tree length  T  (expected substitutions / site)")
    ax.set_ylabel("max achievable conserved score  −log10(p_ceiling)")
    ax.set_title("phyloP conserved-score ceiling is set by tree length\n"
                 "(per-site LRT/CONACC; ceiling = perfectly conserved column)")
    ax.set_ylim(0, 26)
    ax.legend(loc="upper left", fontsize=8.5, framealpha=0.9)
    ax.grid(True, which="both", alpha=0.15)

    # annotate the crossing: where curves clear the N=1e7 bar
    ax.axvline(11, color="black", ls=":", lw=1.0, alpha=0.5)
    ax.text(11.4, 1.0, "T ≈ 11\n(threshold for\nN=10^7)", fontsize=7.5, alpha=0.7)

    fig.tight_layout()
    out_png = os.path.join(HERE, "figures", "analytic_ceiling_vs_treelen.png")
    fig.savefig(out_png, dpi=150)
    print(f"wrote {out_png}")

    out_csv = os.path.join(HERE, "analytic_ceiling_table.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_csv}")
    print()
    for r in rows:
        print(f"  {r['family']:<16} T={r['tree_length']:>6.2f}  "
              f"ceiling −log10p={r['ceiling_neglog10p']:>6.2f}  (p={r['ceiling_p']:.2e})  kappa={r['kappa']}")


def c0_p(neglog10p):
    return float(10 ** (-neglog10p))


if __name__ == "__main__":
    main()
