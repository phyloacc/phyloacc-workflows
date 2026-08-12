#!/usr/bin/env python3
"""
03_plot_sim.py - plot the simulation power sweep (reads sim_power_results.csv).

Two panels:
  A. Sensitivity: % of genuinely-conserved (rho=0.1) sites recovered as conserved
     through real phyloP + the pipeline's BH-FDR, vs total tree length T. Vertical
     markers show where the four real dataset families actually fall.
  B. Mechanism: the max conserved score reached in the simulation vs T, overlaid on
     the analytic ceiling curve (same mammal model) and the genome-wide BH bar - the
     sim tracks the analytic ceiling, and detection switches on only once the ceiling
     clears the bar (~T=8-11 depending on N).

Outputs: figures/sim_power_vs_treelen.png
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
MODEL = os.path.join(REPO, "data/mammals-v3/echolocation-v3/02-neutral-model/phylofit/autosomes/chr1-corrected.mod")

# real families: (label, tree_length, color) - where they actually fall on the T axis
FAMILIES = [
    ("hamster", 0.93, "#4C72B0"),
    ("birds", 1.09, "#55A868"),
    ("turtle", 2.81, "#C44E52"),
    ("mammal", 16.07, "#8172B3"),
]


def load():
    with open(os.path.join(HERE, "sim_power_results.csv")) as f:
        rows = list(csv.DictReader(f))
    T = np.array([float(r["T"]) for r in rows])
    sens = np.array([float(r["sensitivity_pct"]) for r in rows])
    maxc = np.array([float(r["max_conserved_score"]) for r in rows])
    N_bg = int(rows[0]["bg_sites"]) + int(rows[0]["conserved_total"])
    alpha = float(rows[0]["alpha"])
    rho = float(rows[0]["rho_cons"])
    return T, sens, maxc, N_bg, alpha, rho


def main():
    T, sens, maxc, N_sim, alpha, rho = load()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.4))

    # ---- Panel A: sensitivity vs T ----
    axA.plot(T, sens, "-o", color="#333333", lw=2, ms=6, zorder=4)
    axA.set_xscale("log")
    axA.set_xlabel("total neutral tree length  T  (expected substitutions / site)")
    axA.set_ylabel(f"% of genuinely-conserved sites recovered\n(true elements at ρ={rho} of neutral rate)")
    axA.set_title("Simulation: power to recover real conservation\n"
                  "through the real phyloP + pipeline BH-FDR")
    axA.set_ylim(-3, 103)
    axA.grid(True, which="both", alpha=0.15)

    # mark where the real families fall
    for label, Tf, color in FAMILIES:
        axA.axvline(Tf, color=color, ls="--", lw=1.4, alpha=0.8)
        ypos = 88 if label != "mammal" else 40
        axA.text(Tf * 1.03, ypos, f"{label}\nT={Tf}", color=color, fontsize=8,
                 rotation=90, va="top", ha="left")
    axA.axhspan(-3, 3, color="red", alpha=0.05)
    axA.text(1.0, 6, "shallow-tree families → 0% recovered\n(genuine conservation is invisible)",
             fontsize=8.5, color="#993333")

    # ---- Panel B: max score vs T, analytic ceiling, BH bars ----
    md = pp.parse_mod(MODEL)
    root = pp.parse_newick(md["tree_newick"])
    Tbase = pp.tree_length(root)
    Tgrid = np.geomspace(0.8, 26, 50)
    ceil = np.array([pp.ceiling(md, rho_scale=t / Tbase)["neglog10p"] for t in Tgrid])

    axB.plot(Tgrid, ceil, color="#8172B3", lw=1.8, label="analytic ceiling (perfect conservation)")
    axB.plot(T, maxc, "o", color="#333333", ms=6, label="simulation: max conserved score reached")
    for N, style in [(N_sim, ":"), (6_000_000, "--")]:
        bar = -np.log10(alpha / N)
        axB.axhline(bar, color="gray", ls=style, lw=1.1)
        axB.text(0.85, bar + 0.3, f"BH bar, N={N:,} sites (−log10 p > {bar:.1f})",
                 fontsize=7.5, color="gray")
    axB.set_xscale("log")
    axB.set_xlabel("total neutral tree length  T  (expected substitutions / site)")
    axB.set_ylabel("max conserved score  −log10(p)")
    axB.set_title("Why: the conserved-score ceiling must clear the FDR bar\n"
                  "(simulation tracks the analytic ceiling)")
    axB.set_ylim(0, 23)
    axB.legend(loc="upper left", fontsize=8.5)
    axB.grid(True, which="both", alpha=0.15)

    fig.tight_layout()
    out = os.path.join(HERE, "figures", "sim_power_vs_treelen.png")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
