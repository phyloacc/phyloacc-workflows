#!/usr/bin/env python3
"""
02_simulate_power.py - the simulation half (uses AliSim; simulated data was
explicitly approved by the user for this task).

Controlled experiment: hold topology + substitution model fixed (the real
241-taxon mammal neutral model) and vary ONLY the total tree length T by
uniformly rescaling the branches. At each T:

  1. Simulate a NEUTRAL background alignment (branches at scale T) and a genuinely
     CONSERVED block (branches at scale RHO_CONS * T, i.e. evolving RHO_CONS as
     fast) with AliSim, under the model's own GTR parameters. Concatenate them.
  2. Score every site with the REAL phyloP (LRT / CONACC), exactly as the pipeline
     does.
  3. Apply the pipeline's exact Benjamini-Hochberg FDR (re-implemented to match
     utils/adjust_pvals.sh) across all sites.
  4. Record: how many of the truly-conserved sites are recovered (class=conserved
     & FDR-significant) = sensitivity; false positives in the neutral background;
     and the max conserved score reached.

This shows the tree-length power wall biting end to end through the real tool +
the real FDR: on shallow trees even genuinely conserved sites can't be recovered;
past T ~ 10 they can.

Everything here is fully reproducible from a real neutral model; the only
fabricated inputs are the AliSim-simulated alignments (approved for this task).

Outputs:
  sim_power_results.csv
"""

import csv
import os
import re
import subprocess
import sys

import numpy as np

import phylop_power as pp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))

IQTREE = "/n/home07/gthomas/miniconda3/envs/iqtree-env/bin/iqtree3"
PHYLOP = "/n/home07/gthomas/miniconda3/envs/phyloacc-workflows/bin/phyloP"

MODEL = os.path.join(REPO, "data/mammals-v3/echolocation-v3/02-neutral-model/phylofit/autosomes/chr1-corrected.mod")

L_BG = 100_000     # neutral background sites
L_CONS = 5_000     # genuinely conserved sites
RHO_CONS = 0.1     # conserved block evolves at 10% of the neutral rate
ALPHA = 0.05
T_GRID = [1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 24.0]
SEED = 20260725


def gtr_string(md):
    """Build an AliSim GTR{...}+F{...} model string from a PHAST .mod (REV model).

    For a reversible model Q_ij = r_ij * pi_j (i != j), so the symmetric
    exchangeabilities are r_ij = Q_ij / pi_j. AliSim order: AC, AG, AT, CG, CT, GT.
    """
    Q, pi = md["Q"], md["pi"]
    r_ac = Q[0, 1] / pi[1]
    r_ag = Q[0, 2] / pi[2]
    r_at = Q[0, 3] / pi[3]
    r_cg = Q[1, 2] / pi[2]
    r_ct = Q[1, 3] / pi[3]
    r_gt = Q[2, 3] / pi[3]
    rates = "/".join(f"{x:.6f}" for x in [r_ac, r_ag, r_at, r_cg, r_ct, r_gt])
    freqs = "/".join(f"{x:.6f}" for x in pi)
    return f"GTR{{{rates}}}+F{{{freqs}}}"


def scale_newick(newick, s):
    """Multiply every branch length in a newick string by s."""
    return re.sub(r":([0-9eE.+\-]+)", lambda m: f":{float(m.group(1)) * s:.8f}", newick)


def write_scaled_model(md, T, path):
    """Write a PHAST .mod identical to the real one but with the tree scaled to
    total length T (this is the neutral model phyloP scores against)."""
    root = pp.parse_newick(md["tree_newick"])
    s = T / pp.tree_length(root)
    scaled_tree = scale_newick(md["tree_newick"], s)
    with open(MODEL) as f:
        lines = f.readlines()
    with open(path, "w") as o:
        for line in lines:
            if line.startswith("TREE:"):
                o.write(f"TREE: {scaled_tree};\n" if not scaled_tree.endswith(";") else f"TREE: {scaled_tree}\n")
            else:
                o.write(line)


def alisim(tree_path, model_str, length, out_prefix, seed):
    cmd = [IQTREE, "--alisim", out_prefix, "-t", tree_path, "-m", model_str,
           "--length", str(length), "-seed", str(seed), "-af", "fasta", "-redo"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return out_prefix + ".fa"


def read_fasta(path):
    seqs = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def bh_fdr(pvals, alpha):
    """Benjamini-Hochberg, matching utils/adjust_pvals.sh. Returns bool sig array."""
    p = np.asarray(pvals, dtype=float)
    N = len(p)
    order = np.argsort(p, kind="mergesort")
    ranks = np.arange(1, N + 1)
    padj_sorted = p[order] * N / ranks
    padj_sorted = np.minimum.accumulate(padj_sorted[::-1])[::-1]  # enforce monotonicity
    padj = np.empty(N)
    padj[order] = np.minimum(padj_sorted, 1.0)
    return padj < alpha


def run_phylop_scores(mod_path, fasta_path, tmp):
    """Run real phyloP LRT/CONACC and return per-site signed scores (in order)."""
    cmd = [PHYLOP, "--method", "LRT", "--mode", "CONACC", "--wig-scores",
           "-i", "FASTA", mod_path, fasta_path]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    scores = []
    for line in out.splitlines():
        if line.startswith("fixedStep") or line.startswith("track") or not line.strip():
            continue
        scores.append(float(line))
    return np.array(scores)


def main():
    md = pp.parse_mod(MODEL)
    model_str = gtr_string(md)
    tmp = os.path.join(HERE, "sim_tmp")
    os.makedirs(tmp, exist_ok=True)
    print(f"model: {model_str}")
    print(f"{'T':>6} {'max_cons':>9} {'sens_%':>8} {'cons_recov':>11} {'FP_bg':>7}")

    rows = []
    for T in T_GRID:
        # scaled neutral model (what phyloP scores against)
        mod_T = os.path.join(tmp, f"neutral_T{T}.mod")
        write_scaled_model(md, T, mod_T)

        # trees: background at T, conserved block at RHO_CONS*T
        root = pp.parse_newick(md["tree_newick"])
        base = pp.tree_length(root)
        bg_tree = os.path.join(tmp, f"bg_T{T}.nwk")
        cons_tree = os.path.join(tmp, f"cons_T{T}.nwk")
        with open(bg_tree, "w") as f:
            f.write(scale_newick(md["tree_newick"], T / base) + ";\n")
        with open(cons_tree, "w") as f:
            f.write(scale_newick(md["tree_newick"], (RHO_CONS * T) / base) + ";\n")

        bg_fa = alisim(bg_tree, model_str, L_BG, os.path.join(tmp, f"bg_T{T}"), SEED)
        cons_fa = alisim(cons_tree, model_str, L_CONS, os.path.join(tmp, f"cons_T{T}"), SEED + 1)

        # concatenate columns per taxon: [background | conserved]
        bg = read_fasta(bg_fa)
        cons = read_fasta(cons_fa)
        combined = os.path.join(tmp, f"combined_T{T}.fa")
        with open(combined, "w") as o:
            for taxon in bg:
                o.write(f">{taxon}\n{bg[taxon]}{cons[taxon]}\n")

        scores = run_phylop_scores(mod_T, combined, tmp)
        if len(scores) != L_BG + L_CONS:
            print(f"  WARN T={T}: {len(scores)} scores != {L_BG + L_CONS} expected", file=sys.stderr)

        p = np.power(10.0, -np.abs(scores))          # raw p per site (both directions)
        conserved_class = scores > 0
        sig = bh_fdr(p, ALPHA)

        is_cons_block = np.zeros(len(scores), dtype=bool)
        is_cons_block[L_BG:L_BG + L_CONS] = True

        recovered = int(np.sum(sig & conserved_class & is_cons_block))
        false_pos = int(np.sum(sig & conserved_class & ~is_cons_block))
        sens = 100.0 * recovered / L_CONS
        max_cons = float(np.max(scores[conserved_class])) if conserved_class.any() else 0.0

        print(f"{T:>6.1f} {max_cons:>9.2f} {sens:>8.2f} {recovered:>11d} {false_pos:>7d}")
        rows.append({
            "T": T, "max_conserved_score": round(max_cons, 3),
            "sensitivity_pct": round(sens, 3),
            "conserved_recovered": recovered, "conserved_total": L_CONS,
            "false_positives_bg": false_pos, "bg_sites": L_BG,
            "rho_cons": RHO_CONS, "alpha": ALPHA,
        })

    out_csv = os.path.join(HERE, "sim_power_results.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out_csv}")


if __name__ == "__main__":
    main()
