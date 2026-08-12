#!/usr/bin/env python3
"""
build_notebook.py - assemble phylop-clustering-comparison.ipynb (tuned / phastCons-anchored
version) from the committed figures/metrics, and auto-render a self-contained HTML.

Figures are attached as cell outputs so it renders without re-execution; cluster_compare.py
(default-param metrics) and concordance.py (tuning + phastCons + tuned figures) are the
reproducible compute. Re-run after regenerating figures.
"""

import base64
import glob
import os
import subprocess

import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))


def md(src):
    return nbf.v4.new_markdown_cell(src)


def code(src, png=None):
    c = nbf.v4.new_code_cell(src)
    c.execution_count = None
    c.outputs = []
    if png:
        b64 = base64.b64encode(open(os.path.join(HERE, png), "rb").read()).decode()
        c.outputs.append(nbf.v4.new_output("display_data", data={"image/png": b64}, metadata={}))
    return c


cells = []

cells.append(md(r"""# Conserved-region clustering, tuned to phastCons

phyloP scores conservation base-by-base. At **Zoonomia (241-mammal) depth** those per-base calls
are effectively single-base resolution (see [../phylop-tree-length-power/](../phylop-tree-length-power/)).
Turning them into conserved **elements** is a separate clustering step — and at this depth it is the
*same* problem the pipeline's phastCons branch already faces: phastCons's own raw output here is
near-per-base (median 2 bp), and only becomes elements after a downstream merge + length filter. So
the real question isn't "phyloP vs phastCons," it's **how do you segment per-base conservation calls
into elements**, and the natural yardstick is **phastCons**.

This notebook runs four clustering methods on real data, **tunes each to best reproduce phastCons's
conserved elements**, and compares the tuned methods (and phastCons) head to head.

**Data:** a 2 Mb slice — `chr1:26,000,000–28,000,000` (human GRCh38) — of the real Zoonomia 241-mammal
phyloP output: 170,335 FDR-significant conserved sites, shifted to a 0-based mini-chromosome. No
alignment needed for the clustering; the phastCons reference *is* built from the alignment (below).

**Methods** (all take significant-site positions → regions; implemented in `lib/clustering.py`):

| method | what it does |
|---|---|
| **windowed** | fixed-size bins; bins with ≥ k sites are conserved; merge adjacent |
| **hmm** | 2-state online HMM over the per-position conserved/not stream (from the legacy script) |
| **hdbscan** | 1-D density clustering of site positions |
| **gap-merge (current)** | the pipeline's existing method: merge sites within a gap, then length-filter |
"""))

cells.append(code(
    "import lib.clustering as C\n"
    "sites = [(int(s), int(e)) for _, s, e in\n"
    "         (ln.split() for ln in open('data/chr1_26-28Mb.conserved-sites.bed'))]\n"
    "LENGTH = 2_000_000\n"
    "# each returns [(start, end, n_sites), ...]:\n"
    "win = C.windowed(sites, LENGTH, window_bp=20, min_sites_per_window=5)\n"
    "hmm = C.hmm(sites, LENGTH, t1_1=0.99, e1_1=0.5, min_len=20)\n"
    "hdb = C.hdbscan_cluster(sites, min_cluster_size=5)\n"
    "# full runs + figures: cluster_compare.py (defaults) and concordance.py (tuning + phastCons)\n",
    png=None,
))

cells.append(md(r"""## Why tuning is required first

At out-of-the-box parameters the four methods span a **~18× range in coverage** — from 2% to 42% of
the slice — which is almost entirely an artifact of arbitrary defaults, not method quality:

| method | coverage (default) | regions |
|---|---|---|
| hmm | **2.4%** | 1,022 |
| gap-merge (current) | 22.9% | 3,989 |
| windowed | 23.1% | 2,780 |
| hdbscan | **42.0%** | 2,802 |

You can't rank methods that are set to call wildly different amounts of sequence. So we first **fix
each method's parameters against a common target** — phastCons — and only then compare. (Full
default-parameter metrics and figures: `cluster_compare.py` → `cluster_metrics.csv`.)"""))

cells.append(md(r"""## Tuning to phastCons

**Building the reference (and a caveat we had to fix).** phastCons was run on the real 241-species
alignment of chr1:26–28 Mb (extracted from the 557 GB chromosome MAF) with the run's neutral model.
Its raw `--most-conserved` output is *near per-base* — 101k fragments, **median 2 bp** — because the
pipeline's `--rho`-only invocation sets no expected element length **and** this region of the 241-way
alignment is thousands of tiny blocks with unaligned gaps, so phastCons never sustains long runs
(even `--expected-length 45` only reached median 5 bp). The pipeline never uses that raw output as
elements — it **merges fragments within `cnee_ces_merge_gap_bp` (5 bp) and drops anything under
`cnee_min_len_bp` (50 bp)**. Doing exactly that yields a sensible element set — **1,254 elements,
median 101 bp, 9.1% coverage** (≈ Zoonomia's ~10% genome-wide constrained fraction) — which is the
reference used here.

**Objective — bp-level F1 vs the phastCons elements.** Two errors a method can make, in base pairs:

- **false positives (FP)** — bp the method calls conserved that are **not** in a phastCons element
  (over-calls). `precision = TP / (TP + FP)`.
- **false negatives (FN)** — phastCons element bp the method **misses**. `recall = TP / (TP + FN)`.

**F1** combines precision and recall into one score (0–1) that is only high when *both* are — so a
method can't win by over-calling (high recall, low precision) or under-calling. For scale, phastCons
covers **182 kb (9.1%)** of the 2 Mb slice. We grid-searched each method's parameters to maximize F1;
**the HMM stays the simple 2-state model** (we only pick better values for its 5 probabilities).

| method | default → tuned F1 | precision | recall | **false pos** | **false neg** | best params |
|---|---|---|---|---|---|---|
| **hmm** | 0.42 → **0.85** | 0.84 | 0.86 | **31 kb** | 25 kb | t11=.99, e11=.5, min_len=20 |
| gap-merge (current) | 0.56 → **0.82** | 0.74 | 0.93 | 61 kb | **12 kb** | gap=10, ≥10 sites, len 20 |
| windowed | 0.56 → **0.80** | 0.73 | 0.89 | 61 kb | 20 kb | 20 bp, ≥5 sites |
| hdbscan | 0.34 → **0.35** | 0.22 | 0.88 | **569 kb** | 22 kb | mcs=5 (fair grid: eom/leaf × mcs × min_samples) |

**Read the FP/FN columns across:** every method *misses* little (FN 12–25 kb — they all recover
phastCons). The whole difference is **over-calling**: the three good methods add only 30–60 kb of
false positives, while **hdbscan adds 569 kb — ~3× phastCons's entire 182 kb footprint** — which is
exactly why its precision (0.22) and F1 (0.35) collapse despite fine recall."""))

cells.append(code("from IPython.display import Image\nImage('figures/concordance_f1.png')\n",
                  png="figures/concordance_f1.png"))

cells.append(md(r"""## The tuned methods, head to head

At their phastCons-tuned settings, the four methods produce very different element sets — and three of
them land right where phastCons is:

| method | coverage | regions | median element |
|---|---|---|---|
| **phastCons (reference)** | 9.1% | 1,254 | 101 bp |
| hmm | 9.4% | 2,081 | 48 bp |
| gap-merge (current) | 11.5% | 1,898 | 75 bp |
| windowed | 11.2% | 2,903 | 20 bp |
| hdbscan | 36.5% | 6,764 | 57 bp |

**Track view** — a 60 kb window: raw significant sites on top, then each tuned element set. hmm,
gap-merge, and windowed track phastCons (gold); hdbscan (green) is a wall of over-coverage."""))

cells.append(code("Image('figures/tuned_track.png')\n", png="figures/tuned_track.png"))

cells.append(md(r"""**Element-length distributions** (tuned) vs phastCons (dashed). gap-merge sits closest to
phastCons's length profile; windowed pins to its 20 bp bin; hdbscan spreads long."""))

cells.append(code("Image('figures/tuned_sizes.png')\n", png="figures/tuned_sizes.png"))

cells.append(md(r"""**Agreement** — Jaccard of the base pairs each tuned method calls conserved. hmm agrees
most with phastCons (0.74); gap-merge (0.70) and windowed (0.67) close behind; the three agree with
each other ~0.74. hdbscan agrees with nothing (~0.21)."""))

cells.append(code("Image('figures/tuned_jaccard.png')\n", png="figures/tuned_jaccard.png"))

cells.append(md(r"""## Observations

- **All four methods have low false negatives (12–25 kb).** They *all* recover phastCons's elements —
  none of them meaningfully *misses* conservation. So the whole comparison comes down to **false
  positives** (over-calling), which is what the F1 differences track.
- **The simple HMM comes out best (F1 0.85, Jaccard 0.74 with phastCons)** — only **31 kb** of false
  positives, precision 0.84, recall 0.86 — while staying the fast 2-state model
  (`t11=0.99, e11=0.5, min_len=20`). A cheap HMM reproduces phastCons's *elements* well; that's its
  reason for existing.
- **gap-merge and windowed are close** (F1 0.82 / 0.80; ~60 kb false positives). The three
  coherent-region methods agree with phastCons and each other (mutual Jaccard ~0.74) and land at its
  ~9–12% coverage.
- **hdbscan fails on false positives (F1 0.35) — and it's not a tuning artifact.** We gave it a fair
  grid (`eom` *and* `leaf` selection, a range of `min_cluster_size`/`min_samples`); its best is still
  **569 kb of false positives** (3× phastCons's footprint), precision 0.22. Density clustering links
  sites across gaps into over-large regions, and no setting we tried fixes it — it's the wrong tool
  for compact conserved elements.
- **Tuning is decisive.** The HMM went 0.42 → 0.85; the raw default comparison (2% vs 42% coverage)
  was almost pure parameter artifact.

**Caveats.** One 2 Mb region of one chromosome; a single phastCons setting; bp-level (not
element-boundary) matching; parameters tuned on the same region they are scored on. And the F1s are
sensitive to how the phastCons reference is defined (raw fragments vs merged elements) — so trust the
**ranking (hmm ≈ gap-merge ≈ windowed ≫ hdbscan)** more than the exact values. A fuller pass would
tune on one region and evaluate on a held-out one.

**Bottom line.** If a clustering step is wired into `phylop_regions.smk` (gated, large-tree-only), the
**HMM is the natural candidate** — best concordance, simple, fast — with gap-merge/windowed as
close, even simpler alternatives; hdbscan is not worth carrying."""))

cells.append(md(r"""## Reproducibility

Run in the `phyloacc-workflows` env (numpy, matplotlib, hdbscan):

```bash
cd analyses/phylop-clustering-comparison
python cluster_compare.py     # default-param metrics + figures -> cluster_metrics.csv
python concordance.py         # phastCons reference, tuning, tuned figures -> concordance_results.csv
python build_notebook.py      # rebuilds this notebook + HTML
```

| file | what |
|---|---|
| `data/chr1_26-28Mb.conserved-sites.bed` | the 2 Mb slice: real Zoonomia chr1 conserved sites, 0-based |
| `data/phastcons_reference.0based.bed` | phastCons elements (merged/filtered) — the tuning reference |
| `data/phastcons_raw_fragments.0based.bed` | raw phastCons `--most-conserved` fragments (provenance) |
| `lib/clustering.py` (repo root) | windowed / hmm / hdbscan_cluster |
| `cluster_compare.py` | default-parameter comparison |
| `concordance.py` | tuning + phastCons concordance + tuned figures |

**Provenance.** Conserved sites: `data/mammals-v3/echolocation-v3/03-phylop/sites/autosomes/chr1-...conserved.bed`,
window `chr1:26,000,000–28,000,000`, shifted −26,000,000. phastCons reference: the same window's
241-species alignment was extracted from the 557 GB `01-maf-prep/.../chr1.maf` (by byte-range via a
partial `mafutils` index), scored with `phastCons --rho 0.3 --most-conserved` on the run's neutral
model, then merged within 5 bp + filtered ≥50 bp (the pipeline's `cnee_ces_merge_gap_bp` /
`cnee_min_len_bp`) and shifted to the 0-based slice frame."""))

nb = nbf.v4.new_notebook()
nb.cells = cells
nb.metadata = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
               "language_info": {"name": "python"}}
out = os.path.join(HERE, "phylop-clustering-comparison.ipynb")
with open(out, "w") as f:
    nbf.write(nb, f)
print(f"wrote {out} ({len(cells)} cells)")

candidates = ["/n/home07/gthomas/miniconda3/envs/mafutils-dev/bin/jupyter"]
candidates += [j for j in sorted(glob.glob("/n/home07/gthomas/miniconda3/envs/*/bin/jupyter"))
               if j not in candidates]
for jupyter in candidates:
    if not os.path.exists(jupyter):
        continue
    try:
        subprocess.run([jupyter, "nbconvert", "--to", "html", "--embed-images", out],
                       check=True, capture_output=True, text=True, timeout=300)
        print(f"wrote {out[:-6]}.html  (via {jupyter})")
        break
    except Exception:
        continue
else:
    print("WARNING: no working nbconvert found; HTML not regenerated.")
