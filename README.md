# phyloacc-workflows

Snakemake workflows for preparing comparative genomic alignment inputs for
[PhyloAcc](https://github.com/phyloacc/phyloacc-workflows), from a whole-genome
alignment (MAF) through neutral model estimation (phyloFit), site- and
region-level conservation scoring (phyloP, phastCons), and CNEE extraction. An
[align-only mode](#align-only-mode) skips all of that and just aligns a directory of
unaligned elements you already have.

A tutorial exists [on the PhyloAcc website](https://phyloacc.github.io/workflow.html).

## Requirements

- `conda` or `mamba`
- For most datasets, a computing cluster with the SLURM job scheduler. Other clusters
may be supported with minimal effort via Snakemake executor plugins, but remain untested.

## Setup

This repo contains a `phyloacc_workflows` wrapper script that manages its own
conda environment (defined in `envs/environment.yml`) and runs Snakemake.

```bash
./phyloacc_workflows setup      # create/update the conda environment
./phyloacc_workflows check      # verify the environment is ready, read-only
./phyloacc_workflows help       # full usage, including --env-name/--verbose
```

## Configuration

Generate a starter config with the wrapper's `init` subcommand, which strips
`config-template.yaml` down to keys and defaults (required fields left blank for you):

```bash
./phyloacc_workflows init -o my-config.yaml                 # full MAF pipeline
./phyloacc_workflows init --mode phylop -o my-config.yaml   # or: phastcons | both | align
```

`--mode` tailors the output to one workflow fork (emitting only that fork's relevant
options with the `run_*` switches preset); `--mode align` uses the dedicated
`config-template-align.yaml`. Or copy `config-template.yaml` directly - it is the
fully-documented reference covering every option in every fork.

Fill in the required inputs at the top (`output_dir`, `maf`, `tree_file`,
`ref_chromosome_groups`, etc.). `ref_fasta` is only required if `split_strategy` is `ns`
or `fixed_windows` (see below) - the default `num_seqs` strategy needs no reference FASTA
at all. Everything below has working defaults but can be adjusted - see the comments in
the template for what each option controls. Example configs for past runs are in `cfgs/`.

### Resource requirements

The `rule_resources` defaults in `config-template.yaml` were set for our own
datasets and cluster and will not fit every genome/alignment size - review
them rather than trusting them blindly, especially the rules with the
largest per-job footprints:

- `maf_split_chr_by_group` (240 GB by default) and `extract_4d_codons_by_chr`
  (128 GB) scale with whole-alignment/whole-chromosome size, since they load
  a full MAF or codon set at once.
- `run_phylop` (64 GB) and `adjust_pvals` (32 cpus) scale with the number of
  sites/species being processed.
- Most other rules operate on small per-chunk/per-region inputs and default
  to 4 GB / 1 cpu / about an hour, which is unlikely to need much adjustment
  regardless of dataset size.

If a job fails from an out-of-memory (OOM) kill or a timeout, that rule's entry is
what to raise - there's no need to inflate every rule's resources up front.

## Running

Dry-run first to sanity-check the plan before submitting anything for real:

```bash
./phyloacc_workflows run --configfile my-config.yaml -j <N> -e slurm --dryrun
```

Once the job list and target files look right, remove `--dryrun` to actually
execute. `run` defaults to this repo's `Snakefile` automatically, so no `-s`
is needed unless you want a different one.

Which stages run is controlled by `run_phylofit` / `run_phylop` /
`run_phastcons` / `build_cnees` in your config file.

The `run_phylop` stage (per-base phyloP scoring → FDR → conserved/accelerated sites →
region clustering (method-selectable via `phylop_cluster_method`); see the "phyloP-only behavior" section of the config template
for its rules and outputs) has a hard statistical-power limit on shallow trees: it needs
a large total neutral tree length (roughly > ~10 substitutions/site) before any single
conserved site can clear genome-wide FDR. On shallow-tree datasets it returns few/no
conserved sites by construction (not a bug) — prefer phastCons there. A pre-flight gate
(`phylop_power_check`, controlled by `phylop_power_*` config keys) detects this from the
fitted neutral model and stops the phyloP stage early with an explanation rather than
scanning and returning nothing; set `phylop_power_override: true` to run anyway. See
`analyses/phylop-tree-length-power/` for the analysis.

Both conservation stages feed a shared CNEE-building stage (`build_cnees: true`,
`workflow/cnees.smk`): conserved elements from phastCons and/or phyloP are turned into
CNEEs (drop CDS-overlapping elements, length-filter, extract alignments), written under
`05-cnees/{phastcons,phylop}/` — one CNEE set per enabled source.

### Align-only mode

If you already have your elements and only need them aligned, set `align_input_dir`
(instead of `maf`) to a directory of unaligned per-element FASTAs. This bypasses the whole
MAF / phyloFit / phastCons / phyloP / CNEE pipeline and simply aligns each element with
`mafft` (the aligner is configurable via `aligner_commands`), writing per-element aligned
FASTAs plus `manifest.txt` and `skipped.tsv` under `<output_dir>/aligned-elements/`. Input
in `batch<N>/` subdirs is discovered and mirrored in the output; elements are dropped
(with a reason in `skipped.tsv`) below `align_min_seqs` sequences or `align_min_len_bp` bp,
and headers can be trimmed to a chosen field. It runs as a SLURM scatter-gather (or locally
without `-e slurm`). Generate a config with `./phyloacc_workflows init --mode align`;
requires `mafft` in the environment (installed by `setup`).

## Repo layout

- `Snakefile` - entry point; includes the workflow files below based on config toggles
- `workflow/` - the active Snakemake rule files (`phylofit_models.smk`, `phylop_regions.smk`, `phastcons_cnees.smk`, `cnees.smk` - the shared, source-agnostic CNEE-building stage that turns conserved elements from phastCons and/or phyloP into CNEEs - and `align_elements.smk` - the align-only mode); `workflow/legacy/` holds earlier versions of the pipeline not used by the current `Snakefile`
- `lib/` - shared Python helpers used across the workflow files
- `utils/` - standalone scripts/tools invoked by rules
- `envs/environment.yml` - pinned dependencies for the `phyloacc_workflows` environment
- `cfgs/` - example filled-in configs from past analyses

## Running tests (developers)

`tests/` holds three tiers of tests. These are a developer/maintainer concern, not
something needed to run the pipeline itself, so `pytest` is not part of
`envs/environment.yml`. To run them, install pytest into your existing environment and
run it from the repo root:

```bash
mamba install -n phyloacc-workflows pytest    # or: pip install pytest
pytest tests/
```

- **Unit tests** (`test_intervals.py`, `test_parsing.py`) - pure-Python logic in `lib/`
  (interval merging, BED/GFF parsing, Newick tip extraction, etc. - the same code the real
  Snakemake rules import and run). No external tools needed. ~1-2s.
- **DAG/config validation** (`test_dag_validation.py`) - each test invokes a real
  `snakemake -n` subprocess against a minimal config to check that bad config values are
  rejected with the expected error. Needs `snakemake` on `PATH`. ~3-4 min.
- **Integration** (`tests/integration/`) - runs the real pipeline (real
  `mafutils`/`phyloFit`/`phastCons`, not mocked) end-to-end against a small real-data
  fixture; see `tests/integration/README.md` for the fixture's provenance. Needs
  `mafutils`/`phyloFit`/`phastCons` actually on `PATH` - skips cleanly (not a failure) if
  they're missing. ~4-5 min. Includes a soft "silver standard" check that warns (does not
  fail) if a run's summary numbers drift far from a committed reference - see that
  directory's README for how to regenerate it if the fixture or tool versions change.

Running everything together takes roughly 7-8 minutes, dominated by the last two tiers'
real subprocess/tool invocations.

## Releasing (maintainers)

Version bumps are manual, not automated. To cut a release:

1. On a branch, bump `version` (and the `releasedate-*`/`latest-commit-date`
   fields) in `lib/info.yaml`, and `version`/`date-released` in
   `CITATION.cff`. Open a PR and merge it to `main`.
2. From `main`, tag and push:
   ```bash
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```
3. `.github/workflows/release.yml` creates the GitHub Release from that tag
   automatically - it does not check that the tag matches `lib/info.yaml`,
   so a mismatched tag publishes a Release with the wrong version recorded
   in-repo without any warning.

Activate the repo's pre-push hook once per clone to guard against that:

```bash
git config core.hooksPath .githooks
```

With it active, pushing a `vX.Y.Z` tag whose `lib/info.yaml` version doesn't
match is blocked locally with a clear error, before the push happens - this
is the only check in place, so activating it is strongly recommended. It's
optional, dev-only tooling (unrelated to `phyloacc_workflows setup`, which
sets up the environment for *running* the pipeline).
