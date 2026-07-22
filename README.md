# phyloacc-workflows

Snakemake workflows for preparing comparative genomic alignment inputs for
[PhyloAcc](https://github.com/phyloacc/phyloacc-workflows), from a whole-genome
alignment (MAF) through neutral model estimation (phyloFit), site- and
region-level conservation scoring (phyloP, phastCons), and CNEE extraction.

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

Copy `config-template.yaml` and fill in the required inputs at the top
(`output_dir`, `maf`, `tree_file`, `ref_chromosome_groups`, etc.). `ref_fasta` is
only required if `split_strategy` is `ns` or `fixed_windows` (see below) - the
default `num_seqs` strategy needs no reference FASTA at all.
Everything below that section has working defaults but can be adjusted -
see the comments in the template for what each option controls. Example
configs for past runs are in `cfgs/` for reference.

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

## Repo layout

- `Snakefile` - entry point; includes the workflow files below based on config toggles
- `workflow/` - the active Snakemake rule files (`phylofit_models.smk`, `phylop_regions.smk`, `phastcons_cnees.smk`); `workflow/legacy/` holds earlier versions of the pipeline not used by the current `Snakefile`
- `lib/` - shared Python helpers used across the workflow files
- `utils/` - standalone scripts/tools invoked by rules
- `envs/environment.yml` - pinned dependencies for the `phyloacc_workflows` environment
- `cfgs/` - example filled-in configs from past analyses

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
