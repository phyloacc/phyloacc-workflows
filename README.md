# phyloacc-workflows

Snakemake workflows for preparing comparative genomic alignment inputs for
[PhyloAcc](https://github.com/phyloacc/phyloacc-workflows), from a whole-genome
alignment (MAF) through neutral model estimation (phyloFit), site- and
region-level conservation scoring (phyloP, phastCons), and CNEE extraction.

A full tutorial is forthcoming and will be linked here.

## Requirements

- `conda` or `mamba`
- (optional) a Slurm cluster, via `snakemake-executor-plugin-slurm`

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
(`output_dir`, `maf`, `tree_file`, `ref_fasta`, `ref_chromosome_groups`, etc.).
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
