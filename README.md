# phyloacc-workflows

Snakemake workflows for preparing comparative genomic alignment inputs ( *e.g.* for
[PhyloAcc](https://github.com/phyloacc/PhyloAcc)), from a whole-genome
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
conserved sites by construction (not a bug). Prefer phastCons there. A pre-flight gate
(`phylop_power_check`, controlled by `phylop_power_*` config keys) detects this from the
fitted neutral model and stops the phyloP stage early with an explanation rather than
scanning and returning nothing; set `phylop_power_override: true` to run anyway. See
`analyses/phylop-tree-length-power/` for the analysis.

Both conservation stages feed a shared CNEE-building stage (`build_cnees: true`,
`workflow/cnees.smk`): conserved elements from phastCons and/or phyloP are turned into
CNEEs (drop CDS-overlapping elements, length-filter, extract alignments), written under
`05-cnees/{phastcons,phylop}/`, one CNEE set per enabled source.

### Align-only mode

If you already have your elements in FASTA formatand only need them aligned, set `align_input_dir` (instead of
`maf`) and the whole MAF / phyloFit / phastCons / phyloP / CNEE pipeline is bypassed. Each
element is just aligned with `mafft`, run as a SLURM scatter-gather (or locally).

**Required inputs**

- `align_input_dir`: a directory of **unaligned per-element FASTAs** (one multi-FASTA per
  element, one record per sample/species). Files may sit flat in the directory, or one level
  down in `batch<N>/` subdirs (the layout some CNEE extractors write); both are discovered
  automatically, and the output mirrors whichever layout you give it.
- `output_dir` and `tmp_dir`: where results and temporary files go.
- `mafft` on the environment's `PATH` (installed by `setup`; swap in another aligner via
  `aligner_commands`).

**Example config**: generate a starter with `./phyloacc_workflows init --mode align`, or
write one directly.

```yaml
align_input_dir: /path/to/unaligned-elements/   # dir of *.fa / *.fasta / *.fna (or batchN/ subdirs)
output_dir:      /path/to/output/
tmp_dir:         /path/to/tmp/

# optional (defaults shown):
align_min_len_bp:   50     # drop elements whose longest record is < 50 bp (0 = keep all)
align_min_seqs:     2      # drop elements with fewer than this many sequences
align_header_split: " "    # trim each aligned header to one field: ">gGal ce1 chr:1-9" -> ">gGal" (blank = keep as-is)
align_header_field: 0      # which split field to keep (0-based)
align_batch_size:   2000   # elements per SLURM job
```

See `config-template-align.yaml` for every option.

**Run**

```bash
./phyloacc_workflows --env-name <env> run -j <N> -e slurm --configfile my-align.yaml
# drop "-e slurm" to run locally
```

**What it produces** under `<output_dir>/aligned-elements/`:

- `aligned/…`: one aligned FASTA per element, mirroring the input layout (flat, or the same
  `batch<N>/` subdirs)
- `manifest.txt`: the aligned elements, as paths relative to the input directory
- `skipped.tsv`: every dropped element with a reason (`too-short`, `too-few-seqs`)