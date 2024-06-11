#############################################################################
# Snakemake rules to classify windows in a MAF by conservation score and
# alignment depth
# Gregg Thomas, October 2023
#############################################################################

# snakemake -p -s 02_conservation_windows.smk --configfile echolocation-cfg.yaml --profile profiles/slurm_profile/ --use-conda --dryrun --rulegraph | dot -Tpng > dag.png

#############################################################################

import os
import sys

#############################################################################

PROJ_DIR = config["project_dir"];
# The root directory of the project
# All files will be saved relative to this directory

PREFIX = config["prefix"];
# The prefix for all output files

CHROMES = config["chromes"];
# The list of chromosomes to process

ALPHA = str(config["alpha"]);
# The alpha level to consider a site to be conserved

MAF_OUTDIR = config["mafoutdir"];
# The output directory for all MAF-related files

TMPDIR = config["tmp_dir"];
# A temporary directory for intermediate files

## Project and run information
#####

MAF_DIR = config["mafdir"];
MAF_FILE = config["maf"];
MAF_PATH = os.path.join(MAF_DIR, MAF_FILE);
# The MAF file to process

MAF_REF_ID = config["maf_ref_id"];
MAF_CHR_PREFIX = config["maf_chr_prefix"];
MAF_REF_CHR_JOINER = config["maf_ref_chr_joiner"];
MAF_REF_PREFIX = MAF_REF_ID + MAF_REF_CHR_JOINER + MAF_CHR_PREFIX;
# The prefix for the reference chromosome in the MAF file

MAF_SPLIT_CHR_DIR = os.path.join(MAF_DIR, PREFIX + "-mafSplit");
# Info for splitting the MAF file into chromosomes

## MAF 
#####

REF_INDEX = config["ref_genome_index"];

#############################################################################

rule all:
    input:

#############################################################################

# rule run_phylop
## TODO

####################

rule split_maf_chr:
    input:
        maf = MAF_PATH,
        split_bed = REF_CHR_BED
    output:
        outdir = directory(MAF_SPLIT_CHR_DIR),
        split_maf = expand(os.path.join(MAF_SPLIT_CHR_DIR, "{chrome}.00.maf"), chrome=chr_list)
    params:
        prefix = MAF_SPLIT_CHR_DIR
    conda:
        "envs/ucsc-mafsplit.yml"
        # mafSplit needs its own environment right now because the current version (377) requires an outdated
        # openssl library (1.0.0) that doesn't play nice with other tools
    resources:
        partition="intermediate",
        mem="48g",
        time="120:00:00"
    shell:
        """
        mkdir -p {output.outdir}
        mafSplit {input.split_bed} {params.prefix} {input.maf}
        """