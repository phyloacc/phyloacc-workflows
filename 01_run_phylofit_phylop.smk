#############################################################################
# Snakemake rule to extract 4d sites from a MAF
# Gregg Thomas, September 2023
#############################################################################

# snakemake -j 5 -s <snakefile> --configfile <config file> --dryrun

#############################################################################

import os
import sys
import glob
import logging
import traceback

import lib.common as COMMON
from lib.common import spacedOut as SO

from functools import partial

#############################################################################
# System setup

config_flag = config.get("display", False);
version_flag = config.get("version", False);
info_flag = config.get("info", False);
debug = config.get("debug", False);
#debug = True;
# A hacky way to get some custom command line arguments for the pipeline
# These just control preprocessing flags that stop the pipeline early anyways

MAIN, DRY_RUN, OUTPUT_DIR, LOG_DIR, TMPDIR, LOG_LEVEL, LOG_VERBOSITY = COMMON.pipelineSetup(config, sys.argv, version_flag, info_flag, config_flag, debug, workflow);
# Setup the pipeline, including the output directory, log directory, and tmp directory

MLOG = logging.getLogger('META')
# Setup logging if debugging

PIPELINE_DIR = os.path.dirname(workflow.snakefile);
# The directory where the snakemake file is located

ENVS_DIR = os.path.join(PIPELINE_DIR, "envs");
# The directory where the conda environments are located

UTILS_DIR = os.path.join(PIPELINE_DIR, "utils");
# The directory where the utility scripts are located

PREFIX = config["prefix"];
# The prefix for all output files

getRuleResources = partial(COMMON.getResources, config)
# This maps the function to get rule resources from the config file
# so we don't have to pass config each time we call it

#############################################################################
# Input files and output paths

SAMPLE_FILE = config["sample_file"];
sample_file_extension = os.path.splitext(SAMPLE_FILE)[1];
gc_sample_path = SAMPLE_FILE.replace(sample_file_extension, "-gc" + sample_file_extension);
avg_gc_path = SAMPLE_FILE.replace(sample_file_extension, "-avg-gc" + sample_file_extension);

gc_sample_basename = os.path.basename(gc_sample_path)
avg_gc_basename = os.path.basename(avg_gc_path)

GC_SAMPLE_FILE = os.path.join(OUTPUT_DIR, "summary-data", gc_sample_basename);
AVG_GC_FILE = os.path.join(OUTPUT_DIR, "summary-data", avg_gc_basename);
# Sample file setup, along with new GC files

MAF_PATH = config["maf"];
MAF_FILE = os.path.basename(MAF_PATH);
MAF_DIR = os.path.dirname(MAF_PATH);
# MAF info

MAF_INDEX_BLOCK = os.path.join(MAF_DIR, MAF_FILE + ".block.idx");
MAF_INDEX_SCAFF = os.path.join(MAF_DIR, MAF_FILE + ".scaffold.idx");
# This will be created

MAF_REF_ID = config["maf_ref_id"];
MAF_CHR_PREFIX = config["maf_chr_prefix"];
MAF_REF_CHR_JOINER = config["maf_ref_chr_joiner"];
MAF_REF_PREFIX = MAF_REF_ID + MAF_REF_CHR_JOINER + MAF_CHR_PREFIX;
# MAF reference and scaffold id parsing


if config["maf_split_dir"]:
    MAF_SPLIT_DIR = config["maf_split_dir"];
else:
    MAF_SPLIT_DIR = os.path.join(OUTPUT_DIR, "mafs");
# MAF split directory, where the MAF files will be split by chromosome

# REF_CHR_BED_DIR = os.path.join(OUTPUT_DIR, "beds");
# REF_CHR_GFF_DIR = os.path.join(OUTPUT_DIR, "gffs");
# MAF_SPLIT_CHR_DIR = os.path.join(OUTPUT_DIR, PREFIX + "-mafSplit"); # Maybe this should be in a dir specified by the user?
# Various output sub-directories

REF_INDEX = config["ref_genome_index"];
REF_GFF_PATH = config["ref_gff"];
REF_GFF_FILE = os.path.basename(REF_GFF_PATH);
REF_CHROMOSOME_GROUPS = config["ref_chromosome_groups"];
# Reference genome info

#############################################################################
# Basic tree parsing

TREE_FILE = config["tree_file"];
species_tree = open(TREE_FILE, "r").read().strip();
topology = re.sub(r'[)][\d\w<>/.eE_:-]+', ')', species_tree);
topology = re.sub(r':[\d.eE-]+', '', topology);
## Remove the branch lengths and node labels from the input tree string

SPECIES = [ tip_label for tip_label in topology.replace("(","").replace(")","").replace(";","").split(",") ];
## Get the species names from the tree string

#############################################################################
# Other params

SEQ_THRESHOLD_4D = config["filter_threshold_4d"]
# The threshold for filtering 4d sites

ALPHA = str(config["alpha"]);
# The alpha level to consider a site to be conserved

MIN_CLUSTER_SIZES = config["min_cluster_sizes"];
MIN_SAMPLES = config["min_samples"];
# Cluster benchmarking

#############################################################################



# for GROUP in REF_CHROMOSOME_GROUPS:
#     print(GROUP);
#     with open(os.path.join(REF_CHR_BED_DIR, GROUP + ".bed"), "w") as f:
#         for line in open(REF_INDEX, "r"):
#             print(line);
#             chrome, length = line.strip().split("\t")[:2];
#             print(chrome, length);
#             if chrome in REF_CHROMOSOME_GROUPS[GROUP]:
#                 print(chrome);
#                 length = int(length)-1;
#                 f.write(f"{chrome}\t0\t{length}\n");
        
        # for chrome in GROUP:
        #     f.write(f"{chrome}\t0\t{length}\n");
# This gets the chromosome start and end positions from the reference genome index file
# for splitting the MAF

#CHROMES = [MAF_REF + "." + c for c in CHROMES];

flattened_chromosome_groups = [(group, chromosome) for group, chromosome_list in REF_CHROMOSOME_GROUPS.items() for chromosome in chromosome_list];

REF_CHR_GROUPS_LIST, REF_CHROMOSOMES = zip(*flattened_chromosome_groups);
REF_CHR_PATHS = [ REF_CHR_GROUPS_LIST[i] + "/" + REF_CHROMOSOMES[i] for i in range(len(REF_CHROMOSOMES)) ];
# Get two lists of equal length to use for wild cards

####################

if LOG_LEVEL == "debug":
    for i in range(len(REF_CHR_GROUPS_LIST)):
        MLOG.debug(SO(f"CHR GROUP {i}", 15) + f"{REF_CHR_GROUPS_LIST[i]} : {REF_CHROMOSOMES[i]}");
    MLOG.debug("EXITING BEFORE RULES. DEBUG MODE.");
    sys.exit(0);
# Exit before running rules if in debug mode

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

# Helper function to generate combinations
def get_combinations():
    combinations = []
    for chromosome_group, ref_chromosomes in REF_CHROMOSOME_GROUPS.items():
        for ref_chromosome in ref_chromosomes:
            for min_cluster_size in MIN_CLUSTER_SIZES:
                for min_samples in MIN_SAMPLES:
                    combinations.append((chromosome_group, ref_chromosome, min_cluster_size, min_samples))
    return combinations

# Helper function to generate combinations
def generate_combinations():
    combinations = []
    for chromosome_group, ref_chromosomes in REF_CHROMOSOME_GROUPS.items():
        for ref_chromosome in ref_chromosomes:
            for min_cluster_size in MIN_CLUSTER_SIZES:
                for min_samples in MIN_SAMPLES:
                    combinations.append({
                        "chromosome_group": chromosome_group,
                        "ref_chromosome": ref_chromosome,
                        "min_cluster_size": min_cluster_size,
                        "min_samples": min_samples
                    })
    return combinations


# Group A: BED (or manifest) targets
# bed_targets = expand(
#     os.path.join(REF_CHR_BED_DIR, "{group}.bed"),
#     group=list(REF_CHROMOSOME_GROUPS.keys())
# )

# Group B: 4d-codon targets
# codon_targets = expand(
#     os.path.join(OUTPUT_DIR, "01-4d-codons", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-codons.ss"),
#     zip,
#     chromosome_group=REF_CHR_GROUPS_LIST,
#     ref_chromosome=REF_CHROMOSOMES
# )

localrules: all

rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "summary-data", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.conserved-site-counts." + ALPHA + ".tsv"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),
        expand(os.path.join(OUTPUT_DIR, "summary-data", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.accelerated-site-counts." + ALPHA + ".tsv"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),
        
        #expand(os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + ALPHA + ".conserved.bed"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),
        #expand(os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + ALPHA + ".bed"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),
        
        #expand(os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.bed"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),
        #expand(os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.wig"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),

        #expand(os.path.join(OUTPUT_DIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-corrected.mod"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)
        # expand(os.path.join(OUTPUT_DIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.mod"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),

        # gc_sample_file = GC_SAMPLE_FILE,
        # avg_gc_file = AVG_GC_FILE

        # sites_filtered = expand(os.path.join(OUTPUT_DIR, "02-4d-sites", "{chromosome_group}-filtered", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + ".ss"), 
        #     zip,
        #     chromosome_group=REF_CHR_GROUPS_LIST,
        #     ref_chromosome=REF_CHROMOSOMES),
        # summary_file = expand(os.path.join(OUTPUT_DIR, "summary-data", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + "-summary.tsv"),
        #     zip,
        #     chromosome_group=REF_CHR_GROUPS_LIST,
        #     ref_chromosome=REF_CHROMOSOMES),

        #expand(os.path.join(OUTPUT_DIR, "02-4d-sites", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites.ss"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),
        #expand(os.path.join(REF_CHR_BED_DIR, "{chromosome_group}.bed"), chromosome_group=REF_CHR_GROUPS_LIST),
        #expand(os.path.join(OUTPUT_DIR, "gffs", "{chromosome_group}", REF_GFF_FILE.replace(".gff", ".{ref_chromosome}.gff")), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),
        
        #all_codons_output = expand(OUTPUT_DIR, "01-4d-codons", "{ref_chr_path}-4d-codons.ss", ref_chr_path=REF_CHR_PATHS)
        
        #expand(os.path.join(OUTPUT_DIR, "01-4d-codons", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-codons.ss"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),

        #expand(os.path.join(OUTPUT_DIR, "mafs", "{chromosome_group}"), chromosome_group=REF_CHR_GROUPS_LIST),
        #expand(os.path.join(OUTPUT_DIR, "mafs", "{chromosome_group}", "{ref_chromosome}.maf"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),
        #expand(os.path.join(REF_CHR_BED_DIR, "{chromosome_group}.bed"), chromosome_group=REF_CHROMOSOME_GROUPS.keys()),


#############################################################################
# Pipeline rules

rule maf_index:
    input:
        maf = MAF_PATH
    output:
        maf_index_block = MAF_INDEX_BLOCK,
        maf_index_scaff = MAF_INDEX_SCAFF
    params:
        script_path = os.path.join(UTILS_DIR, "maf_index.py")
    log:
        job_log = os.path.join(LOG_DIR, "maf-index.log")
    resources:
        **getRuleResources("maf_index")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                cmd = [ "python", params.script_path, input.maf, output.maf_index_block, output.maf_index_scaff ];
                COMMON.runCommand(cmd, log_stream, log_stream, "maf_index");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise


####################

rule make_group_beds:
    input:
        maf = MAF_PATH
    output:
        chr_group_bed = os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}.bed")
    params:
        ref_chroms = lambda wildcards: REF_CHROMOSOME_GROUPS[wildcards.chromosome_group],
        chr_prefix = MAF_CHR_PREFIX,
        script_path = os.path.join(UTILS_DIR, "make_group_beds.py")
    log:
        job_log = os.path.join(LOG_DIR, "make-maf-split-beds.{chromosome_group}.log")
    resources:
        **getRuleResources("make_group_beds")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                cmd = [ "python", params.script_path,
                        REF_INDEX,
                        params.chr_prefix,
                        output.chr_group_bed,
                        *[str(c) for c in params.ref_chroms] ];

                COMMON.runCommand(cmd, log_stream, log_stream, "make_group_beds", wc=f"{wildcards.chromosome_group}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise
    
    # run:
    #     with open(log.job_log, "w") as log_stream:
    #         try:
    #             with open(os.path.join(REF_CHR_BED_DIR, wildcards.chromosome_group + ".bed"), "w") as f:
    #                 for line in open(REF_INDEX, "r"):
    #                     chrome, length = line.strip().split("\t")[:2];
    #                     if chrome in REF_CHROMOSOME_GROUPS[wildcards.chromosome_group]:
    #                         length = int(length)-1;
    #                         f.write(f"{chrome}\t0\t{length}\t{chrome}\n");
    #         except Exception as e:
    #             traceback.print_exc(file=log_stream)
    #             raise
# This rule generates a bed file for each chromosome group in the reference genome

####################

checkpoint maf_split_chr_by_group:
    input:
        maf = MAF_PATH,
        maf_index_scaff = MAF_INDEX_SCAFF,
        chr_group_bed = os.path.join(OUTPUT_DIR, "beds", "{chromosome_group}.bed")
    output:
        maf_manifest = os.path.join(MAF_SPLIT_DIR, "{chromosome_group}", "manifest.txt")
        #chr_group_dir = directory(os.path.join(OUTPUT_DIR, "mafs", "{chromosome_group}")),
        #chr_maf = os.path.join(OUTPUT_DIR, "mafs", "{chromosome_group}", "{ref_chromosome}.maf")
    params:
        outdir = os.path.join(MAF_SPLIT_DIR, "{chromosome_group}"),
        script_path = os.path.join(UTILS_DIR, "maf_fetch.py"),
        rule_name = "maf_split_by_chr"
    log:
        job_log = os.path.join(LOG_DIR, "maf-split-by-chr.{chromosome_group}.log")
    resources:
        **getRuleResources("maf_split_chr_by_group")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                cmd = [ "python", params.script_path,
                        input.maf,
                        input.maf_index_scaff,
                        input.chr_group_bed,
                        "-o", params.outdir,
                        "-p", str(resources.cpus_per_task),
                        "-m", "scaffold" ];
                        
                COMMON.runCommand(cmd, log_stream, log_stream, params.rule_name, wc=f"{wildcards.chromosome_group}");
                COMMON.writeBedManifest(input.chr_group_bed, output.maf_manifest);
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise

####################

rule ref_gff_split_by_chr:
    input:
        gff = REF_GFF_PATH
    output:
        chromosome_gff = os.path.join(OUTPUT_DIR, "gffs", "{chromosome_group}", REF_GFF_FILE.replace(".gff", ".{ref_chromosome}.gff"))
    params:
        ref_chr = lambda wildcards: wildcards.ref_chromosome,
        prefix = MAF_REF_PREFIX,
        script_path = os.path.join(UTILS_DIR, "ref_gff_split_by_chr.awk")
    log:
        job_log = os.path.join(LOG_DIR, "ref-gff-split-by-chr.{chromosome_group}.{ref_chromosome}.log")
    resources:
        **getRuleResources("ref_gff_split_by_chr")
    run:
        with open(log.job_log, "w") as log_stream, open(output.chromosome_gff, "w") as out_stream:
            try:
                cmd = [ "awk",
                        "-v", f"chr={params.ref_chr}",
                        "-v", f"prefix={params.prefix}",
                        "-f", params.script_path,
                        input.gff ]

                COMMON.runCommand(cmd, log_stream, out_stream, "ref_gff_split_by_chr", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise
    # run:
    #     with open(input.gff, 'r') as fin, open(output.chromosome_gff, 'w') as fout:
    #         for line in fin:
    #             fields = line.split('\t')
    #             if not fields:
    #                 continue

    #             first_field = fields[0];
                
    #             # Condition 1: If first field matches regex "^##[^#]"
    #             # This matches lines that start with "##" but the third character is not "#".
    #             # if re.match(r'^##[^#]', first_field):
    #             #     fout.write(line)
    #             # Condition 2: If the first field starts with "#!".
    #             if first_field.startswith("#!"):
    #                 fout.write(line)
    #             # Condition 3: If the first field equals the specified reference chromosome.
    #             elif first_field == wildcards.ref_chromosome:
    #                 # Prepend the prefix to the entire line and write it.
    #                 fout.write(params.prefix + line)
    #             # Otherwise, the line is skipped (not printed).
    # shell:
    #     """
    #     awk '$1 ~ "^##[^#]"{{print}} $1 ~ "#!"{{print}} $1=="{params.ref_chromosome}"{{print "{params.prefix}"$0}}' {input.gff} > {output.chromosome_gff}
    #     """

####################

# def getChrMAFs(wildcards):
#     # List all files (e.g. *.maf) inside the directory for the given chromosome
#     return sorted(glob.glob(os.path.join(OUTPUT_DIR, "mafs", wildcards.chromosome_group, wildcards.ref_chromosome + ".maf")))[0]

def getChrMAFs(wildcards):
    # Wait for the checkpoint to complete and fetch its output manifest
    #manifest_file = os.path.join(OUTPUT_DIR, "mafs", wildcards.chromosome_group, "manifest.txt")
    # Query checkpoint using just the group (without extra keys):
    cp_output = checkpoints.maf_split_chr_by_group.get(chromosome_group=wildcards.chromosome_group).output
    # Derive the group directory from the marker file's location:
    group_dir = os.path.join(MAF_SPLIT_DIR, wildcards.chromosome_group)
    manifest_file = cp_output.maf_manifest  # This should be at group_dir/manifest.txt
    with open(manifest_file, "r") as mf:
        expected = [line.strip() for line in mf if line.strip()]
    for fname in expected:
        if fname == MAF_CHR_PREFIX + wildcards.ref_chromosome + ".maf":
            return os.path.join(group_dir, fname)
    raise ValueError(f"No MAF file for {wildcards.ref_chromosome} listed in {manifest_file}")

rule extract_4d_codons_by_chr:
    input:
        chromosome_maf = getChrMAFs,
        #chromosome_maf = os.path.join(OUTPUT_DIR, "mafs", "{chromosome_group}", "{ref_chromosome}.maf"),
        chromosome_gff = os.path.join(OUTPUT_DIR, "gffs", "{chromosome_group}", REF_GFF_FILE.replace(".gff", ".{ref_chromosome}.gff")),
        maf_manifest = os.path.join(MAF_SPLIT_DIR, "{chromosome_group}", "manifest.txt")
    output:
        codons_out = os.path.join(OUTPUT_DIR, "01-4d-codons", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-codons.ss"),
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "msa-view-4d-codons-{chromosome_group}-" + MAF_CHR_PREFIX + "{ref_chromosome}.log")
    resources:
        **getRuleResources("extract_4d_codons_by_chr")
    run:
        with open(log.job_log, "w") as log_stream, open(output.codons_out, "w") as out_stream:
            try:  
                cmd = [ "msa_view", input.chromosome_maf, "--4d", "--features", input.chromosome_gff ];

                COMMON.runCommand(cmd, log_stream, out_stream, "msa_view_codons", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise

    # shell:
    #     """
    #     msa_view {input.chromosome_maf} --4d --features {input.chromosome_gff} 2> {log} > {output}
    #     """

####################

rule extract_4d_sites:
    input:
        codons_in = os.path.join(OUTPUT_DIR, "01-4d-codons", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-codons.ss")
    output:
        sites_out = os.path.join(OUTPUT_DIR, "02-4d-sites", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites.ss")
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "msa-view-4d-sites-{chromosome_group}-" + MAF_CHR_PREFIX + "{ref_chromosome}.log")
    resources:
        **getRuleResources("extract_4d_sites")
    run:
        with open(log.job_log, "w") as log_stream, open(output.sites_out, "w") as out_stream:
            try:  
                cmd = [ "msa_view", input.codons_in,
                        "--in-format", "SS",
                        "--out-format", "SS",
                        "--tuple-size", "1" ];

                COMMON.runCommand(cmd, log_stream, out_stream, "msa_view_sites", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise        
    # shell:
    #     """
    #     msa_view {input} --in-format SS --out-format SS --tuple-size 1 2> {log} > {output}
    #     """ 

# ####################

rule filter_4d_sites:
    input:
        sites_in = os.path.join(OUTPUT_DIR, "02-4d-sites", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites.ss")
    output:
        sites_filtered = os.path.join(OUTPUT_DIR, "02-4d-sites", "{chromosome_group}-filtered", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + ".ss"),
        summary_file = os.path.join(OUTPUT_DIR, "summary-data", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + "-summary.tsv")
    params:
        filter_threshold = SEQ_THRESHOLD_4D,
        script_path = os.path.join(PIPELINE_DIR, "utils", "filter_4d_sites.py")
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "filter-4d-sites-{chromosome_group}-" + MAF_CHR_PREFIX + "{ref_chromosome}.log")
    resources:
        **getRuleResources("filter_4d_sites")
    run:
        with open(log.job_log, "w") as log_stream:
            try:  
                cmd = [ "python", params.script_path,
                        input.sites_in,
                        str(params.filter_threshold),
                        output.summary_file,
                        output.sites_filtered ];

                COMMON.runCommand(cmd, log_stream, log_stream, "filter_4d_sites", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise

    # shell:
    #     """
    #     python {params.script_path} {input} {params.filter_threshold} {output.summary_file} {output.ss_aln_4d_filtered}
    #     """

# ####################

rule run_phylofit:
    input:
        sites_filtered = os.path.join(OUTPUT_DIR, "02-4d-sites", "{chromosome_group}-filtered", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + ".ss"),
        tree = TREE_FILE
    output:
        mod_file = os.path.join(OUTPUT_DIR, "03-phylofit", "{chromosome_group}", "uncorrected-mods", MAF_CHR_PREFIX + "{ref_chromosome}.mod")
    params:
        prefix = os.path.join(OUTPUT_DIR, "03-phylofit", "{chromosome_group}", "uncorrected-mods", MAF_CHR_PREFIX + "{ref_chromosome}")
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "phylofit-{chromosome_group}-" + MAF_CHR_PREFIX + "{ref_chromosome}.log")
    resources:
        **getRuleResources("run_phylofit")
    run:
        with open(log.job_log, "w") as log_stream:
            try:  
                cmd = [ "phyloFit", "--tree", input.tree,
                        "--msa-format", "SS",
                        "--out-root", params.prefix,
                        input.sites_filtered ];

                COMMON.runCommand(cmd, log_stream, log_stream, "phylofit", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise
    # shell:
    #     """
    #     phyloFit --tree {input.tree} --msa-format SS --out-root {params.prefix} {input.ss_aln_4d_filtered}
    #     """

# ####################

rule get_gc_content:
    input:
        sample_file = SAMPLE_FILE
    output:
        gc_sample_file = GC_SAMPLE_FILE,
        avg_gc_file = AVG_GC_FILE
    params:
        script_path = os.path.join(PIPELINE_DIR, "utils", "get_gc_content.py")
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "get-gc-content.log")
    resources:
        **getRuleResources("get_gc_content")
    run:
        with open(log.job_log, "w") as log_stream:
            try:  
                cmd = [ "python", params.script_path,
                        input.sample_file,
                        output.gc_sample_file,
                        output.avg_gc_file ];

                COMMON.runCommand(cmd, log_stream, log_stream, "get_gc_content");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise

    # shell:
    #     """
    #     python {params.script_path} {input.sample_file} &> {log}
    #     """

# ####################

rule run_mod_freqs:
    input:
        mod_file = os.path.join(OUTPUT_DIR, "03-phylofit", "{chromosome_group}", "uncorrected-mods", MAF_CHR_PREFIX + "{ref_chromosome}.mod"),
        avg_gc_file = AVG_GC_FILE
    output:
        adj_mod_file = os.path.join(OUTPUT_DIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-corrected.mod")
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-modfreqs.log")
    resources:
        **getRuleResources("run_mod_freqs")
    run:
        with open(log.job_log, "w") as log_stream, open(output.adj_mod_file, "w") as out_stream:
            try:
                with open(input.avg_gc_file, "r") as f:
                    gc_value = f.read().strip()

                cmd = [ "modFreqs", input.mod_file, gc_value ];

                COMMON.runCommand(cmd, log_stream, out_stream, "run_mod_freqs", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise
    # shell:
    #     """
    #     gc=$(cat {input.avg_gc_file})
    #     modFreqs {input.mod_file} $gc > {output.adj_mod_file} 2> {log}
    #     """

# ####################

rule run_phylop:
    input:
        adj_mod_file = os.path.join(OUTPUT_DIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-corrected.mod"),
        chromosome_maf = getChrMAFs
    output:
        phylop_wig_file = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.wig")
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.log")
    resources:
        **getRuleResources("run_phylop")
    run:
        with open(log.job_log, "w") as log_stream, open(output.phylop_wig_file, "w") as out_stream:
            try:
                cmd = [ "phyloP", "--method", "LRT",
                        "--mode", "CONACC",
                        "--wig-scores",
                        "-i", "MAF",
                        input.adj_mod_file,
                        input.chromosome_maf ];

                COMMON.runCommand(cmd, log_stream, out_stream, "run_phylop", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise
    # shell:
    #     """
    #     phyloP --method LRT --mode CONACC --wig-scores -i MAF {input.adj_mod_file} {input.chromosome_maf} > {output.phylop_wig_file} 2> {log}
    #     """

# ####################

rule convert_wig_to_bed:
    input:
        phylop_wig_file = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.wig")
    output:
        phylop_bed_file = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.bed")
    params:
        script_path = os.path.join(PIPELINE_DIR, "utils", "convert_wig_to_bed.awk"),
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-wig-to-bed.log")
    resources:
        **getRuleResources("convert_wig_to_bed")
    run:
        with open(log.job_log, "w") as log_stream, open(output.phylop_bed_file, "w") as out_stream:
            try:
                cmd = [ "awk", "-f", params.script_path, input.phylop_wig_file ];

                COMMON.runCommand(cmd, log_stream, out_stream, "convert_wig_to_bed", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise
    # shell:
    #     """
    #     awk -v OFS="\t" 'BEGIN {{chrom=""; start=0; step=0; counter=0}} /^fixedStep/ {{split($2,a,"="); chrom=a[2]; split($3,b,"="); start=b[2]-1; split($4,c,"="); step=c[2]; counter=0; next}} {{pos=start+counter*step; end=pos+1; score=$1; if (score > 0) print chrom, pos, end, score, "0", 10^-score; else if (score < 0) print chrom, pos, end, score, "2", 10^score; else print chrom, pos, end, score, "1", 10^score; counter++}}' {input.phylop_wig_file} 2> {log} > {output.phylop_bed_file}
    #     """
# Convert phylop wig to bed and add new columns with the conservation status, raw p-value
# HEADERS:
# 1. chromosome
# 2. start position
# 3. end position
# 4. raw score/log(raw p-value)
# 5. conservation status (0 = conserved, 1 = neutral, 2 = accelerated)
# 6. raw p-value

# ####################

rule adjust_pvals:
    input:
        phylop_bed_file = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.bed")
    output:
        phylop_bed_file_fdr = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + ALPHA + ".bed")
    params:
        alpha = ALPHA,
        tmp_dir = TMPDIR,
        script_path = os.path.join(PIPELINE_DIR, "utils", "adjust_pvals.sh")
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.phylop-scores.fdr-" + ALPHA + ".log")
    resources:
        **getRuleResources("adjust_pvals")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                cmd = [ "bash", params.script_path,
                        input.phylop_bed_file,
                        output.phylop_bed_file_fdr,
                        params.tmp_dir,
                        str(resources.cpus_per_task),
                        str(params.alpha) ];

                COMMON.runCommand(cmd, log_stream, log_stream, "adjust_pvals", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise

    # shell:
    #     """
    #     total=$(wc -l {input.phylop_bed_file} | cut -d ' ' -f1)

    #     sort -k6,6g --parallel={resources.cpus} -T {params.tmp_dir} {input.phylop_bed_file} | \
    #     awk -v OFS="\\t" '{{ print $0, NR }}' | \
    #     sort -rn -k7 --parallel={resources.cpus} -T {params.tmp_dir} | \
    #     awk -v OFS="\\t" -v total=${{total}} 'BEGIN{{prev_pn=999}} {{p_adj=(total/$7)*$6; if(p_adj>prev_pn){{p_adj=prev_pn}}; prev_pn=p_adj; if(p_adj>1){{p_adj=1}}; sig=0; if(p_adj<0.05){{sig=1}}; print $0, p_adj, sig}}' | \
    #     sort -k 1,1 -k2,2n --parallel={resources.cpus} -T {params.tmp_dir} |\
    #     cut -f 1,2,3,5,8,9 2> {log} > {output.phylop_bed_file_fdr}
    #     """
# Adjust the raw phylop p-values for multiple testing using the Benjamini-Hochberg method
# STEPS:
# 1. Sort by p-value in ascending order
# 2. Add a column with the row number as the rank of that p-value
# 3. Sort by p-value in descending order so it is easy to check the adjusted p-value of the next highest raw p-value as per the BH method
# 4. Calculate the adjusted p-value using the formula: p_adj = (total number of tests / rank of p-value) * raw p-value
#       If the adjusted p-value is greater than the previous adjusted p-value, set the adjusted p-value to the previous adjusted p-value
#       If the adjusted p-value is greater than 1, set the adjusted p-value to 1
#       Add a column with the adjusted p-value
#       Add a column with the significance status (0 = not significant, 1 = significant)
# 9. Sort by chromosome and then by position numerically
# 10. Remove the row number and raw p-value columns
# 11. Write the output to a new file
# HEADERS:
# 1. chromosome
# 2. start position
# 3. end position
# 4. conservation status (0 = conserved, 1 = neutral, 2 = accelerated)
# 5. adjusted p-value
# 6. significant (0 = not significant, 1 = significant)

# ####################

rule get_conserved_sites:
    input:
        phylop_bed_file_fdr = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + ALPHA + ".bed")
    output:
        conserved_sites_bed = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + ALPHA + ".conserved.bed")
    params:
        script_path = os.path.join(PIPELINE_DIR, "utils", "get_conserved_sites.awk"),
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.phylop-scores.fdr-" + ALPHA + ".conserved.log")
    resources:
        **getRuleResources("get_conserved_sites")
    run:
        with open(log.job_log, "w") as log_stream, open(output.conserved_sites_bed, "w") as out_stream:
            try:
                cmd = [ "awk", "-f", params.script_path,
                        input.phylop_bed_file_fdr ];

                COMMON.runCommand(cmd, log_stream, out_stream, "get_conserved_sites", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise
    # shell:
    #     """
    #     awk -v OFS='\\t' '$4=="0" && $6=="1"{{print $1, $2, $3, $4$6}}' {input.phylop_bed_file_fdr} 2> {log} > {output.conserved_sites_bed}
    #     """

# ####################

rule get_accelerated_sites:
    input:
        phylop_bed_file_fdr = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + ALPHA + ".bed")
    output:
        accelerated_sites_bed = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + ALPHA + ".accelerated.bed")
    params:
        script_path = os.path.join(PIPELINE_DIR, "utils", "get_accelerated_sites.awk"),
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.phylop-scores.fdr-" + ALPHA + ".accelerated.log")
    resources:
        **getRuleResources("get_accelerated_sites")
    run:
        with open(log.job_log, "w") as log_stream, open(output.accelerated_sites_bed, "w") as out_stream:
            try:
                cmd = [ "awk", "-f", params.script_path,
                        input.phylop_bed_file_fdr ];

                COMMON.runCommand(cmd, log_stream, out_stream, "get_accelerated_sites", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise

# ####################

rule site_counts:
    input:
        conserved_sites_bed = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + ALPHA + ".conserved.bed"),
        accelerated_sites_bed = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + ALPHA + ".accelerated.bed")
    output:
        conserved_site_counts = os.path.join(OUTPUT_DIR, "summary-data", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.conserved-site-counts." + ALPHA + ".tsv"),
        accelerated_site_counts = os.path.join(OUTPUT_DIR, "summary-data", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.accelerated-site-counts." + ALPHA + ".tsv")
    # params:
    #     script_path = os.path.join(PIPELINE_DIR, "utils", "conserved_site_counts.awk"),
    log:
        job_log = os.path.join(OUTPUT_DIR, "logs", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.site-counts.log")
    resources:
        **getRuleResources("site_counts")
    run:
        with open(log.job_log, "w") as log_stream, open(output.conserved_site_counts, "w") as out_stream:
            try:
                cnt = sum(1 for _ in open(input.conserved_sites_bed))
                out_stream.write(f"{wildcards.ref_chromosome}\t{cnt}\n")
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise

        with open(log.job_log, "a") as log_stream, open(output.accelerated_site_counts, "w") as out_stream:
            try:
                cnt = sum(1 for _ in open(input.accelerated_sites_bed))
                out_stream.write(f"{wildcards.ref_chromosome}\t{cnt}\n")
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise

        # with open(log.job_log, "w") as log_stream, open(output.conserved_site_counts, "w") as out_stream:
        #     try:
        #         #cmd = [ "awk", "-f", params.script_path,
        #         #        input.conserved_sites_bed ];

        #         cmd = [ "wc", "-l", input.conserved_sites_bed,
        #                 ">", output.conserved_site_counts ];

        #         COMMON.runCommand(cmd, log_stream, out_stream, "conserved_site_counts", wc=f"{wildcards.chromosome_group}-{wildcards.ref_chromosome}");
        #     except Exception as e:
        #         traceback.print_exc(file=log_stream)
        #         raise
    # shell:
    #     """
    #     awk '{{ count[$1]++ }} END {{ for (num in count) print num, count[num] }}' {input.conserved_sites_bed} > {output.conserved_site_counts}
    #     """

# ####################

# rule cluster_conserved_sites:
#     input:
#         conserved_sites_bed = os.path.join(OUTPUT_DIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + ALPHA + ".conserved.bed")
#     output:
#         conserved_site_clusters = os.path.join(OUTPUT_DIR, "04-phylop", "clustering", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop-fdr-" + ALPHA + ".conserved.clusters-{min_cluster_size}-{min_samples}.bed")
#     log:
#         os.path.join(OUTPUT_DIR, "logs", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.phylop-scores.fdr-" + ALPHA + ".conserved.clusters-{min_cluster_size}-{min_samples}.log")
#     params:
#         script_path = os.path.join(PIPELINE_DIR, "scripts", "cluster_conserved_sites.py"),
#         chromosome = "{ref_chromosome}",
#         min_cluster_size = "{min_cluster_size}",
#         min_samples = "{min_samples}"
#     resources:
#         mem = "24g",
#         time = "4:00:00"
#     shell:
#         """
#         python {params.script_path} {input.conserved_sites_bed} {output.conserved_site_clusters} {params.chromosome} {params.min_cluster_size} {params.min_samples} 2> {log}
#         """

#############################################################################
