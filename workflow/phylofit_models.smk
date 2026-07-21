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


def _as_bool(x, default=False):
    if x is None:
        return default
    if isinstance(x, bool):
        return x
    if isinstance(x, (int, float)):
        return bool(x)
    s = str(x).strip().lower()
    if s in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if s in {"0", "false", "f", "no", "n", "off", ""}:
        return False
    raise ValueError(f"Cannot interpret boolean value from '{x}'")

#############################################################################
# System setup

config_flag = config.get("display", False);
version_flag = config.get("version", False);
info_flag = config.get("info", False);
debug = config.get("debug", False);
#debug = True;
# A hacky way to get some custom command line arguments for the pipeline
# These just control preprocessing flags that stop the pipeline early anyways

if config.get("__master_workflow__", False):
    setup = config["__pipeline_setup__"];
    MAIN = setup["MAIN"];
    DRY_RUN = setup["DRY_RUN"];
    OUTPUT_DIR = setup["OUTPUT_DIR"];
    LOG_DIR = setup["LOG_DIR"];
    TMPDIR = setup["TMPDIR"];
    LOG_LEVEL = setup["LOG_LEVEL"];
    LOG_VERBOSITY = setup["LOG_VERBOSITY"];
else:
    MAIN, DRY_RUN, OUTPUT_DIR, LOG_DIR, TMPDIR, LOG_LEVEL, LOG_VERBOSITY = COMMON.pipelineSetup(config, sys.argv, version_flag, info_flag, config_flag, debug, workflow);
# Setup the pipeline once in the master workflow, or locally when this module
# is run as a standalone Snakefile.

MLOG = logging.getLogger('META')
# Setup logging if debugging

MODULE_DIR = os.path.dirname(workflow.snakefile);
PIPELINE_DIR = os.path.dirname(MODULE_DIR) if os.path.basename(MODULE_DIR) == "workflow" else MODULE_DIR;
# Resolve the repository root whether this file is run directly from workflow/
# or included from the top-level Snakefile.

UTILS_DIR = os.path.join(PIPELINE_DIR, "utils");
# The directory where the utility scripts are located

getRuleResources = partial(COMMON.getResources, config)
# This maps the function to get rule resources from the config file
# so we don't have to pass config each time we call it

#############################################################################
# Input files and output paths

USE_GC_CORRECTED_MODELS = _as_bool(
    config.get("use_gc_corrected_models", config.get("apply_gc_correction", True)),
    True,
)
SAMPLE_FILE = str(config.get("sample_file") or "").strip()
GC_SOURCE = "sample_file" if SAMPLE_FILE else "maf"
# Auto-detected, not user-configurable: if sample_file is set, use it (matches
# pre-existing configs unchanged); otherwise compute GC directly from the MAF.
# sample_file has no other purpose in this pipeline, so this is unambiguous, and
# it means blanking sample_file is a complete way to opt into MAF-based GC.
if not USE_GC_CORRECTED_MODELS:
    MLOG.warning(
        "use_gc_corrected_models=false; skipping GC correction. If the GC content "
        "of the 4d sites used to fit the neutral model differs from the genome-wide "
        "GC content, results may be affected."
    )

if GC_SOURCE == "sample_file":
    sample_file_extension = os.path.splitext(SAMPLE_FILE)[1];
    gc_sample_basename = os.path.basename(SAMPLE_FILE.replace(sample_file_extension, "-gc" + sample_file_extension))
    avg_gc_basename = os.path.basename(SAMPLE_FILE.replace(sample_file_extension, "-avg-gc" + sample_file_extension))
else:
    gc_sample_basename = ""
    avg_gc_basename = ""

MAF_PATH = config.get("maf") or "";
if not MAF_PATH:
    raise ValueError("run_phylofit=true requires maf to be set in config.")
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


REFERENCE_DIR = os.path.join(OUTPUT_DIR, "00-reference")
REFERENCE_GFF_DIR = os.path.join(REFERENCE_DIR, "gffs")
GROUP_BEDS_DIR = os.path.join(REFERENCE_DIR, "group-beds")

MAF_PREP_DIR = os.path.join(OUTPUT_DIR, "01-maf-prep")
MAF_SPLIT_DIR = COMMON.getOptionalConfigPath(
    config,
    "maf_split_chr_dir",
    os.path.join(MAF_PREP_DIR, "maf-by-chromosome"),
)

NEUTRAL_MODEL_DIR = os.path.join(OUTPUT_DIR, "02-neutral-model")
CODONS_DIR = os.path.join(NEUTRAL_MODEL_DIR, "4d-codons")
SITES_DIR = os.path.join(NEUTRAL_MODEL_DIR, "4d-sites")
SITES_RAW_DIR = os.path.join(SITES_DIR, "raw")
SITES_FILTERED_DIR = os.path.join(SITES_DIR, "filtered")
NEUTRAL_SUMMARY_DIR = os.path.join(NEUTRAL_MODEL_DIR, "summary")
PHYLOFIT_DIR = COMMON.getOptionalConfigPath(
    config,
    "phylofit_chr_dir",
    os.path.join(NEUTRAL_MODEL_DIR, "phylofit"),
)
PHYLOFIT_UNCORRECTED_DIR = os.path.join(PHYLOFIT_DIR, "{chromosome_group}", "uncorrected-mods")
PHYLOFIT_ACTIVE_MODEL_PATH = (
    os.path.join(PHYLOFIT_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-corrected.mod")
    if USE_GC_CORRECTED_MODELS
    else os.path.join(PHYLOFIT_DIR, "{chromosome_group}", "uncorrected-mods", MAF_CHR_PREFIX + "{ref_chromosome}.mod")
)
config["use_gc_corrected_models"] = USE_GC_CORRECTED_MODELS

GC_SUMMARY_DIR = os.path.join(NEUTRAL_SUMMARY_DIR, "gc")
GC_MAF_OUTPUT_PREFIX = os.path.join(GC_SUMMARY_DIR, MAF_FILE)
GC_SAMPLE_FILE = (
    GC_MAF_OUTPUT_PREFIX + ".gc.csv" if GC_SOURCE == "maf"
    else os.path.join(GC_SUMMARY_DIR, gc_sample_basename)
);
AVG_GC_FILE = (
    GC_MAF_OUTPUT_PREFIX + ".gc.mean.txt" if GC_SOURCE == "maf"
    else os.path.join(GC_SUMMARY_DIR, avg_gc_basename)
);
# Chromosome-split MAF directory used across workflows

# REF_CHR_BED_DIR = os.path.join(OUTPUT_DIR, "beds");
# REF_CHR_GFF_DIR = os.path.join(REFERENCE_DIR, "gffs");
# MAF_SPLIT_CHR_DIR = os.path.join(OUTPUT_DIR, PREFIX + "-mafSplit"); # Maybe this should be in a dir specified by the user?
# Various output sub-directories

REF_FASTA = config.get("ref_fasta") or "";
if not REF_FASTA:
    raise ValueError("run_phylofit=true requires ref_fasta to be set in config.")
REF_INDEX = COMMON.getOptionalConfigPath(
    config,
    "ref_fasta_index",
    COMMON.getOptionalConfigPath(config, "ref_genome_index", REF_FASTA + ".fai"),
);
if os.path.abspath(REF_INDEX) != os.path.abspath(REF_FASTA + ".fai"):
    raise ValueError(
        f"ref_fasta_index must match ref_fasta + '.fai' because samtools faidx writes next to the FASTA. "
        f"Got ref_fasta={REF_FASTA}, ref_fasta_index={REF_INDEX}"
    )
REF_GFF_PATH = config.get("ref_gff") or "";
if not REF_GFF_PATH:
    raise ValueError("run_phylofit=true requires ref_gff to be set in config.")
REF_GFF_FILE = os.path.basename(REF_GFF_PATH);
REF_CHROMOSOME_GROUPS = config["ref_chromosome_groups"];
# Reference genome info

#############################################################################
# Basic tree parsing

TREE_FILE = config.get("tree_file") or "";
if not TREE_FILE:
    raise ValueError("run_phylofit=true requires tree_file to be set in config.")
species_tree = open(TREE_FILE, "r").read().strip();
topology = re.sub(r'[)][\d\w<>/.eE_:-]+', ')', species_tree);
topology = re.sub(r':[\d.eE-]+', '', topology);
## Remove the branch lengths and node labels from the input tree string

SPECIES = [ tip_label for tip_label in topology.replace("(","").replace(")","").replace(";","").split(",") ];
## Get the species names from the tree string

#############################################################################
# Other params

SEQ_THRESHOLD_4D = config.get("filter_threshold_4d", 0.5)
# The threshold for filtering 4d sites

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
# def get_combinations():
#     combinations = []
#     for chromosome_group, ref_chromosomes in REF_CHROMOSOME_GROUPS.items():
#         for ref_chromosome in ref_chromosomes:
#             for min_cluster_size in MIN_CLUSTER_SIZES:
#                 for min_samples in MIN_SAMPLES:
#                     combinations.append((chromosome_group, ref_chromosome, min_cluster_size, min_samples))
#     return combinations

# # Helper function to generate combinations
# def generate_combinations():
#     combinations = []
#     for chromosome_group, ref_chromosomes in REF_CHROMOSOME_GROUPS.items():
#         for ref_chromosome in ref_chromosomes:
#             for min_cluster_size in MIN_CLUSTER_SIZES:
#                 for min_samples in MIN_SAMPLES:
#                     combinations.append({
#                         "chromosome_group": chromosome_group,
#                         "ref_chromosome": ref_chromosome,
#                         "min_cluster_size": min_cluster_size,
#                         "min_samples": min_samples
#                     })
#     return combinations

#############################################################################

PHYLOFIT_MODELS_STANDALONE = not bool(config.get("__master_workflow__", False))

if PHYLOFIT_MODELS_STANDALONE:
    localrules: all

    rule all:
        input:
            expand(PHYLOFIT_ACTIVE_MODEL_PATH, zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)

#############################################################################
# Pipeline rules

if not bool(config.get("__ref_fasta_index_rule_defined__", False)):
    config["__ref_fasta_index_rule_defined__"] = True

    rule ref_fasta_index:
        input:
            ref_fasta = REF_FASTA
        output:
            ref_fasta_index = REF_INDEX
        log:
            job_log = os.path.join(LOG_DIR, "ref_fasta_index", "run.log")
        benchmark:
            os.path.join(LOG_DIR, "benchmarks", "ref_fasta_index", "run.txt")
        resources:
            **getRuleResources("ref_fasta_index")
        run:
            with open(log.job_log, "w") as log_stream:
                try:
                    cmd = ["samtools", "faidx", input.ref_fasta]
                    COMMON.runCommand(cmd, log_stream, log_stream, "ref_fasta_index")
                except Exception:
                    traceback.print_exc(file=log_stream)
                    raise

rule maf_index:
    input:
        maf = MAF_PATH
    output:
        maf_index_block = MAF_INDEX_BLOCK,
        maf_index_scaff = MAF_INDEX_SCAFF
    log:
        job_log = os.path.join(LOG_DIR, "maf_index", "run.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "maf_index", "run.txt")
    resources:
        **getRuleResources("maf_index")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                cmd = [ "mafutils", "index",
                        input.maf, output.maf_index_block, output.maf_index_scaff ];
                COMMON.runCommand(cmd, log_stream, log_stream, "maf_index");
            except Exception as e:
                traceback.print_exc(file=log_stream)
                raise


####################

rule make_group_beds:
    input:
        maf = MAF_PATH,
        ref_fasta_index = REF_INDEX
    output:
        chr_group_bed = os.path.join(GROUP_BEDS_DIR, "{chromosome_group}.bed")
    params:
        ref_chroms = lambda wildcards: REF_CHROMOSOME_GROUPS[wildcards.chromosome_group],
        chr_prefix = MAF_CHR_PREFIX,
        script_path = os.path.join(UTILS_DIR, "make_group_beds.py")
    log:
        job_log = os.path.join(LOG_DIR, "make_group_beds", "{chromosome_group}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "make_group_beds", "{chromosome_group}.txt")
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
        chr_group_bed = os.path.join(GROUP_BEDS_DIR, "{chromosome_group}.bed")
    output:
        maf_manifest = os.path.join(MAF_SPLIT_DIR, "{chromosome_group}", "manifest.txt")
        #chr_group_dir = directory(os.path.join(OUTPUT_DIR, "mafs", "{chromosome_group}")),
        #chr_maf = os.path.join(OUTPUT_DIR, "mafs", "{chromosome_group}", "{ref_chromosome}.maf")
    params:
        outdir = os.path.join(MAF_SPLIT_DIR, "{chromosome_group}"),
        rule_name = "maf_split_by_chr"
    log:
        job_log = os.path.join(LOG_DIR, "maf_split_chr_by_group", "{chromosome_group}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "maf_split_chr_by_group", "{chromosome_group}.txt")
    resources:
        **getRuleResources("maf_split_chr_by_group")
    run:
        with open(log.job_log, "w") as log_stream:
            try:
                cmd = [ "mafutils", "fetch",
                        input.maf,
                        input.chr_group_bed,
                        "-i", input.maf_index_scaff,
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
        chromosome_gff = os.path.join(REFERENCE_GFF_DIR, "{chromosome_group}", REF_GFF_FILE.replace(".gff", ".{ref_chromosome}.gff"))
    params:
        ref_chr = lambda wildcards: wildcards.ref_chromosome,
        prefix = MAF_REF_PREFIX,
        script_path = os.path.join(UTILS_DIR, "ref_gff_split_by_chr.awk")
    log:
        job_log = os.path.join(LOG_DIR, "ref_gff_split_by_chr", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "ref_gff_split_by_chr", "{chromosome_group}", "{ref_chromosome}.txt")
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
        chromosome_gff = os.path.join(REFERENCE_GFF_DIR, "{chromosome_group}", REF_GFF_FILE.replace(".gff", ".{ref_chromosome}.gff")),
        maf_manifest = os.path.join(MAF_SPLIT_DIR, "{chromosome_group}", "manifest.txt")
    output:
        codons_out = os.path.join(CODONS_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-codons.ss"),
    log:
        job_log = os.path.join(LOG_DIR, "extract_4d_codons_by_chr", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "extract_4d_codons_by_chr", "{chromosome_group}", "{ref_chromosome}.txt")
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
        codons_in = os.path.join(CODONS_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-codons.ss")
    output:
        sites_out = os.path.join(SITES_RAW_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites.ss")
    log:
        job_log = os.path.join(LOG_DIR, "extract_4d_sites", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "extract_4d_sites", "{chromosome_group}", "{ref_chromosome}.txt")
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
        sites_in = os.path.join(SITES_RAW_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites.ss")
    output:
        sites_filtered = os.path.join(SITES_FILTERED_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + ".ss"),
        summary_file = os.path.join(NEUTRAL_SUMMARY_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + "-summary.tsv")
    params:
        filter_threshold = SEQ_THRESHOLD_4D,
        script_path = os.path.join(PIPELINE_DIR, "utils", "filter_4d_sites.py")
    log:
        job_log = os.path.join(LOG_DIR, "filter_4d_sites", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "filter_4d_sites", "{chromosome_group}", "{ref_chromosome}.txt")
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
        sites_filtered = os.path.join(SITES_FILTERED_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + ".ss"),
        tree = TREE_FILE
    output:
        mod_file = os.path.join(PHYLOFIT_DIR, "{chromosome_group}", "uncorrected-mods", MAF_CHR_PREFIX + "{ref_chromosome}.mod")
    params:
        prefix = os.path.join(PHYLOFIT_DIR, "{chromosome_group}", "uncorrected-mods", MAF_CHR_PREFIX + "{ref_chromosome}")
    log:
        job_log = os.path.join(LOG_DIR, "run_phylofit", "{chromosome_group}", "{ref_chromosome}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "run_phylofit", "{chromosome_group}", "{ref_chromosome}.txt")
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

if USE_GC_CORRECTED_MODELS and GC_SOURCE == "sample_file":
    # Only meaningful (and only produces valid, non-empty output paths) when GC
    # correction is on and sample_file is the selected source; SAMPLE_FILE/
    # GC_SAMPLE_FILE/AVG_GC_FILE are blank otherwise.

    rule get_gc_content:
        input:
            sample_file = SAMPLE_FILE
        output:
            gc_sample_file = GC_SAMPLE_FILE,
            avg_gc_file = AVG_GC_FILE
        params:
            script_path = os.path.join(UTILS_DIR, "get_gc_content.py"),
            accession_header = config.get("accession_header", "")
        log:
            job_log = os.path.join(LOG_DIR, "get_gc_content", "run.log")
        benchmark:
            os.path.join(LOG_DIR, "benchmarks", "get_gc_content", "run.txt")
        resources:
            **getRuleResources("get_gc_content")
        run:
            with open(log.job_log, "w") as log_stream:
                try:
                    cmd = [ "python", params.script_path,
                            input.sample_file,
                            output.gc_sample_file,
                            output.avg_gc_file ];
                    if params.accession_header:
                        cmd.append(params.accession_header);

                    COMMON.runCommand(cmd, log_stream, log_stream, "get_gc_content");
                except Exception as e:
                    traceback.print_exc(file=log_stream)
                    raise

        # shell:
        #     """
        #     python {params.script_path} {input.sample_file} &> {log}
        #     """

# ####################

if USE_GC_CORRECTED_MODELS and GC_SOURCE == "maf":
    # Computes GC directly from the whole input MAF via `mafutils gc`, rather than
    # looking assemblies up externally via sample_file. Runs once on the whole MAF
    # (not per-chromosome) to match the single genome-wide average the sample_file
    # path already produces - could be split to per-chromosome-group later if
    # there's a reason to want that granularity.

    rule get_gc_content_from_maf:
        input:
            maf = MAF_PATH,
            maf_index_block = MAF_INDEX_BLOCK
        output:
            gc_sample_file = GC_SAMPLE_FILE,
            avg_gc_file = AVG_GC_FILE
        params:
            output_prefix = GC_MAF_OUTPUT_PREFIX
        log:
            job_log = os.path.join(LOG_DIR, "get_gc_content_from_maf", "run.log")
        benchmark:
            os.path.join(LOG_DIR, "benchmarks", "get_gc_content_from_maf", "run.txt")
        resources:
            **getRuleResources("get_gc_content_from_maf")
        run:
            with open(log.job_log, "w") as log_stream:
                try:
                    os.makedirs(os.path.dirname(params.output_prefix), exist_ok=True)
                    cmd = [
                        "mafutils", "gc",
                        input.maf, input.maf_index_block,
                        "-o", params.output_prefix,
                        "-p", str(int(resources.cpus_per_task)),
                    ]
                    COMMON.runCommand(cmd, log_stream, log_stream, "get_gc_content_from_maf");
                except Exception:
                    traceback.print_exc(file=log_stream)
                    raise

# ####################

if USE_GC_CORRECTED_MODELS:
    # Needed regardless of which gc_source populated AVG_GC_FILE above.

    rule run_mod_freqs:
        input:
            mod_file = os.path.join(PHYLOFIT_DIR, "{chromosome_group}", "uncorrected-mods", MAF_CHR_PREFIX + "{ref_chromosome}.mod"),
            avg_gc_file = AVG_GC_FILE
        output:
            adj_mod_file = os.path.join(PHYLOFIT_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-corrected.mod")
        log:
            job_log = os.path.join(LOG_DIR, "run_mod_freqs", "{chromosome_group}", "{ref_chromosome}.log")
        benchmark:
            os.path.join(LOG_DIR, "benchmarks", "run_mod_freqs", "{chromosome_group}", "{ref_chromosome}.txt")
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

#############################################################################
