#############################################################################
# Snakemake rule to extract 4d sites from a MAF
# Gregg Thomas, September 2023
#############################################################################

# snakemake -p -s 01_run_phylofit_4d.smk --configfile echolocation-cfg-2.yaml --profile profiles/slurm_profile/ --use-conda --dryrun --rulegraph | dot -Tpng > dags/01-run-phylofit-4d-dag.png

#############################################################################

import os
import sys

#############################################################################

SAMPLE_FILE = config["sample_file"];
sample_file_extension = os.path.splitext(SAMPLE_FILE)[1];
GC_SAMPLE_FILE = SAMPLE_FILE.replace(sample_file_extension, "-gc" + sample_file_extension);
AVG_GC_FILE = SAMPLE_FILE.replace(sample_file_extension, "-avg-gc" + sample_file_extension);

REF_CHROMOSOME_GROUPS = config["ref_chromosome_groups"];

PROJ_DIR = config["project_dir"];

PREFIX = config["prefix"];
# The prefix for all output files

MAF_PATH = config["maf"];
MAF_FILE = os.path.basename(MAF_PATH);
MAF_DIR = os.path.dirname(MAF_PATH);

MAF_REF_ID = config["maf_ref_id"];
MAF_CHR_PREFIX = config["maf_chr_prefix"];
MAF_REF_CHR_JOINER = config["maf_ref_chr_joiner"];
MAF_REF_PREFIX = MAF_REF_ID + MAF_REF_CHR_JOINER + MAF_CHR_PREFIX;

MAF_SPLIT_CHR_DIR = os.path.join(MAF_DIR, PREFIX + "-mafSplit");

SEQ_THRESHOLD_4D = config["filter_threshold_4d"]

TREE_FILE = config["tree_file"];
species_tree = open(TREE_FILE, "r").read().strip();
topology = re.sub('[)][\d\w<>/.eE_:-]+', ')', species_tree);
topology = re.sub(':[\d.eE-]+', '', topology);
## Remove the branch lengths and node labels from the input tree string

SPECIES = [ tip_label for tip_label in topology.replace("(","").replace(")","").replace(";","").split(",") ];
## Get the species names from the tree string

REF_INDEX = config["ref_genome_index"];

REF_GFF_PATH = config["ref_gff"];
REF_GFF_FILE = os.path.basename(REF_GFF_PATH);

PREFIX = config["prefix"];

OUTDIR = config["outdir"];
os.makedirs(OUTDIR, exist_ok=True);
# Have to manually make the outdir because we parse the ref index into beds first

#############################################################################

REF_CHR_BED_DIR = os.path.join(OUTDIR, "beds");
os.makedirs(REF_CHR_BED_DIR, exist_ok=True);

for GROUP in REF_CHROMOSOME_GROUPS:
    with open(os.path.join(REF_CHR_BED_DIR, GROUP + ".bed"), "w") as f:
        for line in open(REF_INDEX, "r"):
            chrome, length = line.strip().split("\t")[:2];
            if chrome in GROUP:
                length = int(length)-1;
                f.write(f"{chrome}\t0\t{length}\n");
        
        for chrome in GROUP:
            f.write(f"{chrome}\t0\t{length}\n");
# This gets the chromosome start and end positions from the reference genome index file
# for splitting the MAF

#CHROMES = [MAF_REF + "." + c for c in CHROMES];

flattened_chromosome_groups = [(group, chromosome) for group, chromosome_list in REF_CHROMOSOME_GROUPS.items() for chromosome in chromosome_list];
REF_CHR_GROUPS_LIST, REF_CHROMOSOMES = zip(*flattened_chromosome_groups);
# Get two lists of equal length to use for wild cards

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(OUTDIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.wig"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)
        # run_phylop

        #expand(os.path.join(OUTDIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-corrected.mod"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)
        # run_mod_freqs

        #expand(os.path.join(OUTDIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.mod"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)
        # run_phylofit

        #expand(os.path.join(OUTDIR, "02-4d-sites", "{chromosome_group}-filtered", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + ".ss"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)
        # filter_4d_sites

        #expand(os.path.join(OUTDIR, "02-4d-sites", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites.ss"), zip, chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES),
        # extract_4d_sites

        #expand(os.path.join(MAF_SPLIT_CHR_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.00.maf"), chromosome_group=REF_CHR_GROUPS_LIST, ref_chromosome=REF_CHROMOSOMES)
        # split_maf_chr

#############################################################################
# Pipeline rules

rule split_maf_chr:
    input:
        maf = MAF_PATH,
        ref_chr_bed = os.path.join(REF_CHR_BED_DIR, "{chromosome_group}.bed")
    output:
        chromosome_maf = expand(os.path.join(MAF_SPLIT_CHR_DIR, "{{chromosome_group}}", MAF_CHR_PREFIX + "{ref_chromosome}.00.maf"), ref_chromosome=REF_CHROMOSOMES)
    params:
        outdir = os.path.join(MAF_SPLIT_CHR_DIR, "{chromosome_group}"),
        prefix = os.path.join(MAF_SPLIT_CHR_DIR, "{chromosome_group}") + "/"
    conda:
        "envs/ucsc-mafsplit.yml"
        # mafSplit needs its own environment right now because the current version (377) requires an outdated
        # openssl library (1.0.0) that doesn't play nice with other tools
    resources:
        partition="intermediate",
        mem="64g",
        time="120:00:00"
    shell:
        """
        mkdir -p {params.outdir}
        mafSplit {input.ref_chr_bed} {params.prefix} {input.maf}
        """

####################

rule split_ref_gff:
    input:
        gff = REF_GFF_PATH
    output:
        chromosome_gff = os.path.join(OUTDIR, "gffs", "{chromosome_group}", REF_GFF_FILE.replace(".gff", "{ref_chromosome}.gff"))
    params:
        ref_chromosome = "{ref_chromosome}",
        prefix = MAF_REF_PREFIX
    resources:
        mem = "4g",
        time = "1:00:00"
    shell:
        """
        awk '$1 ~ "^##[^#]"{{print}} $1 ~ "#!"{{print}} $1=="{params.ref_chromosome}"{{print "{params.prefix}"$0}}' {input.gff} > {output.chromosome_gff}
        """

####################

rule extract_4d_codons:
    input:
        chromosome_maf = os.path.join(MAF_SPLIT_CHR_DIR, "{chromosome_group}",  MAF_CHR_PREFIX + "{ref_chromosome}.00.maf"),
        chromosome_gff = os.path.join(OUTDIR, "gffs", "{chromosome_group}", REF_GFF_FILE.replace(".gff", "{ref_chromosome}.gff"))
    output:
        os.path.join(OUTDIR, "01-4d-codons", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-codons.ss")
    log:
        os.path.join(OUTDIR, "logs", "msa-view-4d-codons-{chromosome_group}-" + MAF_CHR_PREFIX + "{ref_chromosome}.log")
    resources:
        cpus = 1,
        mem = "256g",
        partition = "intermediate,holy-cow",
        time = "120:00:00"
    shell:
        """
        msa_view {input.chromosome_maf} --4d --features {input.chromosome_gff} 2> {log} > {output}
        """

####################

rule extract_4d_sites:
    input:
        os.path.join(OUTDIR, "01-4d-codons", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-codons.ss")
    output:
        os.path.join(OUTDIR, "02-4d-sites", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites.ss")
    log:
        os.path.join(OUTDIR, "logs", "msa-view-4d-sites-{chromosome_group}-" + MAF_CHR_PREFIX + "{ref_chromosome}.log")
    resources:
        cpus = 1,
        mem = "4g",
        time = "1:00:00"
    shell:
        """
        msa_view {input} --in-format SS --out-format SS --tuple-size 1 2> {log} > {output}
        """

####################

# rule aggregate_4d_sites:
#     input:
#         expand(os.path.join(OUTDIR, "02-4d-sites", "chr{chrome}-4d-sites.ss"), chrome=CHROMES)
#     output:
#         os.path.join(OUTDIR, "02-4d-sites", "all-4d-sites.ss")
#     params:
#         single_file = expand(os.path.join(OUTDIR, "02-4d-sites", "chr{chrome}-4d-sites.ss"), chrome=CHROMES)[1]
#     log:
#         os.path.join(OUTDIR, "logs", "msa-view-aggregate.log")
#     resources:
#         cpus = 1,
#         mem = "4g",
#         time = "1:00:00"
#     shell:
#         """
#         spec_list=$(grep "NAMES = " {params.single_file} | sed 's/NAMES = //g')
#         msa_view --unordered-ss --out-format SS --aggregate $spec_list {input} 2> {log} > {output}
#         """    

####################

rule filter_4d_sites:
    input:
        os.path.join(OUTDIR, "02-4d-sites", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites.ss")
    output:
        ss_aln_4d_filtered = os.path.join(OUTDIR, "02-4d-sites", "{chromosome_group}-filtered", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + ".ss"),
        summary_file = os.path.join(PROJ_DIR, "summary-data", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + "-summary.tsv")
    params:
        filter_threshold = SEQ_THRESHOLD_4D,
        script_path = os.path.join(PROJ_DIR, "scripts", "filter_4d_sites.py")
    resources:
        mem = "4g",
        time = "1:00:00"
    shell:
        """
        python {params.script_path} {input} {params.filter_threshold} {output.summary_file} {output.ss_aln_4d_filtered}
        """

####################

rule run_phylofit:
    input:
        ss_aln_4d_filtered = os.path.join(OUTDIR, "02-4d-sites", "{chromosome_group}-filtered", MAF_CHR_PREFIX + "{ref_chromosome}-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + ".ss"),
        tree = TREE_FILE
    output:
        mod_file = os.path.join(OUTDIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.mod")
    params:
        prefix = os.path.join(OUTDIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}")
    resources:
        mem = "32g",
        time = "48:00:00",
	    partition = "intermediate"
    shell:
        """
        phyloFit --tree {input.tree} --msa-format SS --out-root {params.prefix} {input.ss_aln_4d_filtered}
        """

####################

rule get_gc_content:
    input:
        sample_file = SAMPLE_FILE
    output:
        gc_sample_file = GC_SAMPLE_FILE,
        avg_gc_file = AVG_GC_FILE
    params:
        script_path = os.path.join(PROJ_DIR, "scripts", "get_gc_content.py")
    log:
        os.path.join(OUTDIR, "logs", "get_gc_content.log")
    resources:
        mem = "4g",
        time = "1:00:00"
    shell:
        """
        python {params.script_path} {input.sample_file} &> {log}
        """

####################

rule run_mod_freqs:
    input:
        mod_file = os.path.join(OUTDIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.mod"),
        avg_gc_file = AVG_GC_FILE
    output:
        adj_mod_file = os.path.join(OUTDIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-corrected.mod")
    log:
        os.path.join(OUTDIR, "logs", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-modfreqs.log")
    resources:
        mem = "4g",
        time = "1:00:00"
    shell:
        """
        gc=$(cat {input.avg_gc_file})
        modFreqs {input.mod_file} $gc > {output.adj_mod_file} 2> {log}
        """

####################

rule run_phylop:
    input:
        adj_mod_file = os.path.join(OUTDIR, "03-phylofit", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-corrected.mod"),
        chromosome_maf = os.path.join(MAF_SPLIT_CHR_DIR, "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}.00.maf")
    output:
        wig_file = os.path.join(OUTDIR, "04-phylop", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.wig")
    log:
        os.path.join(OUTDIR, "logs", "{chromosome_group}", MAF_CHR_PREFIX + "{ref_chromosome}-phylop.log")
    resources:
        mem = "128g",
        time = "96:00:00",
        partition = "intermediate"
    shell:
        """
        phyloP --method LRT --mode CONACC --wig-scores -i MAF {input.adj_mod_file} {input.chromosome_maf} > {output.wig_file} 2> {log}
        """

#############################################################################

####################
# STASH
####################

# rule prepare_maffilter:
#     input:
#         maf = MAF_PATH
#     output:
#         maffilter_script = os.path.join(OUTDIR, "03-maffilter", "scripts", "{chrome}-maffilter.script")
#     params:
#         script_path = os.path.join(PROJ_DIR, "scripts", "prep_maffilter.py"),
#         specs = ",".join(SPECIES),
#         ref = MAF_REF_ID,
#         chrome = MAF_CHR_PREFIX + "{chrome}",
#         maffilter_log = os.path.join(OUTDIR, "03-maffilter", "logs", "{chrome}-maffilter.log"),
#         maffilter_out = os.path.join(OUTDIR, "03-maffilter", "{chrome}-maffilter.tsv")        
#     resources:
#         time = "1:00:00",
#         mem = "4g"
#     shell:
#         """
#         python {params.script_path} {input.maf} {params.specs} {params.ref} {params.chrome} {params.maffilter_log} {params.maffilter_out} {output.maffilter_script}
#         """

# ####################

# rule run_maffilter:
#     input:
#         maf = MAF_PATH,
#         maffilter_script = os.path.join(OUTDIR, "03-maffilter", "scripts", "{chrome}-maffilter.script")
#     output:
#         maffilter_out = os.path.join(OUTDIR, "03-maffilter", "{chrome}-maffilter.tsv") 
#     log:
#         os.path.join(OUTDIR, "03-maffilter", "logs", "{chrome}-maffilter.log")
#     resources:
#         mem = "32g",
#         time = "120:00:00",
#         partition = "intermediate"
#     shell:
#         """
#         /usr/bin/time -f '%Uu %Ss %er %MkB %x %C' maffilter param={input.maffilter_script}
#         """
