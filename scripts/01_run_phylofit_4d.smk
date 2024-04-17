#############################################################################
# Snakemake rule to extract 4d sites from a MAF
# Gregg Thomas, September 2023
#############################################################################

# snakemake -p -s 01_run_phylofit_4d.smk --configfile echolocation-cfg.yaml --profile profiles/slurm_profile/ --dryrun --rulegraph | dot -Tpng > dags/01-run-phylofit-4d-dag.png

#############################################################################

import os
import sys

#############################################################################

SAMPLE_FILE = config["sample_file"];
sample_file_extension = os.path.splitext(SAMPLE_FILE)[1];
GC_SAMPLE_FILE = SAMPLE_FILE.replace(sample_file_extension, "-gc" + sample_file_extension);
AVG_GC_FILE = SAMPLE_FILE.replace(sample_file_extension, "-avg-gc" + sample_file_extension);

CHROMES = config["chromes"];

PROJ_DIR = config["project_dir"];

MAF_DIR = config["mafdir"];
MAF_OUTDIR = config["mafoutdir"];
MAF_FILE = config["maf"];
MAF_PATH = os.path.join(MAF_DIR, MAF_FILE);

MAF_REF_ID = config["maf_ref_id"];
MAF_CHR_PREFIX = config["maf_chr_prefix"];
MAF_REF_CHR_JOINER = config["maf_ref_chr_joiner"];
MAF_REF_PREFIX = MAF_REF_ID + MAF_REF_CHR_JOINER + MAF_CHR_PREFIX;

SEQ_THRESHOLD_4D = config["filter_threshold_4d"]

TREE_FILE = config["tree_file"];
species_tree = open(TREE_FILE, "r").read().strip();
topology = re.sub('[)][\d\w<>/.eE_:-]+', ')', species_tree);
topology = re.sub(':[\d.eE-]+', '', topology);
## Remove the branch lengths and node labels from the input tree string

SPECIES = [ tip_label for tip_label in topology.replace("(","").replace(")","").replace(";","").split(",") ];
## Get the species names from the tree string

REF_DIR = config["ref_dir"];
REF_GFF_FILE = config["ref_gff"];
REF_GFF_PATH = os.path.join(REF_DIR, REF_GFF_FILE);

PREFIX = config["prefix"];

#############################################################################

#CHROMES = [MAF_REF + "." + c for c in CHROMES];

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        os.path.join(MAF_OUTDIR, "03-phylop", PREFIX + "-phylop.wig")
        # run_phylop

        #os.path.join(MAF_OUTDIR, "02-phylofit", PREFIX + "-corrected.mod")
        # run_mod_freqs

        # expand(os.path.join(REF_DIR, "chromosome-gffs", REF_GFF_FILE.replace(".gff", "{chrome}.gff")), chrome=CHROMES),
        # # split_gff

        # expand(os.path.join(MAF_OUTDIR, "01-4d-codons", "chr{chrome}-4d-codons.ss"), chrome=CHROMES),
        # # extract_4d_codons

        # expand(os.path.join(MAF_OUTDIR, "01-4d-sites", "chr{chrome}-4d-sites.ss"), chrome=CHROMES),
        # # extract_4d_sites

        # os.path.join(MAF_OUTDIR, "01-4d-sites", "all-4d-sites.ss"),
        # # aggregate_4d_sites

        # os.path.join(MAF_OUTDIR, "02-phylofit", PREFIX + ".mod"),
        # # run_phylofit

        # GC_SAMPLE_FILE,
        # AVG_GC_FILE,
        # get_gc_content



#############################################################################
# Pipeline rules

rule split_ref_gff:
    input:
        gff = REF_GFF_PATH
    output:
        chrome_gff = os.path.join(REF_DIR, "chromosome-gffs", REF_GFF_FILE.replace(".gff", "{chrome}.gff"))
    params:
        chrome = "{chrome}",
        prefix = MAF_REF_PREFIX
    resources:
        mem = "4g",
        time = "1:00:00"
    shell:
        """
        awk '$1 ~ "^##[^#]"{{print}} $1 ~ "#!"{{print}} $1=="{params.chrome}"{{print "{params.prefix}"$0}}' {input.gff} > {output.chrome_gff}
        """

####################

rule extract_4d_codons:
    input:
        maf = MAF_PATH,
        chrome_gff = os.path.join(REF_DIR, "chromosome-gffs", REF_GFF_FILE.replace(".gff", "{chrome}.gff"))
    output:
        os.path.join(MAF_OUTDIR, "01-4d-codons", "chr{chrome}-4d-codons.ss")
    log:
        os.path.join(MAF_OUTDIR, "logs", "msa-view-4d-codons-chr{chrome}.log")
    resources:
        cpus = 1,
        mem = "128g",
        partition = "intermediate,holy-cow",
        time = "96:00:00"
    shell:
        """
        msa_view {input.maf} --4d --features {input.chrome_gff} 2> {log} > {output}
        """

####################

rule extract_4d_sites:
    input:
        os.path.join(MAF_OUTDIR, "01-4d-codons", "chr{chrome}-4d-codons.ss")
    output:
        os.path.join(MAF_OUTDIR, "01-4d-sites", "chr{chrome}-4d-sites.ss")
    log:
        os.path.join(MAF_OUTDIR, "logs", "msa-view-4d-sites-chr{chrome}.log")
    resources:
        cpus = 1,
        mem = "4g",
        time = "1:00:00"
    shell:
        """
        msa_view {input} --in-format SS --out-format SS --tuple-size 1 2> {log} > {output}
        """

####################

rule aggregate_4d_sites:
    input:
        expand(os.path.join(MAF_OUTDIR, "01-4d-sites", "chr{chrome}-4d-sites.ss"), chrome=CHROMES)
    output:
        os.path.join(MAF_OUTDIR, "01-4d-sites", "all-4d-sites.ss")
    params:
        single_file = expand(os.path.join(MAF_OUTDIR, "01-4d-sites", "chr{chrome}-4d-sites.ss"), chrome=CHROMES)[1]
    log:
        os.path.join(MAF_OUTDIR, "logs", "msa-view-aggregate.log")
    resources:
        cpus = 1,
        mem = "4g",
        time = "1:00:00"
    shell:
        """
        spec_list=$(grep "NAMES = " {params.single_file} | sed 's/NAMES = //g')
        msa_view --unordered-ss --out-format SS --aggregate $spec_list {input} 2> {log} > {output}
        """    

####################

rule filter_4d_sites:
    input:
        os.path.join(MAF_OUTDIR, "01-4d-sites", "all-4d-sites.ss")
    output:
        ss_aln_4d_filtered = os.path.join(MAF_OUTDIR, "01-4d-sites", "all-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + ".ss"),
        summary_file = os.path.join(PROJ_DIR, "summary-data", "all-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + "-summary.tsv")
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
        ss_aln_4d_filtered = os.path.join(MAF_OUTDIR, "01-4d-sites", "all-4d-sites-filtered-" + str(SEQ_THRESHOLD_4D) + ".ss"),
        tree = TREE_FILE
    output:
        mod_file = os.path.join(MAF_OUTDIR, "02-phylofit", PREFIX + ".mod")
    params:
        prefix = os.path.join(MAF_OUTDIR, "02-phylofit", PREFIX)
    resources:
        mem = "48g",
        time = "96:00:00",
	    partition = "holy-info"
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
        os.path.join(MAF_OUTDIR, "logs", "get_gc_content.log")
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
        mod_file = os.path.join(MAF_OUTDIR, "02-phylofit", PREFIX + ".mod"),
        avg_gc_file = AVG_GC_FILE
    output:
        adj_mod_file = os.path.join(MAF_OUTDIR, "02-phylofit", PREFIX + "-corrected.mod")
    log:
        os.path.join(MAF_OUTDIR, "logs", "modfreqs.log")
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
        mod_file = os.path.join(MAF_OUTDIR, "02-phylofit", PREFIX + "-corrected.mod"),
        maf = MAF_PATH
    output:
        wig_file = os.path.join(MAF_OUTDIR, "03-phylop", PREFIX + "-phylop.wig")
    log:
        os.path.join(MAF_OUTDIR, "logs", "phylop.log")
    resources:
        mem = "48g",
        time = "72:00:00"
    shell:
        """
        phyloP --method LRT --mode CONACC --wig-scores -i MAF {input.mod_file} {input.maf_file} > {output.wig_file} 2> {log}
        """

#############################################################################

####################
# STASH
####################

# rule prepare_maffilter:
#     input:
#         maf = MAF_PATH
#     output:
#         maffilter_script = os.path.join(MAF_OUTDIR, "03-maffilter", "scripts", "{chrome}-maffilter.script")
#     params:
#         script_path = os.path.join(PROJ_DIR, "scripts", "prep_maffilter.py"),
#         specs = ",".join(SPECIES),
#         ref = MAF_REF_ID,
#         chrome = MAF_CHR_PREFIX + "{chrome}",
#         maffilter_log = os.path.join(MAF_OUTDIR, "03-maffilter", "logs", "{chrome}-maffilter.log"),
#         maffilter_out = os.path.join(MAF_OUTDIR, "03-maffilter", "{chrome}-maffilter.tsv")        
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
#         maffilter_script = os.path.join(MAF_OUTDIR, "03-maffilter", "scripts", "{chrome}-maffilter.script")
#     output:
#         maffilter_out = os.path.join(MAF_OUTDIR, "03-maffilter", "{chrome}-maffilter.tsv") 
#     log:
#         os.path.join(MAF_OUTDIR, "03-maffilter", "logs", "{chrome}-maffilter.log")
#     resources:
#         mem = "32g",
#         time = "120:00:00",
#         partition = "intermediate"
#     shell:
#         """
#         /usr/bin/time -f '%Uu %Ss %er %MkB %x %C' maffilter param={input.maffilter_script}
#         """