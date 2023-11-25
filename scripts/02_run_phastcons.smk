#############################################################################
# Snakemake rule to extract 4d sites from a MAF
# Gregg Thomas, September 2023
#############################################################################

# snakemake -p -s 02_run_phastcons.smk --configfile echolocation-cfg.yaml --profile profiles/slurm_profile/ --dryrun --rulegraph | dot -Tpng > dag.png

#############################################################################

import os
import sys

#############################################################################

CHROMES = config["chromes"];

PROJ_DIR = config["project_dir"];

PREFIX = config["prefix"];

MAF_DIR = config["mafdir"];
MAF_OUTDIR = config["mafoutdir"];
MAF_FILE = config["maf"];
MAF_PATH = os.path.join(MAF_DIR, MAF_FILE);

MAF_REF_ID = config["maf_ref_id"];
MAF_CHR_PREFIX = config["maf_chr_prefix"];
MAF_REF_CHR_JOINER = config["maf_ref_chr_joiner"];
MAF_REF_PREFIX = MAF_REF_ID + MAF_REF_CHR_JOINER + MAF_CHR_PREFIX;

MAF_SPLIT_MB = str(config["maf_split_size_mb"]);

MAF_SPLIT_DIR = os.path.join(MAF_DIR, PREFIX + "-mafSplit-" + MAF_SPLIT_MB + "Mb", "");

REF_DIR = config["ref_dir"];
REF_INDEX_FILE = config["ref_fasta_index"];
REF_INDEX_PATH = os.path.join(REF_DIR, REF_INDEX_FILE);

REF_FA_FILE = config["ref_fasta"];
REF_FA_PATH = os.path.join(REF_DIR, REF_FA_FILE);

REF_GFF_CHR_FILE = config["ref_gff_chr"];
REF_GFF_CHR_PATH = os.path.join(REF_DIR, REF_GFF_CHR_FILE);

REF_SPLIT_BED = config["ref_split_bed"];
REF_SPLIT_BED = os.path.join(REF_DIR, REF_SPLIT_BED);

MOD_FILE = os.path.join(MAF_OUTDIR, "02-phylofit", PREFIX + "-corrected.mod")

chr_counts = {};
chr_list, win_list = [], [];

for line in open(REF_SPLIT_BED):
    cur_chr = line.strip().split("\t")[0];

    if cur_chr not in chr_counts:
        chr_counts[cur_chr] = 0;

    chr_list.append(cur_chr);

    if chr_counts[cur_chr] == 0:
        chr_counts[cur_chr] += 1;
        continue;
    # Not sure why mafSplit starts at 1, but it does

    cur_win = str(chr_counts[cur_chr]);
    while len(cur_win) < 2:
        cur_win = "0" + cur_win;
    win_list.append(cur_win);

    chr_counts[cur_chr] += 1;

# print(chr_counts["chr6"]);
# print(len(chr_counts["chr6"]));
# sys.exit();

#############################################################################

rule all:
    input:
        os.path.join(MAF_OUTDIR, "04-phyloboot-mods", MAF_REF_ID + ".avg.cons.mod"),
        os.path.join(MAF_OUTDIR, "04-phyloboot-mods", MAF_REF_ID + ".avg.noncons.mod")
        # run_phyloboot
        
        #expand(os.path.join(MAF_SPLIT_DIR, "{chrome}.{win}.maf"), zip, chrome=chr_list, win=win_list),
        # split_maf
        

        #expand(os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{chrome}.{win}.cons.mod"), zip, chrome=chr_list, win=win_list),
        #expand(os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{chrome}.{win}.noncons.mod"), zip, chrome=chr_list, win=win_list),     
        # # estimate_rho
#############################################################################

# rule get_split_bed:
#     input:
#         ref_fasta_index = REF_INDEX_PATH
#     output:
#         split_bed = os.path.join(REF_DIR, REF_INDEX_PATH.replace(".fa.fai", ".chr." + MAF_SPLIT_MB + "Mb.bed"))
#     resources:
#         mem = "1g",
#         time = "1:00:00"
#     shell:
#         """
#         bedtool makewindows -g {input.ref_fasta_index} -w {params.window_size} > {output.split_bed}
#         """

####################

rule split_maf:
    input:
        maf = MAF_PATH,
        split_bed = REF_SPLIT_BED
    output:
        outdir = MAF_SPLIT_DIR,
        chrome_maf = expand(os.path.join(MAF_SPLIT_DIR, "{chrome}.{win}.maf"), zip, chrome=chr_list, win=win_list)
    params:
        prefix = MAF_SPLIT_DIR
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
# On holybioinf:
# Splitting 1 files from /n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-sicb/data/00-human-ref-ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.chr.1Mb.bed to /n/holyscratch01/informatics/gwct/241-mammalian-2020v2b-mafSplit-1Mb/
# splitting /n/holyscratch01/informatics/gwct/241-mammalian-2020v2b.maf
# 79436.85user 11923.47system 25:28:11elapsed 99%CPU (0avgtext+0avgdata 7304maxresident)k
# 17477696inputs+20200869936outputs (0major+9728minor)pagefaults 0swaps

####################

checkpoint estimate_rho:
    input:
        mod_file = MOD_FILE,
        chrome_maf = os.path.join(MAF_SPLIT_DIR, "{chrome}.{win}.maf")
    output:
        rho_cons_mod = os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{chrome}.{win}.cons.mod"),
        rho_noncons_mod = os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{chrome}.{win}.noncons.mod")
    params:
        output_prefix = os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{chrome}.{win}"),
        length = 45,
        coverage = 0.3
    log:
        os.path.join(MAF_OUTDIR, "03-phastcons-rho", "logs", MAF_REF_ID + ".{chrome}.{win}.log")
    resources:
        mem = "24g",
        time = "4:00:00",
	    partition = "holy-info,shared"
    shell:
        """
        phastCons --expected-length {params.length} --target-coverage {params.coverage} --no-post-probs --msa-format MAF --estimate-rho {params.output_prefix} {input.chrome_maf} {input.mod_file} &> {log}
        """
    # Many of these will time out for some reason ... need a way to deal with this

####################

# def get_completed_mods(wildcard):
#     # completed_cons_mods = [];
#     # for f in os.listdir(os.path.join(MAF_OUTDIR, "03-phastcons-rho")):
#     #     if f.endswith(".cons.mod"):
#     #         completed_cons_mods.append(os.path.join(MAF_OUTDIR, "03-phastcons-rho", f));
#     # return(completed_cons_mods);

#     completed_cons_mods = [];
#     chrome = wildcard.chrome;
#     win = wildcard.win;

#     cons_mod = os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + "." + chrome + "." + win + ".cons.mod");
#     if os.path.isfile(cons_mod):
#         completed_cons_mods.append(cons_mod);

#     return completed_cons_mods;

def input_for_b(*wildcards):
    print(checkpoints.estimate_rho.get().output);
    return checkpoints.estimate_rho.get().output

rule run_phyloboot:
    input:
        #input_for_b
        rho_cons_mod = expand(os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{{chrome}}.{{win}}.cons.mod"), zip, chrome=chr_list, win=win_list),
        rho_noncons_mod = expand(os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{{chrome}}.{{win}}.noncons.mod"), zip, chrome=chr_list, win=win_list)
    output:
        avg_cons_mod = os.path.join(MAF_OUTDIR, "04-phyloboot-mods", MAF_REF_ID + ".avg.cons.mod"),
        avg_noncons_mod = os.path.join(MAF_OUTDIR, "04-phyloboot-mods", MAF_REF_ID + ".avg.noncons.mod")
    params:
        rho_mod_dir = os.path.join(MAF_OUTDIR, "03-phastcons-rho"),
        cons_mod_list_file = os.path.join(MAF_OUTDIR, "04-phyloboot-mods", MAF_REF_ID + ".cons.mod.list"),
        noncons_mod_list_file = os.path.join(MAF_OUTDIR, "04-phyloboot-mods", MAF_REF_ID + ".noncons.mod.list")
    shell:
        """
        ls {params.rho_mod_dir}/*.cons.mod > {params.cons_mod_list_file}
        phyloBoot --read-mods '*{params.cons_mod_list_file}' --output-average {output.avg_cons_mod}
        ls {params.rho_mod_dir}/*.noncons.mod > {params.noncons_mod_list_file}
        phyloBoot --read-mods '*{params.noncons_mod_list_file}' --output-average {output.avg_noncons_mod}
        """

#############################################################################

# rule all:
#     input:
#         os.path.join(REF_DIR, REF_INDEX_PATH.replace(".fa.fai", ".chr.bed")),
#         # get_chr_bed

#         expand(os.path.join(MAF_SPLIT_DIR, MAF_CHR_PREFIX + "{chrome}.00.maf"), chrome=CHROMES),
#         # split_maf
        
#         expand(os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_PREFIX + "{chrome}.cons.mod"), chrome=CHROMES),
#         expand(os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_PREFIX + "{chrome}.noncons.mod"), chrome=CHROMES)        
#         # estimate_rho

# #############################################################################

# rule get_chr_bed:
#     input:
#         ref_fasta_index = REF_INDEX_PATH
#     output:
#         chr_bed = os.path.join(REF_DIR, REF_INDEX_PATH.replace(".fa.fai", ".chr.bed"))
#     resources:
#         mem = "1g",
#         time = "1:00:00"
#     shell:
#         """
#         awk 'BEGIN{OFS="\t"} {print $1, 0, $2-1}' {input.ref_fasta_index} > {output.chr_bed}
#         """

# ####################

# rule split_maf:
#     input:
#         maf = MAF_PATH,
#         outdir = MAF_SPLIT_DIR,
#         chr_bed = os.path.join(REF_DIR, REF_INDEX_PATH.replace(".fa.fai", ".chr.bed"))
#     output:
#         chrome_maf = os.path.join(MAF_SPLIT_DIR, MAF_CHR_PREFIX + "{chrome}.00.maf")
#     params:
#         script_path = os.path.join(PROJ_DIR, "scripts", "split_maf.sh")
#     resources:
#         time="24:00:00"
#     shell:
#         """
#         source mafSplit  {input.maf} {input.outdir}
#         """

# ####################

# rule estimate_rho:
#     input:
#         mod_file = MOD_FILE,
#         chrome_maf = os.path.join(MAF_SPLIT_DIR, MAF_CHR_PREFIX + "{chrome}.00.maf")
#     output:
#         rho_cons_mod = os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_PREFIX + "{chrome}.cons.mod"),
#         rho_noncons_mod = os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_PREFIX + "{chrome}.noncons.mod")
#     params:
#         output_prefix = os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_PREFIX + "{chrome}"),
#         length = 45,
#         coverage = 0.3
#     log:
#         os.path.join(MAF_OUTDIR, "03-phastcons-rho", "logs", MAF_REF_PREFIX + "{chrome}.log")
#     resources:
#         mem = "48g",
#         time = "96:00:00",
# 	    partition = "holy-info"
#     shell:
#         """
#         phastCons --expected-length {params.length} --target-coverage {params.coverage} --no-post-probs --msa-format MAF --estimate-rho {params.output_prefix} {input.chrome_maf} {input.mod_file} &> {log}
#         """
#
## These rules require the entire maf to be read into memory, which is not feasible for large alignments like zoonomia


