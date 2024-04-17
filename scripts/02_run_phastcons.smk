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
# print(MAF_PATH);
# print(REF_SPLIT_BED);
# sys.exit();

#############################################################################

rule all:
    input:
        expand(os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{chrome}.{win}.cons.mod"), zip, chrome=chr_list, win=win_list),
        expand(os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{chrome}.{win}.noncons.mod"), zip, chrome=chr_list, win=win_list),    
        # estimate_rho

        # expand(os.path.join(MAF_SPLIT_DIR, "{chrome}.{win}.maf"), zip, chrome=chr_list, win=win_list),
        # split_maf

        os.path.join(MAF_OUTDIR, "04-phyloboot-mods", MAF_REF_ID + ".avg.cons.mod"),
        os.path.join(MAF_OUTDIR, "04-phyloboot-mods", MAF_REF_ID + ".avg.noncons.mod"),
        # run_phyloboot
        
        os.path.join(MAF_OUTDIR, "05-phastcons", MAF_REF_ID + ".phastCons.wig"),
        os.path.join(MAF_OUTDIR, "05-phastcons", MAF_REF_ID + ".phastCons.bed")
        # run_phastcons
        

 
        
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

rule split_maf_windows:
    input:
        maf = MAF_PATH,
        split_bed = REF_SPLIT_BED
    output:
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
        mkdir -p {params.prefix}
        mafSplit {input.split_bed} {params.prefix} {input.maf}
        """
# On holybioinf:
# Splitting 1 files from /n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-sicb/data/00-human-ref-ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.chr.1Mb.bed to /n/holyscratch01/informatics/gwct/241-mammalian-2020v2b-mafSplit-1Mb/
# splitting /n/holyscratch01/informatics/gwct/241-mammalian-2020v2b.maf
# 79436.85user 11923.47system 25:28:11elapsed 99%CPU (0avgtext+0avgdata 7304maxresident)k
# 17477696inputs+20200869936outputs (0major+9728minor)pagefaults 0swaps

####################

def get_time(wildcards, attempt):
    if attempt == 1:
        return "4:0:00";
    elif attempt == 2:
        return "0:02:00";

def get_mem(wildcards, attempt):
    if attempt == 1:
        return "24g";
    elif attempt == 2:
        return "1g";
# Changing resources based on attempt
# Attempt 1 tries to run phastCons while attempt 2 simply checks why attempt 1 failed

rule estimate_rho:
    input:
        mod_file = MOD_FILE,
        chrome_maf = os.path.join(MAF_SPLIT_DIR, "{chrome}.{win}.maf")
    output:
        rho_cons_mod = os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{chrome}.{win}.cons.mod"),
        rho_noncons_mod = os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{chrome}.{win}.noncons.mod"),
    params:
        output_prefix = os.path.join(MAF_OUTDIR, "03-phastcons-rho", MAF_REF_ID + ".{chrome}.{win}"),
        length = 45,
        coverage = 0.3,
        job_id_dir = os.path.join(MAF_OUTDIR, "03-phastcons-rho", "job-ids"),
        job_id_file = os.path.join(MAF_OUTDIR, "03-phastcons-rho", "job-ids", MAF_REF_ID + ".{chrome}.{win}.jobid"),
        err_file = os.path.join(MAF_OUTDIR, "03-phastcons-rho", "logs", "errors.log"),
        timeout_file = os.path.join(MAF_OUTDIR, "03-phastcons-rho", "logs", "timeouts.log")
    log:
        os.path.join(MAF_OUTDIR, "03-phastcons-rho", "logs", MAF_REF_ID + ".{chrome}.{win}.log")
    resources:
        mem = lambda wildcards, attempt: get_mem(wildcards, attempt),
        time = lambda wildcards, attempt: get_time(wildcards, attempt),
	    partition = "holy-info,shared",
        cur_attempt = lambda wildcards, attempt: attempt
    retries: 
        1
    run:
        print("CUR ATTEMPT:", resources.cur_attempt);
        if resources.cur_attempt == 1:
            shell(
                    """
                    echo $SLURM_JOB_ID
                    mkdir -p {params.job_id_dir}
                    echo $SLURM_JOB_ID > {params.job_id_file}
                    phastCons --expected-length {params.length} --target-coverage {params.coverage} --no-post-probs --msa-format MAF --estimate-rho {params.output_prefix} {input.chrome_maf} {input.mod_file} &> {log}
                    """
                )
        # Try to run phastCons in attempt 1, but also write the SLURM job id to a file to look up in case of failure
        elif resources.cur_attempt == 2:
            shell(
                    """
                    if [ -f {params.job_id_file} ]; then 
                        echo "Job ID file exists"
                        job_id=$(cat {params.job_id_file})

                        if sacct -j $job_id | grep TIMEOUT; then
                            echo {wildcards.chrome}.{wildcards.win} $job_id >> {params.timeout_file}

                            touch {output.rho_cons_mod}
                            touch {output.rho_noncons_mod}
                        else
                            echo "Job failed with errors. Check log files for {wildcards.chrome}.{wildcards.win} for errors."
                            echo {wildcards.chrome}.{wildcards.win} $job_id >> {params.err_file}
                        fi                        
                    else 
                        echo "Job ID file does not exist. Check log files for {wildcards.chrome}.{wildcards.win} for errors."
                        echo {wildcards.chrome}.{wildcards.win} >> {params.err_file}
                    fi
                    """
                    )
        # If attempt 1 fails, check the job's status in attempt 2. If it timed out, create empty output files. If it failed with errors, write the job ID to a file to look up later
# The estimate_rho rule times out on some jobs, and we want to continue with the pipeline by creating empty output files
# This is accomplished by specifying that this job can be retried once
# During the first attempt, we write the job ID to a file and try to run phastCons, and then, if that fails, check the job's status in the second attempt
#
# real 419894.88
# user 444.03
# sys 134.16


####################

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
        find {params.rho_mod_dir} -name "*.cons.mod" ! -empty > {params.cons_mod_list_file}
        phyloBoot --read-mods '*{params.cons_mod_list_file}' --output-average {output.avg_cons_mod}
        find {params.rho_mod_dir} -name "*.noncons.mod" ! -empty > {params.noncons_mod_list_file}
        phyloBoot --read-mods '*{params.noncons_mod_list_file}' --output-average {output.avg_noncons_mod}
        """

####################

rule run_phastcons:
    input:
        avg_cons_mod = os.path.join(MAF_OUTDIR, "04-phyloboot-mods", MAF_REF_ID + ".avg.cons.mod"),
        avg_noncons_mod = os.path.join(MAF_OUTDIR, "04-phyloboot-mods", MAF_REF_ID + ".avg.noncons.mod"),
        maf = MAF_PATH
    output:
        wig = os.path.join(MAF_OUTDIR, "05-phastcons", MAF_REF_ID + ".phastCons.wig"),
        bed = os.path.join(MAF_OUTDIR, "05-phastcons", MAF_REF_ID + ".phastCons.bed")
    params:
        length = 45,
        coverage = 0.3,
        rho = 0.4
    log:
        os.path.join(MAF_OUTDIR, "05-phastcons", "logs", MAF_REF_ID + ".phastCons.log")
    resources:
        mem = "48g",
        time = "72:00:00",
        partition = "holy-info,shared"
    shell:
        """
        phastCons --expected-length {params.length} --target-coverage {params.coverage} --most-conserved {output.bed} --msa-format MAF {input.maf} {input.avg_cons_mod},{input.avg_noncons_mod} > {output.wig} 2> {log}
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


