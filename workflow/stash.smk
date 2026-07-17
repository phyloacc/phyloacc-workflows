#############################################################################
# Retired rules/helpers from workflow/phastcons_cnees.smk.
#
# This file is NOT included by the Snakefile - nothing here is part of the
# active pipeline, and none of it runs. It's kept only as a reference for
# what these rules used to do and why they stopped being used, in case any
# of the logic is useful again later.
#
# Everything below depends on names (MAF_SPLIT_NS_DIR, RHO_STATS_DIR,
# CONSERVE_DIR, PHYLOFIT_ACTIVE_MODEL_PATH, RHO_MODE, LOG_DIR,
# getRuleResources, rules.global_rho, checkpoints.filter_maf_by_gap, ...)
# that are only defined in workflow/phastcons_cnees.smk. It will not parse
# or run standalone - do not add an `include:` for this file.
#
# What happened: per-chunk phastCons scoring and rho estimation used to be
# their own Snakemake rules (one job per chunk). That's what the two rules
# below implement. At some point the same work got reimplemented inline as a
# Python loop (with its own thread pool) inside run_phastcons_chr (rho
# scoring) and global_rho (rho estimation, only under rho_mode: estimate),
# so nothing in the DAG depends on these rules' outputs anymore - they're
# fully orphaned. See RULES.md for the current, active rule list.
#############################################################################

# Unused even by the rule below it - never called anywhere in the pipeline.
def rho_values_for_chr(wc):
    ckpt = checkpoints.filter_maf_by_gap.get(
        chromosome_group=wc.chromosome_group,
        ref_chromosome=wc.ref_chromosome
    )
    filtered_manifest = ckpt.output.filtered_manifest
    outdir = os.path.dirname(filtered_manifest)

    rho_files = []
    with open(filtered_manifest) as mf:
        for line in mf:
            line = line.strip()
            if not line:
                continue
            chunk = os.path.splitext(os.path.basename(line))[0]
            rho_files.append(
                os.path.join(RHO_STATS_DIR, wc.chromosome_group, wc.ref_chromosome, f"{chunk}.rho.txt")
            )
    return rho_files

rule phastcons_estimate_rho_chunk:
    input:
        maf = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.maf"),
        mod = PHYLOFIT_ACTIVE_MODEL_PATH,
    output:
        stderr_file = os.path.join(
            LOG_DIR,
            "phastcons_estimate_rho_chunk",
            "{chromosome_group}",
            "{ref_chromosome}",
            "{chunk}.stderr.log",
        ),
        rho_value = os.path.join(RHO_STATS_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.rho.txt"),
    params:
        outdir = os.path.join(RHO_STATS_DIR, "{chromosome_group}", "{ref_chromosome}"),
        rho_prefix = os.path.join(RHO_STATS_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.rho", "rho"),
        rule_name = "phastcons_estimate_rho_chunk",
    log:
        job_log = os.path.join(LOG_DIR, "phastcons_estimate_rho_chunk", "{chromosome_group}", "{ref_chromosome}-{chunk}.log"),
    resources:
        **getRuleResources("phastcons_per_chunk")
    run:
        import glob, os, re, shlex, subprocess, traceback, math

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(params.outdir, exist_ok=True)

                # Count unique sequence IDs without invoking a shell.
                uniq_ids = 0
                with open(input.maf) as maf_stream:
                    uniq_ids = len({line.split()[1] for line in maf_stream if line.startswith("s ")})

                if uniq_ids < 2:
                    log_stream.write(f"SKIP: only {uniq_ids} unique sequence id(s) in {input.maf}\n")
                    with open(output.rho_value, "w") as rf:
                        rf.write("nan\n")
                    with open(output.stderr_file, "w") as sf:
                        sf.write(f"SKIP: only {uniq_ids} unique sequence id(s)\n")
                    return

                tmp_prefix = os.path.join(params.outdir, f".{wildcards.chunk}.estimate_rho.tmp")
                cmd = [
                    "phastCons",
                    input.maf,
                    input.mod,
                    "--estimate-rho",
                    tmp_prefix,
                ]

                log_stream.write(f"Running: {' '.join(shlex.quote(c) for c in cmd)}\n")
                log_stream.flush()

                result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
                stderr_text = result.stderr or ""
                with open(output.stderr_file, "w") as sf:
                    sf.write(stderr_text)
                if stderr_text:
                    log_stream.write(stderr_text)
                    log_stream.flush()

                if result.returncode != 0:
                    # Don't crash pipeline; keep stderr, but make other outputs empty
                    log_stream.write(f"FAIL: phastCons exit {result.returncode} on {input.maf}\n")
                    with open(output.rho_value, "w") as rf:
                        rf.write("nan\n")
                    return

                # Parse rho directly from phastCons stderr text.
                rho = float("nan")
                try:
                    m_all = re.findall(r"rho\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", stderr_text)
                    if m_all:
                        rho = float(m_all[-1])
                except Exception:
                    pass

                with open(output.rho_value, "w") as rf:
                    rf.write(f"{rho}\n")

                # Remove temporary files written by --estimate-rho prefix.
                for fp in glob.glob(tmp_prefix + "*"):
                    try:
                        os.remove(fp)
                    except OSError:
                        pass

            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule phastcons_chunk:
    input:
        maf = os.path.join(MAF_SPLIT_NS_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.maf"),
        mod = PHYLOFIT_ACTIVE_MODEL_PATH,
        rho_value = (lambda wc: os.path.join(RHO_STATS_DIR, wc.chromosome_group, wc.ref_chromosome, f"{wc.chunk}.rho.txt")
                     if RHO_MODE == "estimate"
                     else rules.global_rho.output.global_rho),
        global_rho = rules.global_rho.output.global_rho,
    output:
        bed = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.conserved.bed"),
        wig = temp(os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}", "{chunk}.scores.wig")),
        err = os.path.join(
            LOG_DIR,
            "phastcons_chunk",
            "{chromosome_group}",
            "{ref_chromosome}",
            "{chunk}.stderr.log",
        ),
    params:
        outdir = os.path.join(CONSERVE_DIR, "{chromosome_group}", "{ref_chromosome}"),
        rule_name = "phastcons_chunk",
    log:
        job_log = os.path.join(LOG_DIR, "phastcons_chunk", "{chromosome_group}", "{ref_chromosome}-{chunk}.log"),
    resources:
        **getRuleResources("phastcons_per_chunk")
    run:
        import os, shlex, subprocess, traceback, math

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(params.outdir, exist_ok=True)

                # Count unique sequence IDs without invoking a shell.
                uniq_ids = 0
                with open(input.maf) as maf_stream:
                    uniq_ids = len({line.split()[1] for line in maf_stream if line.startswith("s ")})

                if uniq_ids < 2:
                    log_stream.write(f"SKIP: only {uniq_ids} unique sequence id(s) in {input.maf}\n")
                    open(output.bed, "w").close()
                    open(output.wig, "w").close()
                    with open(output.err, "w") as ef:
                        ef.write(f"SKIP: only {uniq_ids} unique sequence id(s)\n")
                    return

                with open(input.global_rho) as gf:
                    global_str = gf.read().strip()
                try:
                    rho_global = float(global_str)
                except Exception:
                    rho_global = float("nan")

                if RHO_MODE == "estimate":
                    with open(input.rho_value) as rf:
                        rho_str = rf.read().strip()
                    try:
                        rho_chunk = float(rho_str)
                    except Exception:
                        rho_chunk = float("nan")
                else:
                    rho_str = global_str
                    rho_chunk = rho_global

                if not (rho_global > 0.0 and math.isfinite(rho_global)):
                    log_stream.write(f"SKIP: invalid global rho '{global_str}'\n")
                    open(output.bed, "w").close()
                    open(output.wig, "w").close()
                    with open(output.err, "w") as ef:
                        ef.write(f"SKIP: invalid global rho '{global_str}'\n")
                    return

                if not (rho_chunk > 0.0 and math.isfinite(rho_chunk)):
                    log_stream.write(f"SKIP: invalid chunk rho '{rho_str}' for {input.maf}\n")
                    open(output.bed, "w").close()
                    open(output.wig, "w").close()
                    with open(output.err, "w") as ef:
                        ef.write(f"SKIP: invalid chunk rho '{rho_str}'\n")
                    return

                if rho_chunk > rho_global:
                    log_stream.write(f"SKIP: rho_chunk {rho_chunk} > global_rho {rho_global}\n")
                    open(output.bed, "w").close()
                    open(output.wig, "w").close()
                    with open(output.err, "w") as ef:
                        ef.write(f"SKIP: rho_chunk {rho_chunk} > global_rho {rho_global}\n")
                    return

                # KEEP post probs (do NOT add --no-post-probs)
                cmd = [
                    "phastCons",
                    input.maf,
                    input.mod,
                    "--rho",
                    str(rho_global),
                    "--most-conserved",
                    output.bed,
                ]

                log_stream.write(f"Running: {' '.join(shlex.quote(c) for c in cmd)} > {output.wig} 2> {output.err}\n")
                log_stream.flush()

                with open(output.wig, "w") as wig_stream, open(output.err, "w") as err_stream:
                    result = subprocess.run(cmd, stdout=wig_stream, stderr=err_stream, text=True)

                if result.returncode != 0:
                    log_stream.write(f"FAIL: phastCons exit {result.returncode} on {input.maf}\n")
                    raise RuntimeError(f"phastCons failed with exit code {result.returncode} on {input.maf}")

            except Exception:
                traceback.print_exc(file=log_stream)
                raise
