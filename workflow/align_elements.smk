#############################################################################
# Align-only entry point.
#
# Bypasses the entire MAF-based prediction pipeline: take a directory of
# UNALIGNED per-element FASTA files (one multi-FASTA per element, one record
# per species) and align each one with an aligner (mafft by default), producing
# a directory of per-element aligned FASTAs (mirroring the input's dir layout) + a
# manifest.txt listing the aligned elements as paths relative to align_input_dir - the
# same output shape as the CNEE terminus in workflow/cnees.smk.
#
# This module is included ONLY when align_input_dir is set (see Snakefile). It
# never coexists with the MAF modules (those are all switched off in align-only
# mode), so nothing here reads maf/ref_gff/tree_file/ref_chromosome_groups.
#
# Scatter-gather: the element list is partitioned into batches of
# align_batch_size elements (round-robin, so batches stay size-balanced); one
# align_batch job aligns each batch, running the aligner over its elements in a
# thread pool sized by cpus_per_task; align_gather concatenates the per-batch
# manifests into the final manifest.
#############################################################################

import os
import glob
import math

from functools import partial

import lib.common as COMMON

_setup = config["__pipeline_setup__"]
OUTPUT_DIR = _setup["OUTPUT_DIR"]
LOG_DIR = _setup["LOG_DIR"]

getRuleResources = partial(COMMON.getResources, config)

#############################################################################
# Config / paths

ALIGN_INPUT_DIR = os.path.abspath(str(config["align_input_dir"]).strip())
if not os.path.isdir(ALIGN_INPUT_DIR):
    raise ValueError(f"align_input_dir '{ALIGN_INPUT_DIR}' is not an existing directory.")

ALIGN_OUTPUT_DIR = COMMON.getOptionalConfigPath(
    config, "align_output_dir", os.path.join(OUTPUT_DIR, "aligned-elements"),
)
ALIGNED_DIR = os.path.join(ALIGN_OUTPUT_DIR, "aligned")

# Which unaligned-element files to pick up. A str is accepted and wrapped.
ALIGN_INPUT_GLOBS = config.get("align_input_glob", ["*.fa", "*.fasta", "*.fna"])
if isinstance(ALIGN_INPUT_GLOBS, str):
    ALIGN_INPUT_GLOBS = [ALIGN_INPUT_GLOBS]

# Also search one level down in subdirectories matching this glob (e.g. the "batchN"
# dirs the CNEE extractor writes, 10k files each). Blank = flat top-level only. This is
# NOT recursive - only the top level and one level of matching subdirs are searched.
ALIGN_BATCH_SUBDIR_GLOB = str(config.get("align_batch_subdir_glob", "batch[0-9]*")).strip()

# Elements with fewer than this many sequences cannot be aligned -> skipped with a
# warning (never fabricated as output). Raise it to also drop sparse elements.
ALIGN_MIN_SEQS = int(config.get("align_min_seqs", 2))
if ALIGN_MIN_SEQS < 1:
    raise ValueError("align_min_seqs must be >= 1.")

# Drop an element whose LONGEST ungapped record is shorter than this (bp); 0 = off.
# Default 50 (matches the pipeline's cnee_min_len_bp). Skipped elements are reported in
# skipped.tsv (reason "too-short"), never aligned.
ALIGN_MIN_LEN_BP = int(config.get("align_min_len_bp", 50))
if ALIGN_MIN_LEN_BP < 0:
    raise ValueError("align_min_len_bp must be >= 0.")

# Optional output-header trimming: if align_header_split is set, each aligned record's
# header is split on that character and field align_header_field (0-based; negatives
# index from the end) is kept - e.g. split " ", field 0 turns ">gGal ce1 coords" into
# ">gGal". Blank = keep headers as-is. Out-of-range field -> keep the full header.
_align_header_split = config.get("align_header_split", None)
ALIGN_HEADER_SPLIT = None if _align_header_split in (None, "") else str(_align_header_split)
ALIGN_HEADER_FIELD = int(config.get("align_header_field", 0))

ALIGN_BATCH_SIZE = int(config.get("align_batch_size", 2000))
if ALIGN_BATCH_SIZE < 1:
    raise ValueError("align_batch_size must be >= 1.")

# Aligner is a config-driven command template so muscle/prank etc. are a data change,
# not a code change. Tokens {input}/{output} are substituted per element; if the
# template contains {output} it is passed as an argument (muscle/prank style),
# otherwise the aligner's stdout is captured to the output file (mafft style).
ALIGNER = str(config.get("aligner", "mafft")).strip().lower()
ALIGNER_COMMANDS = config.get("aligner_commands", {"mafft": ["mafft", "--auto", "{input}"]})
if ALIGNER not in ALIGNER_COMMANDS:
    raise ValueError(
        f"Unknown aligner '{ALIGNER}'. Known aligners: {sorted(ALIGNER_COMMANDS)}. "
        f"Add a command template under aligner_commands to use another."
    )
ALIGNER_TEMPLATE = list(ALIGNER_COMMANDS[ALIGNER])
USES_OUTPUT_TOKEN = any("{output}" in str(tok) for tok in ALIGNER_TEMPLATE)

#############################################################################
# Element discovery + batching (static glob at parse time - inputs are
# user-supplied and exist before the run, so no checkpoint is needed).

# Search the top level plus one level of matching "batchN" subdirs (see
# ALIGN_BATCH_SUBDIR_GLOB). A flat input dir -> the batch globs match nothing; a batched
# tree -> the top-level globs match nothing. Not recursive.
_search_globs = [os.path.join(ALIGN_INPUT_DIR, pat) for pat in ALIGN_INPUT_GLOBS]
if ALIGN_BATCH_SUBDIR_GLOB:
    _search_globs += [
        os.path.join(ALIGN_INPUT_DIR, ALIGN_BATCH_SUBDIR_GLOB, pat) for pat in ALIGN_INPUT_GLOBS
    ]
ELEMENT_FILES = sorted(
    set(f for g in _search_globs for f in glob.glob(g) if os.path.isfile(f))
)

# Each aligned output mirrors its input's path RELATIVE to align_input_dir: a flat input
# dir -> flat output; batched input -> the same batchN subdirs recreated under aligned/.
# Because every input has a unique relative path, output paths (and manifest entries) are
# unique by construction - no basename-collision check is needed.

N_BATCHES = math.ceil(len(ELEMENT_FILES) / ALIGN_BATCH_SIZE) if ELEMENT_FILES else 0
# Round-robin assignment keeps batch sizes balanced even if some elements are large.
BATCHES = {i: ELEMENT_FILES[i::N_BATCHES] for i in range(N_BATCHES)}

# Per-SLURM-batch intermediates live under "run-info/" (deliberately NOT "batches/" - that
# would collide conceptually with the input tool's own batchN dirs; our batches are just how
# we split the work across SLURM jobs). The final manifest/skipped sit at the top level.
BATCH_MANIFEST = os.path.join(ALIGN_OUTPUT_DIR, "run-info", "{batch}.manifest.txt")
BATCH_SKIPPED = os.path.join(ALIGN_OUTPUT_DIR, "run-info", "{batch}.skipped.tsv")
FINAL_MANIFEST = os.path.join(ALIGN_OUTPUT_DIR, "manifest.txt")
FINAL_SKIPPED = os.path.join(ALIGN_OUTPUT_DIR, "skipped.tsv")

# Header for the aggregated skipped table (per-batch sidecars are headerless rows).
SKIPPED_HEADER = "locus\treason\tn_seqs\tmax_len_bp"

# manifest.txt is the sole declared target; skipped.tsv is co-produced by align_gather.
ALL_ALIGN_TARGETS = [FINAL_MANIFEST]

wildcard_constraints:
    batch = r"\d+"

#############################################################################

def _batch_inputs(wildcards):
    return BATCHES[int(wildcards.batch)]

def _trim_headers(path, split, field):
    """Rewrite each FASTA header in `path` in place: split the text after '>' on `split`
    and keep field `field` (0-based; negatives index from the end). If the field is out
    of range for a header, keep that header unchanged. Sequence lines are untouched."""
    out_lines = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                head = line[1:].rstrip("\n")
                parts = head.split(split)
                try:
                    head = parts[field].strip()
                except IndexError:
                    pass  # out-of-range -> keep the full original header
                out_lines.append(">" + head + "\n")
            else:
                out_lines.append(line)
    with open(path, "w") as out:
        out.writelines(out_lines)

rule align_batch:
    input:
        elements = _batch_inputs
    output:
        manifest = BATCH_MANIFEST,
        skipped = BATCH_SKIPPED
    # NOTE: the aligner template / dirs are NOT passed via params - a params string
    # containing "{input}"/"{output}" would be misread by Snakemake as a wildcard to
    # expand. The run block reads the module globals (ALIGNED_DIR, ALIGNER_TEMPLATE,
    # USES_OUTPUT_TOKEN, ALIGN_MIN_SEQS, ALIGNER) directly instead.
    log:
        job_log = os.path.join(LOG_DIR, "align_elements", "batch_{batch}.log")
    benchmark:
        os.path.join(LOG_DIR, "benchmarks", "align_elements", "batch_{batch}.txt")
    resources:
        **getRuleResources("align_batch")
    run:
        import os, subprocess, traceback
        from concurrent.futures import ThreadPoolExecutor, as_completed

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(ALIGNED_DIR, exist_ok=True)
                os.makedirs(os.path.dirname(output.manifest), exist_ok=True)
                elem_log_dir = os.path.join(LOG_DIR, "align_elements", "elements")
                os.makedirs(elem_log_dir, exist_ok=True)

                elements = list(input.elements)
                n = len(elements)
                max_workers = max(1, min(int(resources.cpus_per_task), n if n else 1))
                log_stream.write(
                    f"Batch {wildcards.batch}: aligning {n} element(s) with aligner="
                    f"'{ALIGNER}' parallel_workers={max_workers} "
                    f"(cpus_per_task={int(resources.cpus_per_task)})\n"
                )
                log_stream.flush()

                def measure_element(path):
                    # One pass: number of records + the longest ungapped record length (bp).
                    nrec = 0
                    cur = 0
                    max_len = 0
                    with open(path) as fh:
                        for line in fh:
                            if line.startswith(">"):
                                if cur > max_len:
                                    max_len = cur
                                nrec += 1
                                cur = 0
                            else:
                                cur += len(line.strip())
                    if cur > max_len:
                        max_len = cur
                    return nrec, max_len

                def align_one(infile):
                    # Mirror the input's path RELATIVE to align_input_dir: a flat input dir
                    # gives flat output; a batched input recreates the same batchN subdir
                    # under aligned/ (and the per-element stderr logs). rel is also the
                    # manifest/skipped identity, so it is unique per input by construction.
                    rel = os.path.relpath(infile, ALIGN_INPUT_DIR)
                    outfile = os.path.join(ALIGNED_DIR, rel)
                    err = os.path.join(elem_log_dir, rel + ".stderr.log")
                    os.makedirs(os.path.dirname(err), exist_ok=True)
                    # Results are (rel, status, rc, detail, n_seqs, max_len). n_seqs/max_len
                    # come from measure_element and feed the structured skipped.tsv columns.
                    try:
                        nrec, max_len = measure_element(infile)
                    except OSError as e:
                        with open(err, "w") as ef:
                            ef.write(f"FAIL: could not read input: {e}\n")
                        return rel, "fail-read", 1, f"could not read input: {e}", 0, 0
                    # No records at all -> empty/malformed input -> hard failure.
                    if nrec == 0:
                        with open(err, "w") as ef:
                            ef.write("FAIL: no FASTA records (empty or malformed input)\n")
                        return rel, "fail-empty", 1, "no FASTA records (empty or malformed)", 0, 0
                    # Too few sequences to align -> skip with a warning, no output.
                    if nrec < ALIGN_MIN_SEQS:
                        reason = f"{nrec} sequence(s) < align_min_seqs={ALIGN_MIN_SEQS}"
                        with open(err, "w") as ef:
                            ef.write(f"SKIP: {reason}; cannot align\n")
                        return rel, "skip-fewseqs", 0, reason, nrec, max_len
                    # Element too short (longest record below the bp threshold) -> skip.
                    if ALIGN_MIN_LEN_BP and max_len < ALIGN_MIN_LEN_BP:
                        reason = f"max record {max_len}bp < align_min_len_bp={ALIGN_MIN_LEN_BP}"
                        with open(err, "w") as ef:
                            ef.write(f"SKIP: {reason}\n")
                        return rel, "skip-short", 0, reason, nrec, max_len

                    os.makedirs(os.path.dirname(outfile), exist_ok=True)
                    cmd = [
                        str(tok).replace("{input}", infile).replace("{output}", outfile)
                        for tok in ALIGNER_TEMPLATE
                    ]
                    try:
                        if USES_OUTPUT_TOKEN:
                            with open(err, "w") as ef:
                                proc = subprocess.run(cmd, stderr=ef, text=True)
                        else:
                            with open(outfile, "w") as of, open(err, "w") as ef:
                                proc = subprocess.run(cmd, stdout=of, stderr=ef, text=True)
                        rc = proc.returncode
                    except OSError as e:
                        with open(err, "a") as ef:
                            ef.write(f"FAIL: could not launch aligner: {e}\n")
                        return rel, "fail-exec", 1, f"could not launch aligner: {e}", nrec, max_len

                    if rc != 0:
                        # Don't leave a half-written alignment behind on failure.
                        if not USES_OUTPUT_TOKEN:
                            try:
                                os.remove(outfile)
                            except OSError:
                                pass
                        return rel, "fail-aligner", rc, f"aligner exited {rc}", nrec, max_len

                    # Optional: rewrite the aligned output's headers (e.g. down to species tag).
                    if ALIGN_HEADER_SPLIT is not None:
                        try:
                            _trim_headers(outfile, ALIGN_HEADER_SPLIT, ALIGN_HEADER_FIELD)
                        except OSError as e:
                            with open(err, "a") as ef:
                                ef.write(f"FAIL: could not trim headers: {e}\n")
                            return rel, "fail-trim", 1, f"could not trim headers: {e}", nrec, max_len
                    return rel, "ok", 0, "", nrec, max_len

                # Short reason codes for the skipped.tsv "reason" column.
                SKIP_REASON = {"skip-fewseqs": "too-few-seqs", "skip-short": "too-short"}
                aligned, skipped, failed = [], [], []

                def handle(res):
                    # name = input path relative to align_input_dir
                    name, status, rc, detail, nseq, maxlen = res
                    if status == "ok":
                        aligned.append(name)
                        log_stream.write(f"OK    {name}\n")
                    elif status.startswith("skip"):
                        # Structured row for skipped.tsv: (locus, reason, n_seqs, max_len_bp).
                        skipped.append((name, SKIP_REASON.get(status, status), nseq, maxlen))
                        log_stream.write(f"SKIP  {name} ({detail})\n")
                    else:
                        failed.append((name, status, rc, detail))
                        log_stream.write(f"FAIL  {name} [{status} rc={rc}: {detail}]\n")
                    log_stream.flush()

                if max_workers == 1:
                    for infile in elements:
                        handle(align_one(infile))
                else:
                    with ThreadPoolExecutor(max_workers=max_workers) as pool:
                        futs = {pool.submit(align_one, f): f for f in elements}
                        for fut in as_completed(futs):
                            handle(fut.result())

                if failed:
                    raise RuntimeError(
                        f"{len(failed)} element(s) failed to align in batch "
                        f"{wildcards.batch}; first failure: {failed[0]}"
                    )

                # Manifest lists only the elements this batch actually aligned, as paths
                # relative to align_input_dir (skipped ones go to the sidecar skipped list).
                with open(output.manifest, "w") as out:
                    for name in sorted(aligned):
                        out.write(name + "\n")
                # Per-batch skipped rows (headerless TSV: locus, reason, n_seqs, max_len_bp);
                # align_gather concatenates these and prepends the header.
                with open(output.skipped, "w") as out:
                    for name, reason, nseq, maxlen in sorted(skipped):
                        out.write(f"{name}\t{reason}\t{nseq}\t{maxlen}\n")
                log_stream.write(
                    f"Batch {wildcards.batch} done: {len(aligned)} aligned, "
                    f"{len(skipped)} skipped, {len(failed)} failed\n"
                )
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

####################

rule align_gather:
    input:
        batch_manifests = expand(BATCH_MANIFEST, batch=range(N_BATCHES)),
        batch_skipped = expand(BATCH_SKIPPED, batch=range(N_BATCHES))
    output:
        manifest = FINAL_MANIFEST,
        skipped = FINAL_SKIPPED
    params:
        rule_name = "align_gather"
    log:
        job_log = os.path.join(LOG_DIR, "align_elements", "gather.log")
    resources:
        **getRuleResources("align_gather")
    run:
        import os, traceback

        with open(log.job_log, "w") as log_stream:
            try:
                os.makedirs(os.path.dirname(output.manifest), exist_ok=True)

                # Aligned elements -> manifest.txt (paths relative to align_input_dir).
                names = []
                for bm in sorted(input.batch_manifests):
                    with open(bm) as fh:
                        for line in fh:
                            s = line.strip()
                            if s:
                                names.append(s)
                names = sorted(set(names))
                with open(output.manifest, "w") as out:
                    for nm in names:
                        out.write(nm + "\n")

                # Skipped elements -> skipped.tsv (tab-delimited, header + one row per
                # dropped locus: locus, reason, n_seqs, max_len_bp).
                skips = []
                for sf in sorted(input.batch_skipped):
                    with open(sf) as fh:
                        for line in fh:
                            s = line.rstrip("\n")
                            if s:
                                skips.append(s)
                skips = sorted(set(skips))
                with open(output.skipped, "w") as out:
                    out.write(SKIPPED_HEADER + "\n")
                    for s in skips:
                        out.write(s + "\n")

                log_stream.write(
                    f"Wrote {len(names)} aligned element(s) to {output.manifest} and "
                    f"{len(skips)} skipped element(s) to {output.skipped} "
                    f"from {len(input.batch_manifests)} batch(es)\n"
                )
            except Exception:
                traceback.print_exc(file=log_stream)
                raise

#############################################################################
# Standalone entry (only when NOT run under the master Snakefile).

if not config.get("__master_workflow__", False):
    rule all:
        input:
            ALL_ALIGN_TARGETS
