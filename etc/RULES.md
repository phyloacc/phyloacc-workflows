# Rule reference: phyloFit + phastCons/CNEE workflows

One to two sentence description of every Snakemake rule/checkpoint in
`workflow/phylofit_models.smk` and `workflow/phastcons_cnees.smk`, in file
order. Written by reading each rule's actual `input`/`output`/`run` body,
not inferred from its name.

## `workflow/phylofit_models.smk`

- **`all`** *(only defined when this file runs standalone, not under the master workflow)* - Aggregator with no logic of its own; depends on the final phyloFit `.mod` file (GC-corrected or not, per `use_gc_corrected_models`) for every reference chromosome.

- **`maf_index`** - Runs `mafutils index` on the whole-alignment MAF, producing a block index and a scaffold index used for random access into it.

- **`make_group_beds`** - Runs `make_group_beds.py` against the whole-genome MAF's block index (from `maf_index`) to write a BED file listing the chromosome coordinates belonging to one `chromosome_group`. Chromosome lengths are derived from the block index (`max(ref_start + ref_len)` per scaffold) rather than a reference FASTA index, so this rule has no `ref_fasta` dependency regardless of `split_strategy`.

- **`maf_split_chr_by_group`** *(checkpoint)* - Runs `mafutils fetch` in scaffold mode to pull per-chromosome MAF blocks for a group out of the whole-alignment MAF, and writes a manifest listing the per-chromosome MAF files produced. It's a checkpoint because downstream rules resolve which per-chromosome file corresponds to a given chromosome by reading this manifest at runtime.

- **`ref_gff_split_by_chr`** - Runs the `ref_gff_split_by_chr.awk` script to extract one chromosome's annotation lines from the full reference GFF and prefix them to match the MAF's naming.

- **`extract_4d_codons_by_chr`** - Runs `msa_view --4d --features <gff>` on a chromosome's MAF to pull out fourfold-degenerate codon columns per the GFF annotation.

- **`extract_4d_sites`** - Runs `msa_view --tuple-size 1` on the 4d-codon alignment to collapse it down to single-site 4d positions.

- **`filter_4d_sites`** - Runs `filter_4d_sites.py` to drop sites/sequences not meeting the `filter_threshold_4d` presence threshold, producing a filtered alignment plus a summary TSV.

- **`run_phylofit`** - Runs `phyloFit --tree <tree> --msa-format SS` on the filtered 4d sites to fit a neutral substitution model for that chromosome (the uncorrected `.mod` file).

- **`get_gc_content`** - Runs `get_gc_content.py` against the sample sheet to compute per-sample and average GC content, used by the GC-correction step below.

- **`run_mod_freqs`** - Runs `modFreqs <mod_file> <gc_value>` to rescale the uncorrected phyloFit model's background base frequencies to the target GC content, producing the GC-corrected `.mod` file. This corrected model is what downstream phastCons rules actually use whenever `use_gc_corrected_models` (the default) is true.

## `workflow/phastcons_cnees.smk`

- **`all`** *(only defined when this file runs standalone)* - Aggregator with no logic of its own; depends on the conserved-region BED, and (if enabled) CNEE BED and CNEE-alignment-manifest targets, for every reference chromosome.

- **`ref_fasta_index`** *(guarded; `split_strategy: ns` or `fixed_windows` only - not scheduled under the default `num_seqs` strategy)* - Runs `samtools faidx` on the reference FASTA.

- **`ref_fasta_dict`** *(guarded; `split_strategy: ns` only - `fixed_windows` only needs the `.fai`, not this)* - Runs `picard CreateSequenceDictionary` on the reference FASTA, producing the `.dict` file Picard tools require.

- **`picard_scatter_by_ns`** *(`split_strategy: ns` only)* - Runs `picard ScatterIntervalsByNs --OUTPUT_TYPE ACGT` to list every non-N interval genome-wide, using `min_Ns_to_split_by` as the N-run threshold.

- **`ns_to_bed`** - Pulls one chromosome's intervals out of the genome-wide Picard interval file, writing `chrom:start-end` lines for just that chromosome.

- **`filter_ns_bed_minlen`** - Drops any interval shorter than `min_keep_region_len` from that per-chromosome interval file (`awk`).

- **`ns_minlen_to_bed3`** - Converts the length-filtered `chrom:start-end` lines to standard 0-based BED3 (`awk`) - this is the chunk-boundary BED used to split the chromosome's MAF by N-runs.

- **`fixed_windows_bed`** *(`split_strategy: fixed_windows` only)* - Pure-Python: tiles a chromosome (using its `.fai` length) into fixed-size, optionally-overlapping windows of `window_size_bp`/`window_overlap_bp`, as an alternative to N-run-based chunking.

- **`maf_index_chr`** - Runs `mafutils index` on a chromosome-level MAF, producing a block index and scaffold index for fast random access.

- **`num_seqs_chunk_bed_chr`** *(`split_strategy: num_seqs` only, the default)* - Pure-Python: reads the chromosome's block index (from `maf_index_chr`) and finds runs of blocks with `num_seqs <= num_seqs_max_for_gap` (few aligned species) at least `num_seqs_min_gap_bp` long, treating those as gaps; the complement becomes the chunk-boundary BED, with any resulting chunk shorter than `num_seqs_min_keep_region_len` dropped. Needs no reference FASTA - splits based on cross-species alignment coverage instead of reference-genome N-runs.

- **`maf_split_chunks`** - Runs `mafutils fetch -m block` to pull MAF sub-blocks out for each chunk-boundary BED interval (whichever `split_strategy` is active - `ns`, `num_seqs`, or `fixed_windows`), and writes a manifest of the chunk MAF files produced.

- **`filter_maf_by_gap`** *(checkpoint)* - Computes each chunk's non-reference gap fraction and keeps only chunks at or below `max_gap_pct`, writing a filtered manifest. A checkpoint because downstream per-chunk rules re-resolve their chunk list from this output.

- **`phastcons_estimate_rho_chunk`** - Runs `phastCons --estimate-rho` on one filtered chunk and parses the estimated rho out of its stderr. **Orphaned**: nothing in the file references this rule's output (`rules.phastcons_estimate_rho_chunk...` does not appear anywhere), and the same per-chunk estimation is done again independently, inline, inside `global_rho` below - Snakemake never actually schedules this rule.

- **`global_rho`** - In `estimate` mode: loops over every filtered chunk, itself running `phastCons --estimate-rho` per chunk (duplicating what `phastcons_estimate_rho_chunk` does, without depending on it), then computes mean/median/p90 across chunks and picks one per `global_rho_stat` (default p90) as the chromosome-wide rho. In `fixed` mode: skips estimation entirely and just uses the configured `fixed_rho` constant. Either way, writes the chosen chromosome-wide rho value to its output file.

- **`run_phastcons_chr`** - Loops over every filtered chunk for a chromosome (in a thread pool) and runs `phastCons --rho <global_rho> --most-conserved` on each, producing a per-chunk conserved-regions BED and scores WIG; only emits a `chunks.done` sentinel as its formal Snakemake output. Always applies the one chromosome-wide rho value to every chunk, in both rho modes.

- **`phastcons_concat_chr`** - Waits on `chunks.done`, globs all per-chunk `*.conserved.bed` files phastCons wrote, sorts and concatenates them into one chromosome-level conserved-regions BED. If `cleanup_chunk_intermediates` is set, deletes the per-chunk BED/WIG files and the chunked-MAF directory afterward to reclaim disk space.

- **`extract_cds_bed_chr`** - Pure Python: filters the reference GFF down to `CDS` rows on one chromosome and converts them to BED coordinates.

- **`cnees_from_conserved_chr`** - Pure Python: merges nearby conserved intervals within `cnee_ces_merge_gap_bp`, merges the CDS intervals, and drops any conserved element that overlaps a CDS at all (no partial/flanking fragments kept) - i.e. turns "conserved regions" into "conserved, non-coding regions" (CNEE candidates).

- **`cnees_to_bed4_chr`** - Drops any CNEE candidate shorter than `cnee_min_len_bp` and assigns each survivor a sequential ID (`{chromosome}.cnee{n:07d}`), writing a BED4 (chrom/start/end/id) of the final CNEE set.

- **`cnee_alignments_chr`** *(renamed from `cnee_mafs_chr`, which no longer accurately described it once FASTA output was added)* - Runs `mafutils fetch -m block -b id`, keyed by each CNEE's BED4 ID, to extract one alignment file per conserved element from the chromosome MAF - FASTA (with header style and expected-species list) when `cnee_output_format: fasta`, raw MAF blocks otherwise. In FASTA mode, also deletes any per-element file containing two sequences from the same species. Writes a manifest of the surviving per-element files.

- **`phastcons_chunk`** - Runs `phastCons --rho <rho> --most-conserved` on a single chunk, with its own skip logic for invalid/oversized rho. **Orphaned**: its output path is identical to what `run_phastcons_chr` already writes per-chunk internally, nothing references `rules.phastcons_chunk...` anywhere, and it isn't part of any target list - Snakemake never actually schedules this rule either.

## Notes for whoever reads this next

- **Two orphaned rules exist**: `phastcons_estimate_rho_chunk` and `phastcons_chunk` are fully defined but never reachable from any target - their real work happens inline inside `global_rho` and `run_phastcons_chr` instead. They still get read/parsed by Snakemake (and their names are borrowed as the shared `rule_resources` key for the real per-chunk-ish rules), but never scheduled as jobs. Worth deleting or wiring back in - a separate decision from the `cnee_mafs_chr` -> `cnee_alignments_chr` rename made above.
- `PHYLOFIT_ACTIVE_MODEL_PATH` (the neutral model every phastCons rule consumes) is referenced throughout `phastcons_cnees.smk` but defined in `phylofit_models.smk` - only visible when both files are included together by the master `Snakefile`.
- `ref_fasta_index`/`ref_fasta_dict` are only defined in `phastcons_cnees.smk` now (an identical `ref_fasta_index` copy used to exist in `phylofit_models.smk` too, guarded against duplicate registration - removed once `make_group_beds` stopped needing it). `ref_fasta` itself is only required at all when `split_strategy` is `ns` or `fixed_windows`; the default `num_seqs` strategy and `make_group_beds`'s block-index-based chromosome lengths need no reference FASTA anywhere in the pipeline.
