# Tier-3 integration test fixture

`test_pipeline_e2e.py` runs the real pipeline (real `mafutils`/`phyloFit`/`phastCons`
calls, not mocked) end-to-end against a small real-data fixture in `data/`. Nothing in
the fixture is fabricated - every file is a genuine excerpt of real data already used
elsewhere in this repo, extracted/shifted, never invented.

## Provenance

Source files:

- MAF: `data/hamsters/workflow-tests/maf/cricetid-15spec.Mus_musculus.nodupes.maf`
  (real hamster whole-genome alignment, already relied on throughout this repo's own
  workflow tests)
- GFF: the real `ref_gff` used in those same workflow tests
- Tree: `/n/holylfs05/LABS/informatics/Everyone/support/20251114-raftrey-hamster-wga/data/tree/cricetid-15spec.tre`
  (real 15-species tree)

Window: reference chromosome `CM000994.3`, original genome coordinates
`3,030,000-3,290,000` (260,000 bp, 0-based half-open). Chosen to contain, in one
contiguous window:

- a real assembly-gap region (~3,030,000-3,040,000, `num_seqs=1` - no other species
  aligned there), so the pipeline's `num_seqs`-based chunk-splitting logic
  (`lib.intervals.complement_gaps`) has genuine structure to find. Without this, the
  window is uniformly well-covered and `mafutils fetch`'s block-mode extraction refuses
  to run (its internal QC heuristic rejects an extraction spanning ~100% of a scaffold -
  an edge case that never arises on a real, large chromosome).
- a real CDS-overlapping conserved element (~3,286,244-3,287,231 in original
  coordinates, overlapping a real CDS at 3,286,244-3,287,191) - a genuine drop case for
  the CDS-overlap logic ([[project_cnee_cds_overlap_question]]).
- ~25 other real conserved elements with zero CDS overlap - genuine keep cases.

### Extraction steps

1. `mafutils fetch -m scaffold` (or block, against a single-window bed3) extracted the
   `CM000994.3:3,030,000-3,290,000` window, all 15 species, from the source MAF.
2. The reference species' (`Mus_musculus.CM000994.3`) `s`-line `start` was shifted by
   `-3,030,000` and `srcSize` set to `260000`, so the extracted window behaves as a
   self-contained mini-chromosome starting at 0. This is necessary because
   `mafutils fetch` preserves original absolute genome coordinates in its output, but
   `make_group_beds.py` derives chromosome extent as `0` to `max(ref_start + ref_len)` -
   for a mid-chromosome window that produces a bogus, huge "chromosome" spanning from
   position 0 instead of the real ~260kb window. Only the reference species' coordinates
   were shifted; all other species' lines (their own genomes' unrelated coordinate
   spaces) and all aligned sequence/gap content are untouched.
3. The real GFF's rows for `CM000994.3` falling in the same window were extracted and
   shifted by the same `-3,030,000` offset, so CDS coordinates line up with the shifted
   MAF. No feature/gene content was altered, only start/end coordinates.
4. The tree file was copied byte-for-byte, unmodified (already real, already tiny, all
   15 species match the MAF as-is).
5. `mafutils index` was rerun on the shifted MAF to regenerate `.block.idx`/`.scaffold.idx`.

### Files in `data/`

- `cricetid-window.maf` - the shifted, extracted MAF (see header comment for the exact
  offset)
- `cricetid-window.maf.block.idx`, `cricetid-window.maf.scaffold.idx` - real indices
  rebuilt via `mafutils index`
- `cricetid-window.gff` - the shifted, extracted GFF rows
- `cricetid-15spec.tre` - the real, unmodified tree file
- `config.yaml` - fixture pipeline config (see its own header comment)
- `reference_summary.json` - the "silver standard" reference: real summary numbers
  (`ces_raw`, `ces_merged`, `ces_dropped_cds_overlap`, `cnees_after_cds_drop`) captured
  from an actual successful run of this exact fixture through the real pipeline, via
  `capture_reference.py`. `test_pipeline_e2e.py`'s
  `test_silver_standard_no_dramatic_drift` compares each future run's own fresh numbers
  against this reference and issues a `warnings.warn(...)` (not a failure) if any value
  drifts by more than 50% relative change - real tool output isn't guaranteed stable
  across `phastCons`/`mafutils` versions (evidenced directly this session by a real
  `mafutils` index-format change), so this is a soft drift check, not a hard snapshot
  comparison.

## Regenerating the fixture

If the source MAF/GFF/tree change, or a different window is needed, redo the extraction
steps above, then reseed the reference summary:

```
conda activate phyloacc-workflows2   # or any env with mafutils/phyloFit/phastCons on PATH
python3 tests/integration/capture_reference.py
```

This overwrites `reference_summary.json` with a fresh real run's numbers. Only run this
deliberately - it changes what future test runs are compared against.
