#!/usr/bin/env bash
set -euo pipefail

# Translate a UCSC-style bed's chrom column to GenBank accessions, using the
# same sequence-report header-lookup idiom as
# 04_replace_bed_seqids_with_genbank.sh in the support data repo. Restricts
# output to only the 21 chromosomes this pipeline actually has predictions
# for (drops _random/chrUn_/unplaced scaffolds). Handles both a plain BED3
# (chrom in field 1, e.g. the ATAC peaks bed) and mm39-35way.bed's own
# UCSC-table format (leading "bin" column, chrom in field 2) via the
# CHROM_FIELD argument.
#
# Usage: translate_ucsc_chroms.sh <ucsc_bed> <sequence_report.tsv> <output_bed> <target_accessions_file> <chrom_field>
# <target_accessions_file>: one GenBank accession per line (our 21 real chromosomes) -
# rows mapping to anything else (unplaced scaffolds, etc.) are dropped.
# <chrom_field>: 1-based column index of the chrom field in <ucsc_bed> (1 for a
# plain BED3+, 2 for mm39-35way.bed's bin-prefixed table format).

UCSC_BED="$1"
SEQUENCE_REPORT="$2"
OUTPUT_BED="$3"
TARGET_ACCESSIONS="$4"
CHROM_FIELD="${5:-1}"

awk -v chrom_field="$CHROM_FIELD" '
BEGIN { FS = OFS = "\t" }
FNR == 1 && ARGIND == 1 {
  for (i = 1; i <= NF; i++) header[$i] = i
  ucsc_col = header["UCSC style name"]
  genbank_col = header["GenBank seq accession"]
  next
}
ARGIND == 1 {
  ucsc = $ucsc_col
  genbank = $genbank_col
  if (ucsc == "" || genbank == "") next
  map[ucsc] = genbank
  next
}
ARGIND == 2 { targets[$1] = 1; next }
ARGIND == 3 {
  if ($1 ~ /^#/) next
  chrom = $chrom_field
  if (!(chrom in map)) next
  genbank = map[chrom]
  if (!(genbank in targets)) next
  out = genbank
  for (i = chrom_field + 1; i <= NF; i++) out = out OFS $i
  print out
}
' "$SEQUENCE_REPORT" "$TARGET_ACCESSIONS" "$UCSC_BED" > "$OUTPUT_BED"

echo "Wrote $(wc -l < "$OUTPUT_BED") rows to $OUTPUT_BED"
