#!/usr/bin/env bash
# A script to apply Benjamini-Hochberg FDR correction to a BED file with phyloP scores
# Usage: fdr_adjust.sh input.bed output.bed.tmp /tmp 4 0.05

set -euo pipefail  # Exit on error, unset variables, and fail in pipelines

# === INPUT PARAMETERS ===
PH="$1"      # Input BED file with phyloP p-values in column 6
OUT="$2"     # Output BED file with adjusted p-values and significance
TMP="$3"     # Temporary directory for sorting
CPU="$4"     # Number of CPU threads to use in sorting
ALPHA="$5"   # Significance threshold (e.g., 0.05)

# === COUNT TOTAL NUMBER OF LINES (used for BH procedure) ===
TOTAL=$(wc -l < "$PH")

# === SORT BY P-VALUE ASCENDING (column 6) ===
# This prepares for ranking p-values from smallest to largest
LC_ALL=C sort -k6,6g --parallel="$CPU" -T "$TMP" -S 2G "$PH" |

# === ADD RANK (line number) AS NEW COLUMN ===
# Outputs: original columns + rank (column 7)
awk -v OFS="\t" '{ print $0, NR }' |

# === SORT BY RANK DESCENDING (column 7) ===
# This is required by BH procedure to walk backwards when adjusting
LC_ALL=C sort -rn -k7 --parallel="$CPU" -T "$TMP" -S 2G |

# === APPLY BENJAMINI-HOCHBERG FDR ADJUSTMENT ===
# Compute adjusted p-value: (total / rank) * p
# Apply monotonic correction: adjusted p-values must be non-increasing
# Clamp adjusted p-values at 1 if necessary
# Determine significance (1 if p_adj < alpha, else 0)
# Outputs: original columns, p_adj (col 8), significance (col 9)
awk -v OFS="\t" -v total="$TOTAL" -v alpha="$ALPHA" '
  BEGIN { prev_pn = 999 } 
  {
    p_adj = (total / $7) * $6
    if (p_adj > prev_pn) { p_adj = prev_pn }
    prev_pn = p_adj
    if (p_adj > 1) { p_adj = 1 }
    sig = (p_adj < alpha) ? 1 : 0
    print $0, p_adj, sig
  }
' |

# === SORT OUTPUT BACK TO BED ORDER (chrom, start) ===
# BED format requires sorted output, typically by chrom and start
LC_ALL=C sort -k1,1 -k2,2n --parallel="$CPU" -T "$TMP" -S 2G |

# === OUTPUT FINAL FIELDS ===
# Keep columns: chrom (1), start (2), end (3), score (5), p_adj (8), sig (9)
cut -f1,2,3,5,8,9 > "$OUT"