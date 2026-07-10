#!/usr/bin/awk -f
# --------------------------------------------------------------------------
# convert_wig_to_bed.awk
#
# Converts a fixedStep WIG file into a BED-like format with calculated
# adjusted values based on the input WIG score (typically -log10(p-values)).
# # Supports WIG files with multiple fixedStep blocks, possible gaps, etc.
#
# Fields in output:
#   chrom, start, end, raw score=log(raw p-value), conservation status (0 = conserved, 1 = neutral, 2 = accelerated), raw p_value
#
# Usage example:
#   awk -f wig_to_bed.awk input.wig > output.bed
#
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# BEGIN block runs before any lines are processed
# Set output field separator to tabs (required for BED format)
BEGIN {
    OFS = "\t"          # Use tab-separated values
    chrom = ""          # Current chromosome
    start = 0           # Start position from fixedStep line
    step = 0            # Step size from fixedStep line
    counter = 0         # Number of data lines processed under current track
}

# --------------------------------------------------------------------------
# If a variableStep line is ever found, print error and exit.
/^variableStep/ {
    print "ERROR: variableStep WIG lines detected. This script only supports fixedStep WIG blocks." > "/dev/stderr"
    exit 1
}

# --------------------------------------------------------------------------
# Detect and parse a fixedStep header line
# Example: fixedStep chrom=chr1 start=1000 step=5
/^fixedStep/ {
    # Loop through each field on the fixedStep line
    for (i = 1; i <= NF; i++) {
        split($i, a, "=")
        if (a[1] == "chrom") {
            chrom = a[2]            # Extract chromosome
        } else if (a[1] == "start") {
            start = a[2] - 1        # Convert to 0-based (BED standard)
        } else if (a[1] == "step") {
            step = a[2]             # Step between positions
        }
    }

    # Reset the counter for the following data block
    counter = 0
    next  # Skip processing the fixedStep line as a data line
}

# Ignore comments, track lines, blank lines
/^(track|#|$)/ { next }

# --------------------------------------------------------------------------
# Process each numeric score line following a fixedStep header
# Each score corresponds to a single genomic position
{
    # Compute 0-based genomic coordinate from start + offset
    pos = start + counter * step
    end = pos + 1

    # Read conservation score (assumed to be a -log10(p-value) format)
    score = $1

    # Assign a strand-like code:
    #   "0" for positive scores
    #   "2" for negative scores
    #   "1" for zero
    if (score > 0) {
        cons_score = "0"
        pval = 10^-score
    } else if (score < 0) {
        cons_score = "2"
        pval = 10^score  # note: score is negative
    } else {
        cons_score = "1"
        pval = 1  # 10^0
    }

    # Output fields in BED-like format + p-value
    # Format: chrom, start, end, raw score=log(raw p-value), conservation status (0 = conserved, 1 = neutral, 2 = accelerated), raw p_value
    print chrom, pos, end, score, cons_score, pval

    # Move to next offset for current block
    counter++
}