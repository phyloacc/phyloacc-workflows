#!/usr/bin/awk -f
# Usage:
#   awk -v chr=CHR_NAME -v prefix=PREFIX -f ref_gff_split_by_chr.awk input.gff > output.gff
#
# Retains lines beginning with "#!" unconditionally.
# Keeps lines where $1 == chr and prepends prefix.

BEGIN {
    FS = OFS = "\t"
}

{
    if ($0 ~ /^#!/) {
        print
    }
    else if ($1 == chr) {
        print prefix $0
    }
    # lines not matching chr are skipped
}