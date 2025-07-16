#!/usr/bin/awk -f
# --------------------------------------------------------------------
# count_sites.awk
#
# Counts the number of times each unique value appears in a given column.
# Specifically, it counts occurrences in the first column (e.g., chromosome ID).
#
# Usage:
#   awk -f count_sites.awk input.bed > output.counts
#
# Output:
#   <value>  <count>
# --------------------------------------------------------------------

{
    count[$1]++   # Count occurrences of the value in column 1
}

END {
    if (length(count) == 0) {
        print 0
    } else {
        for (val in count) {
            print val, count[val]
        }
    }
}