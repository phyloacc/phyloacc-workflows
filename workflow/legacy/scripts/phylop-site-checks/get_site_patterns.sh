#!/bin/bash

while read -r chr start end score; do
    site=$((start+1)) # convert BED to 1-based for msa_view
    msa_view "../../data/birds/test2/mafs/group1/${chr}.maf" --start "$site" --end "$site" --out-format SS | \
        awk 'NF==3 && $1 ~ /^[0-9]+$/ {print $0}'
done < bird-chr10-max-scores.bed > bird-chr10-patterns-of-max-sites.txt