#!/bin/bash

# Output combined file name
out=bird-chr10-max-score-sites.maf

# Get the header from the first file
echo "##maf version=1 scoring=roast" > "$out"

# Now append all the MAF blocks (skip header in each file)
for f in bird-chr10-max-score-mafs/*.maf; do
    awk '!/^##/' "$f" >> "$out"
done
