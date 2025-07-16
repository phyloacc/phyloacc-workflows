#!/usr/bin/awk -f
# Extract conserved sites from BED+FDR table
# Keep entries where conservation_code == "0" and significance == "1"

BEGIN { OFS = "\t" }

$4 == "0" && $6 == "1" {
    print $1, $2, $3, $4 $6
}