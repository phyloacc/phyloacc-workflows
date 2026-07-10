#!/usr/bin/awk -f
# Extract accelerated sites from BED+FDR table
# Keep entries where conservation_code == "2" and significance == "1"

BEGIN { OFS = "\t" }

$4 == "2" && $6 == "1" {
    print $1, $2, $3, $4 $6
}