#############################################################################
# Using HDBSCAN to cluster conserved sites
#
# Gregg Thomas, October 2024
#############################################################################

import sys
import pandas as pd
import numpy as np
import hdbscan

#############################################################################

CONSERVED_SITES_FILE, OUTPUT_BED_FILE, CHROMOSOME, MIN_CLUSTER_SIZE, MIN_SAMPLES = sys.argv[1:];

#conserved_sites_file = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01b-zoonomia-aln/04-phylop/all-chromosomes/chr19-phylop-fdr-0.05.conserved.bed";
#conserved_sites_file = "test-conserved-sites.bed";

print("Reading file");
conserved_df = pd.read_csv(CONSERVED_SITES_FILE, sep='\t', header=None);
# Read the BED file with conserved sites

print("Extracting sites");
conserved_sites_array = conserved_df.iloc[:, 1].values;  # Modify the column index if necessary
# Extract start positions; assuming single sites share their start with end

print("Reshaping array");
conserved_sites = conserved_sites_array.reshape(-1, 1);
# Convert to a 2D NumPy array as required for HDBSCAN

print("Initializing clustering");
clusterer = hdbscan.HDBSCAN(min_cluster_size=MIN_CLUSTER_SIZE, min_samples=MIN_SAMPLES);  # Adjust parameters based on your dataset
# Initialize and run HDBSCAN
# min_cluster_size: the smallest size (number of points) that any cluster can have. It determines the minimum number of points required to form a cluster.
# min_samples: the number of samples in a neighborhood for a point to be considered a core point. This includes the point itself.

print("Clustering");
labels = clusterer.fit_predict(conserved_sites);

# Output the cluster labels
print("Cluster labels:", labels)

clusters_ranges = {};
# Initialize a dictionary to store cluster ranges


for idx, label in enumerate(labels):
    if label == -1:
        continue;  # Skip noise points
    if label not in clusters_ranges:
        clusters_ranges[label] = {'start': float('inf'), 'end': float('-inf')};
    # Update the min/max extents of the cluster
    clusters_ranges[label]['start'] = min(clusters_ranges[label]['start'], conserved_sites_array[idx]);
    clusters_ranges[label]['end'] = max(clusters_ranges[label]['end'], conserved_sites_array[idx]);
# Iterate over indices and cluster labels

bed_output = [];
for cluster_id, bounds in clusters_ranges.items():
    start = bounds['start'];
    end = bounds['end'];
    bed_output.append([CHROMOSOME, start, end, f'cluster_{cluster_id}']);
# Iterate over the cluster ranges and store them in a BED format

bed_df_output = pd.DataFrame(bed_output, columns=['chrom', 'chromStart', 'chromEnd', 'cluster_id']);
# Convert list to DataFrame

bed_df_output.to_csv(OUTPUT_BED_FILE, sep='\t', header=False, index=False);
print(f"Clusters saved in BED format to {output_bed_file}");
# Write to a BED file

#############################################################################