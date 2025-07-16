#############################################################################
# Given 4d sites parsed from a maf by msa_view in ss format, retain only
# sites that are present in > than some threshold of species
#
# Notes on SS format:
# LENGTH = number of sites (sum of all site occurences)
# NTUPLES = number of site patterns
#
# Gregg Thomas, September 2023
#############################################################################

import sys

ss_file, threshold, summary_file, out_file = sys.argv[1:];
threshold = float(threshold);

#ss_file = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-sicb/data/01-zoonomia-aln/01-4d-sites/all-4d-sites.ss";
#threshold = 0.5;
#summary_file = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-sicb/summary-data/all-4f-sites-filtered-" + str(threshold) + "-summary.tsv";
#out_file = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-sicb/data/01-zoonomia-aln/01-4d-sites/all-4d-sites-filtered-" + str(threshold) + ".ss";
# Input/output parameters

###############

summary_stat_cats = [ "num.seqs", "aln.len", "num.sites", "avg.site.count", "max.site.count", "num.singletons", "total.missing.cols", "avg.missing.per.site" ];
filter_cats = ["pre.filter", "post.filter"];
summary_stats = { stat_cat : { filter_cat : 0 for filter_cat in filter_cats } for stat_cat in summary_stat_cats }
# Summary stats to collect both before and after filtering, will be written to summary tsv file

###############

headers = [];
sites = [];
# Output headers and filtered sites

header = True;
# A flag while reading the input file

for line in open(ss_file):
    line = line.strip();
    # Remove the newline from the input line
    
    if header:
        if line.startswith("NSEQS = "):
            num_seqs = int(line.split(" = ")[1]);
            summary_stats["num.seqs"]["pre.filter"] = num_seqs;
            # Get the number of sequences from the input file            
            
            spec_missing_counts = [0] * num_seqs;
            # A list of the number of times in a site a species is missing, the index of the list
            # corresponds to the species
            # Initialize all values as 0 here         

            max_missing = round(num_seqs * threshold);
            # For filtering, this is the maximum number of missing sequences (*) allowed in a given site

            headers.append(line);
            # Add the NSEQS header line to the output headers list

        elif line.startswith("LENGTH = "):
            summary_stats["aln.len"]["pre.filter"] = int(line.split(" = ")[1]);
            # Get the pre-filter alignment length

            headers.append("UPDATE LENGTH");
            # Instead of adding this line to the output headers, add a temp string so we can update it after filtering

        elif line.startswith("NTUPLES = "):
            summary_stats["num.sites"]["pre.filter"] = int(line.split(" = ")[1]);
            # Get the pre-filter number of sites

            headers.append("UPDATE NTUPLES");
            # Instead of adding this line to the output headers, add a temp string so we can update it after filtering

        else:
            headers.append(line);
            # Add all other lines to the output header list as they are

        if line == "":
            header = False;
            continue;
        # If the line is blank, that indicates the end of the header
    ## End header block

    ###############

    if not header:
        line = line.split("\t");
        site_count = int(line[2]);
        # Split the line for each site by tabs and get the number of times this site occurs

        num_missing = line[1].count("*") * site_count;
        summary_stats["total.missing.cols"]["pre.filter"] += num_missing;
        # Count the number of missing sequences in this site

        if site_count == 1:
            summary_stats["num.singletons"]["pre.filter"] += 1;
        # If the site count is 1, add it to the pre-filter singletons

        if site_count > summary_stats["max.site.count"]["pre.filter"]:
            summary_stats["max.site.count"]["pre.filter"] = site_count;
        # If the site count exceeds the previous maximum count, set the pre-filter max site count to the current site count

        ##########

        if num_missing <= max_missing:
        # If the number of missing sequences in this site pattern is lower than the number required by the threshold for filtering,
        # enter the output site block -- this site will be written to output and others will not

            summary_stats["total.missing.cols"]["post.filter"] += num_missing;
            # Add the number of missing sequences in this site pattern to the post-filter count of missing columns

            if site_count == 1:
                summary_stats["num.singletons"]["post.filter"] += 1;
            # If the site count is 1, add it to the post-filter singletons

            if site_count > summary_stats["max.site.count"]["post.filter"]:
                summary_stats["max.site.count"]["post.filter"] = site_count;
            # If the site count exceeds the previous maximum count, set the post-filter max site count to the current site count

            for col_index in range(len(line[1])):
                if line[1][col_index] == "*":
                    spec_missing_counts[col_index] += 1;
            # For each sequence in the alignment, check if it is missing in the current site pattern and add 1 to the
            # missing counts list if so

            summary_stats["aln.len"]["post.filter"] += site_count;
            # Add the number of times this site occurs to the post-filter alignment lengths

            site = [str(len(sites)), line[1], line[2]];
            sites.append("\t".join(site));
            # Compile the output site string
        ## End filter block
    ## End site block
## End line loop

###############

summary_stats["num.sites"]["post.filter"] = len(sites);
# Get the number of sites retaind after filtering

for spec in spec_missing_counts:
    if spec != summary_stats["num.sites"]["post.filter"]:
        summary_stats["num.seqs"]["post.filter"] += 1;
# Count how many sequences are present in at least one site pattern after filtering (i.e. check if a sequence is missing in all the site patterns now)

assert summary_stats["num.seqs"]["pre.filter"] == summary_stats["num.seqs"]["post.filter"], "\n LOST SOME SEQUENCES DURING FILTERING";
# If any of the sequences is missing in the filtered alignment, print an error message here
# TODO: A better type of filtering?

###############

with open(out_file, "w") as out_stream:
    for header in headers:
        if header == "UPDATE LENGTH":
            out_stream.write("LENGTH = " + str(summary_stats["aln.len"]["post.filter"]) + "\n");
        elif header == "UPDATE NTUPLES":
            out_stream.write("NTUPLES = " + str(summary_stats["num.sites"]["post.filter"]) + "\n");
        else:
            out_stream.write(header + "\n");
        # Output the header lines and update those that need updating
    ## End header output

    for site in sites:
        out_stream.write(site + "\n");
    ## End site output
## Close output file

###############        

summary_stats["avg.site.count"]["pre.filter"] = summary_stats["aln.len"]["pre.filter"] / summary_stats["num.sites"]["pre.filter"];
summary_stats["avg.site.count"]["post.filter"] = summary_stats["aln.len"]["post.filter"] / summary_stats["num.sites"]["post.filter"];
# Calculate the average site count

summary_stats["avg.missing.per.site"]["pre.filter"] = summary_stats["total.missing.cols"]["pre.filter"] / summary_stats["aln.len"]["pre.filter"];
summary_stats["avg.missing.per.site"]["post.filter"] = summary_stats["total.missing.cols"]["post.filter"] / summary_stats["aln.len"]["post.filter"];
# Calculate the average number of missing sequences per site

###############

with open(summary_file, "w") as out_stream:
    out_stream.write("# input file: " + ss_file + "\n");
    out_stream.write("# missing threshold: " + str(threshold) + "\n");
    out_stream.write("# output file: " + out_file + "\n");
    # Some log info for the summary file

    summary_headers = [ "var", "filter.cat", "value" ];
    out_stream.write("\t".join(summary_headers) + "\n");
    # Output headers for the summary file

    for filter_cat in filter_cats:
        for summary_stat in summary_stat_cats:
            outline = "\t".join([ str(col) for col in [ summary_stat, filter_cat, summary_stats[summary_stat][filter_cat] ] ]);
            out_stream.write(outline + "\n");
    # Write the summary statistics
## Close summary output file

#############################################################################