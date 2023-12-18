#############################################################################
# This script prepares a script for maffilter to count nucleotide frequencies
# for every species in a MAF file.
#
# Gregg Thomas, October 2023
#############################################################################

import sys

######################
# maffilter template
######################

maffilter_template = """input.file={maf_file}
input.file.compression=none
input.format=Maf
output.log={log_file}
maf.filter=\\
	SelectChr(\\
		ref_species={ref_species},\\
		chromosome={cur_chrome},\\
	),\\
	Subset(\\
		species=({species_list}),\\
		remove_duplicates=yes,\\
	),\\
	SequenceStatistics(\\
		statistics=(\\
			BlockLength,\\
{block_counts}
		),\\
		ref_species={ref_species},\\
		file={out_file}\\
	)
"""

#############################################################################

maf_file, species_list_str, ref_species, cur_chrome, maffilter_log_file, maffilter_out_file, maffilter_script_file = sys.argv[1:];


species_list = species_list_str.split(",");
block_counts = [];
for species in species_list:
    block_counts.append("\t\t\tBlockCounts(species={species}, suffix=.{species}),\\".format(species=species));
block_counts_str = "\n".join(block_counts);
block_counts_str = block_counts_str[:-2] + "\\";

with open(maffilter_script_file, "w") as script_out_stream:
    script_out_stream.write(maffilter_template.format(maf_file=maf_file, log_file=maffilter_log_file, ref_species=ref_species, cur_chrome=cur_chrome, species_list=species_list_str, block_counts=block_counts_str, out_file=maffilter_out_file));

#############################################################################