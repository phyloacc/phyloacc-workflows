import csv
import gzip

def load_mapping(mapping_file):
    """
    Load mapping from sequence_report.tsv.
    
    Reads a TSV file with headers including 'GenBank seq accession'
    and 'RefSeq seq accession' and returns a dictionary mapping RefSeq IDs to GenBank IDs.
    """
    mapping = {}
    with open(mapping_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            refseq = row["RefSeq seq accession"].strip()
            genbank = row["GenBank seq accession"].strip()
            if refseq and genbank:
                mapping[refseq] = genbank
    return mapping

def update_gff(gff_in, gff_out, mapping):
    """
    Updates the scaffold identifiers in the GFF file and inserts a conversion documentation line
    at the end of the header block.
    
    This function reads a gzipped GFF input file (if the filename ends with '.gz') in text mode,
    buffers all header lines (lines starting with '#') and writes them first after updating RefSeq IDs.
    Then it writes an extra header line documenting the conversion before processing the remaining lines.
    """
    # Open the input file; use gzip if the file is gzipped.
    if gff_in.endswith('.gz'):
        infile = gzip.open(gff_in, 'rt')
    else:
        infile = open(gff_in, 'r')
    
    # Buffer for header lines.
    header_lines = []
    header_written = False

    with infile, open(gff_out, 'w') as outfile:
        for line in infile:
            # If still in the header block.
            if line.startswith("#") and not header_written:
                # Replace any occurrence of a RefSeq ID with the GenBank ID in the header.
                for refseq, genbank in mapping.items():
                    if refseq in line:
                        line = line.replace(refseq, genbank)
                header_lines.append(line)
            else:
                # First non-header line: write out all buffered header lines first,
                # then the conversion documentation, and mark header as written.
                if not header_written:
                    for header in header_lines:
                        outfile.write(header)
                    outfile.write("# Conversion: Scaffold IDs replaced from RefSeq to GenBank using sequence_report.tsv\n")
                    header_written = True

                # Process non-header lines: update the first column (scaffold ID) if it matches.
                parts = line.rstrip().split("\t")
                if parts and parts[0] in mapping:
                    parts[0] = mapping[parts[0]]
                outfile.write("\t".join(parts) + "\n")
        
        # In case the file contained only header lines, write the conversion line.
        if not header_written:
            for header in header_lines:
                outfile.write(header)
            outfile.write("# Conversion: Scaffold IDs replaced from RefSeq to GenBank using sequence_report.tsv\n")

def main():
    # Filenames for mapping table and GFF files.
    mapping_file = "../summary-data/turtles/sequence_report.tsv"  # Mapping table with your sequence information.
    gff_in = "/n/holylfs05/LABS/informatics/Users/gthomas/turtles/genomes/Dcoriacea/GCA_009764565.3_rDerCor1.pri.v3/GCF_009764565.2_rDerCor1.pri.v3_genomic.gff.gz"                # Input RefSeq GFF file (gzipped).
    gff_out = "../data/turtles/GCF_009764565.2_rDerCor1.pri.v3_genomic.GENBANK_IDS.gff"     # Output file (uncompressed).
    
    # Load mapping from the TSV file.
    mapping = load_mapping(mapping_file)
    
    # Process the GFF file to update scaffold IDs and append the conversion header at the end of the header block.
    update_gff(gff_in, gff_out, mapping)
    
    print("GFF file updated successfully. Output written to:", gff_out)

if __name__ == "__main__":
    main()


