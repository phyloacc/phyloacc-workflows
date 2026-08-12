#!/usr/bin/awk -f
# Usage:
#   awk -v match_chrom=GFF_NAME -v out_seqname=MAF_SEQNAME -f ref_gff_split_by_chr.awk in.gff > out.gff
#
# Keeps rows whose chromosome (column 1) equals match_chrom (the chromosome name as it
# appears in the GFF), and rewrites column 1 to out_seqname (the sequence name as it
# appears in the MAF, e.g. "Homo_sapiens.chr1") so msa_view --features lines up with the
# alignment. Retains "#!" header lines unconditionally.

BEGIN {
    FS = OFS = "\t"
}

{
    if ($0 ~ /^#!/) {
        print
    }
    else if ($1 == match_chrom) {
        $1 = out_seqname
        print
    }
    # rows for other chromosomes are skipped
}
