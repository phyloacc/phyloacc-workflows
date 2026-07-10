import os
import sys
from collections import Counter

def parse_maf_block(block_lines):
    seqs = []
    for line in block_lines:
        if line.startswith('s '):
            fields = line.strip().split()
            seqs.append(fields[6])
    return seqs

def count_site_patterns_maf_file(filepath):
    patterns = Counter()
    with open(filepath) as f:
        block = []
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            if line.startswith('a'):
                if block:
                    seqs = parse_maf_block(block)
                    if len(seqs) > 1:
                        aln_len = len(seqs[0])
                        for col in range(aln_len):
                            pattern = ''.join(seq[col].upper() for seq in seqs)
                            patterns[pattern] += 1
                    block = []
                continue
            if line.startswith('s '):
                block.append(line)
        # process any remaining block
        if block:
            seqs = parse_maf_block(block)
            if len(seqs) > 1:
                aln_len = len(seqs[0])
                for col in range(aln_len):
                    pattern = ''.join(seq[col].upper() for seq in seqs)
                    patterns[pattern] += 1
    return patterns

def main(maf_dir, out_fn="site_pattern_counts.tsv"):
    all_patterns = Counter()
    num_files = 0
    for fn in os.listdir(maf_dir):
        if not fn.endswith('.maf'):
            continue
        num_files += 1
        if num_files % 100 == 0:
            print("Processed {} files...".format(num_files))
        fp = os.path.join(maf_dir, fn)
        #print("Processing {}".format(fp))
        patterns = count_site_patterns_maf_file(fp)
        all_patterns.update(patterns)
    with open(out_fn, "w") as out:
        out.write("pattern\tcount\n")
        for pattern, count in sorted(all_patterns.items(), key=lambda x: -x[1]):
            out.write("{}\t{}\n".format(pattern, count))
    print("Wrote results to {}".format(out_fn))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python count_maf_site_patterns_all.py <maf_directory>")
        sys.exit(1)
    main(sys.argv[1])