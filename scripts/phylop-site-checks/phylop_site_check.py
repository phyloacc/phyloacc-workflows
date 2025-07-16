import re

wig_file = "../../data/birds/test2/04-phylop/group1/chr10-phylop.wig"
bed_file = "bird-chr10-max-scores.bed"

# --- PASS 1: Find the max score ---
max_score = None

with open(wig_file) as f:
    for line in f:
        line = line.strip()
        # Skip header or empty/comment
        if line == "" or line.startswith(("fixedStep", "#")):
            continue
        try:
            score = float(line)
            if (max_score is None) or (score > max_score):
                max_score = score
        except ValueError:
            continue

print("Max score is", max_score)

# --- PASS 2: Extract sites with max score as BED ---
with open(wig_file) as f, open(bed_file, "w") as out:
    chrom = None
    start = None
    step = None
    pos = None
    for line in f:
        line = line.strip()
        if line.startswith("fixedStep"):
            m_chrom = re.search(r"chrom=(\S+)", line)
            m_start = re.search(r"start=(\d+)", line)
            m_step  = re.search(r"step=(\d+)", line)
            chrom = m_chrom.group(1)
            start = int(m_start.group(1))
            step  = int(m_step.group(1))
            pos = start - 1  # 0-based BED
            continue
        if line == "" or line.startswith("#"):
            continue
        try:
            score = float(line)
        except ValueError:
            continue
        if score == max_score:
            # BED format: chrom, start, end, score
            print(f"{chrom}\t{pos}\t{pos+1}\t{score}", file=out)
        pos += step