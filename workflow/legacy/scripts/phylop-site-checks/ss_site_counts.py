import sys

filename = "../../data/birds/test2/mafs/group1/chr10.ss"  # Change as needed

total = 0

with open(filename) as f:
    for line in f:
        # Skip header lines (all caps or field lines)
        if line[0].isalpha():
            continue
        parts = line.strip().split()
        if len(parts) != 3:
            continue
        idx, alleles, count = parts
        present = [base for base in alleles if base not in ['-', '*', 'N']]
        # print(line);
        # sys.exit();
        if present and len(present) == 45 and len(set(present)) == 1:
            print(line)
            total += int(count)

print(f"Total sites where all present alleles are identical: {total}")