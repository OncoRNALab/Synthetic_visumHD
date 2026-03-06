#!/usr/bin/env python3
# Paired R1/R2 FASTQ from read assignments and reference FASTA. Deterministic with seed.

import argparse
import random
import sys

BASES = "ACGT"
RC = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(s):
    return s.translate(RC)[::-1]

def load_fasta(path):
    name, seq = None, []
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(seq)
                name = line[1:].split()[0].strip()
                seq = []
            else:
                seq.append(line.rstrip())
        if name:
            out[name] = "".join(seq)
    return out

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--assignments", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--len-r1", type=int, default=75)
    p.add_argument("--len-r2", type=int, default=100)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--use-art", action="store_true")
    p.add_argument("--out-r1", default="sim_R1.fastq")
    p.add_argument("--out-r2", default="sim_R2.fastq")
    args = p.parse_args()
    random.seed(args.seed)
    ref = load_fasta(args.reference)
    if not ref:
        sys.stderr.write("No sequences in reference\n")
        return 1
    rows = []
    with open(args.assignments) as f:
        next(f)
        for line in f:
            t = line.rstrip().split("\t")
            if len(t) >= 3:
                rows.append((t[0], t[1], t[2]))
    L1, L2 = args.len_r1, args.len_r2
    with open(args.out_r1, "w") as r1, open(args.out_r2, "w") as r2:
        for read_id, _, tid in rows:
            seq = ref.get(tid) or list(ref.values())[0]
            if len(seq) < L1 + L2:
                seq += "".join(random.choices(BASES, k=L1 + L2 - len(seq) + 50))
            start = random.randint(0, max(0, len(seq) - L1 - L2))
            s1 = seq[start : start + L1]
            s2 = revcomp(seq[start + L1 : start + L1 + L2])
            if len(s2) < L2:
                s2 += "".join(random.choices(BASES, k=L2 - len(s2)))
            q1, q2 = "F" * L1, "F" * L2
            r1.write("@%s\n%s\n+\n%s\n" % (read_id, s1, q1[:len(s1)]))
            r2.write("@%s\n%s\n+\n%s\n" % (read_id, s2, q2[:len(s2)]))
    return 0

if __name__ == "__main__":
    sys.exit(main())
