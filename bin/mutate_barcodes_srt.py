#!/usr/bin/env python3
"""
Streaming barcode mutation: substitution, insertion, deletion on first barcode_length
bases of R1. Outputs TSV with read_id, original_barcode, mutated_barcode, mutation_events.
Deterministic given --seed.
"""
import argparse
import random
import sys

BASES = "ACGT"


def mutate_region(seq, sub_rate, ins_rate, del_rate, rng):
    out = []
    desc = []
    i = 0
    while i < len(seq):
        if rng.random() < del_rate:
            desc.append("del:%d:%s" % (i, seq[i]))
            i += 1
            continue
        if rng.random() < ins_rate:
            b = rng.choice(BASES)
            out.append(b)
            desc.append("ins:%d:%s" % (i, b))
        if i < len(seq):
            b = seq[i]
            if rng.random() < sub_rate:
                new_b = rng.choice(BASES)
                out.append(new_b)
                desc.append("sub:%d:%s>%s" % (i, b, new_b))
            else:
                out.append(b)
            i += 1
    return "".join(out), desc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--barcode-length", type=int, default=36)
    ap.add_argument("--sub-rate", type=float, default=0.001)
    ap.add_argument("--ins-rate", type=float, default=0.0005)
    ap.add_argument("--del-rate", type=float, default=0.0005)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-r1", default="mutated_R1.fastq")
    ap.add_argument("--out-r2", default="mutated_R2.fastq")
    ap.add_argument("--out-mutations", default="barcode_mutations.tsv")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    L = args.barcode_length

    with open(args.r1) as r1_in, open(args.r2) as r2_in, \
         open(args.out_r1, "w") as r1_out, open(args.out_r2, "w") as r2_out, \
         open(args.out_mutations, "w") as mut_out:
        mut_out.write("read_id\toriginal_barcode\tmutated_barcode\tmutation_events\n")
        while True:
            h1 = r1_in.readline()
            if not h1:
                break
            s1 = r1_in.readline().rstrip()
            _ = r1_in.readline()
            q1 = r1_in.readline().rstrip()
            h2 = r2_in.readline()
            s2 = r2_in.readline()
            _ = r2_in.readline()
            q2 = r2_in.readline()

            read_id = h1[1:].split()[0].strip()
            orig_bc = s1[:L]
            rest_seq = s1[L:]
            rest_qual = q1[L:] if len(q1) >= L else ""

            mut_bc, desc_list = mutate_region(orig_bc, args.sub_rate, args.ins_rate, args.del_rate, rng)
            events_str = ";".join(desc_list) if desc_list else "none"

            new_s1 = mut_bc + rest_seq
            new_q1 = ("F" * len(mut_bc))[:len(mut_bc)] + rest_qual
            if len(new_q1) > len(new_s1):
                new_q1 = new_q1[:len(new_s1)]
            elif len(new_q1) < len(new_s1):
                new_q1 += "F" * (len(new_s1) - len(new_q1))

            r1_out.write(h1)
            r1_out.write(new_s1 + "\n")
            r1_out.write("+\n")
            r1_out.write(new_q1 + "\n")
            r2_out.write(h2)
            r2_out.write(s2)
            r2_out.write("+\n")
            r2_out.write(q2)
            mut_out.write("%s\t%s\t%s\t%s\n" % (read_id, orig_bc, mut_bc, events_str))

    return 0


if __name__ == "__main__":
    sys.exit(main())
