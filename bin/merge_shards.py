#!/usr/bin/env python3
"""
Concatenate per-shard FASTQ and TSV files into final outputs.

Usage
-----
merge_shards.py \\
    --r1-list       "exp_R1_0.fastq exp_R1_1.fastq ..." \\
    --r2-list       "exp_R2_0.fastq ..." \\
    --mutations-list "barcode_mutations_0.tsv ..." \\
    --mapping-list   "read_to_spot_0.tsv ..." \\
    --out-r1        simulated_R1.fastq.gz \\
    --out-r2        simulated_R2.fastq.gz \\
    --out-mutations barcode_mutations.tsv \\
    --out-mapping   read_to_spot.tsv

All list args accept either space-separated strings or can be repeated.
"""
import argparse
import gzip
import sys


def cat_fastq_gz(file_list, out_gz):
    with gzip.open(out_gz, "wt") as out:
        for path in file_list:
            with open(path) as f:
                for line in f:
                    out.write(line)


def cat_tsv(file_list, out_path):
    """Concatenate TSV files, keeping header from first file only."""
    wrote_header = False
    with open(out_path, "w") as out:
        for path in file_list:
            with open(path) as f:
                header = f.readline()
                if not wrote_header:
                    out.write(header)
                    wrote_header = True
                for line in f:
                    out.write(line)


def parse_list(values):
    """Flatten a list of potentially space-separated strings into file paths."""
    result = []
    for v in values:
        result.extend(v.split())
    return result


def main():
    ap = argparse.ArgumentParser(description="Merge per-shard simulation outputs")
    ap.add_argument("--r1-list",        nargs="+", required=True)
    ap.add_argument("--r2-list",        nargs="+", required=True)
    ap.add_argument("--mutations-list", nargs="+", required=True)
    ap.add_argument("--mapping-list",   nargs="+", required=True)
    ap.add_argument("--out-r1",         default="simulated_R1.fastq.gz")
    ap.add_argument("--out-r2",         default="simulated_R2.fastq.gz")
    ap.add_argument("--out-mutations",  default="barcode_mutations.tsv")
    ap.add_argument("--out-mapping",    default="read_to_spot.tsv")
    args = ap.parse_args()

    r1_files  = parse_list(args.r1_list)
    r2_files  = parse_list(args.r2_list)
    mut_files = parse_list(args.mutations_list)
    map_files = parse_list(args.mapping_list)

    sys.stderr.write("Merging %d shards\n" % len(r1_files))

    cat_fastq_gz(r1_files, args.out_r1)
    cat_fastq_gz(r2_files, args.out_r2)
    cat_tsv(mut_files, args.out_mutations)
    cat_tsv(map_files, args.out_mapping)

    sys.stderr.write("Done: %s  %s  %s  %s\n"
                     % (args.out_r1, args.out_r2,
                        args.out_mutations, args.out_mapping))
    return 0


if __name__ == "__main__":
    sys.exit(main())
