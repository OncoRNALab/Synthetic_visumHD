#!/usr/bin/env python3
"""
Split molecules.tsv into N roughly equal shards.
Outputs: shard_0.tsv, shard_1.tsv, ..., shard_{N-1}.tsv
Each shard keeps the full header row.
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser(description="Shard molecules.tsv for parallel processing")
    ap.add_argument("--molecules", required=True, help="molecules.tsv from prepare_transcript_abundance")
    ap.add_argument("--n-shards",  type=int, required=True)
    ap.add_argument("--prefix",    default="shard_", help="Output file prefix")
    args = ap.parse_args()

    if args.n_shards < 1:
        sys.stderr.write("n-shards must be >= 1\n")
        return 1

    with open(args.molecules) as f:
        header = f.readline()
        rows   = f.readlines()

    n   = len(rows)
    sz  = max(1, (n + args.n_shards - 1) // args.n_shards)

    n_written = 0
    for i in range(args.n_shards):
        chunk = rows[i * sz : (i + 1) * sz]
        if not chunk:
            break
        out_path = "%s%d.tsv" % (args.prefix, i)
        with open(out_path, "w") as fout:
            fout.write(header)
            fout.writelines(chunk)
        n_written += 1

    sys.stderr.write("Wrote %d shards from %d molecules\n" % (n_written, n))
    return 0


if __name__ == "__main__":
    sys.exit(main())
