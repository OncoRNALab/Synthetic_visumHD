#!/usr/bin/env python3
# Write summary_metrics.json from spot_coordinates, barcode_mutations, read_to_spot. CWD = process work dir.

import argparse
import json
import os


def count_lines(path):
    if not os.path.isfile(path):
        return 0
    with open(path) as f:
        return sum(1 for _ in f) - 1  # exclude header


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--grid", default="spot_coordinates.tsv")
    p.add_argument("--mutations", default="barcode_mutations.tsv")
    p.add_argument("--mapping", default="read_to_spot.tsv")
    p.add_argument("--reads-requested", type=int, default=0)
    p.add_argument("--out", default="summary_metrics.json")
    args = p.parse_args()

    n_spots = max(0, count_lines(args.grid))
    n_mut = max(0, count_lines(args.mutations))
    n_reads = max(0, count_lines(args.mapping))

    metrics = {
        "reads_total": args.reads_requested,
        "reads_simulated": n_reads,
        "spots": n_spots,
        "barcode_mutations_recorded": n_mut,
        "output": {
            "R1": "simulated_R1.fastq.gz",
            "R2": "simulated_R2.fastq.gz",
            "spot_coordinates": "spot_coordinates.tsv",
            "read_to_spot": "read_to_spot.tsv",
            "barcode_mutations": "barcode_mutations.tsv",
        },
    }
    with open(args.out, "w") as f:
        json.dump(metrics, f, indent=2)
    return 0


if __name__ == "__main__":
    exit(main())
