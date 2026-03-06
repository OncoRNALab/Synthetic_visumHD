#!/usr/bin/env python3
"""
Convert read_to_spot.tsv → CSV with columns barcode,x,y (one row per unique barcode).
Barcode is the true (unmutated) sequence from spot_coordinates.tsv.
50M reads map to 100k unique spots/barcodes; output has one row per spot.

Join key: spot_id (col 1 in read_to_spot.tsv, col 0 in spot_coordinates.tsv)

Usage:
    python3 read_to_spot_csv.py \
        --read_to_spot  <path/to/read_to_spot.tsv> \
        --spot_coords   <path/to/spot_coordinates.tsv> \
        --out           <path/to/output.csv>
"""

import argparse
import os
import sys

WRITE_BUF = 64 * 1024 * 1024   # 64 MB write buffer
LOG_EVERY  = 10_000_000


def build_spot_map(coords_path: str) -> dict:
    """Load spot_coordinates.tsv → {spot_id: 'true_barcode,x,y'}"""
    spot_map = {}
    with open(coords_path, "rb") as fh:
        header = fh.readline().decode().rstrip("\n").split("\t")
        try:
            i_id  = header.index("spot_id")
            i_bc  = header.index("true_barcode")
            i_x   = header.index("x")
            i_y   = header.index("y")
        except ValueError as e:
            sys.exit(f"[ERROR] spot_coordinates header missing column: {e}")

        for raw in fh:
            cols = raw.decode().rstrip("\n").split("\t")
            spot_map[cols[i_id]] = f"{cols[i_bc]},{cols[i_x]},{cols[i_y]}"

    print(f"  Loaded {len(spot_map):,} spots from {os.path.basename(coords_path)}")
    return spot_map


def convert(read_to_spot_path: str, spot_map: dict, out_path: str):
    seen = set()   # spot_ids already written (one row per unique barcode)
    n_written = 0
    n_missing = 0
    n_skipped = 0  # reads whose spot already written

    with open(read_to_spot_path, "rb") as fin, \
         open(out_path, "w", buffering=WRITE_BUF) as fout:

        header = fin.readline().decode().rstrip("\n").split("\t")
        try:
            i_spot = header.index("spot_id")
        except ValueError:
            sys.exit("[ERROR] read_to_spot.tsv has no 'spot_id' column")

        fout.write("barcode,x,y\n")

        for raw in fin:
            cols = raw.decode().rstrip("\n").split("\t")
            spot_id = cols[i_spot]
            if spot_id in seen:
                n_skipped += 1
                continue
            entry = spot_map.get(spot_id)
            if entry is None:
                n_missing += 1
                continue
            seen.add(spot_id)
            fout.write(entry)
            fout.write("\n")
            n_written += 1
            if n_written % 10_000 == 0:
                print(f"  {n_written:,} unique barcodes written ...")

    print(f"\nDone.")
    print(f"  Unique barcodes written : {n_written:,}")
    print(f"  Reads skipped (dup)     : {n_skipped:,}")
    print(f"  Missing (no coords)     : {n_missing:,}")
    print(f"  Output                  : {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--read_to_spot", required=True,
                        help="Path to read_to_spot.tsv")
    parser.add_argument("--spot_coords", required=True,
                        help="Path to spot_coordinates.tsv (must have true_barcode column)")
    parser.add_argument("--out", required=True,
                        help="Output CSV path")
    args = parser.parse_args()

    for p in (args.read_to_spot, args.spot_coords):
        if not os.path.exists(p):
            sys.exit(f"[ERROR] File not found: {p}")

    print(f"Building spot map ...")
    spot_map = build_spot_map(args.spot_coords)

    print(f"Converting {os.path.basename(args.read_to_spot)} ...")
    convert(args.read_to_spot, spot_map, args.out)


if __name__ == "__main__":
    main()
