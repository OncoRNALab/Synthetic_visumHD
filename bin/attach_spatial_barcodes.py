#!/usr/bin/env python3
"""
Insert true spatial barcode and UMI at the start of R1; preserve pairing.

UMIs come from the molecules shard (column 'umi') – this ensures that all PCR
duplicates of the same molecule later carry identical UMI sequences.

R1 layout after attachment:
    [barcode_length nt barcode][umi_length nt UMI][cDNA suffix]

Outputs
-------
barcoded_R1.fastq
barcoded_R2.fastq  (pass-through)
mol_read_map.tsv   : molecule_id, spot_id, barcode, umi
"""
import argparse
import sys

BASES = "ACGT"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1",             required=True)
    ap.add_argument("--r2",             required=True)
    ap.add_argument("--molecules",      required=True,
                    help="molecules shard TSV (molecule_id, spot_id, "
                         "transcript_id, umi, family_size, frag_start)")
    ap.add_argument("--coordinates",    required=True,
                    help="spot_coordinates.tsv (spot_id, barcode, x, y)")
    ap.add_argument("--barcode-length", type=int, default=36)
    ap.add_argument("--umi-length",     type=int, default=10)
    ap.add_argument("--read-length-r1", type=int, default=50)
    ap.add_argument("--out-r1",         default="barcoded_R1.fastq")
    ap.add_argument("--out-r2",         default="barcoded_R2.fastq")
    ap.add_argument("--out-mapping",    default="mol_read_map.tsv")
    args = ap.parse_args()

    # spot_id → barcode
    spot_barcode = {}
    with open(args.coordinates) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) >= 2:
                spot_barcode[parts[0]] = parts[1]

    if not spot_barcode:
        sys.stderr.write("No spots loaded from coordinates file\n")
        return 1
    fallback_spot    = next(iter(spot_barcode))
    fallback_barcode = spot_barcode[fallback_spot]

    # molecule_id → (spot_id, umi)
    mol_info = {}
    with open(args.molecules) as f:
        hdr = next(f).rstrip().split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for line in f:
            t = line.rstrip().split("\t")
            if len(t) < 4:
                continue
            mol_id  = t[col["molecule_id"]]
            spot_id = t[col["spot_id"]]
            umi     = t[col["umi"]]
            mol_info[mol_id] = (spot_id, umi)

    bc_len    = args.barcode_length
    umi_len   = args.umi_length
    prefix_l  = bc_len + umi_len
    suffix_l  = max(0, args.read_length_r1 - prefix_l)

    with open(args.r1) as r1_in, open(args.r2) as r2_in, \
         open(args.out_r1, "w") as r1_out, open(args.out_r2, "w") as r2_out, \
         open(args.out_mapping, "w") as map_out:

        map_out.write("molecule_id\tspot_id\tbarcode\tumi\n")

        while True:
            h1 = r1_in.readline()
            if not h1:
                break
            s1 = r1_in.readline().rstrip()
            r1_in.readline()          # '+'
            q1 = r1_in.readline().rstrip()
            h2 = r2_in.readline()
            s2 = r2_in.readline()
            r2_in.readline()
            q2 = r2_in.readline()

            mol_id  = h1[1:].split()[0].strip()
            spot_id, umi = mol_info.get(mol_id, (fallback_spot, "A" * umi_len))
            barcode = spot_barcode.get(spot_id, fallback_barcode)

            suffix_seq  = s1[prefix_l : prefix_l + suffix_l] if len(s1) > prefix_l else ""
            suffix_qual = q1[prefix_l : prefix_l + suffix_l] if len(q1) > prefix_l else ""

            new_s1 = barcode + umi + suffix_seq
            new_q1 = "F" * bc_len + "F" * umi_len + suffix_qual
            # equalise lengths
            if len(new_q1) > len(new_s1):
                new_q1 = new_q1[:len(new_s1)]
            elif len(new_q1) < len(new_s1):
                new_q1 += "F" * (len(new_s1) - len(new_q1))

            r1_out.write(h1)
            r1_out.write(new_s1 + "\n+\n" + new_q1 + "\n")
            r2_out.write(h2 + s2 + "+\n" + q2)
            map_out.write("%s\t%s\t%s\t%s\n" % (mol_id, spot_id, barcode, umi))

    return 0


if __name__ == "__main__":
    sys.exit(main())
