#!/usr/bin/env python3
"""
Expand molecule-level FASTQ into read-level FASTQ by replicating each record
family_size times (PCR duplicates) and applying independent barcode errors to
every copy.

Input
-----
  --r1            barcoded_R1.fastq  (one record per molecule, barcode already attached)
  --r2            barcoded_R2.fastq
  --molecules     molecules shard TSV  (columns: molecule_id, family_size, ...)
  --barcode-length  length of barcode prefix in R1
  --sub-rate / --ins-rate / --del-rate  per-base error probabilities
  --seed

Output
------
  exp_R1.fastq, exp_R2.fastq  – one record per sequenced read
  barcode_mutations.tsv       – read_id, molecule_id, original_barcode,
                                mutated_barcode, mutation_events
  read_to_spot.tsv            – read_id, spot_id, transcript_id, umi,
                                molecule_id
"""
import argparse
import random
import sys

BASES = "ACGT"


def mutate_region(seq, sub_rate, ins_rate, del_rate, rng):
    out   = []
    desc  = []
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
    ap = argparse.ArgumentParser(
        description="Expand molecule FASTQs to read-level with PCR duplication"
    )
    ap.add_argument("--r1",             required=True)
    ap.add_argument("--r2",             required=True)
    ap.add_argument("--molecules",      required=True)
    ap.add_argument("--coordinates",    default=None,
                    help="spot_coordinates.tsv produced by run_sccigen_sim.R. "
                         "When provided, the true_barcode column is used as "
                         "original_barcode in barcode_mutations.tsv so that "
                         "ground truth is always the perfect whitelist barcode, "
                         "not the synthesis-mutated one.")
    ap.add_argument("--barcode-length", type=int,   default=36)
    # Sequencing-error rates (applied independently per read, models base-calling
    # errors after PCR).  These should be much lower than synthesis rates.
    ap.add_argument("--seq-sub-rate",   type=float, default=0.001,
                    help="Per-base substitution rate during sequencing (default 0.001)")
    ap.add_argument("--seq-ins-rate",   type=float, default=0.0001,
                    help="Per-base insertion rate during sequencing (default 0.0001)")
    ap.add_argument("--seq-del-rate",   type=float, default=0.0001,
                    help="Per-base deletion rate during sequencing (default 0.0001)")
    ap.add_argument("--seed",           type=int,   default=42)
    ap.add_argument("--shard-id",       default=None,
                    help="Shard identifier (e.g. 0, 1, ...) prepended to read IDs "
                         "as 'read_s{shard_id}_{n}' to ensure global uniqueness "
                         "across shards. If omitted, IDs are 'read_{n}'.")
    ap.add_argument("--out-r1",         default="exp_R1.fastq")
    ap.add_argument("--out-r2",         default="exp_R2.fastq")
    ap.add_argument("--out-mutations",  default="barcode_mutations.tsv")
    ap.add_argument("--out-mapping",    default="read_to_spot.tsv")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    L   = args.barcode_length
    shard_prefix = ("s%s_" % args.shard_id) if args.shard_id is not None else ""

    # Load molecule metadata: molecule_id → (spot_id, transcript_id, umi, family_size)
    mol_meta = {}
    with open(args.molecules) as f:
        hdr = next(f).rstrip().split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for line in f:
            t = line.rstrip().split("\t")
            if not t or len(t) < 4:
                continue
            mid = t[col["molecule_id"]]
            mol_meta[mid] = (
                t[col["spot_id"]],
                t[col["transcript_id"]],
                t[col["umi"]],
                int(t[col.get("family_size", 4)]) if "family_size" in col else 1,
            )

    # Load true (pre-synthesis-error) barcodes per spot from spot_coordinates.tsv.
    # Column layout: spot_id, barcode (synthesis-mutated), x, y, true_barcode.
    # Falls back to using s1[:L] (synthesis-mutated barcode from FASTQ) when the
    # file is absent or the true_barcode column is missing.
    spot_true_bc = {}
    if args.coordinates:
        try:
            with open(args.coordinates) as f:
                hdr = next(f).rstrip().split("\t")
                if "true_barcode" in hdr:
                    idx_sid = hdr.index("spot_id")
                    idx_tbc = hdr.index("true_barcode")
                    for line in f:
                        t = line.rstrip().split("\t")
                        if len(t) > idx_tbc:
                            spot_true_bc[t[idx_sid]] = t[idx_tbc]
                    sys.stderr.write(
                        "Loaded true_barcode for %d spots from %s\n"
                        % (len(spot_true_bc), args.coordinates))
                else:
                    sys.stderr.write(
                        "WARNING: true_barcode column not found in %s; "
                        "using synthesis-mutated barcode as original_barcode\n"
                        % args.coordinates)
        except OSError as e:
            sys.stderr.write("WARNING: cannot open coordinates file: %s\n" % e)

    read_counter = 0

    with open(args.r1) as r1_in, \
         open(args.r2) as r2_in, \
         open(args.out_r1, "w") as r1_out, \
         open(args.out_r2, "w") as r2_out, \
         open(args.out_mutations, "w") as mut_out, \
         open(args.out_mapping,   "w") as map_out:

        mut_out.write("read_id\tmolecule_id\toriginal_barcode\t"
                      "synthesis_barcode\tmutated_barcode\tmutation_events\n")
        map_out.write("read_id\tspot_id\ttranscript_id\tumi\tmolecule_id\n")

        while True:
            h1 = r1_in.readline()
            if not h1:
                break
            s1 = r1_in.readline().rstrip()
            r1_in.readline()
            q1 = r1_in.readline().rstrip()
            h2 = r2_in.readline()
            s2 = r2_in.readline()
            r2_in.readline()
            q2 = r2_in.readline()

            mol_id = h1[1:].split()[0].strip()
            spot_id, tid, umi, family_size = mol_meta.get(
                mol_id, ("unknown", "unknown", "N" * 10, 1))

            # synthesis_bc: barcode as printed on the array (from FASTQ, already
            #               carries synthesis errors applied in run_sccigen_sim.R)
            # orig_bc     : true whitelist barcode (ground truth for evaluation)
            synthesis_bc = s1[:L]
            orig_bc      = spot_true_bc.get(spot_id, synthesis_bc)
            rest_seq     = s1[L:]
            rest_qual    = q1[L:] if len(q1) >= L else ""

            for dup in range(family_size):
                read_id = "read_%s%d" % (shard_prefix, read_counter)
                read_counter += 1

                # Apply sequencing errors independently to each read
                mut_bc, desc_list = mutate_region(
                    synthesis_bc,
                    args.seq_sub_rate, args.seq_ins_rate, args.seq_del_rate, rng)
                events_str = ";".join(desc_list) if desc_list else "none"

                new_s1 = mut_bc + rest_seq
                new_q1 = "F" * len(mut_bc) + rest_qual
                if len(new_q1) > len(new_s1):
                    new_q1 = new_q1[:len(new_s1)]
                elif len(new_q1) < len(new_s1):
                    new_q1 += "F" * (len(new_s1) - len(new_q1))

                r1_out.write("@%s\n%s\n+\n%s\n" % (read_id, new_s1, new_q1))
                r2_out.write("@%s\n%s+\n%s" % (read_id, s2, q2))
                # original_barcode = perfect whitelist barcode (ground truth)
                # synthesis_barcode = after array synthesis errors (shared by spot)
                # mutated_barcode   = after sequencing errors (per-read)
                mut_out.write("%s\t%s\t%s\t%s\t%s\t%s\n"
                              % (read_id, mol_id, orig_bc,
                                 synthesis_bc, mut_bc, events_str))
                map_out.write("%s\t%s\t%s\t%s\t%s\n"
                              % (read_id, spot_id, tid, umi, mol_id))

    sys.stderr.write("Expanded to %d reads\n" % read_counter)
    return 0


if __name__ == "__main__":
    sys.exit(main())
