#!/usr/bin/env python3
"""
Simulate paired-end cDNA FASTQ from a molecules shard.

Input
-----
molecules shard TSV columns:
    molecule_id, spot_id, transcript_id, umi, family_size, frag_start

One FASTQ record is written per molecule (family_size handled later by
expand_pcr_duplicates.py).  The fragment position is taken from frag_start
(deterministic per molecule) so PCR duplicates will later share the same
cDNA sequence.

Output
------
mol_R1.fastq, mol_R2.fastq  – one record per molecule (not per read)
"""
import argparse
import collections
import hashlib
import os
import random
import shutil
import subprocess
import sys
import tempfile

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
                raw  = line[1:].split()[0].strip()
                # GENCODE headers are pipe-separated; keep only the transcript ID
                name = raw.split('|')[0] if '|' in raw else raw
                seq = []
            else:
                seq.append(line.rstrip())
        if name:
            out[name] = "".join(seq)
    return out


def frag_start_for_seq(mol_id, seq_len, L1, L2):
    """
    Recompute a deterministic fragment start using the ACTUAL transcript length.
    This prevents the clamping artifact that occurs when the precomputed
    frag_start (based on a hardcoded transcript_len=10000) exceeds the real
    sequence length, which would force all molecules of short transcripts to
    the same start position and produce millions of identical reads.
    """
    need      = L1 + L2
    max_start = max(0, seq_len - need)
    if max_start == 0:
        return 0
    h = int(hashlib.md5(mol_id.encode()).hexdigest(), 16)
    return h % (max_start + 1)


def get_seq_slice(seq, mol_id, L1, L2, rng):
    """Extract R1/R2 sequence from a transcript.
    Fragment start is recomputed from the actual sequence length."""
    need = L1 + L2
    if len(seq) < need:
        pad = "".join(rng.choices(BASES, k=need - len(seq) + 50))
        seq = seq + pad
    start = frag_start_for_seq(mol_id, len(seq), L1, L2)
    s1 = seq[start : start + L1]
    # Make R2 match the sense (transcript) strand like the library design.
    s2 = seq[start + L1 : start + L1 + L2]
    if len(s2) < L2:
        s2 += "".join(rng.choices(BASES, k=L2 - len(s2)))
    return s1, s2


def python_simulate(rows, ref, L1, L2, seed, out_r1, out_r2):
    """rows: list of (mol_id, spot_id, tid, umi, family_size, frag_start)"""
    rng = random.Random(seed)
    fallback_seq = list(ref.values())[0] if ref else ""
    with open(out_r1, "w") as r1, open(out_r2, "w") as r2:
        for mol_id, _spot, tid, _umi, _fs, _frag_start in rows:
            seq = ref.get(tid) or fallback_seq
            s1, s2 = get_seq_slice(seq, mol_id, L1, L2, rng)
            q1 = "F" * len(s1)
            q2 = "F" * len(s2)
            r1.write("@%s\n%s\n+\n%s\n" % (mol_id, s1, q1))
            r2.write("@%s\n%s\n+\n%s\n" % (mol_id, s2, q2))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--molecules",  required=True, help="molecules shard TSV")
    ap.add_argument("--reference",  required=True, help="Ensembl cDNA FASTA")
    ap.add_argument("--use-art",    action="store_true")
    ap.add_argument("--seed",       type=int, default=42)
    ap.add_argument("--len-r1",     type=int, default=50)
    ap.add_argument("--len-r2",     type=int, default=100)
    ap.add_argument("--out-r1",     default="mol_R1.fastq")
    ap.add_argument("--out-r2",     default="mol_R2.fastq")
    args = ap.parse_args()

    ref = load_fasta(args.reference)
    if not ref:
        sys.stderr.write("No sequences in reference\n")
        return 1

    rows = []
    with open(args.molecules) as f:
        hdr = next(f).rstrip().split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for line in f:
            t = line.rstrip().split("\t")
            if len(t) < 4:
                continue
            mol_id    = t[col["molecule_id"]]
            spot_id   = t[col["spot_id"]]
            tid       = t[col["transcript_id"]]
            umi       = t[col["umi"]]
            family_sz = int(t[col.get("family_size", 4)]) if "family_size" in col else 1
            frag_st   = int(t[col.get("frag_start",  5)]) if "frag_start"  in col else 0
            rows.append((mol_id, spot_id, tid, umi, family_sz, frag_st))

    if not rows:
        sys.stderr.write("No molecules in shard, writing empty FASTQ\n")
        open(args.out_r1, "w").close()
        open(args.out_r2, "w").close()
        return 0

    art_bin = shutil.which("art_illumina")
    if not args.use_art or not art_bin:
        return python_simulate(rows, ref, args.len_r1, args.len_r2,
                               args.seed, args.out_r1, args.out_r2)

    # ----------------------------------------------------------------
    # ART path: generate one read per molecule using its frag_start
    # ----------------------------------------------------------------
    rng = random.Random(args.seed)
    fallback_seq = list(ref.values())[0] if ref else ""

    by_transcript = collections.defaultdict(list)
    for row in rows:
        by_transcript[row[2]].append(row)

    with open(args.out_r1, "w") as r1out, open(args.out_r2, "w") as r2out:
        with tempfile.TemporaryDirectory(prefix="art_") as tmp:
            for tid, mol_list in by_transcript.items():
                seq = ref.get(tid) or fallback_seq
                need = args.len_r1 + args.len_r2
                if len(seq) < need:
                    seq += "".join(rng.choices(BASES, k=need - len(seq) + 100))

                tfa = os.path.join(tmp, "t.fa")
                with open(tfa, "w") as f:
                    f.write(">%s\n%s\n" % (tid, seq))

                prefix = os.path.join(tmp, "out")
                cmd = [
                    art_bin, "-i", tfa,
                    "-l", str(args.len_r1),
                    "-c", str(len(mol_list)),
                    "-p", "-m", "200", "-s", "10",
                    "-o", prefix,
                    "--rndSeed", str(args.seed + abs(hash(tid)) % 100000),
                ]
                ret = subprocess.run(cmd, capture_output=True, cwd=tmp)
                if ret.returncode != 0 or not os.path.isfile(prefix + "1.fq"):
                    for row in mol_list:
                        s1, s2 = get_seq_slice(seq, row[0],
                                               args.len_r1, args.len_r2, rng)
                        r1out.write("@%s\n%s\n+\n%s\n" % (row[0], s1, "F"*len(s1)))
                        r2out.write("@%s\n%s\n+\n%s\n" % (row[0], s2, "F"*len(s2)))
                    continue

                with open(prefix + "1.fq") as r1f, open(prefix + "2.fq") as r2f:
                    for row in mol_list:
                        h1 = r1f.readline(); s1 = r1f.readline()
                        r1f.readline();      q1 = r1f.readline()
                        h2 = r2f.readline(); s2 = r2f.readline()
                        r2f.readline();      q2 = r2f.readline()
                        if not h1:
                            break
                        r1out.write("@%s\n%s+\n%s" % (row[0], s1, q1))
                        r2out.write("@%s\n%s+\n%s" % (row[0], s2, q2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
