#!/usr/bin/env python3
"""
Convert expression matrix into a molecule table with UMI and PCR family sizes,
and produce transcript-level abundance counts.

Each row in molecules.tsv represents one unique cDNA molecule that will produce
one or more sequencing reads (duplicates) depending on family_size.

Outputs
-------
molecules.tsv        : molecule_id, spot_id, transcript_id, umi, family_size,
                       frag_start  (deterministic fragment start offset)
transcript_abundance.tsv : transcript_id, total_molecules
"""
import argparse
import hashlib
import random
import re
import sys


BASES = "ACGT"


# ---------------------------------------------------------------------------
# Gene → transcript mapping from transcript FASTA (GENCODE or Ensembl)
# ---------------------------------------------------------------------------

def build_gene_transcript_map(fasta_path):
    """
    Parse transcript FASTA headers and return two dicts:
      symbol_map : gene_symbol               -> list[transcript_id]
      geneid_map : ensembl_gene_id (no ver.) -> list[transcript_id]

    Supports two formats detected automatically:

    GENCODE (pipe-separated):
      >ENSMUST00000070533.5|ENSMUSG00000051951.6|...|Xkr4-201|Xkr4|3634|protein_coding|
      fields: [0]=transcript_id  [1]=gene_id  [5]=gene_symbol

    Ensembl (space-separated key:value):
      >ENSMUST00000103301.3 cdna ... gene:ENSMUSG00000076500.3 ... gene_symbol:Gm20730 ...
    """
    symbol_map = {}
    geneid_map = {}
    sym_re  = re.compile(r'gene_symbol:(\S+)')
    gene_re = re.compile(r'\bgene:(\S+)')

    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            header = line[1:].rstrip()

            if '|' in header:
                # ── GENCODE pipe-separated ────────────────────────────────
                fields = header.split('|')
                tid = fields[0]
                gid = fields[1].split('.')[0] if len(fields) > 1 else ''
                sym = fields[5]              if len(fields) > 5 else ''
            else:
                # ── Ensembl space-separated ───────────────────────────────
                tid   = header.split()[0]
                sym_m = sym_re.search(header)
                gen_m = gene_re.search(header)
                sym   = sym_m.group(1)              if sym_m else ''
                gid   = gen_m.group(1).split('.')[0] if gen_m else ''

            if sym:
                symbol_map.setdefault(sym, []).append(tid)
            if gid:
                geneid_map.setdefault(gid, []).append(tid)

    return symbol_map, geneid_map


def resolve_transcript(gene_name, symbol_map, geneid_map, rng):
    candidates = symbol_map.get(gene_name)
    if not candidates:
        base = gene_name.split('.')[0]
        candidates = geneid_map.get(base)
    if candidates:
        return rng.choice(candidates)
    return gene_name


# ---------------------------------------------------------------------------
# PCR family-size sampling
# ---------------------------------------------------------------------------

def sample_family_size(dist, geom_p, max_fam, rng):
    """Return an integer family size >= 1."""
    if dist == "const":
        return 1
    # geometric: P(X=k) = (1-p)^(k-1) * p  →  mean = 1/p
    size = 1
    while size < max_fam and rng.random() > geom_p:
        size += 1
    return size


def random_umi(length, rng):
    return "".join(rng.choices(BASES, k=length))


def deterministic_frag_start(mol_id, transcript_len, read_len):
    """Stable fragment start derived from molecule ID hash, not from rng state."""
    h = int(hashlib.md5(mol_id.encode()).hexdigest(), 16)
    max_start = max(0, transcript_len - read_len)
    return h % (max_start + 1) if max_start > 0 else 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Prepare molecule table and transcript abundances"
    )
    ap.add_argument("--matrix",          required=True,
                    help="expression_matrix.tsv (spots x genes)")
    ap.add_argument("--coordinates",     required=True,
                    help="spot_coordinates.tsv")
    ap.add_argument("--total-reads",     type=int, required=True,
                    help="Target total sequenced reads (after PCR expansion)")
    ap.add_argument("--seed",            type=int, default=42)
    ap.add_argument("--reference-fasta", default=None)
    # PCR model
    ap.add_argument("--pcr-enable",      action="store_true",
                    help="Assign variable PCR family sizes")
    ap.add_argument("--pcr-family-dist", default="geometric",
                    choices=["geometric", "const"])
    ap.add_argument("--pcr-geom-p",      type=float, default=0.5)
    ap.add_argument("--pcr-max-family",  type=int,   default=20)
    # UMI
    ap.add_argument("--umi-length",      type=int, default=10)
    # read lengths (for frag_start bound only)
    ap.add_argument("--read-length-r2",  type=int, default=100)
    # outputs
    ap.add_argument("--out-molecules",   default="molecules.tsv")
    ap.add_argument("--out-abundance",   default="transcript_abundance.tsv")
    args = ap.parse_args()

    rng = random.Random(args.seed)

    # Build gene → transcript map
    symbol_map, geneid_map = {}, {}
    if args.reference_fasta:
        sys.stderr.write("Building gene→transcript map from %s ...\n"
                         % args.reference_fasta)
        symbol_map, geneid_map = build_gene_transcript_map(args.reference_fasta)
        sys.stderr.write("  %d gene symbols, %d gene IDs mapped\n"
                         % (len(symbol_map), len(geneid_map)))

    # Read expression matrix → pool of (spot_id, gene_name) with multiplicity
    with open(args.matrix) as f:
        header   = next(f).rstrip().split("\t")
        gene_ids = [h for h in header[1:] if h]
        pool = []
        for line in f:
            parts = line.rstrip().split("\t")
            if not parts:
                continue
            spot_id = parts[0]
            for i, g in enumerate(gene_ids):
                try:
                    c = int(parts[i + 1]) if i + 1 < len(parts) else 0
                except ValueError:
                    c = 0
                for _ in range(c):
                    pool.append((spot_id, g))

    if not pool:
        with open(args.coordinates) as f:
            next(f)
            spots = [line.split("\t")[0] for line in f if line.strip()]
        spots    = spots or ["spot_1"]
        gene_ids = gene_ids or ["gene_0"]
        pool     = [(spots[0], gene_ids[0])]

    total_pool = len(pool)

    # ------------------------------------------------------------------
    # Determine how many unique molecules we need.
    # With PCR enabled we over-sample molecules then trim so that the
    # sum of family sizes equals exactly total_reads.
    # ------------------------------------------------------------------
    target_reads = args.total_reads

    if not args.pcr_enable:
        # Each molecule = 1 read; sample exactly target_reads molecules
        n_molecules = target_reads
        family_sizes_fixed = [1] * n_molecules
        indices = [rng.randint(0, total_pool - 1) for _ in range(n_molecules)]
    else:
        # Draw molecules until cumulative reads >= target_reads, then trim last
        mean_fam = 1.0 / args.pcr_geom_p if args.pcr_family_dist == "geometric" else 1.0
        n_estimate = max(1, int(target_reads / mean_fam * 1.05))

        indices      = []
        family_sizes = []
        cum_reads    = 0

        # Draw in batches
        batch = max(1000, n_estimate)
        while cum_reads < target_reads:
            batch_idx = [rng.randint(0, total_pool - 1) for _ in range(batch)]
            for idx in batch_idx:
                fs = sample_family_size(args.pcr_family_dist,
                                        args.pcr_geom_p,
                                        args.pcr_max_family, rng)
                remaining = target_reads - cum_reads
                if fs > remaining:
                    fs = remaining
                indices.append(idx)
                family_sizes.append(fs)
                cum_reads += fs
                if cum_reads >= target_reads:
                    break
            if cum_reads >= target_reads:
                break

        family_sizes_fixed = family_sizes

    n_molecules = len(indices)
    sys.stderr.write("Molecules sampled: %d  (target reads: %d)\n"
                     % (n_molecules, target_reads))

    # Diagnostic: count how many distinct gene names will resolve to transcripts
    if args.reference_fasta and pool:
        unique_genes = set(g for _, g in pool)
        resolved_genes = sum(
            1 for g in unique_genes
            if symbol_map.get(g) or geneid_map.get(g.split('.')[0])
        )
        sys.stderr.write(
            "Gene resolution: %d / %d unique genes found in reference "
            "(%.1f%%) — unresolved genes use fallback sequence\n"
            % (resolved_genes, len(unique_genes),
               100.0 * resolved_genes / max(len(unique_genes), 1))
        )

    # ------------------------------------------------------------------
    # Write molecules.tsv and transcript_abundance.tsv
    # ------------------------------------------------------------------
    transcript_counts = {}
    matched = 0

    with open(args.out_molecules, "w") as mol_out:
        mol_out.write("molecule_id\tspot_id\ttranscript_id\tumi\t"
                      "family_size\tfrag_start\n")
        for i, idx in enumerate(indices):
            spot_id, gene_name = pool[idx]
            tid = resolve_transcript(gene_name, symbol_map, geneid_map, rng)
            if tid != gene_name:
                matched += 1
            mol_id    = "mol_%d" % i
            umi       = random_umi(args.umi_length, rng)
            fs        = family_sizes_fixed[i]
            frag_st   = deterministic_frag_start(
                            mol_id, 10000, args.read_length_r2)
            mol_out.write("%s\t%s\t%s\t%s\t%d\t%d\n"
                          % (mol_id, spot_id, tid, umi, fs, frag_st))
            transcript_counts[tid] = transcript_counts.get(tid, 0) + fs

    if args.reference_fasta:
        sys.stderr.write("Transcript resolution: %d / %d molecules mapped "
                         "(%.1f%%)\n"
                         % (matched, n_molecules,
                            100.0 * matched / max(n_molecules, 1)))

    with open(args.out_abundance, "w") as ab_out:
        ab_out.write("transcript_id\ttotal_reads\n")
        for tid, count in sorted(transcript_counts.items()):
            ab_out.write("%s\t%d\n" % (tid, count))

    return 0


if __name__ == "__main__":
    sys.exit(main())
