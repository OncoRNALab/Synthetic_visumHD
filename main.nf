#!/usr/bin/env nextflow
/*
 * Biologically realistic Visium-HD-like spatial FASTQ simulation
 *
 * Workflow stages
 * ---------------
 * 1. RUN_SRT_SIM               – spatial expression via R (sCCIgen / NB fallback)
 * 2. PREPARE_TRANSCRIPT_ABUNDANCE – build molecule table (UMI, PCR family sizes)
 * 3. SHARD_MOLECULES            – split molecule table for parallel execution
 * 4. Per-shard (fan-out, joined on shard_id):
 *    a. SIMULATE_CDNA_READS_SRT – one cDNA FASTQ record per molecule
 *    b. ATTACH_SPATIAL_BARCODES – prepend barcode + UMI to R1
 *    c. EXPAND_PCR_DUPLICATES   – replicate by family_size, apply barcode errors
 * 5. MERGE_SHARDS               – concatenate shard outputs → final FASTQs + TSVs
 * 6. FINALIZE_OUTPUTS           – publish files + write metrics JSON
 */
nextflow.enable.dsl = 2

include { RUN_SRT_SIM }                  from './modules/run_srt_sim.nf'
include { PREPARE_TRANSCRIPT_ABUNDANCE } from './modules/prepare_transcript_abundance.nf'
include { SHARD_MOLECULES }              from './modules/shard_molecules.nf'
include { SIMULATE_CDNA_READS_SRT }      from './modules/simulate_cdna_reads_srt.nf'
include { ATTACH_SPATIAL_BARCODES }      from './modules/attach_spatial_barcodes.nf'
include { EXPAND_PCR_DUPLICATES }        from './modules/expand_pcr_duplicates.nf'
include { MERGE_SHARDS }                 from './modules/merge_shards.nf'
include { FINALIZE_OUTPUTS }             from './modules/finalize_outputs.nf'

def check_params() {
    if (params.total_reads <= 0)
        exit 1, "total_reads must be positive (got ${params.total_reads})"
    if (params.seed == null)
        exit 1, "seed must be set"
    def ref = file(params.reference_fasta)
    if (!ref.exists())
        exit 1, "reference_fasta not found: ${params.reference_fasta}"
    if (params.synth_sub_rate < 0 || params.synth_ins_rate < 0 || params.synth_del_rate < 0)
        exit 1, "Synthesis error rates (synth_*_rate) must be >= 0"
    if (params.seq_sub_rate < 0 || params.seq_ins_rate < 0 || params.seq_del_rate < 0)
        exit 1, "Sequencing error rates (seq_*_rate) must be >= 0"
    if (params.sim_shards < 1)
        exit 1, "sim_shards must be >= 1"
}

check_params()

workflow {

    // ------------------------------------------------------------------
    // 1. Spatial expression simulation
    // ------------------------------------------------------------------
    def feature_h5_val     = params.feature_h5      ?: ''
    def required_genes_val = params.required_genes  ?: ''
    RUN_SRT_SIM(
        params.seed,
        params.n_genes,
        params.n_spots,
        feature_h5_val,
        required_genes_val,
        params.synth_sub_rate,
        params.synth_ins_rate,
        params.synth_del_rate
    )

    // ------------------------------------------------------------------
    // 2. Molecule table: UMIs, PCR family sizes, transcript assignments
    // ------------------------------------------------------------------
    PREPARE_TRANSCRIPT_ABUNDANCE(
        RUN_SRT_SIM.out.expression_matrix,
        RUN_SRT_SIM.out.spot_coordinates,
        params.total_reads,
        params.seed,
        file(params.reference_fasta),
        params.pcr_enable,
        params.pcr_family_dist,
        params.pcr_geom_p,
        params.pcr_max_family,
        params.umi_length,
        params.read_length_r2
    )

    // ------------------------------------------------------------------
    // 3. Shard the molecule table
    //    Each element: [ shard_id (Int), shard_file (Path) ]
    // ------------------------------------------------------------------
    SHARD_MOLECULES(
        PREPARE_TRANSCRIPT_ABUNDANCE.out.molecules,
        params.sim_shards
    )

    shards_ch = SHARD_MOLECULES.out.shards
        .flatten()
        .map { f -> [ f.baseName.replaceAll(/^shard_/, '').toInteger(), f ] }

    // ------------------------------------------------------------------
    // 4a. Simulate cDNA reads – one record per molecule
    //     Input  : [ shard_id, shard_file ]
    //     Output : [ shard_id, mol_R1, mol_R2 ]
    // ------------------------------------------------------------------
    SIMULATE_CDNA_READS_SRT(
        shards_ch,
        file(params.reference_fasta),
        params.use_art,
        params.seed,
        params.read_length_r1,
        params.read_length_r2
    )

    // ------------------------------------------------------------------
    // 4b. Attach barcodes + UMIs
    //     Join simulation outputs with the shard file on shard_id
    //     Input  : [ shard_id, mol_R1, mol_R2, shard_file ]
    //     Output : [ shard_id, barcoded_R1, barcoded_R2 ]
    // ------------------------------------------------------------------
    attach_in = SIMULATE_CDNA_READS_SRT.out.sim_out
        .join( shards_ch )
        // sim_out: [id, r1, r2]  join shards: [id, shard] → [id, r1, r2, shard]
        .map { id, r1, r2, shard -> [ id, r1, r2, shard ] }

    ATTACH_SPATIAL_BARCODES(
        attach_in,
        RUN_SRT_SIM.out.spot_coordinates,
        params.barcode_length,
        params.umi_length,
        params.read_length_r1
    )

    // ------------------------------------------------------------------
    // 4c. PCR expansion + barcode mutation
    //     Join attach outputs with shard file on shard_id
    //     Input  : [ shard_id, barcoded_R1, barcoded_R2, shard_file ]
    //     Output : exp_R1, exp_R2, barcode_mutations, read_to_spot
    // ------------------------------------------------------------------
    expand_in = ATTACH_SPATIAL_BARCODES.out.attach_out
        .join( shards_ch )
        // attach_out: [id, r1, r2]  join shards: [id, shard] → [id, r1, r2, shard]
        .map { id, r1, r2, shard -> [ id, r1, r2, shard ] }

    EXPAND_PCR_DUPLICATES(
        expand_in,
        RUN_SRT_SIM.out.spot_coordinates,
        params.barcode_length,
        params.seq_sub_rate,
        params.seq_ins_rate,
        params.seq_del_rate,
        params.seed
    )

    // ------------------------------------------------------------------
    // 5. Merge all shard outputs into final files
    // ------------------------------------------------------------------
    MERGE_SHARDS(
        EXPAND_PCR_DUPLICATES.out.exp_r1.collect(),
        EXPAND_PCR_DUPLICATES.out.exp_r2.collect(),
        EXPAND_PCR_DUPLICATES.out.barcode_mutations.collect(),
        EXPAND_PCR_DUPLICATES.out.read_to_spot.collect()
    )

    // ------------------------------------------------------------------
    // 6. Publish results + write metrics JSON
    // ------------------------------------------------------------------
    FINALIZE_OUTPUTS(
        MERGE_SHARDS.out.r1_gz,
        MERGE_SHARDS.out.r2_gz,
        MERGE_SHARDS.out.barcode_mutations,
        RUN_SRT_SIM.out.spot_coordinates,
        MERGE_SHARDS.out.read_to_spot,
        params.total_reads,
        params.out_dir
    )
}
