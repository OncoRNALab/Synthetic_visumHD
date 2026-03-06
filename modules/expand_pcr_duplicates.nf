process EXPAND_PCR_DUPLICATES {
    input:
        tuple val(shard_id), path(barcoded_r1), path(barcoded_r2), path(molecules_shard)
        path  spot_coordinates   // spot_coordinates.tsv with true_barcode column
        val   barcode_length
        val   seq_sub_rate       // sequencing error: substitution rate per base
        val   seq_ins_rate       // sequencing error: insertion rate per base
        val   seq_del_rate       // sequencing error: deletion rate per base
        val   seed

    output:
        path "exp_R1.fastq",          emit: exp_r1
        path "exp_R2.fastq",          emit: exp_r2
        path "barcode_mutations.tsv", emit: barcode_mutations
        path "read_to_spot.tsv",      emit: read_to_spot

    script:
        """
        expand_pcr_duplicates.py \\
            --r1              ${barcoded_r1} \\
            --r2              ${barcoded_r2} \\
            --molecules       ${molecules_shard} \\
            --coordinates     ${spot_coordinates} \\
            --barcode-length  ${barcode_length} \\
            --seq-sub-rate    ${seq_sub_rate} \\
            --seq-ins-rate    ${seq_ins_rate} \\
            --seq-del-rate    ${seq_del_rate} \\
            --seed            ${seed} \\
            --shard-id        ${shard_id} \\
            --out-r1          exp_R1.fastq \\
            --out-r2          exp_R2.fastq \\
            --out-mutations   barcode_mutations.tsv \\
            --out-mapping     read_to_spot.tsv
        """
}
