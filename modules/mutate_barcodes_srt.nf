process MUTATE_BARCODES_SRT {
    input:
        path barcoded_r1
        path barcoded_r2
        val barcode_length
        val sub_rate
        val ins_rate
        val del_rate
        val seed

    output:
        path "mutated_R1.fastq", emit: mutated_r1
        path "mutated_R2.fastq", emit: mutated_r2
        path "barcode_mutations.tsv", emit: barcode_mutations

    script:
        """
        mutate_barcodes_srt.py \\
            --r1 barcoded_R1.fastq \\
            --r2 barcoded_R2.fastq \\
            --barcode-length ${barcode_length} \\
            --sub-rate ${sub_rate} \\
            --ins-rate ${ins_rate} \\
            --del-rate ${del_rate} \\
            --seed ${seed} \\
            --out-r1 mutated_R1.fastq \\
            --out-r2 mutated_R2.fastq \\
            --out-mutations barcode_mutations.tsv
        """
}
