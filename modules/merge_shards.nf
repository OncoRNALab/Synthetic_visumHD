process MERGE_SHARDS {
    input:
        path r1_files,        stageAs: "exp_R1_*.fastq"
        path r2_files,        stageAs: "exp_R2_*.fastq"
        path mutation_files,  stageAs: "barcode_mutations_*.tsv"
        path mapping_files,   stageAs: "read_to_spot_*.tsv"

    output:
        path "simulated_R1.fastq.gz", emit: r1_gz
        path "simulated_R2.fastq.gz", emit: r2_gz
        path "barcode_mutations.tsv", emit: barcode_mutations
        path "read_to_spot.tsv",      emit: read_to_spot

    script:
        """
        merge_shards.py \\
            --r1-list        ${r1_files} \\
            --r2-list        ${r2_files} \\
            --mutations-list ${mutation_files} \\
            --mapping-list   ${mapping_files} \\
            --out-r1         simulated_R1.fastq.gz \\
            --out-r2         simulated_R2.fastq.gz \\
            --out-mutations  barcode_mutations.tsv \\
            --out-mapping    read_to_spot.tsv
        """
}
