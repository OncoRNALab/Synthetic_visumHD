process SIMULATE_CDNA_READS_SRT {
    input:
        tuple val(shard_id), path(molecules_shard)
        path  reference_fasta
        val   use_art
        val   seed
        val   len_r1
        val   len_r2

    output:
        tuple val(shard_id), path("mol_R1.fastq"), path("mol_R2.fastq"), emit: sim_out

    script:
        def use_art_flag = use_art ? '--use-art' : ''
        """
        run_art_simulation.py \\
            --molecules ${molecules_shard} \\
            --reference ${reference_fasta} \\
            ${use_art_flag} \\
            --seed   ${seed} \\
            --len-r1 ${len_r1} \\
            --len-r2 ${len_r2} \\
            --out-r1 mol_R1.fastq \\
            --out-r2 mol_R2.fastq
        """
}
