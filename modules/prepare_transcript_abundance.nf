process PREPARE_TRANSCRIPT_ABUNDANCE {
    input:
        path expression_matrix
        path spot_coordinates
        val  total_reads
        val  seed
        path reference_fasta
        val  pcr_enable
        val  pcr_family_dist
        val  pcr_geom_p
        val  pcr_max_family
        val  umi_length
        val  read_length_r2

    output:
        path "molecules.tsv",            emit: molecules
        path "transcript_abundance.tsv", emit: transcript_abundance

    script:
        def pcr_flag = pcr_enable ? '--pcr-enable' : ''
        """
        prepare_transcript_abundance.py \\
            --matrix            ${expression_matrix} \\
            --coordinates       ${spot_coordinates} \\
            --total-reads       ${total_reads} \\
            --seed              ${seed} \\
            --reference-fasta   ${reference_fasta} \\
            ${pcr_flag} \\
            --pcr-family-dist   ${pcr_family_dist} \\
            --pcr-geom-p        ${pcr_geom_p} \\
            --pcr-max-family    ${pcr_max_family} \\
            --umi-length        ${umi_length} \\
            --read-length-r2    ${read_length_r2} \\
            --out-molecules     molecules.tsv \\
            --out-abundance     transcript_abundance.tsv
        """
}
