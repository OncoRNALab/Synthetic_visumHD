process ATTACH_SPATIAL_BARCODES {
    input:
        tuple val(shard_id), path(mol_r1), path(mol_r2), path(molecules_shard)
        path  spot_coordinates
        val   barcode_length
        val   umi_length
        val   read_length_r1

    output:
        tuple val(shard_id), path("barcoded_R1.fastq"), path("barcoded_R2.fastq"), emit: attach_out

    script:
        """
        attach_spatial_barcodes.py \\
            --r1             ${mol_r1} \\
            --r2             ${mol_r2} \\
            --molecules      ${molecules_shard} \\
            --coordinates    ${spot_coordinates} \\
            --barcode-length ${barcode_length} \\
            --umi-length     ${umi_length} \\
            --read-length-r1 ${read_length_r1} \\
            --out-r1         barcoded_R1.fastq \\
            --out-r2         barcoded_R2.fastq \\
            --out-mapping    mol_read_map.tsv
        """
}
