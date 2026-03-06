process FINALIZE_OUTPUTS {
    publishDir = [path: "${params.out_dir}", mode: 'copy']

    input:
        path r1_gz
        path r2_gz
        path barcode_mutations
        path spot_coordinates
        path read_to_spot
        val  total_reads
        val  out_dir

    output:
        path "simulated_R1.fastq.gz",  emit: r1_gz
        path "simulated_R2.fastq.gz",  emit: r2_gz
        path "barcode_mutations.tsv",  emit: barcode_mutations
        path "spot_coordinates.tsv",   emit: spot_coordinates
        path "read_to_spot.tsv",       emit: read_to_spot
        path "barcodes_coords.csv",    emit: barcodes_coords
        path "summary_metrics.json",   emit: summary_metrics

    script:
        """
        # FASTQs and TSVs are already in working dir under their final names
        # (Nextflow stages them); just symlink/copy if names differ
        [ "${r1_gz}"             = "simulated_R1.fastq.gz"    ] || cp ${r1_gz}             simulated_R1.fastq.gz
        [ "${r2_gz}"             = "simulated_R2.fastq.gz"    ] || cp ${r2_gz}             simulated_R2.fastq.gz
        [ "${barcode_mutations}" = "barcode_mutations.tsv"    ] || cp ${barcode_mutations} barcode_mutations.tsv
        [ "${spot_coordinates}"  = "spot_coordinates.tsv"     ] || cp ${spot_coordinates}  spot_coordinates.tsv
        [ "${read_to_spot}"      = "read_to_spot.tsv"         ] || cp ${read_to_spot}      read_to_spot.tsv

        # Build barcode whitelist CSV for Panorama-seq (one row per unique spot:
        # true_barcode,x,y — deduplicated from read_to_spot.tsv + spot_coordinates.tsv)
        read_to_spot_csv.py \\
            --read_to_spot  read_to_spot.tsv \\
            --spot_coords   spot_coordinates.tsv \\
            --out           barcodes_coords.csv

        write_metrics_json.py \\
            --grid      spot_coordinates.tsv \\
            --mutations barcode_mutations.tsv \\
            --mapping   read_to_spot.tsv \\
            --reads-requested ${total_reads} \\
            --out       summary_metrics.json
        """
}
