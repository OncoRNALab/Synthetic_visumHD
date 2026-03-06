process RUN_SRT_SIM {
    input:
        val seed
        val n_genes
        val n_spots
        val feature_h5        // empty string when not provided
        val required_genes    // comma-separated gene symbols to force-include, or ''
        val synth_sub_rate    // per-base substitution rate during array synthesis
        val synth_ins_rate    // per-base insertion rate during array synthesis
        val synth_del_rate    // per-base deletion rate during array synthesis

    output:
        path "expression_matrix.tsv", emit: expression_matrix
        path "spot_coordinates.tsv",  emit: spot_coordinates

    script:
        def h5_arg  = (feature_h5     && feature_h5     != '') ? "--feature-h5=${feature_h5}"         : ""
        def req_arg = (required_genes && required_genes != '') ? "--required-genes=${required_genes}"  : ""
        """
        Rscript ${projectDir}/scripts/run_sccigen_sim.R \\
            --seed=${seed} \\
            --n-genes=${n_genes} \\
            --n-spots=${n_spots} \\
            ${h5_arg} \\
            ${req_arg} \\
            --synth-sub-rate=${synth_sub_rate} \\
            --synth-ins-rate=${synth_ins_rate} \\
            --synth-del-rate=${synth_del_rate} \\
            --output-matrix=expression_matrix.tsv \\
            --output-coords=spot_coordinates.tsv
        """
}
