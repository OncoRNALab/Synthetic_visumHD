process SHARD_MOLECULES {
    input:
        path molecules
        val  n_shards

    output:
        path "shard_*.tsv", emit: shards

    script:
        """
        shard_molecules.py \\
            --molecules ${molecules} \\
            --n-shards  ${n_shards} \\
            --prefix    shard_
        """
}
