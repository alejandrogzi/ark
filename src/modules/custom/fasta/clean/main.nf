/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTA_CLEAN — Uncompress a (gzipped) genome FASTA and keep only the first word of each
    header. NCBI/Ensembl headers carry descriptions ('>chr12  AC:CM000674.2 ...'); minimap2
    cuts them but xloci keys sequences by the whole line and then misses every chromosome.
    IUPAC ambiguity codes (GRCh38 has a few) become N/n, case kept: xloci 0.0.6 panics on them.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTA_CLEAN {
    tag "$fasta"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data' :
        'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8' }"

    input:
    tuple val(meta), path(fasta, stageAs: 'input/*')

    output:
    tuple val(meta), path("$clean"), emit: fasta
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    clean = fasta.name - '.gz' // INFO: same name as a plain gunzip, so downstream file names do not change
    """
    zcat -f $fasta \\
        | awk '/^>/ { print \$1; next } { gsub(/[RYKMSWBDHV]/, "N"); gsub(/[rykmswbdhv]/, "n"); print }' \\
        > $clean

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version | head -n1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    clean = fasta.name - '.gz'
    """
    touch $clean

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: "stub"
    END_VERSIONS
    """
}
