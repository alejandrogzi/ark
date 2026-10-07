/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTX_ORIENT — Orient a primer-free, partly reversed FASTA/FASTQ by its tails: reads with
    only a 5' polyT are reverse-complemented; nothing is dropped or renamed.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTX_ORIENT {
    tag "$meta.id"
    label 'process_low'

    container 'ghcr.io/alejandrogzi/isotools:v0.0.45'

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("*.oriented.fast*"), emit: reads
    path "*.orient.tsv"                      , emit: report
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ext = (fastx.name =~ /\.fast[aq](\.gz)?$/)[0][0] // INFO: same format and compression out
    """
    iso-fastx orient \\
        $args \\
        --fastx $fastx \\
        --output ${meta.id}.oriented${ext} \\
        --prefix ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-fastx: \$( iso-fastx --version | sed 's/iso-fastx //g' )
    END_VERSIONS
    """

    stub:
    def ext = (fastx.name =~ /\.fast[aq](\.gz)?$/)[0][0]
    """
    cp $fastx ${meta.id}.oriented${ext}
    touch ${meta.id}.orient.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-fastx: stub
    END_VERSIONS
    """
}
