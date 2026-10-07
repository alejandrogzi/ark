/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTX_INSPECT — Classify a FASTA/FASTQ from its first reads: subreads, ccs (primers),
    fl (tails), flnc, mixed (partly reversed), clustered, ambiguous or empty.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTX_INSPECT {
    tag "$meta.id"
    label 'process_single'

    container 'ghcr.io/alejandrogzi/isotools:v0.0.45'

    input:
    tuple val(meta), path(fastx)
    path primers // INFO: optional user primers ([] when absent)

    output:
    tuple val(meta), path(fastx), path("*.inspect.tsv"), emit: inspected
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def user_primers = primers ? "--primers ${primers}" : ''
    """
    iso-fastx inspect \\
        $args \\
        --fastx $fastx \\
        --prefix ${meta.id} \\
        $user_primers

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-fastx: \$( iso-fastx --version | sed 's/iso-fastx //g' )
    END_VERSIONS
    """

    stub:
    """
    printf 'file\\treads\\tstate\\tkit\\nX\\t0\\tfl\\tnone\\n' > ${meta.id}.inspect.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-fastx: stub
    END_VERSIONS
    """
}
