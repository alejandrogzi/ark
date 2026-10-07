/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHAIN_COLLAPSE — Collapse segmented reads of one [sample, chr, class] into transcript
    models by intron chain; models keep a real read line and carry #CN<support>.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CHAIN_COLLAPSE {
    tag "$meta.id:$meta.chr:$meta.class"
    label 'process_medium'

    container 'ghcr.io/alejandrogzi/isox-rs:v2.1.0'

    input:
    tuple val(meta), path(beds, stageAs: 'beds/*')
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.models.bed")                         , emit: models
    tuple val(meta), path("*.support.tsv"), path("*.counts.tsv")    , emit: support
    tuple val(meta), path("*.members.tsv.gz"), path("*.excluded.bed"), emit: members
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chr}.${meta.class}"
    def queries = (beds instanceof List ? beds : [beds]).join(',')
    """
    collapse chain \\
        $args \\
        --bed $queries \\
        --ref $reference \\
        --prefix $prefix \\
        --threads ${task.cpus}

    gzip -f ${prefix}.members.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collapse: \$( collapse --version | sed 's/collapse //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chr}.${meta.class}"
    """
    touch ${prefix}.models.bed ${prefix}.support.tsv ${prefix}.counts.tsv ${prefix}.excluded.bed
    echo | gzip > ${prefix}.members.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collapse: stub
    END_VERSIONS
    """
}
