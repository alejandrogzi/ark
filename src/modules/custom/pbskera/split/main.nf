/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PBSKERA_SPLIT — Split BAM files using skera
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PBSKERA_SPLIT {
    tag "$meta.id"
    label 'process_extreme'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbskera:1.4.0--hdfd78af_0' :
        'biocontainers/pbskera:1.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path primers

    output:
    tuple val(meta), path("*skera.bam")               , emit: bam
    tuple val(meta), path("*.skera.bam.pbi")          , emit: pbi
    tuple val(meta), path("*.non_passing.bam")        , emit: non_passing_bam
    tuple val(meta), path("*.found_adapters.csv.gz")  , emit: found_adapters
    tuple val(meta), path("*.summary.csv")            , emit: summary
    tuple val(meta), path("*.summary.json")           , emit: summary_json
    tuple val(meta), path("*.ligations.csv")          , emit: ligations
    tuple val(meta), path("*.read_lengths.csv")       , emit: read_lengths
    path  "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
   """
    skera \\
        split \\
        $bam \\
        $primers \\
        ${prefix}.skera.bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skera: \$( skera --version | head -n 1 | sed 's/skera split //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.skera.bam
    touch ${prefix}.non_passing.bam
    touch *summary.csv
    touch *ligations.csv
    touch *read_lengths.csv
    touch *found_adapters.csv.gz
    touch *summary.json
    touch *skera.bam.pbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skera: \$( skera --version | head -n 1 | sed 's/skera split //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """
}
