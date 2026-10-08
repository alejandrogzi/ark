process SPLICEAI_DERIVE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/isox-rs:latest' }"

    input:
    tuple val(meta), path(genome)
    tuple val(meta1), path(annotation)
    tuple val(meta2), path(spliceai)

    output:
    tuple val(meta), path("*.derived.tsv"), emit: scores
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    splicing \\
        $args \\
        -t ${task.cpus} \\
        --bigwig-dir ${spliceai} \\
        --sequence ${genome} \\
        --regions ${annotation} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splicing: \$( splicing --version | sed 's/splicing //g' )
    END_VERSIONS
    """

    stub:
    """
    touch *.derived.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splicing: \$( splicing --version | sed 's/splicing //g' )
    END_VERSIONS
    """
}
