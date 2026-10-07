process ISOTOOLS_TRUNCATION_DETECTOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/isotools:v0.0.45' }"

    input:
    tuple val(meta), path(bed), path(_), path(_)

    output:
    tuple val(meta), path("*.tsv")       , optional: true, emit: descriptor
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    iso-utr \\
        $args \\
        --ref $bed \\
        --query $bed \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        -O cds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-utr: \$( iso-utr --version | sed 's/iso-utr //g' )
    END_VERSIONS
    """

    stub:
    """
    touch *.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-utr: \$( iso-utr --version | sed 's/iso-utr //g' )
    END_VERSIONS
    """
}
