process COLLAPSE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/isox-rs:latest' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.collapsed.bed"), emit: collapsed
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def chr = meta.chr ? ".${meta.chr}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}${chr}"
    """
    collapse run \\
        $args \\
        --bed ${bed} \\
        --extend

    sort -k1,1 -k2,2n collapsed/collapsed.bed > ${prefix}.collapsed.bed
    rm -rf collapsed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collapse: \$( collapse --version | sed 's/collapse //g' )
    END_VERSIONS
    """

    stub:
    """
    touch *.collapsed.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collapse: \$( collapse --version | sed 's/collapse //g' )
    END_VERSIONS
    """
}
