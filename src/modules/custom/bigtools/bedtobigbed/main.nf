process BEDTOBIGBED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bigtools:0.5.6--hc1c3326_1':
        'biocontainers/bigtools:0.5.6--hc1c3326_1' }"

    input:
    tuple val(meta), path(bed)
    path chrom_sizes
    path autosql

    output:
    tuple val(meta), path("*.bb"), optional: true, emit: bigbed
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def as = autosql ? "--autosql $autosql" : ''
    // bigtools 0.5.6 binary-searches chromosome boundaries once a BED is
    // >= 200 MB. A short chromosome between two long ones is skipped, then
    // the writer reports "File is not sorted" and panics. -p no reads in order.
    """
    bigtools bedtobigbed \\
        -p no \\
        $args \\
        $as \\
        $bed \\
        $chrom_sizes \\
        ${prefix}.bb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigtools: \$(bigtools --version | sed -e "s/bigtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigtools: \$(bigtools --version | sed -e "s/bigtools v//g")
    END_VERSIONS
    """
}
