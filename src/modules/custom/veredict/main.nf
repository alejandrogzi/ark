process VEREDICT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'ghcr.io/alejandrogzi/isox-py:v2.1.0' }"

    input:
    tuple val(meta), path(bed), path(_), path(orfs), path(retentions), path(intraprimming), path(truncations)

    output:
    tuple val(meta), path("veredict/*.pass.bed")          , optional: true, emit: pass
    tuple val(meta), path("veredict/*.trash.bed")         , optional: true, emit: trash
    tuple val(meta), path("veredict/*.retentions.bed")    , optional: true, emit: retentions
    tuple val(meta), path("veredict/*.truncations.bed")   , optional: true, emit: truncations
    tuple val(meta), path("veredict/*.intraprimming.bed") , optional: true, emit: intraprimming
    tuple val(meta), path("veredict/*.rt.bed")            , optional: true, emit: rt
    tuple val(meta), path("veredict/*.artifacts.bed")     , optional: true, emit: artifacts
    env(ADDITIONAL_BED_COLUMNS)                                           , emit: additional_bed_columns
    path "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    veredict \\
        $args \\
        --reads $bed \\
        --orfs $orfs \\
        --retentions $retentions \\
        --intraprimming $intraprimming \\
        --truncations $truncations \\
        --outdir veredict \\
        --prefix $prefix

    ADDITIONAL_BED_COLUMNS=\$(head -1 veredict/${prefix}.pass.bed | awk -F'\t' '{print NF}')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        veredict: \$(veredict --version 2>&1 | sed 's/veredict //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch veredict/${prefix}.pass.bed 
    touch veredict/${prefix}.trash.bed
    touch veredict/${prefix}.retentions.bed
    touch veredict/${prefix}.truncations.bed
    touch veredict/${prefix}.intraprimming.bed
    touch veredict/${prefix}.rt.bed
    touch veredict/${prefix}.artifacts.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        veredict: \$(veredict --version 2>&1 | sed 's/veredict //')
    END_VERSIONS
    """
}
