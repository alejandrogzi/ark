process ISOTOOLS_SEGMENT {
    tag "$meta.id chunk $meta.chunk"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/isotools:latest' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.hq.bam")                    , optional: true, emit: hq_bam
    tuple val(meta), path("*.hq.bed")                    , optional: true, emit: hq_bed
    tuple val(meta), path("*.lq.bam")                    , optional: true, emit: lq_bam
    tuple val(meta), path("*.lq.bed")                    , optional: true, emit: lq_bed
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def singleton = meta.singleton ? "--singleton" : ""
    def keep_temp = task.ext.keep_temp ?: false
    def chunk     = meta.chunk.toInteger()
    """
    iso-segment \\
        $args \\
        --bam $bam \\
        --prefix ${prefix} \\
        --batch ${chunk} \\
        $singleton

    # INFO: -X writes one hq / lq BED per chromosome; both are renamed per chunk (and singleton class) so that
    # INFO: published files of different chunks do not overwrite each other. lq is removed unless kept.
    suffix=".${chunk}"
    if [ "${singleton}" ]; then suffix=".${chunk}.singleton"; fi
    for kind in hq lq; do
        for f in *${prefix}.\$kind.bed; do
            [ -f "\$f" ] || continue
            if [ \$kind = lq ] && [ "${keep_temp}" = false ]; then
                rm "\$f"
            else
                mv "\$f" "\${f%.\$kind.bed}\${suffix}.\$kind.bed"
            fi
        done
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-segment: \$( iso-segment --version | sed 's/iso-segment //g' )
    END_VERSIONS
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    touch *${prefix}*.hq.bam
    touch *${prefix}*.hq.bed
    touch *${prefix}*.lq.bam
    touch *${prefix}*.lq.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-segment: \$( iso-segment --version | sed 's/iso-segment //g' )
    END_VERSIONS
    """
}
