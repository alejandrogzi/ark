process ISOTOOLS_CLASSIFY_INTRON {
    tag "$meta.id"
    label 'process_high_low_memory_fast'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/isotools:latest' }"

    input:
    tuple val(meta), path(reads), path(intronic)
    tuple val(meta1), path(genome)
    tuple val(meta2), path(annotation)
    tuple val(meta3), path(repeats)
    tuple val(meta4), path(bigwigs)

    output:
    tuple val(meta), path("*.tsv")      , optional: true, emit: tsv
    tuple val(meta), path("*.introns_track.bed"), optional: true, emit: track
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def spliceai  = bigwigs ? "--bigwig $bigwigs" : ''
    def repeats   = repeats ? "--repeats $repeats" : ''
    def iic       = intronic && intronic.size() > 0 ? "--iic $intronic" : ''
    """
    iso-classify intron \\
        --input $reads \\
        --sequence $genome \\
        --reference $annotation \\
        --prefix ${prefix} \\
        $spliceai \\
        $repeats \\
        $iic \\
        --outdir . \\
        $args

    # INFO: isotools names the --intron-track file '<prefix>.introns_track,bed' (comma typo upstream)
    if [ -s "${prefix}.introns_track,bed" ]; then
        mv "${prefix}.introns_track,bed" ${prefix}.introns_track.bed
    else
        rm -f "${prefix}.introns_track,bed"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-classify: \$( iso-classify --version | sed 's/iso-classify //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.introns_track.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-classify: \$( iso-classify --version | sed 's/iso-classify //g' )
    END_VERSIONS
    """
}
