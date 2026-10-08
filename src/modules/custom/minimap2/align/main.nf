process MINIMAP2_ALIGN {
    tag "$meta.id chunk $meta.chunk"
    label 'process_medium_fast'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:b411340b52d82a9c276d87c7a3dcffc880be762f-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:b411340b52d82a9c276d87c7a3dcffc880be762f-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta1), path(reference)
    tuple val(meta2), path(splice_scores)
    tuple val(meta3), path(junc_bed)

    output:
    // INFO: minimap2 pipes into samtools sort, so no SAM is written and nothing downstream deletes this task's outputs (-resume keeps it cached)
    tuple val(meta), path("*.bam")                       , emit: bam
    tuple val(meta), path("*.bai")                       , emit: bai
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def singleton = meta.singleton ? ".singleton" : ""
    def bam = "${meta.id}.${meta.chunk}${singleton}.bam"
    def spsc = splice_scores ? "--spsc=${splice_scores}" : ''
    def junc = task.ext.use_junc_bed ? "--junc-bed ${junc_bed}" : ''
    """
    minimap2 \\
        $args \\
        $spsc \\
        $junc \\
        -t $task.cpus \\
        ${reference} \\
        ${reads} \\
        | samtools sort -@ ${task.cpus} -o $bam -

    samtools index -@ ${task.cpus} $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chunk}${singleton}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
