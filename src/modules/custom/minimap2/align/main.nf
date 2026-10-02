process MINIMAP2_ALIGN {
    tag "$meta.id chunk $meta.chunk"
    label 'process_medium_fast'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.31--h118bc1c_0' :
        'biocontainers/minimap2:2.31--h118bc1c_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta1), path(reference)
    tuple val(meta2), path(splice_scores)
    tuple val(meta3), path(junc_bed)

    output:
    // WARN: not optional. SAMTOOLS_BAM deletes the SAM; on -resume a missing required output reruns this task instead of silently emitting nothing
    tuple val(meta), path("*.sam")                       , emit: sam
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def keep_temp = task.ext.keep_temp ?: false
    
    def singleton = meta.singleton ? ".singleton" : ""
    def sam = "${meta.id}.${meta.chunk}${singleton}.sam"
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
        -o $sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chunk}${singleton}"
    """
    touch ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
