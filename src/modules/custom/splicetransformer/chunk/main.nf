// INFO: ported from hillerlab/core modules/nextflow/splicetransformer/chunk (container pinned, stub fixed)
// INFO: splits the genome into tmp.<contig>.chunk.<n>.fa(.gz) shards with --flank-size flanks on both ends
process SPLICETRANSFORMER_CHUNK {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/hillerlab/splicetransformer:0.0.5' }"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("chunks/*.fa"),    optional: true, emit: fasta
    tuple val(meta), path("chunks/*.fa.gz"), optional: true, emit: fasta_gz
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    splicetransformer chunk \\
        $args \\
        -t ${task.cpus} \\
        --sequence ${genome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splicetransformer chunk: \$( splicetransformer chunk --version | sed 's/splicetransformer-chunk //g' )
    END_VERSIONS
    """

    stub:
    def suffix = task.ext.args?.contains('--gzip') ? 'fa.gz' : 'fa'
    """
    mkdir -p chunks
    touch chunks/tmp.${meta.id}.chunk.0.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splicetransformer chunk: stub
    END_VERSIONS
    """
}
