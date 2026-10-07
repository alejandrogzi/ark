process SAMTOOLS_FASTA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // INFO: FASTA, not FASTQ: nothing downstream reads qualities, so they would only double the bytes.
    // INFO: -0 takes the unpaired reads (all PacBio reads); .gz makes it multi-threaded BGZF.
    """
    samtools fasta -@ ${task.cpus} -0 ${prefix}.fasta.gz ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | sed 's/samtools //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version | head -n 1 | sed 's/samtools //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """
}
