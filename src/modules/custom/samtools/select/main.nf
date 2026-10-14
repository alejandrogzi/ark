/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SAMTOOLS_SELECT_READS — Keep (ext.mode 'keep') or drop (ext.mode 'drop') every record
    (primary, supplementary, secondary) of the reads named in a list, then sort and index.
    The sort is needed: iso-cigar BAMs are not strictly coordinate-sorted (edited starts).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SAMTOOLS_SELECT_READS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(bam), path(bai), path(names)

    output:
    // INFO: optional: a 'keep' selection with no records emits nothing (an empty BAM would still be segmented)
    tuple val(meta), path("*.${suffix}.bam"), path("*.${suffix}.bam.bai"), optional: true, emit: bam
    path "versions.yml",                                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode = task.ext.mode ?: 'keep'
    if (!(mode in ['keep', 'drop'])) {
        error "SAMTOOLS_SELECT_READS: ext.mode must be 'keep' or 'drop', got '${mode}'"
    }
    suffix = mode == 'keep' ? 'selected' : 'filtered'
    def out = "${task.ext.prefix ?: bam.baseName}.${suffix}.bam"
    def select = mode == 'keep'
        ? "-N ${names} -o selected.unsorted.bam"
        : "-N ${names} -U selected.unsorted.bam -o /dev/null"
    """
    samtools view -@ ${task.cpus} -b ${select} ${bam}

    if [ "${mode}" = "keep" ] && [ "\$(samtools view -c selected.unsorted.bam)" -eq 0 ]; then
        rm selected.unsorted.bam
    else
        samtools sort -@ ${task.cpus} -o ${out} selected.unsorted.bam
        rm selected.unsorted.bam
        samtools index -@ ${task.cpus} ${out}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    suffix = (task.ext.mode ?: 'keep') == 'keep' ? 'selected' : 'filtered'
    """
    touch ${task.ext.prefix ?: bam.baseName}.${suffix}.bam
    touch ${task.ext.prefix ?: bam.baseName}.${suffix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
