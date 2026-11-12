/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GAWK_POOL — Pool the reads of every sample of one [chr, class] for a pooled reconstruction:
    one BED per sample with read names <sample>.<read>, and the same rows in one file.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GAWK_POOL {
    tag "$meta.id:$meta.chr:$meta.class"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0':
        'biocontainers/gawk:5.3.0' }"

    input:
    // INFO: meta.samples holds the sample ids in the order of beds; every sample's file has the same
    // INFO: name (fusions.free.bed), so each one is staged in a numbered directory
    tuple val(meta), path(beds, stageAs: 'in??/*')

    output:
    tuple val(meta), path("samples/*.bed"), emit: tagged // INFO: one BED per sample, named <sample>.bed
    tuple val(meta), path("*.pooled.bed") , emit: pooled // INFO: the same rows in one file, in sample order
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chr}.${meta.class}"
    def files = beds instanceof List ? beds : [beds]
    if (files.size() != meta.samples.size()) {
        error "ERROR: GAWK_POOL got ${files.size()} files for ${meta.samples.size()} samples (${meta.samples.join(', ')})"
    }
    // INFO: iso-segment names reads per sample chunk (R<n>@<chunk>_<chr>), so two samples share read names;
    // INFO: the sample prefix keeps them apart and stays in the name of the model a read represents
    def tag = [meta.samples, files].transpose().collect { sample, bed ->
        "gawk -v s='${sample}' 'BEGIN { FS = OFS = \"\\t\" } { \$4 = s \".\" \$4; print }' ${bed} > 'samples/${sample}.bed'"
    }.join('\n    ')
    def tagged = meta.samples.collect { sample -> "'samples/${sample}.bed'" }.join(' ')
    """
    mkdir samples
    ${tag}
    cat ${tagged} > ${prefix}.pooled.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | sed 's/^.*awk version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.chr}.${meta.class}"
    """
    mkdir samples
    touch ${meta.samples.collect { sample -> "'samples/${sample}.bed'" }.join(' ')} ${prefix}.pooled.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: stub
    END_VERSIONS
    """
}
