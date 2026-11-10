/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RESCUE_FLNC_FROM_UNPREDICTED_ORFS — models that went into ORF calling minus models with a
    row in any output (VEREDICT buckets, NMD) -> <prefix>.unpredicted.bed. With the reference,
    an unpredicted model that carries a reference CDS intact takes it and is rescued to pass.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RESCUE_FLNC_FROM_UNPREDICTED_ORFS {
    tag "$meta.id"
    label 'process_single'

    // INFO: rescue_flnc.py (src/bin) is stdlib-only; any python3 image works
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(models, stageAs: 'in/models.bed'), path(predicted, stageAs: 'in/predicted.bed'), path(outputs, stageAs: 'out/*')
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("rescue/*.unpredicted.bed")           , emit: unpredicted
    tuple val(meta), path("rescue/*.rescued_unpredicted.bed")   , emit: rescued
    tuple val(meta), path("rescue/*.rescued_unpredicted.tsv")   , emit: audit
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def xorf      = predicted ? "--predicted ${predicted}" : ''
    def went_out  = outputs ? '--out out/*' : '' // INFO: one staged file arrives as a path, not a list
    """
    rescue_flnc.py unpredicted \\
        --models ${models} \\
        $xorf \\
        $went_out \\
        --reference ${reference} \\
        --prefix ${prefix} \\
        --outdir rescue \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rescue_flnc: \$(rescue_flnc.py --version | sed 's/rescue_flnc.py //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p rescue
    touch rescue/${prefix}.unpredicted.bed rescue/${prefix}.rescued_unpredicted.bed
    printf 'id\\toriginal_id\\torigin\\tdestination\\trule\\treference\\n' > rescue/${prefix}.rescued_unpredicted.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rescue_flnc: stub
    END_VERSIONS
    """
}
