/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RESCUE_FLNC_FROM_ARTIFACTS — rewrites bucket files without the rows it moves: the seven
    VEREDICT buckets + NMD of one sample@chromosome, or (as RESCUE_FLNC_FROM_ORPHANS) the orphan
    finder's hq/scraps of one sample. rt/artifact rows whose RT/ARTIFACT introns are all annotated
    are re-bucketed by their other flaws; rows of the rescue categories with exactly a reference
    CDS go to pass (hq) tagged #RS<origin>.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RESCUE_FLNC_FROM_ARTIFACTS {
    tag "$meta.id"
    label 'process_single'

    // INFO: rescue_flnc.py (src/bin) is stdlib-only; any python3 image works
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(buckets, stageAs: 'veredict/*'), path(introns, stageAs: 'in/introns.tsv'), path(nmd, stageAs: 'in/rescue.nmd.bed')
    tuple val(meta2), path(reference)

    output:
    // INFO: each input bucket comes back under its own label; pass is `hq` for the orphan finder
    tuple val(meta), path("rescue/*.{pass,hq}.bed")         , emit: pass
    tuple val(meta), path("rescue/*.trash.bed")             , optional: true, emit: trash
    tuple val(meta), path("rescue/*.retentions.bed")        , optional: true, emit: retentions
    tuple val(meta), path("rescue/*.truncations.bed")       , optional: true, emit: truncations
    tuple val(meta), path("rescue/*.intraprimming.bed")     , optional: true, emit: intraprimming
    tuple val(meta), path("rescue/*.rt.bed")                , optional: true, emit: rt
    tuple val(meta), path("rescue/*.artifacts.bed")         , optional: true, emit: artifacts
    tuple val(meta), path("rescue/*.nmd.bed")               , optional: true, emit: nmd
    tuple val(meta), path("rescue/*.scraps.bed")            , optional: true, emit: scraps
    tuple val(meta), path("rescue/*.rescued_artifacts.tsv") , emit: audit
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def beds   = [ buckets ? 'veredict/*.bed' : '', nmd ? 'in/rescue.nmd.bed' : '' ].findAll { it }.join(' ')
    def table  = introns ? '--introns in/introns.tsv' : ''
    """
    rescue_flnc.py artifacts \\
        --buckets ${beds} \\
        $table \\
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
    for bucket in pass trash retentions truncations intraprimming rt artifacts nmd; do
        touch rescue/${prefix}.\$bucket.bed
    done
    printf 'id\\toriginal_id\\torigin\\tdestination\\trule\\treference\\n' > rescue/${prefix}.rescued_artifacts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rescue_flnc: stub
    END_VERSIONS
    """
}
