process ISOTOOLS_ORPHAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/isotools:latest' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta1), path(reference)
    tuple val(meta2), path(splice_scores)

    output:
    tuple val(meta), path("*/*.hq.bed")             , optional: true, emit: hq
    tuple val(meta_scraps), path("*/*.scraps.bed")  , optional: true, emit: scraps
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def scores    = splice_scores ? "--splicing-scores $splice_scores" : ''

    meta_scraps = meta.clone()
    meta_scraps.id = meta.id + '.scraps'
    """
    # INFO: specific for isopipe -> input is bed12+n
    cut -f1-12 ${bed} > tmp.bed

    iso-orphan \\
        $args \\
        $scores \\
        --ref $reference \\
        --query tmp.bed \\
        --all \\
        --threads ${task.cpus} \\
        --prefix ${prefix} 

    # INFO: rehydrate bed12+n
    if [ ! -s orphans/${prefix}.hq.bed ]; then
        rm orphans/${prefix}.hq.bed
    else
        grep -f \\
        <(cut -f4 orphans/${prefix}.hq.bed) ${bed} \\
        > tmp.hq.bed && \\
        sort -k1,1 -k2,2n tmp.hq.bed >| orphans/${prefix}.hq.bed
    fi

    if [ ! -s orphans/${prefix}.scraps.bed ]; then
        rm orphans/${prefix}.scraps.bed
    else
        grep -f \\
        <(cut -f4 orphans/${prefix}.scraps.bed) ${bed} \\
        > tmp.scraps.bed && \\
        sort -k1,1 -k2,2n tmp.scraps.bed >| orphans/${prefix}.scraps.bed
    fi

    rm tmp*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-orphan: \$( iso-orphan --version | sed 's/iso-orphan //g' )
    END_VERSIONS
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    touch orphans/${prefix}.hq.bed
    touch orphans/${prefix}.scraps.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-orphan: \$( iso-orphan --version | sed 's/iso-orphan //g' )
    END_VERSIONS
    """
}
