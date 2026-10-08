process XLOCI_INTRON {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/xloci:latest' }"

    input:
    tuple val(_), path(genome)
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fa")      , optional: true, emit: fasta
    tuple val(meta), path("*.tsv")     , optional: true, emit: tsv
    path  "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    xloci \\
        $args \\
        -f intron \\
        -o . \\
        -s $genome \\
        -r $reads \\
        --unmask \\
        -t $task.cpus \\
        --prefix ${prefix}

    # INFO: xloci writes one row per read and intron (17M rows for 59k introns on a giraffe chromosome);
    # INFO: the name is the intron's coordinates, and intronIC only needs each intron once
    for tsv in *.tsv; do
        [ -e "\$tsv" ] || continue
        awk -F'\t' '!seen[\$1]++' "\$tsv" > "\$tsv.unique" && mv "\$tsv.unique" "\$tsv"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xloci: \$( xloci --version | head -n 1 | sed 's/xloci //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xloci: \$( xloci --version | head -n 1 | sed 's/xloci //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """
}
