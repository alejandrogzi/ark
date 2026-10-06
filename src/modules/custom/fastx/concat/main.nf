/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTX_CONCAT — Byte-concatenate staged read files into one pooled file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTX_CONCAT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data' :
        'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("${meta.outfile}"), emit: reads
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def outfile = meta.outfile
    """
    # INFO: channel order is not guaranteed, so sort here for a deterministic pooled file
    # INFO: (plain cat is format-safe for .gz: concatenated gzip members stay valid)
    cat \$(printf '%s\\n' ${files} | sort -u | tr '\\n' ' ') > ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -n 1 | sed 's/cat (GNU coreutils) //')
        sort: \$(sort --version 2>&1 | head -n 1 | sed 's/sort (GNU coreutils) //')
    END_VERSIONS
    """

    stub:
    def outfile = meta.outfile
    """
    touch ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -n 1 | sed 's/cat (GNU coreutils) //')
        sort: \$(sort --version 2>&1 | head -n 1 | sed 's/sort (GNU coreutils) //')
    END_VERSIONS
    """
}
