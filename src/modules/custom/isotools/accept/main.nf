/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ISOTOOLS_ALIGN_ACCEPT — Decide which realigned fragment reads replace their pass-1
    records (iso-align accept mode, primary records only). A read is accepted when:
      1. AS(realigned) > AS(pass-1), strictly (pass-1 BAM: iso-cigar drops AS on edited records);
      2. realigned identity >= --min-identity (98, iso-segment formula; no NM = reject);
      3. realigned clip <= --max-clip (20 bp) on both ends.
    Emits the accepted names (sorted) and a per-read report (AS, identity, clips, decision).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ISOTOOLS_ALIGN_ACCEPT {
    tag "$meta.id"
    label 'custom_process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/isotools:latest' }"

    input:
    tuple val(meta), path(pass1_bam), path(pass1_bai), path(realigned_bam), path(realigned_bai)

    output:
    tuple val(meta), path("*.accepted.txt"), emit: names
    tuple val(meta), path("*.accept.tsv"),   emit: report
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # INFO: accept mode (--realigned-bam) exists from iso-align 0.0.5; fail clearly on an older isotools image
    help=\$(iso-align --help 2>&1 || true)
    if [[ "\$help" != *--realigned-bam* ]]; then
        echo "ERROR: fragment acceptance needs iso-align >= 0.0.5 (found: \$(iso-align --version))" >&2
        exit 1
    fi

    iso-align \\
        $args \\
        --threads ${task.cpus} \\
        --bam $pass1_bam \\
        --realigned-bam $realigned_bam \\
        --output ${prefix}.accepted.txt \\
        --report ${prefix}.accept.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-align: \$( iso-align --version | sed 's/iso-align //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.accepted.txt
    touch ${prefix}.accept.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iso-align: \$( iso-align --version | sed 's/iso-align //g' )
    END_VERSIONS
    """
}
