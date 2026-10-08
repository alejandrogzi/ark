/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTX_RENAME — Give every read a PacBio CCS name (<movie>/<n>/ccs) before lima.
    lima 26.2.1 hangs on FASTA/FASTQ whose names are not movie/zmw/ccs (SRA renames reads to
    SRRxxx.N), even with --per-read. The movie is the sanitized sample id, so pooled samples
    never share read names.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTX_RENAME {
    tag "$meta.id"
    label 'process_single'

    container 'biocontainers/lima:26.2.1--h9ee0642_0' // INFO: the lima image already ships awk and gzip

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("*.ccs_names.fast*"), emit: reads
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def movie = meta.sample_id.replaceAll(/[^A-Za-z0-9]/, '_')
    def ext = (fastx.name =~ /\.fast[aq]/)[0]
    """
    zcat -f $fastx \\
        | awk -v m=$movie 'NR == 1 { fq = substr(\$0, 1, 1) == "@" }
            fq && NR % 4 == 1 { printf "@%s/%d/ccs\\n", m, ++n; next }
            !fq && /^>/ { printf ">%s/%d/ccs\\n", m, ++n; next }
            { print }' \\
        | gzip > ${meta.id}.ccs_names${ext}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$( awk --version 2>&1 | head -n 1 )
    END_VERSIONS
    """

    stub:
    def ext = (fastx.name =~ /\.fast[aq]/)[0]
    """
    cp $fastx ${meta.id}.ccs_names${ext}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: stub
    END_VERSIONS
    """
}
